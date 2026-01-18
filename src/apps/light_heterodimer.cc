#include "basic/types.hh"
#include "scene/Residue.hh"
#include "scene/Conformation.hh"
#include "scene/Pose.hh"
#include "basic/macros.hh"
#include "utils/hash_util.hh"
#include "sampling/FragMover.hh"
#include "sampling/mcmc_protocols.hh"
#include "scoring/ScoreFunction.hh"
#include "utils/utils.hh"
#include "basic/assert.hh"
#include "utils/random_util.hh"
#include "basic/assert.hh"
#include "utils/dssp_util.hh"
#include "utils/math_util.hh"

#include "metric/LigandNeighbors.hh"
#include "metric/SidechainNeighbors.hh"

#include "apps/args.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <ctime>
#include <random>
#include <algorithm>
#include <iomanip>
#include <chrono>

#include <parallel_hashmap/phmap.h>
#include <memory>

#include <chrono>
typedef std::chrono::high_resolution_clock Clock;

using namespace scene;
using namespace scoring;
using namespace sampling;
using namespace utils;

// for monomers only

int main(int argc, char *argv[])
{
    Options opts;
    opts.parse_args(argc, argv);

    std::cout << "Load and prepare the pdbs: " << std::endl
              << "    ==>  Chain1 len: " << opts.chain1_len << std::endl
              << "    ==>  Chain1 pdb: " << opts.chain1_pdb << std::endl
              << "    ==>  Chain2 len: " << opts.chain2_len << std::endl
              << "    ==>  Chain2 pdb: " << opts.chain2_pdb << std::endl;

    Pose pose(opts.chain1_len);
    PoseOP target_pose = std::make_shared<Pose>(opts.chain2_len);

    pose.load_pdb(opts.chain1_pdb, 1, false, true, true);
    target_pose->load_pdb(opts.chain2_pdb, 1, true, true, true);
    std::string ss1 = utils::get_dssp_from_pose(pose);
    std::string ss2 = utils::get_dssp_from_pose(*target_pose);
    std::cout << "The secondary structure of chain1 is automatically determined to be:" << std::endl
              << "    ==> " << ss1 << std::endl;
    std::cout << "The secondary structure of chain2 is automatically determined to be:" << std::endl
              << "    ==> " << ss2 << std::endl;
    pose.set_dssp(ss1);
    target_pose->set_dssp(ss2);

    // no need to recenter the pose
    // the target has to be fixed 
    // as the hash table stores the relative position of the rif residues
    // but to make sure this is OK, it should not matter
    if (true) {
        Vec target_center = target_pose->conformation().center_vec();
        EigenXform new_root = target_pose->conformation().root();
        new_root.translation() = new_root.translation() - target_center;
        target_pose->conformation().set_root(new_root);
        target_pose->conformation().update_coordinates();
    }
    pose.set_target_pose(target_pose);


    std::cout << "Prepare the score function ...\n";
    std::cout << "Use secondary informtion for rpx scoring: " << opts.use_ss << std::endl;
    ScoreFunction sfxn;
    sfxn.regist_method("rpx",ScoreMethodOP(new RpxScoreTargetMethod(opts.rpx_db, opts.rpx_cart_resl, opts.rpx_ang_resl, opts.use_ss)));
    sfxn.regist_method("clash",ScoreMethodOP(new VoxelClashScoreMethod(*target_pose, "BB_CB", 1.5f*opts.CB_swelling_factor,0.25f,false,"PROTEIN")));
    sfxn.get_method("clash")->set_score_type(TARGET_CONTEXT_CLASH);
    sfxn.get_method("clash")->set_weight(opts.context_clash_weight);
    // std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash"))->visualize_voxel_grid("clash_grid.pdb");
    pose.energy_manager().add_energy_onebody(TARGET_CONTEXT_CLASH);

    std::cout << "Loading privileged motif hash tables " << std::endl;
    sfxn.regist_method("privileged_motif", ScoreMethodOP(new PrivilegedMotifScoreMethod(opts.privileged_motif_hash_tables, 
                                                                                        opts.privileged_motif_hash_cart_resl,
                                                                                        opts.privileged_motif_hash_ang_resl)));
    sfxn.get_method("privileged_motif")->set_weight(opts.privileged_motif_weight);
    sfxn.get_method("privileged_motif")->set_score_type(PRIVILEGED_MOTIF);
    pose.energy_manager().add_energy_onebody(PRIVILEGED_MOTIF);
    // the AzoF rif is relative to the stub of the 1st residue:
    std::static_pointer_cast<PrivilegedMotifScoreMethod>(sfxn.get_method("privileged_motif"))->set_relative_pos(target_pose->stub(1));


    metric::SidechainNeighbors sidechain_neighbors;

    auto t1 = Clock::now();

    Size total_len = opts.chain1_len + opts.chain2_len;

    std::vector<EigenXform> dumped_poses;
    Size ref_resi = pose.size()/2;

    Real rg=15;

    for(Size itry = 1; itry <= opts.nstruct; ++itry) {
        
        std::cout << "mcmc ... " << itry << std::endl;
        // re-random the root (this function will update the coordinates automatically)
        pose.random_root(true, true, true, true);
        // not perturbing the jumps between chains
        bool success = sampling::mcmc_light_heterodimer(pose, sfxn, sidechain_neighbors, opts);

        if( !success ) continue;

        bool is_redundant(false);
        EigenXform ref_stub = pose.stub(ref_resi);
        EigenXform ref_stub_inverse = ref_stub.inverse();
        for(Size i_pose=0; i_pose<dumped_poses.size(); ++i_pose) {
            Real rmsd = utils::xform_magnitude(ref_stub_inverse * dumped_poses[i_pose], rg); // rg equals 15 for 65 aa scaffolds
            if(rmsd<opts.redundancy_rmsd_cutoff) {
                is_redundant = true;
                break;
            }
        }
        if(is_redundant) {
            std::cout << "Find a redundant dock with rmsd_cutoff" << opts.redundancy_rmsd_cutoff << std::endl;
            continue;
        } else {
            dumped_poses.push_back(ref_stub);
        }

        Real total_sc = sfxn.score(pose) / pose.size();

        std::stringstream name;

        if(opts.output_dir != "")
            name << opts.output_dir << "/";
        // prefix, itry, time
        if(opts.output_prefix != "") name << opts.output_prefix << "_";
        
        name << std::chrono::high_resolution_clock::now().time_since_epoch().count();
        name << "_TotalSc_" << total_sc;

        
        Real score_interchain = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->score(pose)+std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("clash"))->score(pose)) / total_len;
        name << "_InterChain_" << score_interchain;

        Real privileged_motif_score = opts.privileged_motif_weight * std::static_pointer_cast<scoring::PrivilegedMotifScoreMethod>(sfxn.get_method("privileged_motif"))->score(pose);
        name << "_PIM_" << privileged_motif_score; 
        


        name << ".pdb"<<(opts.gzip?".gz":"");

        std::cout << "Dump model: " << name.str() << std::endl;

        pose.dump_pdb(name.str(), false, opts.gzip, false);

    }

    auto t2 = Clock::now();    
    std::cout << "Time usage: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count()
              << " nanoseconds" << std::endl;

    return 0;
}

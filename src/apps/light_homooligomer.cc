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

    std::cout << "Protein length: " << opts.len << std::endl
              << "Num Repeats: "    << opts.num_repeats << std::endl
              << "Symmetry: "       << opts.symmetry << std::endl;


    std::cout << "Prepare the pose ...\n";
    Pose pose(opts.len, opts.num_repeats, opts.symmetry);
    pose.set_dssp(opts.ss);
    // when the ss is set, but seq is empty; the seq will be automatically set
    // when the scaffold is loaded, the sequence will be overwrite
    pose.set_sequence(opts.seq);

    // load the scaffold pdb
    assert(opts.motif_pdb != "");
    {
        std::cout << "Loading motif " << opts.motif_pdb << " into position " << opts.motif_pos << std::endl;
        pose.load_pdb(opts.motif_pdb, opts.motif_pos, false, true, true);
        pose.set_root_index(opts.motif_pos);
        pose.update_coordinates();
    }
    // prepare the scoring function
    std::cout << "Prepare the score function ...\n";
    std::cout << "Use secondary informtion for rpx scoring: " << opts.use_ss << std::endl;
    ScoreFunction sfxn;
    sfxn.regist_method("clash",ScoreMethodOP(new ClashScoreMethod()));
    // inter chain pair only
    sfxn.get_method("clash")->set_energy_type(TWO_BODY_INTER_CHAIN);
    std::static_pointer_cast<ClashScoreMethod>(sfxn.get_method("clash"))->set_CB_swelling_factor(opts.CB_swelling_factor);

    sfxn.regist_method("rpx",ScoreMethodOP(new RpxScoreMethod(opts.rpx_db, opts.rpx_cart_resl, opts.rpx_ang_resl, opts.use_ss)));
    // inter chain pair only
    sfxn.get_method("rpx")->set_energy_type(TWO_BODY_INTER_CHAIN);
    if(opts.designable_residues != "") {
        sfxn.get_method("rpx")->set_designable_res(opts.designable_residues);
    }
    //
    // privileged intarface motif score method
    sfxn.regist_method("privileged_interface_motif", ScoreMethodOP(new PrivilegedInterfaceMotifScoreMethod(opts.privileged_interface_motif_hash_table, opts.privileged_interface_motif_hash_cart_resl, opts.privileged_interface_motif_hash_ang_resl)));
    std::cout << "Loading the privileged interface motif hash table " << opts.privileged_interface_motif_hash_table << std::endl;
    sfxn.get_method("privileged_interface_motif")->set_weight(opts.privileged_interface_motif_weight);
    if(opts.designable_residues != "") {
        sfxn.get_method("privileged_interface_motif")->set_designable_res(opts.designable_residues);
    }

    metric::SidechainNeighbors sidechain_neighbors;

    auto t1 = Clock::now();

    for(Size itry = 1; itry <= opts.nstruct; ++itry) {
        
        std::cout << "mcmc ... " << itry << std::endl;
        // re-random the root (this function will update the coordinates automatically)
        pose.random_root(true, true, true, true);
        // not perturbing the jumps between chains
        bool success = sampling::mcmc_light_assembly(pose, sfxn, sidechain_neighbors, opts);

        if( !success ) continue;

        Real scs = sidechain_neighbors.compute(pose,false,true);
        Real total_sc = sfxn.score(pose) / pose.size();

        std::stringstream name;

        if(opts.output_dir != "")
            name << opts.output_dir << "/";
        // prefix, itry, time
        if(opts.output_prefix != "") name << opts.output_prefix << "_";
        
        name << std::chrono::high_resolution_clock::now().time_since_epoch().count();
        name << "_SN_" << std::setprecision(3) << scs;
        name << "_TotalSc_" << total_sc;

        
        Real score_interchain =(std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1) ) / pose.size();
        name << "_ScoreInterChain_" << score_interchain;

        Real privileged_interface_motif_score = opts.privileged_interface_motif_weight * std::static_pointer_cast<PrivilegedInterfaceMotifScoreMethod>(sfxn.get_method("privileged_interface_motif"))->score(pose);
        name << "_PIM_" << privileged_interface_motif_score; 
        


        name << ".pdb"<<(opts.gzip?".gz":"");

        std::cout << "Dump model: " << name.str() << std::endl;

        pose.dump_pdb(name.str(),false,opts.gzip);

    }

    auto t2 = Clock::now();    
    std::cout << "Time usage: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count()
              << " nanoseconds" << std::endl;

    return 0;
}

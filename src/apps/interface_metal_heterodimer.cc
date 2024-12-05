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
#include "utils/dssp.hh"

#include "scene/InterfaceMetal.hh"
#include "scene/InterfaceZinc.hh"
#include "scene/InterfaceCalcium.hh"

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

    std::cout << "Code refactored at 2024-01-25, please do test first!!!" << std::endl;
    exit(0);

    std::cout << "Prepare the metal configuration" << std::endl
              << "    ==>  Metal config1: " << opts.interface_metal_config1 << std::endl
              << "    ==>  Metal config2: " << opts.interface_metal_config2 << std::endl
              << "    ==>  Metal Radius: " << opts.interface_metal_radius << std::endl
              << "    ==>  Required number of Metals: " << opts.num_interface_metals << std::endl
              << "    ==>  Metal score cutoff: " << opts.interface_metal_score_cutoff << std::endl
              << "    ==>  Metal distance optimization weight: " << opts.interface_metal_distance_optimization_weight << std::endl
              << "    ==>  Excluded metal clusters: " << opts.exclude_interface_metal_names << std::endl;

    InterfaceMetalOP interface_metal_p;
    if(opts.interface_metal_type == "zinc") {
        interface_metal_p = std::make_shared<InterfaceZinc>(InterfaceZinc());
    } else if (opts.interface_metal_type == "calcium") {
        interface_metal_p = std::make_shared<InterfaceCalcium>(InterfaceCalcium());
    } else {
        std::cout << "Allowed metal types include zinc and calcium. If you want to work on a new metal, let's talk first." << std::endl; 
    }
    interface_metal_p->load_metal_configs(opts.interface_metal_config1, opts.interface_metal_config2);
    if(opts.exclude_interface_metal_names != "") {
        interface_metal_p->exclude_metals(opts.exclude_interface_metal_names);
    }
    interface_metal_p->set_metal_radius(opts.interface_metal_radius);

    std::cout << "Load and prepare the pdbs: " << std::endl
              << "    ==>  Chain1 len: " << opts.interface_metal_chain1_len << std::endl
              << "    ==>  Chain1 pdb: " << opts.interface_metal_chain1_pdb << std::endl
              << "    ==>  Chain2 len: " << opts.interface_metal_chain2_len << std::endl
              << "    ==>  Chain2 pdb: " << opts.interface_metal_chain2_pdb << std::endl;

    Pose pose(opts.interface_metal_chain1_len);
    PoseOP target_pose = std::make_shared<Pose>(opts.interface_metal_chain2_len);

    pose.load_pdb(opts.interface_metal_chain1_pdb, 1, false, true, true);
    target_pose->load_pdb(opts.interface_metal_chain2_pdb, 1, true, true, true);
    std::string ss1 = utils::get_dssp_from_pose(pose);
    std::string ss2 = utils::get_dssp_from_pose(*target_pose);
    std::cout << "The secondary structure of chain1 is automatically determined to be:" << std::endl
              << "    ==> " << ss1 << std::endl;
    std::cout << "The secondary structure of chain2 is automatically determined to be:" << std::endl
              << "    ==> " << ss2 << std::endl;
    pose.set_dssp(ss1);
    target_pose->set_dssp(ss2);

    Vec target_center = target_pose->conformation().center_vec();
    EigenXform new_root = target_pose->conformation().root();
    new_root.translation() = new_root.translation() - target_center;
    target_pose->conformation().set_root(new_root);
    target_pose->conformation().update_coordinates();
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
    std::cout << "Done!" << std::endl;

    // pose.dump_pdb("test.pdb", false,opts.gzip,false);

    auto t1 = Clock::now();

    Size total_len = opts.interface_metal_chain1_len + opts.interface_metal_chain2_len;

    for(Size itry = 1; itry <= opts.nstruct; ++itry) {
        
        std::cout << "mcmc ... " << itry << std::endl;
        // not perturbing the jumps between chains
        bool success = sampling::interface_metal_heterodimer_mcmc(pose, interface_metal_p, sfxn, opts);

        if( !success ) continue;

        Real rpx_sc = std::static_pointer_cast<RpxScoreTargetMethod>(sfxn.get_method("rpx"))->score(pose);
        Real clash_sc = std::static_pointer_cast<VoxelClashScoreMethod>(sfxn.get_method("clash"))->score(pose);
        Real interface_metal_sc = interface_metal_p->topN_pair_distance(opts.num_interface_metals) * opts.interface_metal_distance_optimization_weight;

        Real total_sc = (rpx_sc + clash_sc)/total_len + interface_metal_sc;

        std::vector<std::string> interface_metal_pairs = interface_metal_p->report_topN_pairs(opts.num_interface_metals);
        std::string interface_metal_str = interface_metal_pairs[0];
        for(Size i=1; i<opts.num_interface_metals;++i) interface_metal_str += "_" + interface_metal_pairs[i];

        std::stringstream name;

        if(opts.output_dir != "")
            name << opts.output_dir << "/";
        // prefix, itry, time
        if(opts.output_prefix != "") name << opts.output_prefix << "_";
        
        name << std::chrono::high_resolution_clock::now().time_since_epoch().count();
        name << "_TotalScore_" << total_sc;
        name << "_RPX_" << rpx_sc / total_len;
        name << "_Clash_" << clash_sc / total_len;
        name << "_MetalScore_" << interface_metal_sc;
        name << "_MetalDist_" << interface_metal_sc/opts.interface_metal_distance_optimization_weight;
        name << "__" << interface_metal_str;
        name << ".pdb"<<(opts.gzip?".gz":"");

        std::cout << "Dump model: " << name.str() << std::endl;
        
        pose.dump_pdb(name.str(),false,opts.gzip,false);

    }

    auto t2 = Clock::now();    
    std::cout << "Time usage: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count()
              << " nanoseconds" << std::endl;

    return 0;
}

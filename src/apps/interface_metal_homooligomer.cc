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
#include "scene/InterfaceNickle.hh"

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
using namespace metric;

// for monomers only

int main(int argc, char *argv[])
{
    Options opts;
    opts.parse_args(argc, argv);

    std::cout << "General Configurations" << std::endl
              << "    ==> Scaffold pdb: " << opts.interface_metal_scaffold_pdb << std::endl
              << "    ==> Protein len: " << opts.len << std::endl
              << "    ==> Symmetry: " << opts.symmetry << std::endl
              << "    ==> Dssp: " << opts.ss << std::endl
              << "    ==> rpx_path: " << opts.rpx_db << std::endl
              << "    ==> use_ss: " << opts.use_ss << std::endl
              << "    ==> Output prefix: " << opts.output_prefix << std::endl;

    assert(opts.symmetry != "C1" && opts.num_repeats == 1);

    std::cout << "Prepare the pose ...\n";
    Pose pose(opts.len, opts.num_repeats, opts.symmetry);
    pose.set_sequence(opts.seq);;
    
    assert(opts.interface_metal_scaffold_pdb != "");
    {
        std::cout << "Loading scaffold pdb file " << opts.interface_metal_scaffold_pdb <<  std::endl;
        pose.load_pdb(opts.interface_metal_scaffold_pdb, 1, false, true, true);
        pose.set_root_index(1);
    }
    pose.set_dssp(opts.ss);

    std::cout << "Metal Configurations" << std::endl
              << "    ==>  Metal type: " << opts.interface_metal_type << std::endl
              << "    ==>  Config1: " << opts.interface_metal_config1 << std::endl
              << "    ==>  Config2: " << (opts.interface_metal_config2 == ""? opts.interface_metal_config1:opts.interface_metal_config2) << std::endl
              << "    ==>  Metal Radius: " << opts.interface_metal_radius << std::endl
              << "    ==>  Required number of metals: " << opts.num_interface_metals << std::endl
              << "    ==>  Interface metal score cutoff: " << opts.interface_metal_score_cutoff << std::endl
              << "    ==>  Interface metal distance optimization weight: " << opts.interface_metal_distance_optimization_weight << std::endl
              << "    ==>  Excluded interface metal clusters: " << opts.exclude_interface_metal_names << std::endl;

    InterfaceMetalOP interface_metal_p;
    InterfaceMetalOP interface_metal_p2;
    if(opts.interface_metal_type == "zinc") {
        interface_metal_p = std::make_shared<InterfaceZinc>(InterfaceZinc());
        interface_metal_p2 = std::make_shared<InterfaceZinc>(InterfaceZinc());
    } else if (opts.interface_metal_type == "calcium") {
        interface_metal_p = std::make_shared<InterfaceCalcium>(InterfaceCalcium());
        interface_metal_p2 = std::make_shared<InterfaceCalcium>(InterfaceCalcium());
    } else if (opts.interface_metal_type == "nickle") {
        interface_metal_p = std::make_shared<InterfaceNickle>(InterfaceNickle());
        interface_metal_p2 = std::make_shared<InterfaceNickle>(InterfaceNickle());
    } else {
        std::cout << "Allowed metal types include zinc and calcium. If you want to work on a new metal, let's talk first." << std::endl; 
    }

    interface_metal_p->load_metal_configs(opts.interface_metal_config1, opts.interface_metal_config2 == ""? opts.interface_metal_config1:opts.interface_metal_config2);
    if(opts.exclude_interface_metal_names != "") {
        interface_metal_p->exclude_metals(opts.exclude_interface_metal_names);
    }
    interface_metal_p->set_metal_radius(opts.interface_metal_radius);


    if(opts.symmetry.at(0) == 'H') {
        interface_metal_p2->load_metal_configs(opts.interface_metal_config1, opts.interface_metal_config2 == ""? opts.interface_metal_config1:opts.interface_metal_config2);
        if(opts.exclude_interface_metal_names != "") {
            interface_metal_p2->exclude_metals(opts.exclude_interface_metal_names);
        }
        interface_metal_p2->set_metal_radius(opts.metal_radius);
    }


    std::cout << "Prepare the score function ...\n";
    std::cout << "Use secondary informtion for rpx scoring: " << opts.use_ss << std::endl;
    ScoreFunction sfxn;
    sfxn.regist_method("clash",ScoreMethodOP(new ClashScoreMethod()));
    sfxn.get_method("clash")->set_energy_type(TWO_BODY_INTER_CHAIN);
    std::static_pointer_cast<ClashScoreMethod>(sfxn.get_method("clash"))->set_CB_swelling_factor(opts.CB_swelling_factor);
    sfxn.regist_method("rpx",ScoreMethodOP(new RpxScoreMethod(opts.rpx_db, opts.rpx_cart_resl, opts.rpx_ang_resl, opts.use_ss)));
    // inter chain pair only
    sfxn.get_method("rpx")->set_energy_type(TWO_BODY_INTER_CHAIN);
    std::cout << "Done!" << std::endl;

    SidechainNeighbors sidechain_neighbors;

    auto t1 = Clock::now();

    const Size protein_len = pose.size();

    Size metal_chain_idx(-1);
    Size metal2_chain_idx(-1);

    for(Size itry = 1; itry <= opts.nstruct; ++itry) {
        
        std::cout << "mcmc ... " << itry << std::endl;
        // not perturbing the jumps between chains
        bool success;
        if(opts.symmetry.at(0) != 'H') {
            success = sampling::interface_metal_homooligomer_mcmc(pose, interface_metal_p, sfxn, sidechain_neighbors, opts);
        } else {
            success = sampling::interface_metal_fiber_mcmc(pose, interface_metal_p, interface_metal_p2, metal_chain_idx, metal2_chain_idx, sfxn, sidechain_neighbors, opts);
        }

        if( !success ) continue;

        Real rpx_sc = std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1);
        Real clash_sc = std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1);
        Real interface_metal_sc = interface_metal_p->topN_pair_distance(opts.num_interface_metals) * opts.interface_metal_distance_optimization_weight;
        std::vector<std::string> interface_metal_pairs = interface_metal_p->report_topN_pairs(opts.num_interface_metals);
        std::string interface_metal_str = interface_metal_pairs[0];
        for(Size i=1; i<opts.num_interface_metals;++i) interface_metal_str += "_" + interface_metal_pairs[i];

        // report the second
        Real interface_metal2_sc(0.0);
        std::string interface_metal2_str("");
        if(opts.symmetry.at(0) == 'H') {
            interface_metal2_sc = interface_metal_p2->topN_pair_distance(opts.num_interface_metals) * opts.interface_metal_distance_optimization_weight;
            std::vector<std::string> interface_metal2_pairs = interface_metal_p2->report_topN_pairs(opts.num_interface_metals);
            interface_metal2_str = interface_metal2_pairs[0];
            for(Size i=1; i<opts.num_interface_metals;++i) interface_metal2_str += "_" + interface_metal2_pairs[i];
        }

        Real sidechain_neighbor_v = sidechain_neighbors.compute(pose,false,true);
        Real total_sc = (rpx_sc + clash_sc) / protein_len + interface_metal_sc + interface_metal2_sc;

        std::stringstream name;

        if(opts.output_dir != "")
            name << opts.output_dir << "/";
        // prefix, itry, time
        if(opts.output_prefix != "") name << opts.output_prefix << "_";
        
        name << std::chrono::high_resolution_clock::now().time_since_epoch().count();
        name << "_SN_" << std::setprecision(3) << sidechain_neighbor_v;
        name << "_TotalSc_" << total_sc;
        name << "_RPX_" << rpx_sc / protein_len;
        name << "_Clash_" << clash_sc / protein_len;
        
        // report the first contacting chain
        if(opts.symmetry.at(0) == 'H') {
            name << "_MetalChIdx_" << metal_chain_idx;
        }
        name << "_MetalDist_" << interface_metal_sc/opts.interface_metal_distance_optimization_weight;
        name << "_" << interface_metal_str;

        if(opts.symmetry.at(0) == 'H') {
            name << "_Metal2ChIdx_" << metal2_chain_idx;
            name << "_Metal2Dist_" << interface_metal2_sc/opts.interface_metal_distance_optimization_weight;
            name << "_" << interface_metal2_str;
        }
        
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

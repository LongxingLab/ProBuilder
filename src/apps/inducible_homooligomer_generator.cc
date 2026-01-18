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
#include "utils/ligand_util.hh"

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
#include <cstring>

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

    assert(opts.motif_pdb != "" && opts.rif_table != "" && opts.context_pdb != "");
    assert(opts.ligand_name != "");

    std::cout << "Protein length: " << opts.len << std::endl
              << "Num Repeats: "    << opts.num_repeats << std::endl
              << "Symmetry: "       << opts.symmetry << std::endl
              << "Ligand: "         << opts.ligand_name << std::endl;

    if(opts.ligand_name != "AMA" || opts.ligand_name != "HQA") {
        std::cout << "Currently, ligand-neighbors metric only supports AMA and HQA," << std::endl;
        std::cout << "if you want to work on another ligand, modify the code first." << std::endl;
    }

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
    //
    // load the rif table and context pdb into voxel array
    assert(opts.rif_table != "" && opts.context_pdb != "");

    Size num_ligands(1);
    std::cout << "Loading rif table " << opts.rif_table << std::endl;
    sfxn.regist_method("target_rif",ScoreMethodOP(new RifScoreMethod(opts.rif_table, opts.rif_cart_resl, opts.rif_ang_resl)));
    sfxn.get_method("target_rif")->set_weight(opts.rif_weight);
    sfxn.get_method("target_rif")->set_score_type(TARGET_RIF);
    pose.energy_manager().add_energy_onebody(TARGET_RIF);

    std::cout << "Loading context pdb " << opts.context_pdb << " for clash check" << std::endl;
    sfxn.regist_method("target_context_clash",ScoreMethodOP(new VoxelClashScoreMethod(opts.context_pdb,"ALL",opts.context_clash_radius,0.25F,false,"MOLECULE")));
    sfxn.get_method("target_context_clash")->set_score_type(TARGET_CONTEXT_CLASH);
    sfxn.get_method("target_context_clash")->set_weight(opts.context_clash_weight);
    pose.energy_manager().add_energy_onebody(TARGET_CONTEXT_CLASH);


    EigenXform extra_ligand_xform(EigenXform::Identity());
    std::string extra_xform_str;
    if(opts.extra_ligand != "") {
        num_ligands = 2;

        if( opts.extra_ligand == "RANDOM") {
            extra_ligand_xform = utils::random_place_AMA(extra_xform_str);
        } else {
            std::vector<std::string> splited_str;
            std::string delimiters = "_";
            std::string::size_type lastPos = opts.extra_ligand.find_first_not_of(delimiters,0);
            std::string::size_type pos = opts.extra_ligand.find_first_of(delimiters, lastPos);
            while( std::string::npos != pos || std::string::npos != lastPos ) {
                splited_str.push_back(opts.extra_ligand.substr(lastPos, pos-lastPos));
                lastPos = opts.extra_ligand.find_first_not_of(delimiters, pos);
                pos = opts.extra_ligand.find_first_of(delimiters, lastPos);
            }

            assert(splited_str.size() == 3); // flip, translate, rotate
            Real ligand_flip = stof(splited_str[0]);
            Real ligand_translate = stof(splited_str[1]);
            Real ligand_rotate = stof(splited_str[2]);

            EigenXform xform_flip(EigenXform::Identity());
            if( ligand_flip == -1 ) {
                xform_flip.rotate( AngleAxis( 3.14159265359, Vec(1.0,0.0,0.0) ) );
            }
            EigenXform xform_translate(EigenXform::Identity());
            xform_translate.translation() = Vec(0.0, 0.0, ligand_translate);
            EigenXform xform_rotate(EigenXform::Identity());
            xform_rotate.rotate( AngleAxis( ligand_rotate /180.0*3.14159265359, Vec(0.0,0.0,1.0) ) );

            extra_ligand_xform = xform_rotate * xform_translate * xform_flip;


            std::stringstream tmp_xform_stream;
            tmp_xform_stream << "REMARK PDBinfo-LABEL:    1 LIGAND_XFORM:";
            tmp_xform_stream << extra_ligand_xform(0,0) << "_" << extra_ligand_xform(0,1) << "_" << extra_ligand_xform(0,2) << "_" << extra_ligand_xform(0,3) << "_"
                      << extra_ligand_xform(1,0) << "_" << extra_ligand_xform(1,1) << "_" << extra_ligand_xform(1,2) << "_" << extra_ligand_xform(1,3) << "_"
                      << extra_ligand_xform(2,0) << "_" << extra_ligand_xform(2,1) << "_" << extra_ligand_xform(2,2) << "_" << extra_ligand_xform(2,3);

            extra_xform_str = tmp_xform_stream.str();
        }

        sfxn.regist_method("target_rif2",ScoreMethodOP(new RifScoreMethod(opts.rif_table, opts.rif_cart_resl, opts.rif_ang_resl)));
        sfxn.get_method("target_rif2")->set_weight(opts.rif_weight);
        sfxn.get_method("target_rif2")->set_score_type(TARGET_RIF2);
        std::static_pointer_cast<RifScoreMethod>(sfxn.get_method("target_rif2"))->set_relative_pos(extra_ligand_xform);
        pose.energy_manager().add_energy_onebody(TARGET_RIF2);

        sfxn.regist_method("target_context_clash2",ScoreMethodOP(new VoxelClashScoreMethod(opts.context_pdb,"ALL",opts.context_clash_radius,0.25F,false,"MOLECULE")));
        sfxn.get_method("target_context_clash2")->set_score_type(TARGET_CONTEXT_CLASH2);
        std::static_pointer_cast<VoxelClashScoreMethod>(sfxn.get_method("target_context_clash2"))->set_relative_pos(extra_ligand_xform);
        sfxn.get_method("target_context_clash2")->set_weight(opts.context_clash_weight);
        pose.energy_manager().add_energy_onebody(TARGET_CONTEXT_CLASH2);
    }


    auto t1 = Clock::now();

    for(Size itry = 1; itry <= opts.nstruct; ++itry) {


        // random_place_AMA
        if(itry % 20 == 0 && opts.extra_ligand == "RANDOM" ) {
            extra_ligand_xform = utils::random_place_AMA(extra_xform_str);
            std::static_pointer_cast<RifScoreMethod>(sfxn.get_method("target_rif2"))->set_relative_pos(extra_ligand_xform);
            std::static_pointer_cast<VoxelClashScoreMethod>(sfxn.get_method("target_context_clash2"))->set_relative_pos(extra_ligand_xform);
            std::cout << "Randomly place the second AMA to a new position ...";
        }
        
        std::cout << "mcmc ... " << itry << std::endl;
        // re-random the root (this function will update the coordinates automatically)
        pose.random_root(true, true, true, true);
        // not perturbing the jumps between chains
        bool success = sampling::mcmc_jump(pose, sfxn, opts, extra_ligand_xform);

        if( !success ) continue;

        Real scs = utils::sidechain_neighbors(pose,false,true);
        Real total_sc = sfxn.score(pose) / pose.size();

        std::stringstream name;

        if(opts.output_dir != "")
            name << opts.output_dir << "/";
        // prefix, itry, time
        if(opts.output_prefix != "") name << opts.output_prefix << "_";
        
        name << std::chrono::high_resolution_clock::now().time_since_epoch().count();
        name << "_SN_" << std::setprecision(3) << scs;
        name << "_TS_" << total_sc;

        
        Real score_interchain =(std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1) ) / pose.size();
        name << "_SIC_" << score_interchain;
        Real ligand1_rif_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score(pose);
        Real ligand1_clash_score = opts.context_clash_weight * std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash"))->score(pose);
        name << "_Lig1RIF_" << ligand1_rif_score << "_Lig1C_" << ligand1_clash_score;
        
        // which ligand
        Real ligand_neighbors(0.0);
        if(opts.ligand_name == "AMA") {
            ligand_neighbors = utils::AMA_neighbors(pose);
        }  else if(opts.ligand_name == "HQA") {
            ligand_neighbors = utils::HQA_neighbors(pose);
        } else {
            std::cout << "Holy shit!!! I am died." << std::endl;
            exit(0);
        }

        name << "_Lig1N_" << ligand_neighbors;
        if(opts.extra_ligand != "") {
            Real ligand2_rif_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif2"))->score(pose);
            Real ligand2_clash_score = opts.context_clash_weight * std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash2"))->score(pose);
            name << "_Lig2RIF_" << ligand2_rif_score << "_Lig2C_" << ligand2_clash_score;
            
            // which ligand
            Real ligand2_neighbors(0.0);
            if(opts.ligand_name == "AMA") {
                ligand2_neighbors = utils::AMA_neighbors(pose, extra_ligand_xform);
            }  else if(opts.ligand_name == "HQA") {
                ligand2_neighbors = utils::HQA_neighbors(pose, extra_ligand_xform);
            } else {
                std::cout << "Holy shit!!! I am died." << std::endl;
                exit(0);
            }
            name << "_Lig2N_" << ligand2_neighbors;
        }

        name << ".pdb"<<(opts.gzip?".gz":"");

        std::cout << "Dump model: " << name.str() << std::endl;

        if(opts.extra_ligand == "") {
            pose.dump_pdb(name.str(),false,opts.gzip);
        }
        else {
            pose.dump_pdb(name.str(),false,opts.gzip, false, false, extra_xform_str);
        }

    }

    auto t2 = Clock::now();    
    std::cout << "Time usage: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count()
              << " nanoseconds" << std::endl;

    return 0;
}

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
    if(!opts.hallucinate) {
        pose.set_dssp(opts.ss);
        pose.set_sequence(opts.seq);
    }
    std::cout << "Prepare the score function ...\n";
    std::cout << "Use secondary informtion for rpx scoring: " << opts.use_ss << std::endl;
    ScoreFunction sfxn;
    sfxn.regist_method("clash",ScoreMethodOP(new ClashScoreMethod()));
    std::static_pointer_cast<ClashScoreMethod>(sfxn.get_method("clash"))->set_CB_swelling_factor(opts.CB_swelling_factor);
    if(opts.num_repeats>1){
        sfxn.regist_method("repeat",ScoreMethodOP(new RepeatScoreMethod(std::pair<Real,Real>(opts.repeat_min_radius,opts.repeat_max_radius),std::pair<Real,Real>(opts.repeat_min_rise,opts.repeat_max_rise),std::pair<Real,Real>(opts.repeat_min_omega,opts.repeat_max_omega))));
        sfxn.get_method("repeat")->set_weight(opts.repeat_weight);
    }
    ScoreFunction sfxn_refine;
    if(opts.rpx_db_refine != "") {
        sfxn_refine.regist_method("clash",ScoreMethodOP(new ClashScoreMethod()));
        sfxn_refine.regist_method("rpx",ScoreMethodOP(new RpxScoreMethod(opts.rpx_db_refine, opts.rpx_cart_resl_refine, opts.rpx_ang_resl_refine, opts.use_ss)));
    }

    if(opts.disulfide_hash_table != "") {
        sfxn.regist_method("disulfide", ScoreMethodOP(new DisulfideScoreMethod(opts.disulfide_hash_table, opts.disulfide_min_sequence_separation, opts.disulfide_hash_cart_resl, opts.disulfide_hash_ang_resl)));
        pose.energy_manager().add_energy_twobody_intra(DISULFIDE,opts.num_repeats);
        sfxn.get_method("disulfide")->set_weight(opts.disulfide_weight);
        std::static_pointer_cast<DisulfideScoreMethod>(sfxn.get_method("disulfide"))->set_scnb_burial_cutoff(opts.gold_burial_cutoff);
    }

    sfxn.regist_method("rpx",ScoreMethodOP(new RpxScoreMethod(opts.rpx_db, opts.rpx_cart_resl, opts.rpx_ang_resl, opts.use_ss)));
    std::cout << "Load the fragment_library: " << opts.frag_path << std::endl;
    FragMapMover frag_mover(opts.len, 7);
    frag_mover.load_fragments(opts.frag_path);


    auto t1 = Clock::now();

    std::string hallucinate_dssp_prefix = "";

    for(Size itry = 1; itry <= opts.nstruct; ++itry) {
        
        if(opts.hallucinate && (itry-1) % 20 == 0){
            std::string sequence = "";
            std::string dssp, dssp_prefix;
            // std::string dssp = random_helix_dssp(opts.min_helix_len,opts.max_helix_len,opts.min_loop_len,opts.max_loop_len,opts.max_len,opts.segment_num);
            bool dssp_pass;

            if(opts.num_repeats == 1)
                dssp_pass = random_bundle_dssp(opts.len, opts.num_helix, opts.min_helix_len, opts.max_helix_len, opts.min_loop_len, opts.max_loop_len, dssp_prefix, dssp);
            else
                dssp_pass = random_repeat_bundle_dssp(opts.len, opts.num_repeats, opts.num_helix, opts.min_helix_len, opts.max_helix_len, opts.min_loop_len, opts.max_loop_len, dssp_prefix, dssp);

            std::cout << "dssp : " << dssp  << std::endl;
            for( Size ii=0; ii<dssp.length(); ++ii ) {
                if(dssp.at(ii) == 'L')
                    sequence += "G";
                else
                    sequence += "V";
            }
            std::cout << "sequence: " << sequence << std::endl;
            if (dssp_pass || itry == 1 ) {
                pose.set_dssp(dssp);
                pose.set_sequence(sequence);

                hallucinate_dssp_prefix = dssp_prefix;
            }
        }


        pose.reset_coords();
        
        std::cout << "mcmc ... " << itry << std::endl;
        // not perturbing the jumps between chains
        // bool success = opts.peptide?sampling::mcmc_repeat_peptide(pose, sfxn, frag_mover, opts):sampling::mcmc(pose, sfxn, frag_mover, opts);
        bool success = sampling::mcmc(pose, sfxn, frag_mover, opts);

        if( !success ) continue;

        if(opts.rpx_db_refine != "") {
            std::cout << "refine the structure with a different rpx table ..." << std::endl;
            sampling::mcmc_refine(pose, sfxn_refine, frag_mover, opts);
        }
        
        // for correct repeat rpx score for 2 repeat
        
        Real scs = utils::sidechain_neighbors(pose);
        Real total_sc = sfxn.score(pose) / pose.size();

        std::stringstream name;

        if(opts.output_dir != "")
            name << opts.output_dir << "/";
        // prefix, itry, time
        if(opts.output_prefix != "") name << opts.output_prefix << "_";
        if( hallucinate_dssp_prefix != "" ) name << hallucinate_dssp_prefix << "_";
        
        name << std::chrono::high_resolution_clock::now().time_since_epoch().count();
        name << "_Sc_" << std::setprecision(3) << scs;
        name << "_TotalScore_" << total_sc;
        if(pose.num_chains()>1)name << "_Cscore_" << std::static_pointer_cast<RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,pose.symmetry().at(0)=='C'?-1:pose.num_chains()/2)/pose.size();
        if(pose.symmetry().at(0)=='D') {
            name << "_Dscore_" << std::static_pointer_cast<RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,pose.num_chains()/2+1,-1)/pose.size();
        }
        name << "_Clash_" << std::static_pointer_cast<ClashScoreMethod>(sfxn.get_method("clash"))->score(pose)/pose.size();

        if(opts.disulfide_hash_table != "") {
            Real disulfide_score = std::static_pointer_cast<DisulfideScoreMethod>(sfxn.get_method("disulfide"))->score(pose);
            name << "_DisulfideNum_" << -disulfide_score; 
        }

        if(opts.num_repeats>1){
            // std::cout<<std::static_pointer_cast<RepeatScoreMethod>(sfxn.get_method("repeat"))->score(pose)<<std::endl;
            Real radius,rise,omega;
            Vec tmpvec;
            get_repeat_parameters_from_coords(pose,rise,radius,omega,tmpvec,tmpvec);
            Real repeat_score = opts.repeat_weight * std::static_pointer_cast<RepeatScoreMethod>(sfxn.get_method("repeat"))->score(pose);
	    omega = omega/3.1415926*180;
            name << "_RepeatScore_" << repeat_score;
            name  << "_Radius_"<<radius<<"_Omega_"<<abs(omega)<<"_Rise_"<<abs(rise);
        }
        name << ".pdb"<<(opts.gzip?".gz":"");

        std::cout << "Dump model: " << name.str() << std::endl;
        
        std::string original_seq = "";
        if(opts.dump_polyG) {
            original_seq = pose.sequence();
            std::string polyG_str(pose.size(), 'G');
            pose.set_sequence(polyG_str);
        }
        if(opts.dump_trajactory&&opts.peptide){
            for(auto & i :pose.get_trajctory()){
                i.set_target_pose(nullptr);
            }
        }
        pose.dump_pdb(name.str(),false,opts.gzip,opts.dump_trajactory);
        pose.clear_trajctory();
        if(opts.dump_polyG) {
            // reverse the sequence back
            pose.set_sequence(original_seq);
        }

    }

    auto t2 = Clock::now();    
    std::cout << "Time usage: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count()
              << " nanoseconds" << std::endl;

    return 0;
}

#include "sampling/mcmc_protocols.hh"

#include "utils/utils.hh"
#include "basic/types.hh"
#include "scene/Pose.hh"
#include "sampling/FragMover.hh"
#include "scoring/ScoreFunction.hh"
#include "utils/random_util.hh"
#include "scene/InterfaceMetal.hh"
#include "scene/InterfaceZinc.hh"
#include "scene/InterfaceCalcium.hh"
#include "utils/ligand_util.hh"

#include <cmath>
#include <ctime>
#include <iostream>
#include<iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>

namespace sampling {

using namespace basic;


// assemble oligomers for metal coordination from partially satisfied metals
bool interface_metal_heterodimer_mcmc(scene::Pose & pose, scene::InterfaceMetalOP & interface_metal_p, scoring::ScoreFunction & sfxn, Options const & opts)
{
    std::uniform_real_distribution<Real> runif(0,1);

    Size total_len = opts.interface_metal_chain1_len + opts.interface_metal_chain2_len;


    Size root_index = utils::random_int(1,pose.size());
    //std::cout << "Root Index: " << root_index << std::endl;
    pose.set_root_index(root_index);
    pose.update_coordinates();

    // random root
    pose.random_root(true, true, true, true);


    interface_metal_p->update_metal_position(pose.stub(1), pose.get_target_pose()->stub(1));

    // take a snapshoot right before sampling
    //pose.root_snapshot();
    // score the pose first before taking the snapshot
    // sfxn.score(pose);
    sfxn.score(pose);
    
    pose.snapshot();
    interface_metal_p->snapshot();
    //pose.snapshot(1);

    // global parameters
    Size multi_cool_cycles = opts.mcmc_outer_cycles;
    Size total_steps = opts.mcmc_inner_cycles==-1 ? total_len * 10 : opts.mcmc_inner_cycles;
    Size M = total_steps / 10;


    Real best_score = 9e9;
    Real best_interface_metal_score = 9e9;
    Real best_rpx_score = 9e9;
    
    Real decrease_factor = 1.4;
    Real increase_factor = 1.1;

    Real num_interface_metals = opts.num_interface_metals;
    Real interface_metal_weight = opts.interface_metal_distance_optimization_weight;

    Real pre_rpx_sc = sfxn.score(pose);
    
    Real T = 0.5;
    for( Size ii_outer=0; ii_outer<multi_cool_cycles; ++ii_outer) {

        // local parameters for MCMC
        // Real T = 3.0 / (ii_outer+1); // multicycle cool good ???? perturb the pose out of local minima?? I don't know.
        T *= 2.0;

        Real perturb_ang_mag(30);
        Real perturb_trans_mag(5);

        for(Size ii_inner=0; ii_inner<total_steps; ++ii_inner) {
            

            pose.perturb_root(perturb_ang_mag, perturb_trans_mag);
            interface_metal_p->update_chain1_metal_position(pose.stub(1));


            Real rpx_sc = sfxn.score(pose) / total_len;
            Real interface_metal_sc = interface_metal_p->topN_pair_distance(opts.num_interface_metals) * interface_metal_weight;

            Real cur_score = rpx_sc + interface_metal_sc;
            Real delta = cur_score - best_score;
            bool pass = pass_metropolis(T, delta, runif(utils::rng()));

            if( pass ) {
            //    std::cout << "Accepted: " << cur_score << std::endl;
                pre_rpx_sc = rpx_sc;
                best_score = cur_score;
                best_rpx_score = rpx_sc;
                best_interface_metal_score = interface_metal_sc;
                pose.snapshot(true);
                interface_metal_p->snapshot(true);
            } else {
                rpx_sc = pre_rpx_sc;
                pose.rollback(true);
                interface_metal_p->rollback(true);
            }

            
            if(rpx_sc > 3000 || rpx_sc == 0.0) { // clash too much or no contact

                Size num_try(100);
                bool perturb_root_pass = false;
                while(--num_try){
                    pose.random_root(true, true, true, true);
                    rpx_sc = sfxn.score(pose) / total_len;
                    if(rpx_sc <= 3000 && rpx_sc != 0.0){
                        pre_rpx_sc = rpx_sc;
                        best_rpx_score = rpx_sc;
                        interface_metal_p->update_chain1_metal_position(pose.stub(1));
                        best_interface_metal_score = interface_metal_p->topN_pair_distance(opts.num_interface_metals) * interface_metal_weight;
                        best_score = rpx_sc + interface_metal_sc;
                        pose.snapshot(true);
                        interface_metal_p->snapshot(true);
                        perturb_root_pass = true;
                        break;
                    }
                }
                if( !perturb_root_pass ) {
                    pose.rollback(true);
                }  
            }

            
            if( ii_inner % M == 0 && ii_inner != 0 ) {
                T /= decrease_factor;
                perturb_ang_mag /= decrease_factor;
                perturb_ang_mag = perturb_ang_mag>5?perturb_ang_mag:5;
                perturb_trans_mag /= decrease_factor;
                perturb_trans_mag = perturb_trans_mag>0.5?perturb_trans_mag:0.5;
            }

            if(opts.verbose_level!=-1 && ii_inner%opts.verbose_level==0) {
                std::cout << std::setprecision(3) 
                          << "Outer: " << ii_outer 
                          << "; Inner: " << ii_inner 
                          << "; Total score: " << cur_score 
                          << "; RPX score: " << rpx_sc
                          << "; Metal score: " << interface_metal_sc 
                          << "; Metal distance: " << interface_metal_sc/interface_metal_weight 
                          <<  std::endl;
            }
            // after each move, check satisfication and return if good
            bool pass_requirements(true);
            // if(best_score > opts.total_score_cutoff) pass_requirements = false;
            if(pass_requirements && best_rpx_score > opts.inter_chain_score_cutoff) pass_requirements = false;
            if(pass_requirements && best_interface_metal_score > opts.interface_metal_score_cutoff) pass_requirements=false;

            if(pass_requirements) {
                std::cout << "Solution found! Eearly stop." << std::endl;
                return true;
            }

        } // inner loop
    } // outer loop

    std::cout << "No solution found! Too bad." << std::endl;
    return false;
}


// assemble oligomers for metal coordination from partially satisfied metals
bool interface_metal_homooligomer_mcmc(scene::Pose & pose, scene::InterfaceMetalOP & interface_metal_p, scoring::ScoreFunction & sfxn, metric::SidechainNeighbors & sidechain_neighbors, Options const & opts)
{
    std::uniform_real_distribution<Real> runif(0,1);

    Size protein_len = opts.len;


    Size root_index = utils::random_int(-10,10) + Size(protein_len/2);
    //std::cout << "Root Index: " << root_index << std::endl;
    pose.set_root_index(root_index);
    pose.update_coordinates();

    // random root
    pose.random_root(true, true, true, true);


    interface_metal_p->update_metal_position(pose.stub(1, 1), pose.stub(1, 2));

    // take a snapshoot right before sampling
    //pose.root_snapshot();
    // score the pose first before taking the snapshot
    // sfxn.score(pose);
    sfxn.score(pose);
    
    pose.snapshot();
    interface_metal_p->snapshot();
    //pose.snapshot(1);

    // global parameters
    Size multi_cool_cycles = opts.mcmc_outer_cycles;
    Size total_steps = opts.mcmc_inner_cycles==-1 ? protein_len * 10 : opts.mcmc_inner_cycles;
    Size M = total_steps / 10;


    Real best_score = 9e9;
    Real best_interface_metal_score = 9e9;
    Real best_rpx_score = 9e9;
    
    Real decrease_factor = 1.4;
    Real increase_factor = 1.1;

    Real num_interface_metals = opts.num_interface_metals;
    Real interface_metal_weight = opts.interface_metal_distance_optimization_weight;

    
    Real T = 0.5;
    for( Size ii_outer=0; ii_outer<multi_cool_cycles; ++ii_outer) {

        // local parameters for MCMC
        // Real T = 3.0 / (ii_outer+1); // multicycle cool good ???? perturb the pose out of local minima?? I don't know.
        T *= 2.0;

        Real perturb_ang_mag(30);
        Real perturb_trans_mag(8);

        for(Size ii_inner=0; ii_inner<total_steps; ++ii_inner) {
            

            pose.perturb_root(perturb_ang_mag, perturb_trans_mag);
            interface_metal_p->update_metal_position(pose.stub(1, 1), pose.stub(1, 2));


            Real rpx_sc = sfxn.score(pose)/protein_len;
            Real interface_metal_sc = interface_metal_p->topN_pair_distance(opts.num_interface_metals) * interface_metal_weight;


            Real cur_score = rpx_sc + interface_metal_sc;
            Real delta = cur_score - best_score;
            bool pass = pass_metropolis(T, delta, runif(utils::rng()));

            if( pass ) {
            //    std::cout << "Accepted: " << cur_score << std::endl;
                best_score = cur_score;
                best_rpx_score = rpx_sc;
                best_interface_metal_score = interface_metal_sc;
                pose.snapshot(true);
                interface_metal_p->snapshot(false);
            } else {
                pose.rollback(true);
                interface_metal_p->rollback(false);
            }


            // reroot

            if(best_rpx_score>opts.inter_chain_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor) {
                bool perturb_root_pass = false;
                Size num_try(100);
                while(--num_try){
                    pose.random_root(true, true, true, true);

                    rpx_sc = sfxn.score(pose)/protein_len;
                    
                    if(rpx_sc<=opts.inter_chain_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor){
                        interface_metal_p->update_metal_position(pose.stub(1, 1), pose.stub(1, 2));
                        interface_metal_sc = interface_metal_p->topN_pair_distance(opts.num_interface_metals) * interface_metal_weight;
                        best_rpx_score = rpx_sc;
                        best_interface_metal_score = interface_metal_sc;
                        best_score = best_rpx_score + best_interface_metal_score;
                        pose.snapshot(true);
                        interface_metal_p->snapshot(false);
                        perturb_root_pass = true;
                        break;
                    }  
                }
                if( !perturb_root_pass ) {
                    pose.rollback(true);
                }
            } else {

                pose.perturb_root(perturb_ang_mag, perturb_trans_mag);
                interface_metal_p->update_metal_position(pose.stub(1, 1), pose.stub(1, 2));

                rpx_sc = sfxn.score(pose)/protein_len;
                interface_metal_sc = interface_metal_p->topN_pair_distance(opts.num_interface_metals) * interface_metal_weight;
                cur_score = rpx_sc + interface_metal_sc;

                delta = cur_score - best_score;
                pass = pass_metropolis(T, delta, runif(utils::rng()));

                if( pass ) {
                //    std::cout << "Accepted: " << cur_score << std::endl;
                    best_score = cur_score;
                    best_rpx_score = rpx_sc;
                    best_interface_metal_score = interface_metal_sc;
                    pose.snapshot(true);
                    interface_metal_p->snapshot(false);
                } else {
                    pose.rollback(true);
                    interface_metal_p->rollback(false);
                }

            }

            
            if( ii_inner % M == 0 && ii_inner != 0 ) {
                T /= decrease_factor;
                perturb_ang_mag /= decrease_factor;
                perturb_ang_mag = perturb_ang_mag>5?perturb_ang_mag:5;
                perturb_trans_mag /= decrease_factor;
                perturb_trans_mag = perturb_trans_mag>1.0?perturb_trans_mag:1.0;
            }

            if(opts.verbose_level!=-1 && ii_inner%opts.verbose_level==0) {
                std::cout << std::setprecision(3) 
                          << "Outer: " << std::left << std::setw(10) << ii_outer 
                          << "Inner: " << std::left << std::setw(10) << ii_inner 
                          << "Total sc: " << std::left << std::setw(10) << best_score 
                          << "RPX sc: " << std::left << std::setw(10) << best_rpx_score
                          << "Metal sc: " << std::left << std::setw(10) << best_interface_metal_score
                          << "Metal Dist: " << std::left << std::setw(10) << best_interface_metal_score / interface_metal_weight
                          << "Sc neig: " << std::left << std::setw(10) << sidechain_neighbors.compute(pose,false,true) 
                          <<  std::endl;
            }
            // after each move, check satisfication and return if good
            bool pass_requirements(true);
            // if(best_score > opts.total_score_cutoff) pass_requirements = false;
            if(pass_requirements && best_rpx_score > opts.inter_chain_score_cutoff) pass_requirements = false;
            if(pass_requirements && best_interface_metal_score > opts.interface_metal_score_cutoff) pass_requirements=false;
            if(pass_requirements && sidechain_neighbors.compute(pose,false,true) < opts.sidechain_neighbor_cutoff) pass_requirements=false;


            if(pass_requirements) {
                std::cout << "Solution found! Eearly stop." << std::endl;
                return true;
            }

        } // inner loop
    } // outer loop

    std::cout << "No solution found! Too bad." << std::endl;
    return false;
}



// I assume the center of the metal(s) is at the origin. There is no point in putting the metal randomly far from the origin,
// right? Of course, you can put the metals anywhere you want, and I can make the center Vec as a variable. But WHY????
bool metal_coordination_mcmc(scene::Pose & pose, scoring::ScoreFunction & sfxn, sampling::FragMover & frag_mover, metric::SidechainNeighbors & sidechain_neighbors,
                                metric::LigandNeighbors & ligand_neighbors, Options const & opts)
{
    std::uniform_real_distribution<Real> runif(0,1);

    Size protein_len = pose.size();
    std::string dssp = opts.fake_ss;
    if(dssp == ""){
        dssp = pose.dssp();
    }
    Size frag_len = frag_mover.frag_length();

    assert(frag_mover.protein_length() == protein_len);

    Size root_index = utils::random_int(-12,12) + Size(protein_len/2);
    //std::cout << "Root Index: " << root_index << std::endl;
    pose.set_root_index(root_index);
    pose.update_coordinates();

    // take a snapshoot right before sampling
    //pose.root_snapshot();
    // score the pose first before taking the snapshot
    // sfxn.score(pose);
    sfxn.score(pose);
    
    pose.snapshot();
    //pose.snapshot(1);

    // global parameters
    Size multi_cool_cycles = opts.mcmc_outer_cycles;
    Real best_score = sfxn.score(pose)/pose.size();
    Real best_score_inter_chain = 99999.0;
    Real best_metal_coordination_score = 99999.0;
    Real best_rif_score = 99999.0;
    Real best_privileged_motif_score = 99999.0;
    
    Size total_steps = opts.mcmc_inner_cycles==-1 ? protein_len * 10 : opts.mcmc_inner_cycles;
    Real decrease_factor = 1.4;
    Real increase_factor = 1.1;
    Size M = total_steps / 10;

    // for C/D-symmetry
    Real C_symmetry_score = -9e9;
    Real D_symmetry_score = -9e9;

    // inter chain pair scores
    // get the symmetry
    std::string symmetry = pose.symmetry();
    Size num_chains = pose.num_chains();
    Size num_chains_C_symmetry = symmetry.at(0)=='C'?num_chains:num_chains/2;

    std::vector<Size> metal_coordination_res;
    Size metal_token;

    Real T = 0.75;

    if(opts.temperature != -1.0) {
        T = opts.temperature/2.0;
    }

    Size mcmc_failed_attempt_counter(0);
    const Size max_failed_attempts_before_change_root(M/10);

    const Size maximum_allowed_rejects = M * 15;
    Size accumulated_rejects = 0;

    for( Size ii_outer=0; ii_outer<multi_cool_cycles; ++ii_outer) {

        // local parameters for MCMC
        // Real T = 3.0 / (ii_outer+1); // multicycle cool good ???? perturb the pose out of local minima?? I don't know.
        T *= 2.0;
        Real loop_weight(10.0);

        Real change_root_prob(0.1);

        // initialize
        frag_mover.bias_loop_sampling_by_dssp(dssp, loop_weight, 1);

        for(Size ii_inner=0; ii_inner<total_steps; ++ii_inner) {
            //std::cout << "Step: " << ii << std::endl;
            Size ires = frag_mover.pick_position(true);

            // dynamicly change root index
            // this works pretty well
            // need to test one large protein to see the speed boost
            // if( protein_len - ires + 1 > ires + frag_len ) {
            //     pose.set_root_index(std::min(ires+frag_len+1, protein_len), false);
            // } else {
            //     pose.set_root_index(1, false);
            // }

            // do fragment insertion
            frag_mover.apply(pose,ires);


            Real cur_score = sfxn.score(pose) / protein_len;
            Real cur_C_symmetry_score =0;
            Real cur_D_symmetry_score =0;
            Real cur_score_interchain = 0;
            Real metal_coordination_score = 0.0;
            Real delta = cur_score - best_score;
            bool pass = pass_metropolis(T, delta, runif(utils::rng()));

            if( pass ) {
            //    std::cout << "Accepted: " << cur_score << std::endl;
                best_score = cur_score;
                if( num_chains>1 ) {
                    best_score_inter_chain = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1)) / protein_len;
                }
                if(opts.rif_table != "") {
                    best_rif_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score(pose);
                }
                if(opts.privileged_motif_hash_tables.size() != 0) {
                    best_privileged_motif_score = opts.privileged_motif_weight * std::static_pointer_cast<scoring::PrivilegedMotifScoreMethod>(sfxn.get_method("privileged_motif"))->score(pose);
                }
                best_metal_coordination_score = opts.metal_coordination_score_weight * std::static_pointer_cast<scoring::MetalCoordinationScoreMethod>(sfxn.get_method("metal"))->score(pose);
                pose.snapshot();

                mcmc_failed_attempt_counter = 0;
                accumulated_rejects = 0;
            } else {
                pose.rollback();
                mcmc_failed_attempt_counter += 1;
                accumulated_rejects +=1;
            }

            // test code
            // Real s1 = sfxn.score_intra_chain(pose)/protein_len;
            // Real s2 = sfxn.score_inter_chain(pose)/protein_len;

            // if(s2>0.0) {
            //     std::stringstream name;
            //     name << "Clash_" << ii_inner << ".pdb";
            //     std::cout << "I found a clash at iteration " << ii_inner << ", with the s2 score " << s2 << std::endl;
            //     pose.dump_pdb(name.str());
            // }
            
            // cur_score = sfxn.score_inter_chain(pose) / protein_len;
            // std::cout << "Inter Chain Score: " << cur_score << std::endl;

            // random chain perturb

            if (num_chains > 1) {
                // total score interchain
                cur_score_interchain = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1)) / protein_len;
            
                // total C symmetry score
                if(num_chains_C_symmetry>1){
                    cur_C_symmetry_score = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,num_chains_C_symmetry)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,num_chains_C_symmetry) ) / protein_len;
                }
                // total D symmetry score
                if(symmetry.at(0)=='D'){
                    cur_D_symmetry_score  = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,num_chains_C_symmetry+1,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,num_chains_C_symmetry+1,-1) ) / protein_len;    
                }
            }

            metal_coordination_score = opts.metal_coordination_score_weight * std::static_pointer_cast<scoring::MetalCoordinationScoreMethod>(sfxn.get_method("metal"))->score_with_tokens(pose, metal_coordination_res, metal_token);
            if( mcmc_failed_attempt_counter > max_failed_attempts_before_change_root && metal_coordination_score < 0 && runif(utils::rng()) < change_root_prob) {
                pose.set_root_index(metal_coordination_res[utils::random_int(0, metal_coordination_res.size()-1)]);
                pose.snapshot(true);

                mcmc_failed_attempt_counter = 0;
            }


            if( metal_coordination_score >= 0
                || (num_chains > 1 && cur_score_interchain>=opts.inter_chain_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor) 
                || (num_chains_C_symmetry>1 && cur_C_symmetry_score>=opts.C_symmetry_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor) 
                || (symmetry.at(0) == 'D' && cur_D_symmetry_score>=opts.D_symmetry_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor) ) {

                Size num_try(100);
                bool random_root_pass = false;
                while(--num_try) {
                    Size root_index = utils::random_int(1,protein_len);
                    //std::cout << "Root Index: " << root_index << std::endl;
                    pose.set_root_index(root_index);
                    pose.random_root(true, true, true, true);

                    if (num_chains > 1) {
                        // total score interchain
                        cur_score_interchain = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1)) / protein_len;
                    
                        // total C symmetry score
                        if(num_chains_C_symmetry>1){
                            cur_C_symmetry_score = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,num_chains_C_symmetry)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,num_chains_C_symmetry) ) / protein_len;
                        }
                        // total D symmetry score
                        if(symmetry.at(0)=='D'){
                            cur_D_symmetry_score  = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,num_chains_C_symmetry+1,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,num_chains_C_symmetry+1,-1) ) / protein_len;    
                        }
                    }
                    metal_coordination_score = opts.metal_coordination_score_weight * std::static_pointer_cast<scoring::MetalCoordinationScoreMethod>(sfxn.get_method("metal"))->score_with_tokens(pose, metal_coordination_res, metal_token);


                    if( (num_chains == 1 || cur_score_interchain<opts.inter_chain_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor) 
                        && (num_chains_C_symmetry == 1 || cur_C_symmetry_score<opts.C_symmetry_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor) 
                        && (symmetry.at(0) != 'D' || cur_D_symmetry_score<opts.D_symmetry_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor)
                        && metal_coordination_score < 0 ) {

                        pose.set_root_index(metal_coordination_res[utils::random_int(0, metal_coordination_res.size()-1)]);

                        best_score = sfxn.score(pose) / protein_len;
                        if( num_chains>1 ) {
                            best_score_inter_chain = cur_score_interchain;
                        }
                        if(opts.rif_table != "") {
                            best_rif_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score(pose);;
                        }
                        if(opts.privileged_motif_hash_tables.size() != 0) {
                            best_privileged_motif_score = opts.privileged_motif_weight * std::static_pointer_cast<scoring::PrivilegedMotifScoreMethod>(sfxn.get_method("privileged_motif"))->score(pose);
                        }
                        best_metal_coordination_score = metal_coordination_score;
                        pose.snapshot(true);
                        random_root_pass = true;
                        mcmc_failed_attempt_counter = 0;
                        break;
                    }
                }
                if( !random_root_pass ) {
                    pose.rollback(true);
                }
            }

            // do fragment insertion again
            ires = frag_mover.pick_position(true);
            frag_mover.apply(pose,ires);


            cur_score = sfxn.score(pose) / protein_len;
            delta = cur_score - best_score;
            pass = pass_metropolis(T, delta, runif(utils::rng()));

            if( pass ) {
            //    std::cout << "Accepted: " << cur_score << std::endl;
                best_score = cur_score;
                if( num_chains>1 ) {
                    best_score_inter_chain = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1)) / protein_len;
                }
                if(opts.rif_table != "") {
                    best_rif_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score(pose);;
                }
                if(opts.privileged_motif_hash_tables.size() != 0) {
                    best_privileged_motif_score = opts.privileged_motif_weight * std::static_pointer_cast<scoring::PrivilegedMotifScoreMethod>(sfxn.get_method("privileged_motif"))->score(pose);
                }
                best_metal_coordination_score = opts.metal_coordination_score_weight * std::static_pointer_cast<scoring::MetalCoordinationScoreMethod>(sfxn.get_method("metal"))->score(pose);
                pose.snapshot();
                mcmc_failed_attempt_counter = 0;
                accumulated_rejects = 0;
            } else {
                pose.rollback();
                mcmc_failed_attempt_counter += 1;
                accumulated_rejects += 1;
            }
                

            if(opts.verbose_level!=-1 && ii_inner%opts.verbose_level==0) {

                Real tmp_cur_score = sfxn.score(pose) / protein_len;

                Real tmp_cur_rpx = std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->score(pose)/pose.size();
                Real tmp_cur_clash = std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->score(pose)/pose.size();

                Real tmp_cur_sidechain_neighbors = sidechain_neighbors.compute(pose);

                Real tmp_cur_target_clash = opts.context_clash_weight * std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash"))->score(pose);

                Real tmp_cur_metal_coordination_score = opts.metal_coordination_score_weight * std::static_pointer_cast<scoring::MetalCoordinationScoreMethod>(sfxn.get_method("metal"))->score(pose);

                std::cout <<std::setprecision(3) 
                                                  << "Outer: " << std::left << std::setw(4) << ii_outer 
                                                  << "Inner: " << std::left << std::setw(7) << ii_inner 
                                                  << "TotalSc: " << std::left << std::setw(10) << tmp_cur_score 
                                                  << "RPX: " << std::left << std::setw(10) << tmp_cur_rpx 
                                                  << "Clash: " << std::left << std::setw(10) << tmp_cur_clash;
                                                  if(opts.rif_table!="") {
                                                    Real tmp_cur_rif_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score(pose);
                                                    std::cout << std::setprecision(3) << "Rif: " << std::left << std::setw(10) << tmp_cur_rif_score;
                                                  }
                                                  if(opts.privileged_motif_hash_tables.size()!=0) {
                                                    Real tmp_cur_privileged_motif_score = opts.privileged_motif_weight * std::static_pointer_cast<scoring::PrivilegedMotifScoreMethod>(sfxn.get_method("privileged_motif"))->score(pose);
                                                    std::cout << std::setprecision(3) << "PvMotif: " << std::left << std::setw(10) << tmp_cur_privileged_motif_score;
                                                  }
                                                  if(opts.ligand_pdb != "") {
                                                    Real tmp_cur_ligand_neighbors = ligand_neighbors.compute(pose);
                                                    std::cout << std::setprecision(3) << "LN: " <<std::left << std::setw(10) << tmp_cur_ligand_neighbors;
                                                  }
                                                  std::cout << std::setprecision(3) << "TargetClash: " << std::left << std::setw(10) << tmp_cur_target_clash
                                                  << "SC: " << std::left << std::setw(10) << tmp_cur_sidechain_neighbors
                                                  << "MetalCoordination: " << std::left << std::setw(10) << tmp_cur_metal_coordination_score
                                                  <<  std::endl;
                
            }


            // cur_score = sfxn.score_intra_chain(pose) / protein_len;
            // std::cout << "The inter chain score: " << cur_score << std::endl;
            // cur_score = sfxn.score_inter_chain(pose) / protein_len;
            // std::cout << "The inter chain score: " << cur_score << std::endl;
            
            // if(false) {
            //     std::cout << best_score << std::endl;
            // }

            if( ii_inner % M == 0 && ii_inner != 0 ) {
                T /= decrease_factor;
                loop_weight /= decrease_factor;
                loop_weight = loop_weight>=1?loop_weight:1;
                frag_mover.bias_loop_sampling_by_dssp(dssp, loop_weight, 1);
            }

            if(opts.early_stop && accumulated_rejects > maximum_allowed_rejects) {
                std::cout << "No acceptions after " << accumulated_rejects << " moves. Stop!" << std::endl;
                return false;
            }

            // after each move, check satisfication and return if good
            bool pass_requirements(true);
            if(best_score > opts.total_score_cutoff) pass_requirements = false;
            if(pass_requirements && num_chains>1 && best_score_inter_chain > opts.inter_chain_score_cutoff) pass_requirements=false;
            if(pass_requirements && best_metal_coordination_score > opts.metal_coordination_score_cutoff) pass_requirements=false;
            if(pass_requirements && opts.rif_table != "" && best_rif_score > opts.rif_score_cutoff) { pass_requirements=false; }
            if(pass_requirements && opts.privileged_motif_hash_tables.size()!=0 && best_privileged_motif_score > opts.privileged_motif_score_cutoff) {pass_requirements=false; }
            if(pass_requirements && opts.ligand_pdb!="" && opts.ligand_neighbors_cutoff != -1 && ligand_neighbors.compute(pose) < opts.ligand_neighbors_cutoff) pass_requirements=false;
            if(pass_requirements && sidechain_neighbors.compute(pose) < opts.sidechain_neighbor_cutoff) pass_requirements=false;

            if(pass_requirements) {
                std::cout << "Solution found! Eearly stop." << std::endl;

                Real tmp_cur_score = sfxn.score(pose) / protein_len;

                Real tmp_cur_rpx = std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->score(pose)/pose.size();
                Real tmp_cur_clash = std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->score(pose)/pose.size();

                Real tmp_cur_sidechain_neighbors = sidechain_neighbors.compute(pose);

                Real tmp_cur_target_clash = opts.context_clash_weight * std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash"))->score(pose);

                Real tmp_cur_metal_coordination_score = opts.metal_coordination_score_weight * std::static_pointer_cast<scoring::MetalCoordinationScoreMethod>(sfxn.get_method("metal"))->score(pose);

                std::cout <<std::setprecision(3)  << "Success: TotalSc: " << std::left << std::setw(10) << tmp_cur_score 
                                                  << "RPX: " << std::left << std::setw(10) << tmp_cur_rpx 
                                                  << "Clash: " << std::left << std::setw(10) << tmp_cur_clash;
                                                  if(opts.rif_table!="") {
                                                    Real tmp_cur_rif_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score(pose);
                                                    std::cout << std::setprecision(3) << "Rif: " << std::left << std::setw(10) << tmp_cur_rif_score;
                                                  }
                                                  if(opts.privileged_motif_hash_tables.size()!=0) {
                                                    Real tmp_cur_privileged_motif_score = opts.privileged_motif_weight * std::static_pointer_cast<scoring::PrivilegedMotifScoreMethod>(sfxn.get_method("privileged_motif"))->score(pose);
                                                    std::cout << std::setprecision(3) << "PvMotif: " << std::left << std::setw(10) << tmp_cur_privileged_motif_score;
                                                  }
                                                  if(opts.ligand_pdb != "") {
                                                    Real tmp_cur_ligand_neighbors = ligand_neighbors.compute(pose);
                                                    std::cout << std::setprecision(3) << "LN: " <<std::left << std::setw(10) << tmp_cur_ligand_neighbors;
                                                  }
                                                  std::cout <<std::setprecision(3) << "TargetClash: " << std::left << std::setw(10) << tmp_cur_target_clash
                                                  << "SC: " << std::left << std::setw(10) << tmp_cur_sidechain_neighbors
                                                  << "MetalCoordination: " << std::left << std::setw(10) << tmp_cur_metal_coordination_score
                                                  <<  std::endl;

                return true;
            }

        } // inner loop
    } // outer loop

    std::cout << "No solution found! Too bad." << std::endl;

    Real tmp_cur_score = sfxn.score(pose) / protein_len;

    Real tmp_cur_rpx = std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->score(pose)/pose.size();
    Real tmp_cur_clash = std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->score(pose)/pose.size();

    Real tmp_cur_sidechain_neighbors = sidechain_neighbors.compute(pose);

    Real tmp_cur_target_clash = opts.context_clash_weight * std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash"))->score(pose);

    Real tmp_cur_metal_coordination_score = opts.metal_coordination_score_weight * std::static_pointer_cast<scoring::MetalCoordinationScoreMethod>(sfxn.get_method("metal"))->score(pose);


    std::cout <<std::setprecision(3)  << "Failed: TotalSc: " << std::left << std::setw(10) << tmp_cur_score 
                                      << "RPX: " << std::left << std::setw(10) << tmp_cur_rpx 
                                      << "Clash: " << std::left << std::setw(10) << tmp_cur_clash;
                                      if(opts.rif_table!="") {
                                        Real tmp_cur_rif_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score(pose);
                                        std::cout << std::setprecision(3) << "Rif: " << std::left << std::setw(10) << tmp_cur_rif_score;
                                      }
                                      if(opts.privileged_motif_hash_tables.size()!=0) {
                                        Real tmp_cur_privileged_motif_score = opts.privileged_motif_weight * std::static_pointer_cast<scoring::PrivilegedMotifScoreMethod>(sfxn.get_method("privileged_motif"))->score(pose);
                                        std::cout << std::setprecision(3) << "PvMotif: " << std::left << std::setw(10) << tmp_cur_privileged_motif_score;
                                      }
                                      if(opts.ligand_pdb != "") {
                                        Real tmp_cur_ligand_neighbors = ligand_neighbors.compute(pose);
                                        std::cout << std::setprecision(3) << "LN: " <<std::left << std::setw(10) << tmp_cur_ligand_neighbors;
                                      }
                                      std::cout <<std::setprecision(3) << "TargetClash: " << std::left << std::setw(10) << tmp_cur_target_clash
                                      << "SC: " << std::left << std::setw(10) << tmp_cur_sidechain_neighbors
                                      << "MetalCoordination: " << std::left << std::setw(10) << tmp_cur_metal_coordination_score
                                      <<  std::endl;

    return false;
}


// fiber mcmc

bool interface_metal_fiber_mcmc(scene::Pose & pose, scene::InterfaceMetalOP & interface_metal_p1, scene::InterfaceMetalOP & interface_metal_p2, Size & metal1_chain_idx, Size & metal2_chain_idx, scoring::ScoreFunction & sfxn, metric::SidechainNeighbors & sidechain_neighbors, Options const & opts)
{
    std::uniform_real_distribution<Real> runif(0,1);

    Size protein_len = opts.len;


    Size root_index = utils::random_int(-30,30) + Size(protein_len/2);
    //std::cout << "Root Index: " << root_index << std::endl;
    pose.set_root_index(root_index);
    pose.update_coordinates();

    // random root
    pose.random_root(true, true, true, true);
    pose.random_jump(opts.fiber_rotation_angle_expectation, opts.fiber_raise_dist_expectation);

    // which chain??
    metal1_chain_idx = std::stoi(opts.symmetry.substr(1));
    // Size tmp_ii = utils::random_int(1,3);
    // zn2_chain_idx = (tmp_ii==1)?2:(tmp_ii==2?zn1_chain_idx-1:zn1_chain_idx+1);
    metal2_chain_idx = 2;

    interface_metal_p1->update_metal_position(pose.stub(1, 1), pose.stub(1, metal1_chain_idx));
    interface_metal_p2->update_metal_position(pose.stub(1, 1), pose.stub(1, metal2_chain_idx));

    // take a snapshoot right before sampling
    //pose.root_snapshot();
    // score the pose first before taking the snapshot
    // sfxn.score(pose);
    sfxn.score(pose);
    
    pose.snapshot();
    interface_metal_p1->snapshot();
    interface_metal_p2->snapshot();
    //pose.snapshot(1);

    // global parameters
    Size multi_cool_cycles = opts.mcmc_outer_cycles;
    Size total_steps = opts.mcmc_inner_cycles==-1 ? protein_len * 10 : opts.mcmc_inner_cycles;
    Size M = total_steps / 10;


    Real best_score = 9e9;
    Real best_interface_metal1_score = 9e9;
    Real best_interface_metal2_score = 9e9;
    Real best_rpx_score = 9e9;
    
    Real decrease_factor = 1.4;
    Real increase_factor = 1.1;

    Real num_interface_metals = opts.num_interface_metals;
    Real interface_metal_weight = opts.interface_metal_distance_optimization_weight;

    
    Real T = 0.5;
    for( Size ii_outer=0; ii_outer<multi_cool_cycles; ++ii_outer) {

        // local parameters for MCMC
        // Real T = 3.0 / (ii_outer+1); // multicycle cool good ???? perturb the pose out of local minima?? I don't know.
        T *= 2.0;

        Real perturb_ang_mag(30);
        Real perturb_trans_mag(8);

        for(Size ii_inner=0; ii_inner<total_steps; ++ii_inner) {
            

            pose.perturb_root(perturb_ang_mag, perturb_trans_mag);
            pose.perturb_jump(perturb_ang_mag, perturb_trans_mag);
            interface_metal_p1->update_metal_position(pose.stub(1, 1), pose.stub(1, metal1_chain_idx));
            interface_metal_p2->update_metal_position(pose.stub(1, 1), pose.stub(1, metal2_chain_idx));


            Real rpx_sc = sfxn.score(pose)/protein_len;
            Real interface_metal1_sc = interface_metal_p1->topN_pair_distance(opts.num_interface_metals) * interface_metal_weight;
            Real interface_metal2_sc = interface_metal_p2->topN_pair_distance(opts.num_interface_metals) * interface_metal_weight;



            


            Real cur_score = rpx_sc + interface_metal1_sc + interface_metal2_sc;
            Real delta = cur_score - best_score;
            bool pass = pass_metropolis(T, delta, runif(utils::rng()));

            if( pass ) {
            //    std::cout << "Accepted: " << cur_score << std::endl;
                best_score = cur_score;
                best_rpx_score = rpx_sc;
                best_interface_metal1_score = interface_metal1_sc;
                best_interface_metal2_score = interface_metal2_sc;
                pose.snapshot(true);
                interface_metal_p1->snapshot(false);
                interface_metal_p2->snapshot(false);
            } else {
                pose.rollback(true);
                interface_metal_p1->rollback(false);
                interface_metal_p2->rollback(false);
            }


            // reroot
            Real chain1_rpx_sc = std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,metal1_chain_idx,metal1_chain_idx);
            Real chain2_rpx_sc = std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,metal2_chain_idx,metal2_chain_idx);

            if(best_rpx_score>opts.inter_chain_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor ||
                chain1_rpx_sc == 0 ||
                chain2_rpx_sc == 0 ) {
                bool perturb_root_pass = false;

                Size num_try(100);
                while(--num_try){
                    pose.random_root(true, true, true, true);
                    pose.random_jump(opts.fiber_rotation_angle_expectation, opts.fiber_raise_dist_expectation);

                    rpx_sc = sfxn.score(pose)/protein_len;
                    chain1_rpx_sc = std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,metal1_chain_idx,metal1_chain_idx);
                    chain2_rpx_sc = std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,metal2_chain_idx,metal2_chain_idx);
                    
                    if(rpx_sc<=opts.inter_chain_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor &&
                        chain1_rpx_sc != 0 &&
                        chain2_rpx_sc != 0 ){
                        interface_metal_p1->update_metal_position(pose.stub(1, 1), pose.stub(1, metal1_chain_idx));
                        interface_metal1_sc = interface_metal_p1->topN_pair_distance(opts.num_interface_metals) * interface_metal_weight;
                        interface_metal_p2->update_metal_position(pose.stub(1, 1), pose.stub(1, metal2_chain_idx));
                        interface_metal2_sc = interface_metal_p2->topN_pair_distance(opts.num_interface_metals) * interface_metal_weight;
                        best_rpx_score = rpx_sc;
                        best_interface_metal1_score = interface_metal1_sc;
                        best_interface_metal2_score = interface_metal2_sc;
                        best_score = best_rpx_score + best_interface_metal1_score + best_interface_metal2_score;
                        pose.snapshot(true);
                        interface_metal_p1->snapshot(false);
                        interface_metal_p2->snapshot(false);
                        perturb_root_pass = true;
                        break;
                    }  
                }
                if( !perturb_root_pass ) {
                    pose.rollback(true);
                }
            } else {

                pose.perturb_root(perturb_ang_mag, perturb_trans_mag);
                interface_metal_p1->update_metal_position(pose.stub(1, 1), pose.stub(1, metal1_chain_idx));
                interface_metal_p2->update_metal_position(pose.stub(1, 1), pose.stub(1, metal2_chain_idx));

                rpx_sc = sfxn.score(pose)/protein_len;
                interface_metal1_sc = interface_metal_p1->topN_pair_distance(opts.num_interface_metals) * interface_metal_weight;
                interface_metal2_sc = interface_metal_p2->topN_pair_distance(opts.num_interface_metals) * interface_metal_weight;
                cur_score = rpx_sc + interface_metal1_sc + interface_metal2_sc;

                delta = cur_score - best_score;
                pass = pass_metropolis(T, delta, runif(utils::rng()));

                if( pass ) {
                //    std::cout << "Accepted: " << cur_score << std::endl;
                    best_score = cur_score;
                    best_rpx_score = rpx_sc;
                    best_interface_metal1_score = interface_metal1_sc;
                    best_interface_metal2_score = interface_metal2_sc;
                    pose.snapshot(true);
                    interface_metal_p1->snapshot(false);
                    interface_metal_p2->snapshot(false);
                } else {
                    pose.rollback(true);
                    interface_metal_p1->rollback(false);
                    interface_metal_p2->rollback(false);
                }

            }

            
            if( ii_inner % M == 0 && ii_inner != 0 ) {
                T /= decrease_factor;
                perturb_ang_mag /= decrease_factor;
                perturb_ang_mag = perturb_ang_mag>5?perturb_ang_mag:5;
                perturb_trans_mag /= decrease_factor;
                perturb_trans_mag = perturb_trans_mag>1.0?perturb_trans_mag:1.0;
            }

            if(opts.verbose_level!=-1 && ii_inner%opts.verbose_level==0) {
                std::cout << std::setprecision(3) 
                          << "Outer: " << std::left << std::setw(10) << ii_outer 
                          << "Inner: " << std::left << std::setw(10) << ii_inner 
                          << "Total sc: " << std::left << std::setw(10) << best_score 
                          << "RPX sc: " << std::left << std::setw(10) << best_rpx_score
                          << "Metal1 Sc: " << std::left << std::setw(10) << best_interface_metal1_score
                          << "Metal1 Dist: " << std::left << std::setw(10) << best_interface_metal1_score / interface_metal_weight
                          << "Metal2 sc: " << std::left << std::setw(10) << best_interface_metal2_score
                          << "Metal2 Dist: " << std::left << std::setw(10) << best_interface_metal2_score / interface_metal_weight
                          << "Sc Neigh: " << std::left << std::setw(10) << sidechain_neighbors.compute(pose,false,true) 
                          <<  std::endl;
            }
            // after each move, check satisfication and return if good
            bool pass_requirements(true);
            if(best_score > opts.total_score_cutoff) pass_requirements = false;
            if(pass_requirements && best_rpx_score > opts.inter_chain_score_cutoff) pass_requirements = false;
            if(pass_requirements && best_interface_metal1_score > opts.interface_metal_score_cutoff) pass_requirements=false;
            if(pass_requirements && best_interface_metal2_score > opts.interface_metal_score_cutoff) pass_requirements=false;
            if(pass_requirements && 
                std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,metal1_chain_idx,metal1_chain_idx) == 0) {
                pass_requirements = false;
            }
            if(pass_requirements && 
                std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,metal2_chain_idx,metal2_chain_idx) == 0) {
                pass_requirements = false;
            }
            if(pass_requirements && sidechain_neighbors.compute(pose,false,true) < opts.sidechain_neighbor_cutoff) pass_requirements=false;


            if(pass_requirements) {
                std::cout << "Solution found! Eearly stop." << std::endl;
                return true;
            }

        } // inner loop
    } // outer loop

    std::cout << "No solution found! Too bad." << std::endl;
    return false;
}

}

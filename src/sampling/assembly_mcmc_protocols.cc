#include "sampling/mcmc_protocols.hh"

#include "utils/utils.hh"
#include "basic/types.hh"
#include "scene/Pose.hh"
#include "sampling/FragMover.hh"
#include "scoring/ScoreFunction.hh"
#include "utils/random_util.hh"

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


bool mcmc_pH_assembly(scene::Pose & pose, scoring::ScoreFunction & sfxn, Options const & opts)
{
    std::uniform_real_distribution<Real> runif(0,1);

    Size protein_len = pose.size();

    Size root_index = utils::random_int(-10,10) + Size(protein_len/2);
    //std::cout << "Root Index: " << root_index << std::endl;
    pose.set_root_index(root_index);
    pose.update_coordinates();

    // take a snapshoot right before sampling
    //pose.root_snapshot();
    // score the pose first before taking the snapshot
    sfxn.score(pose);
    pose.snapshot();
    //pose.snapshot(1);

    // global parameters
    Size multi_cool_cycles = opts.mcmc_outer_cycles;
    Real best_score = 9999999.0;
    Real best_score_inter_chain = 99999.0;
    Real best_privileged_pair_score = 99999.0;
    
    Size total_steps = opts.mcmc_inner_cycles==-1 ? protein_len * 10 : opts.mcmc_inner_cycles;
    Real decrease_factor = 1.4;
    Size M = total_steps / 10;

    Real T = 1.5;

    for( Size ii_outer=1; ii_outer<=multi_cool_cycles; ++ii_outer) {

        // local parameters for MCMC
        T *= 2.0; // multicycle cool good ???? perturb the pose out of local minima?? I don't know.

        Real perturb_ang_mag(30);
        Real perturb_trans_mag(10);

        for(Size ii_inner=0; ii_inner<total_steps; ++ii_inner) {

            pose.perturb_root(perturb_ang_mag, perturb_trans_mag);

            Real cur_score = sfxn.score(pose) / protein_len;
            Real cur_score_interchain(0.0);
            
            // Real z_cen     = std::fabs(pose_center(pose)[2]);

            //
            Real delta = cur_score - best_score;

            bool pass = pass_metropolis(T, delta, runif(utils::rng()));

            if( pass ) {
            //    std::cout << "Accepted: " << cur_score << std::endl;
                best_score = cur_score;
                best_score_inter_chain = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1)) / protein_len;
                best_privileged_pair_score = opts.privileged_pair_weight * std::static_pointer_cast<scoring::PrivilegedPairScoreMethod>(sfxn.get_method("privileged_pair"))->score(pose);
                pose.snapshot(true);
            } else {
                pose.rollback(true);
            }

            cur_score_interchain = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1) ) / protein_len;
            

            if (cur_score_interchain>=opts.inter_chain_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor) {
                //
                bool perturb_root_pass = false;
                Size num_try(100);
                while(--num_try){
                    pose.random_root(true, true, true, true);

                    // check again
                    cur_score_interchain =(std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1) ) / protein_len;
                    if (cur_score_interchain<opts.inter_chain_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor) {
                        best_score = sfxn.score(pose) / protein_len;
                        best_score_inter_chain = cur_score_interchain;
                        best_privileged_pair_score = opts.privileged_pair_weight * std::static_pointer_cast<scoring::PrivilegedPairScoreMethod>(sfxn.get_method("privileged_pair"))->score(pose);
                        pose.snapshot(true);
                        perturb_root_pass = true;
                        break;
                    }
                }
                if( !perturb_root_pass ) {
                    pose.rollback(true);
                }
            } else {

                pose.perturb_root(perturb_ang_mag, perturb_trans_mag);
                cur_score = sfxn.score(pose) / protein_len;
                delta = cur_score - best_score;
                pass = pass_metropolis(T, delta, runif(utils::rng()));

                if( pass ) {
                //    std::cout << "Accepted: " << cur_score << std::endl;
                    best_score = cur_score;
                    best_score_inter_chain = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1)) / protein_len;
                    best_privileged_pair_score = opts.privileged_pair_weight * std::static_pointer_cast<scoring::PrivilegedPairScoreMethod>(sfxn.get_method("privileged_pair"))->score(pose);
                    pose.snapshot(true);
                } else {
                    pose.rollback(true);
                }
            }


            if( ii_inner % M == 0 && ii_inner != 0 ) {
                T /= decrease_factor;
                //change_root_prob /= decrease_factor;
                perturb_ang_mag /= decrease_factor;
                perturb_ang_mag = perturb_ang_mag>5?perturb_ang_mag:5;
                perturb_trans_mag /= decrease_factor;
                perturb_trans_mag = perturb_trans_mag>0.5?perturb_trans_mag:0.5;
            }
            // after each move, check satisfication and return if good
            bool pass_requirements(true);
            if(pass_requirements && best_score_inter_chain > opts.inter_chain_score_cutoff) { pass_requirements=false; }
            if(pass_requirements && best_privileged_pair_score > opts.privileged_pair_score_cutoff) { pass_requirements=false; }
            if(pass_requirements && utils::sidechain_neighbors(pose,false,true) < opts.sidechain_neighbor_cutoff) { pass_requirements=false; }

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
bool mcmc_light_heterodimer(scene::Pose & pose, scoring::ScoreFunction & sfxn, metric::SidechainNeighbors & sidechain_neighbors, Options const & opts)
{
    std::uniform_real_distribution<Real> runif(0,1);

    Size total_len = opts.chain1_len + opts.chain2_len;

    // randomly select a position as a root
    Size root_index = utils::random_int(1,pose.size());
    //std::cout << "Root Index: " << root_index << std::endl;
    pose.set_root_index(root_index);
    pose.update_coordinates();

    // random root
    pose.random_root(true, true, true, true);

    sfxn.score(pose);    
    pose.snapshot();

    // global parameters
    Size multi_cool_cycles = opts.mcmc_outer_cycles;
    Size total_steps = opts.mcmc_inner_cycles==-1 ? total_len * 10 : opts.mcmc_inner_cycles;
    Size M = total_steps / 10;


    Real best_score = 9e9;
    Real best_privileged_motif_score = 9e9;
    Real best_inter_chain_score = 9e9;
    
    Real decrease_factor = 1.4;
    Real increase_factor = 1.1;
    
    Real T = 0.5;
    for( Size ii_outer=0; ii_outer<multi_cool_cycles; ++ii_outer) {

        // local parameters for MCMC
        // Real T = 3.0 / (ii_outer+1); // multicycle cool good ???? perturb the pose out of local minima?? I don't know.
        T *= 2.0;

        Real perturb_ang_mag(30);
        Real perturb_trans_mag(5);

        for(Size ii_inner=0; ii_inner<total_steps; ++ii_inner) {

            pose.perturb_root(perturb_ang_mag, perturb_trans_mag);

            Real cur_score = sfxn.score(pose) / total_len;

            Real delta = cur_score - best_score;
            bool pass = pass_metropolis(T, delta, runif(utils::rng()));

            if( pass ) {
            //    std::cout << "Accepted: " << cur_score << std::endl;
                best_score = cur_score;
                best_privileged_motif_score = opts.privileged_motif_weight * std::static_pointer_cast<scoring::PrivilegedMotifScoreMethod>(sfxn.get_method("privileged_motif"))->score(pose);
                best_inter_chain_score = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->score(pose)+std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("clash"))->score(pose)) / total_len;
                pose.snapshot(true);
            } else {
                pose.rollback(true);
            }

            
            if(best_score > 50 || best_score == 0.0) { // clash too much or no contact

                Size num_try(100);
                bool perturb_root_pass = false;
                while(--num_try){
                    pose.random_root(true, true, true, true);
                    cur_score = sfxn.score(pose) / total_len;
                    if(cur_score <= 50 && cur_score != 0.0){
                        best_score = cur_score;
                        best_privileged_motif_score = opts.privileged_motif_weight * std::static_pointer_cast<scoring::PrivilegedMotifScoreMethod>(sfxn.get_method("privileged_motif"))->score(pose);
                        best_inter_chain_score = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->score(pose)+std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("clash"))->score(pose)) / total_len;
                        pose.snapshot(true);
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
                          << "; Total score: " << best_score 
                          << "; RPX score: " << best_inter_chain_score 
                          << "; Motif score: " << best_privileged_motif_score 
                          <<  std::endl;
            }
            // after each move, check satisfication and return if good
            bool pass_requirements(true);
            if(best_score > opts.total_score_cutoff) pass_requirements = false;
            if(pass_requirements && best_inter_chain_score > opts.inter_chain_score_cutoff) pass_requirements = false;
            if(pass_requirements && best_privileged_motif_score > opts.privileged_motif_score_cutoff) pass_requirements=false;

            if(pass_requirements) {
                std::cout << "Solution found! Eearly stop." << std::endl;
                return true;
            }

        } // inner loop
    } // outer loop

    std::cout << "No solution found! Too bad." << std::endl;
    return false;
}

bool mcmc_light_assembly(scene::Pose & pose, scoring::ScoreFunction & sfxn, metric::SidechainNeighbors & sidechain_neighbors, Options const & opts)
{
    std::uniform_real_distribution<Real> runif(0,1);

    Size protein_len = pose.size();

    Size root_index = utils::random_int(-10,10) + Size(protein_len/2);
    //std::cout << "Root Index: " << root_index << std::endl;
    pose.set_root_index(root_index);
    pose.update_coordinates();

    // take a snapshoot right before sampling
    //pose.root_snapshot();
    // score the pose first before taking the snapshot
    sfxn.score(pose);
    pose.snapshot();
    //pose.snapshot(1);

    // global parameters
    Size multi_cool_cycles = opts.mcmc_outer_cycles;
    Real best_score = 9999999.0;
    Real best_score_inter_chain = 99999.0;
    Real best_privileged_interface_motif_score = 99999.0;
    
    Size total_steps = opts.mcmc_inner_cycles==-1 ? protein_len * 10 : opts.mcmc_inner_cycles;
    Real decrease_factor = 1.4;
    Size M = total_steps / 10;

    Size mcmc_failed_attempt_exit_counter(0);
    const Size max_failed_attempts_before_quit(total_steps*2);

    Real T = 1.5;

    for( Size ii_outer=1; ii_outer<=multi_cool_cycles; ++ii_outer) {

        // local parameters for MCMC
        T *= 2.0; // multicycle cool good ???? perturb the pose out of local minima?? I don't know.

        Real perturb_ang_mag(30);
        Real perturb_trans_mag(10);

        for(Size ii_inner=0; ii_inner<total_steps; ++ii_inner) {

            pose.perturb_root(perturb_ang_mag, perturb_trans_mag);

            Real cur_score = sfxn.score(pose) / protein_len;
            Real cur_score_interchain(0.0);
            
            // Real z_cen     = std::fabs(pose_center(pose)[2]);

            //
            Real delta = cur_score - best_score;

            bool pass = pass_metropolis(T, delta, runif(utils::rng()));

            if( pass ) {
            //    std::cout << "Accepted: " << cur_score << std::endl;
                best_score = cur_score;
                best_score_inter_chain = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1)) / protein_len;
                best_privileged_interface_motif_score = opts.privileged_interface_motif_weight * std::static_pointer_cast<scoring::PrivilegedInterfaceMotifScoreMethod>(sfxn.get_method("privileged_interface_motif"))->score(pose);
                pose.snapshot(true);
                mcmc_failed_attempt_exit_counter = 0;
            } else {
                pose.rollback(true);
                ++mcmc_failed_attempt_exit_counter;
            }

            cur_score_interchain = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1) ) / protein_len;
            

            if (cur_score_interchain>=opts.inter_chain_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor) {
                //
                bool perturb_root_pass = false;
                Size num_try(100);
                while(--num_try){
                    pose.random_root(true, true, true, true);

                    // check again
                    cur_score_interchain =(std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1) ) / protein_len;
                    if (cur_score_interchain<opts.inter_chain_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor) {
                        best_score = sfxn.score(pose) / protein_len;
                        best_score_inter_chain = cur_score_interchain;
                        best_privileged_interface_motif_score = opts.privileged_interface_motif_weight * std::static_pointer_cast<scoring::PrivilegedInterfaceMotifScoreMethod>(sfxn.get_method("privileged_interface_motif"))->score(pose);
                        pose.snapshot(true);
                        perturb_root_pass = true;
                        mcmc_failed_attempt_exit_counter=0;
                        break;
                    }
                }
                if( !perturb_root_pass ) {
                    pose.rollback(true);
                }
            } else {

                pose.perturb_root(perturb_ang_mag, perturb_trans_mag);
                cur_score = sfxn.score(pose) / protein_len;
                delta = cur_score - best_score;
                pass = pass_metropolis(T, delta, runif(utils::rng()));

                if( pass ) {
                //    std::cout << "Accepted: " << cur_score << std::endl;
                    best_score = cur_score;
                    best_score_inter_chain = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1)) / protein_len;
                    best_privileged_interface_motif_score = opts.privileged_interface_motif_weight * std::static_pointer_cast<scoring::PrivilegedInterfaceMotifScoreMethod>(sfxn.get_method("privileged_interface_motif"))->score(pose);
                    pose.snapshot(true);
                    mcmc_failed_attempt_exit_counter=0;

                } else {
                    pose.rollback(true);
                    ++mcmc_failed_attempt_exit_counter;
                }
            }

            if( opts.early_stop && mcmc_failed_attempt_exit_counter > max_failed_attempts_before_quit ) {
                std::cout << "Got trapped in local minima, STOP!!!" << std::endl;
                return false;
            }

            if(opts.verbose_level!=-1 && ii_inner%opts.verbose_level==0) {

                Real tmp_score = sfxn.score(pose) / protein_len;

                Real tmp_rpx = std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)/pose.size();
                Real tmp_clash = std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1)/pose.size();

                Real tmp_sidechain_neighbors = sidechain_neighbors.compute(pose, false, true);

                Real tmp_privileged_interface_motif_score = opts.privileged_interface_motif_weight * std::static_pointer_cast<scoring::PrivilegedInterfaceMotifScoreMethod>(sfxn.get_method("privileged_interface_motif"))->score(pose);

                std::cout <<std::setprecision(4) 
                                                  << "Accepted <===> Outer: " << std::left << std::setw(5) << ii_outer 
                                                  << "Inner: " << std::left << std::setw(10) << ii_inner 
                                                  << "TotalSc: " << std::left << std::setw(10) << tmp_score 
                                                  << "RPX: " << std::left << std::setw(10) << tmp_rpx 
                                                  << "Clash: " << std::left << std::setw(10) << tmp_clash 
                                                  << "SCN: " << std::left << std::setw(10) << tmp_sidechain_neighbors
                                                  << "PIM: " << std::left << std::setw(10) << tmp_privileged_interface_motif_score
                                                  << "Counter: " << std::left << std::setw(10) << mcmc_failed_attempt_exit_counter
                                                  <<  std::endl;
                    
            }

            if( ii_inner % M == 0 && ii_inner != 0 ) {
                T /= decrease_factor;
                //change_root_prob /= decrease_factor;
                perturb_ang_mag /= decrease_factor;
                perturb_ang_mag = perturb_ang_mag>5?perturb_ang_mag:5;
                perturb_trans_mag /= decrease_factor;
                perturb_trans_mag = perturb_trans_mag>0.5?perturb_trans_mag:0.5;
            }
            // after each move, check satisfication and return if good
            bool pass_requirements(true);
            if(pass_requirements && best_score_inter_chain > opts.inter_chain_score_cutoff) { pass_requirements=false; }
            if(pass_requirements && best_privileged_interface_motif_score > opts.privileged_interface_motif_score_cutoff) { pass_requirements=false; }
            if(pass_requirements && sidechain_neighbors.compute(pose,false,true) < opts.sidechain_neighbor_cutoff) { pass_requirements=false; }

            if(pass_requirements) {
                std::cout << "Solution found! Eearly stop." << std::endl;

                Real tmp_score = sfxn.score(pose) / protein_len;

                Real tmp_rpx = std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)/pose.size();
                Real tmp_clash = std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1)/pose.size();

                Real tmp_sidechain_neighbors = sidechain_neighbors.compute(pose, false, true);

                Real tmp_privileged_interface_motif_score = opts.privileged_interface_motif_weight * std::static_pointer_cast<scoring::PrivilegedInterfaceMotifScoreMethod>(sfxn.get_method("privileged_interface_motif"))->score(pose);

                std::cout <<std::setprecision(4) 
                          << "Success: TotalSc: " << std::left << std::setw(10) << tmp_score 
                          << "RPX: " << std::left << std::setw(10) << tmp_rpx 
                          << "Clash: " << std::left << std::setw(10) << tmp_clash 
                          << "SCN: " << std::left << std::setw(10) << tmp_sidechain_neighbors
                          << "PIM: " << std::left << std::setw(10) << tmp_privileged_interface_motif_score
                          <<  std::endl;

                return true;
            }
        } // inner loop
    } // outer loop
    std::cout << "No solution found! Too bad." << std::endl;

    Real tmp_score = sfxn.score(pose) / protein_len;

    Real tmp_rpx = std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)/pose.size();
    Real tmp_clash = std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1)/pose.size();

    Real tmp_sidechain_neighbors = sidechain_neighbors.compute(pose, false, true);

    Real tmp_privileged_interface_motif_score = opts.privileged_interface_motif_weight * std::static_pointer_cast<scoring::PrivilegedInterfaceMotifScoreMethod>(sfxn.get_method("privileged_interface_motif"))->score(pose);

    std::cout <<std::setprecision(4) 
              << "Failed: TotalSc: " << std::left << std::setw(10) << tmp_score 
              << "RPX: " << std::left << std::setw(10) << tmp_rpx 
              << "Clash: " << std::left << std::setw(10) << tmp_clash 
              << "SCN: " << std::left << std::setw(10) << tmp_sidechain_neighbors
              << "PIM: " << std::left << std::setw(10) << tmp_privileged_interface_motif_score
              <<  std::endl;

    return false;
}


}

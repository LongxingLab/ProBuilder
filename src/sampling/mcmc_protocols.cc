#include "sampling/mcmc_protocols.hh"

#include "utils/utils.hh"
#include "basic/types.hh"
#include "scene/Pose.hh"
#include "sampling/FragMover.hh"
#include "scoring/ScoreFunction.hh"
#include "utils/random_util.hh"
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

bool mcmc(scene::Pose & pose, scoring::ScoreFunction & sfxn, sampling::FragMover & frag_mover, Options const & opts)
{
    const bool perturb_root = opts.symmetry != "C1";
    if(opts.dump_trajactory)pose.append_trajctory();
    std::uniform_real_distribution<Real> runif(0,1);

    Size protein_len = pose.size();
    std::string dssp = opts.fake_ss;
    if(dssp == ""){
        dssp = pose.dssp();
    }
    Size frag_len = frag_mover.frag_length();

    assert(frag_mover.protein_length() == protein_len);

    Size root_index = utils::random_int(-8,8) + Size(protein_len/2);
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
    // Size multi_cool_cycles = 3.0;
    Size multi_cool_cycles = opts.mcmc_outer_cycles;
    Real best_score = 99999.0;
    Real best_score_inter_chain = 99999.0;
    
    // Size total_steps = protein_len * 10;
    Size total_steps = opts.mcmc_inner_cycles==-1 ? protein_len * 10 : opts.mcmc_inner_cycles;
    Real decrease_factor = 1.4;
    Real increase_factor = 1.1;
    // Size M = protein_len;
    Size M = total_steps / 10;

    // for C/D-symmetry
    Real C_symmetry_score = -9e9;
    Real D_symmetry_score = -9e9;

    // inter chain pair scores
    // get the symmetry
    std::string symmetry = pose.symmetry();
    Size num_chains = pose.num_chains();
    Size num_chains_C_symmetry = symmetry.at(0)=='C'?num_chains:num_chains/2;

    Real T = 1.5;

    for( Size ii_outer=1; ii_outer<=multi_cool_cycles; ++ii_outer) {

        // local parameters for MCMC
        // Real T = 3.0 / (ii_outer+1); // multicycle cool good ???? perturb the pose out of local minima?? I don't know.
        T *= 2.0;
        Real loop_weight(10.0);

        Real perturb_ang_mag(30);
        Real perturb_trans_mag(10);

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
            if (perturb_root) {
                cur_score_interchain = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1)) / protein_len;
            }
            Real delta = cur_score - best_score;
            bool pass = pass_metropolis(T, delta, runif(utils::rng()));

            if( pass ) {
            //    std::cout << "Accepted: " << cur_score << std::endl;
                best_score = cur_score;
                best_score_inter_chain = cur_score_interchain;
                pose.snapshot();
            } else {
                pose.rollback();
            }
            if(opts.dump_trajactory&&pass)pose.append_trajctory();
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

            if(opts.verbose_level!=-1 && ii_inner%opts.verbose_level==0) {


                Real cur_sidechain_neighbors = utils::sidechain_neighbors(pose);

                std::cout << std::setprecision(3) << "Protein Folding Stage: " << ii_outer
                                                  << " Outer: " << ii_outer 
                                                  << "; Total score: " << best_score 
                                                  << "; Sidechain Neighbor: " << cur_sidechain_neighbors;
	    	if(opts.num_repeats>1) {
			std::cout << std::setprecision(3) << "; Repeat Score: " << opts.repeat_weight * std::static_pointer_cast<scoring::RepeatScoreMethod>(sfxn.get_method("repeat"))->score(pose);
	    	}

		std::cout << std::endl;
            }

            if( perturb_root ) {

                bool perturb_root_pass = false;
                
                bool changed=false;
                if(num_chains_C_symmetry>1){
                    cur_C_symmetry_score = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,num_chains_C_symmetry)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,num_chains_C_symmetry) ) / protein_len;
                }
                    
                if(symmetry.at(0)=='D'){
                    cur_D_symmetry_score  = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,num_chains_C_symmetry+1,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,num_chains_C_symmetry+1,-1) ) / protein_len;    
                }
                    
                cur_score_interchain = cur_C_symmetry_score+cur_D_symmetry_score;
                if(cur_score_interchain>=opts.inter_chain_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor || cur_C_symmetry_score>=opts.C_symmetry_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor || cur_D_symmetry_score>=opts.D_symmetry_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor  ) {
                    // pose.root_rollback();
                    if(cur_score_interchain>=opts.inter_chain_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor) {
                        Size num_try(100);
                        perturb_root_pass = false;
                        while(--num_try){
                            pose.random_root(true, true, true, true);
                            cur_score_interchain =(std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1) ) / protein_len;
                            if(cur_score_interchain<opts.inter_chain_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor){
                                changed=true;
                                pose.snapshot(true);
                                best_score_inter_chain = cur_score_interchain;
                                best_score = sfxn.score(pose) / protein_len;
                                pose.snapshot(true);
                                perturb_root_pass = true;
                                break;
                            }
                        }
                        if( !perturb_root_pass ) {
                            pose.rollback(true);
                        }    
                    }
                    cur_C_symmetry_score = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,num_chains_C_symmetry)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,num_chains_C_symmetry) ) / protein_len;
                    if(cur_C_symmetry_score>=opts.C_symmetry_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor) {
                        Size num_try(100);
                        perturb_root_pass = false;
                        while(--num_try){
                            pose.random_root(false, true, true, false);
                            cur_C_symmetry_score  = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,num_chains_C_symmetry)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,num_chains_C_symmetry) ) / protein_len;
                            if(cur_C_symmetry_score<opts.C_symmetry_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor){
                                changed=true;
                                pose.snapshot(true);
                                best_score_inter_chain =(std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1) ) / protein_len;
                                best_score = sfxn.score(pose) / protein_len;
                                pose.snapshot(true);
                                perturb_root_pass = true;
                                break;
                            }
                        }
                        if( !perturb_root_pass ) {
                            pose.rollback(true);
                        }   
                    }


                    if(symmetry.at(0)=='D')
                        cur_D_symmetry_score  = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,num_chains_C_symmetry+1,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,num_chains_C_symmetry+1,-1) ) / protein_len;    
                    if(symmetry.at(0)=='D' && cur_D_symmetry_score>=opts.D_symmetry_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor) {
                        Size num_try(100);
                        perturb_root_pass = false;
                        while(--num_try){
                            pose.random_root(false, false, false, true);
                            cur_D_symmetry_score  = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,num_chains_C_symmetry+1,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,num_chains_C_symmetry+1,-1) ) / protein_len;    
                            if(cur_D_symmetry_score<opts.D_symmetry_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor){
                                changed=true;
                                pose.snapshot(true);
                                best_score_inter_chain =(std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1) ) / protein_len;
                                best_score = sfxn.score(pose) / protein_len;
                                pose.snapshot(true);
                                perturb_root_pass = true;
                                break;
                            }
                        }
                        if( !perturb_root_pass ) {
                            pose.rollback(true);
                        }   
                    }

                } else {
                    pose.perturb_root(perturb_ang_mag, perturb_trans_mag);
                    cur_score_interchain =(std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1) ) / protein_len;
                    delta = cur_score_interchain - best_score_inter_chain;
                    pass = pass_metropolis(T, delta, runif(utils::rng()));
                    
                    if( pass ) {
                        changed = true;
                        //std::cout << "Accepted: " << delta << std::endl;
                        best_score_inter_chain = cur_score_interchain;
                        best_score = sfxn.score(pose) / protein_len;
                        pose.snapshot(true);
                    } else {
                        pose.rollback(true);
                    }

                }
                if(opts.dump_trajactory&&changed)pose.append_trajctory();
            }
                

            // do fragment insertion again
            ires = frag_mover.pick_position(true);
            frag_mover.apply(pose,ires);


            cur_score = sfxn.score(pose) / protein_len;
            if (perturb_root) {
                cur_score_interchain = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1)) / protein_len;
            }
            delta = cur_score - best_score;
            pass = pass_metropolis(T, delta, runif(utils::rng()));

            if( pass ) {
            //    std::cout << "Accepted: " << cur_score << std::endl;
                best_score = cur_score;
                best_score_inter_chain = cur_score_interchain;
                pose.snapshot();
            } else {
                pose.rollback();
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
                //change_root_prob /= decrease_factor;
                perturb_ang_mag /= decrease_factor;
                perturb_ang_mag = perturb_ang_mag>5?perturb_ang_mag:5;
                perturb_trans_mag /= decrease_factor;
                perturb_trans_mag = perturb_trans_mag>0.5?perturb_trans_mag:0.5;
                loop_weight /= decrease_factor;
                loop_weight = loop_weight>=1?loop_weight:1;
                frag_mover.bias_loop_sampling_by_dssp(dssp, loop_weight, 1);
            }
            
            // after each move, check satisfication and return if good
            bool pass_requirements(true);
            if(best_score > opts.total_score_cutoff) pass_requirements = false;
            if(pass_requirements && perturb_root && best_score_inter_chain > opts.inter_chain_score_cutoff) pass_requirements=false;
	    if(pass_requirements && opts.num_repeats>1) {
		    if( opts.repeat_weight * std::static_pointer_cast<scoring::RepeatScoreMethod>(sfxn.get_method("repeat"))->score(pose) > opts.repeat_score_cutoff ) pass_requirements=false;
	    }
            if(pass_requirements && utils::sidechain_neighbors(pose) < opts.sidechain_neighbor_cutoff) pass_requirements=false;
            
            // num disulfide
            if(pass_requirements && opts.disulfide_requirement_num > 0) {
                Real disulfide_score = std::static_pointer_cast<scoring::DisulfideScoreMethod>(sfxn.get_method("disulfide"))->score(pose);
                if(disulfide_score > -opts.disulfide_requirement_num) pass_requirements=false; // -5 is the disulfide score multiplier
            }

            if(pass_requirements) {
                std::cout << "Solution found! Eearly stop." << std::endl;
                return true;
            }

        } // inner loop
    } // outer loop

    std::cout << "No solution found! Too bad." << std::endl;
    return false;
}



bool mcmc_motif_pair(scene::Pose & pose, scoring::ScoreFunction & sfxn, sampling::FragMover & frag_mover, metric::SidechainNeighbors & sidechain_neighbors, Options const & opts)
{
    const bool perturb_root = opts.symmetry != "C1";
    std::uniform_real_distribution<Real> runif(0,1);

    Size protein_len = pose.size();
    std::string dssp = opts.fake_ss;
    if(dssp == ""){
        dssp = pose.dssp();
    }
    Size frag_len = frag_mover.frag_length();

    assert(frag_mover.protein_length() == protein_len);

    Size root_index = utils::random_int(-Size(protein_len/4),Size(protein_len/4)) + Size(protein_len/2);
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
    // Size multi_cool_cycles = 3.0;
    Size multi_cool_cycles = opts.mcmc_outer_cycles;
    Real best_score = 99999.0;
    Real best_score_inter_chain = 99999.0;
    Real best_motif_pair_dist = 999999.0;
    
    // Size total_steps = protein_len * 10;
    Size total_steps = opts.mcmc_inner_cycles==-1 ? protein_len * 10 : opts.mcmc_inner_cycles;
    Real decrease_factor = 1.4;
    Real increase_factor = 1.1;
    // Size M = protein_len;
    Size M = total_steps / 10;

    const Size maximum_allowed_rejects = M * 25;
    Size accumulated_rejects = 0;

    // for C/D-symmetry
    Real C_symmetry_score = -9e9;
    Real D_symmetry_score = -9e9;

    // inter chain pair scores
    // get the symmetry
    std::string symmetry = pose.symmetry();
    Size num_chains = pose.num_chains();
    Size num_chains_C_symmetry = symmetry.at(0)=='C'?num_chains:num_chains/2;



    bool pass_this_stage = false;

    // Folding Stage

    sfxn.get_method("motif_pair")->set_weight(0.0);

    Real T = 0.25;
    for( Size ii_outer=1; ii_outer<=multi_cool_cycles; ++ii_outer) {

        // local parameters for MCMC
        // Real T = 3.0 / (ii_outer+1); // multicycle cool good ???? perturb the pose out of local minima?? I don't know.
        T *= 2.0;
        Real loop_weight(10.0);

        Real perturb_ang_mag(30);
        Real perturb_trans_mag(10);

        // initialize
        frag_mover.bias_loop_sampling_by_dssp(dssp, loop_weight, 1);

        for(Size ii_inner=0; ii_inner<total_steps; ++ii_inner) {
            //std::cout << "Step: " << ii << std::endl;
            Size ires = frag_mover.pick_position(true);

            frag_mover.apply(pose,ires);

            Real cur_score = sfxn.score(pose) / protein_len;
            Real cur_C_symmetry_score =0;
            Real cur_D_symmetry_score =0;
            Real cur_score_interchain = 0;
            if (perturb_root) {
                cur_score_interchain = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1)) / protein_len;
            }
            Real delta = cur_score - best_score;
            bool pass = pass_metropolis(T, delta, runif(utils::rng()));

            if( pass ) {
            //    std::cout << "Accepted: " << cur_score << std::endl;
                best_score = cur_score;
                best_score_inter_chain = cur_score_interchain;
                pose.snapshot();
                accumulated_rejects = 0;
            } else {
                pose.rollback();
                accumulated_rejects += 1;
            }


            // perturb root
            // only perturb root when there exists clashes
            // this might save some time
            if( perturb_root && best_score_inter_chain > 1.0) {

                pose.perturb_root(perturb_ang_mag, perturb_trans_mag);
                cur_score = sfxn.score(pose) / protein_len;
                cur_score_interchain =(std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1) ) / protein_len;
                delta = cur_score - best_score;
                pass = pass_metropolis(T, delta, runif(utils::rng()));
                
                if( pass ) {
                    best_score_inter_chain = cur_score_interchain;
                    best_score = cur_score;
                    pose.snapshot(true);
                    accumulated_rejects = 0;
                } else {
                    pose.rollback(true);
                    accumulated_rejects += 1;
                }

            }


            if(opts.early_stop && accumulated_rejects > maximum_allowed_rejects) {
                std::cout << "No acceptions after " << accumulated_rejects << " moves. Stop!" << std::endl;
                return false;
            }


            if(opts.verbose_level!=-1 && ii_inner%opts.verbose_level==0) {


                Real cur_sidechain_neighbors = sidechain_neighbors.compute(pose);
                Real cur_motif_pair_dist = opts.motif_pair_weight * std::static_pointer_cast<scoring::MotifPairRelativePosScoreMethod>(sfxn.get_method("motif_pair"))->score(pose);

                std::cout << std::setprecision(3) << "Folding Stage  <===>  Outer: " << std::left << std::setw(3) << ii_outer
                                                  << "Inner: " << std::left << std::setw(7) << ii_inner 
                                                  << "Total Sc: " << std::left << std::setw(10) << best_score 
                                                  << "InterChain Sc: " << std::left << std::setw(10) << best_score_inter_chain
                                                  << "SN: " << std::left << std::setw(10) << cur_sidechain_neighbors
                                                  << "Motif Dist: " << std::setprecision(5) << std::left << std::setw(10) << cur_motif_pair_dist;

                std::cout << std::endl;
            }

            
            // modify the weights

            if( ii_inner % M == 0 && ii_inner != 0 ) {
                T /= decrease_factor;
                perturb_ang_mag /= decrease_factor;
                perturb_ang_mag = perturb_ang_mag>5?perturb_ang_mag:5;
                perturb_trans_mag /= decrease_factor;
                perturb_trans_mag = perturb_trans_mag>0.5?perturb_trans_mag:0.5;
                loop_weight /= decrease_factor;
                loop_weight = loop_weight>=1?loop_weight:1;
                frag_mover.bias_loop_sampling_by_dssp(dssp, loop_weight, 1);
            }
            
            // after each move, check satisfication and return if good
            bool pass_requirements(true);
            // 0.8 * total_score ????
            if(best_score > opts.total_score_cutoff*0.75) pass_requirements = false;
            if(pass_requirements && sidechain_neighbors.compute(pose) < opts.sidechain_neighbor_cutoff*0.75) pass_requirements=false;

            if(pass_requirements) {
                std::cout << "Folding Stage Solution found! Eearly stop." << std::endl;
                pass_this_stage = true;
                break;
            }
        } // inner loop
        if(pass_this_stage) {
            break;
        }
    } // outer loop
    if(!pass_this_stage) {
        std::cout << "No Solution Found at Folding Stage!!" << std::endl;
        return false;
    }


    // Assembling Stage
    pass_this_stage = false;

    sfxn.get_method("motif_pair")->set_weight(opts.motif_pair_weight);
    best_score = sfxn.score(pose) / protein_len;
    if (perturb_root) {
        best_score_inter_chain = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1)) / protein_len;
    }
    best_motif_pair_dist = opts.motif_pair_weight * std::static_pointer_cast<scoring::MotifPairRelativePosScoreMethod>(sfxn.get_method("motif_pair"))->score(pose);

    T = 0.1;
    for( Size ii_outer=1; ii_outer<=2; ++ii_outer) {

        T *= 2.0;

        Real perturb_ang_mag(30);
        Real perturb_trans_mag(10);

        for(Size ii_inner=0; ii_inner<2000; ++ii_inner) {

            // perturb root
            // only perturb root when there exists clashes
            // this might save some time
            if( perturb_root ) {

                pose.perturb_root(perturb_ang_mag, perturb_trans_mag);
                Real cur_score = sfxn.score(pose) / protein_len;
                Real cur_score_interchain =(std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1) ) / protein_len;
                Real delta = cur_score - best_score;
                bool pass = pass_metropolis(T, delta, runif(utils::rng()));
                
                if( pass ) {
                    best_score_inter_chain = cur_score_interchain;
                    best_score = cur_score;
                    pose.snapshot(true);
                    accumulated_rejects = 0;
                } else {
                    pose.rollback(true);
                    accumulated_rejects += 1;
                }

            }


            if(opts.verbose_level!=-1 && ii_inner%opts.verbose_level==0) {


                Real cur_sidechain_neighbors = sidechain_neighbors.compute(pose);

                std::cout << std::setprecision(3) << "Assembling Stage  <===>  Outer: " << std::left << std::setw(3) << ii_outer
                                                  << "Inner: " << std::left << std::setw(7) << ii_inner 
                                                  << "Total Sc: " << std::left << std::setw(10) << best_score 
                                                  << "InterChain Sc: " << std::left << std::setw(10) << best_score_inter_chain
                                                  << "SN: " << std::left << std::setw(10) << cur_sidechain_neighbors
                                                  << "Motif Dist: " << std::setprecision(5) << std::left << std::setw(10) << best_motif_pair_dist;

                std::cout << std::endl;
            }

            
            // modify the weights

            if( ii_inner % 200 == 0 && ii_inner != 0 ) {
                T /= decrease_factor;
                perturb_ang_mag /= decrease_factor;
                perturb_ang_mag = perturb_ang_mag>5?perturb_ang_mag:5;
                perturb_trans_mag /= decrease_factor;
                perturb_trans_mag = perturb_trans_mag>0.5?perturb_trans_mag:0.5;
            }
        } // inner loop
    } // outer loop
    



    // co optimization stage
    T = 0.1;
    for( Size ii_outer=1; ii_outer<=multi_cool_cycles; ++ii_outer) {

        // local parameters for MCMC
        // Real T = 3.0 / (ii_outer+1); // multicycle cool good ???? perturb the pose out of local minima?? I don't know.
        T *= 2.0;
        Real loop_weight(10.0);

        Real perturb_ang_mag(30);
        Real perturb_trans_mag(10);

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
            if (perturb_root) {
                cur_score_interchain = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1)) / protein_len;
            }
            Real delta = cur_score - best_score;
            bool pass = pass_metropolis(T, delta, runif(utils::rng()));

            if( pass ) {
            //    std::cout << "Accepted: " << cur_score << std::endl;
                best_score = cur_score;
                best_score_inter_chain = cur_score_interchain;
                best_motif_pair_dist = opts.motif_pair_weight * std::static_pointer_cast<scoring::MotifPairRelativePosScoreMethod>(sfxn.get_method("motif_pair"))->score(pose);
                pose.snapshot();
                accumulated_rejects = 0;
            } else {
                pose.rollback();
                accumulated_rejects += 1;
            }


            // perturb root
            if( perturb_root ) {

                pose.perturb_root(perturb_ang_mag, perturb_trans_mag);
                cur_score = sfxn.score(pose) / protein_len;
                cur_score_interchain =(std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1) ) / protein_len;
                delta = cur_score - best_score;
                pass = pass_metropolis(T, delta, runif(utils::rng()));
                
                if( pass ) {
                    best_score_inter_chain = cur_score_interchain;
                    best_score = cur_score;
                    best_motif_pair_dist = opts.motif_pair_weight * std::static_pointer_cast<scoring::MotifPairRelativePosScoreMethod>(sfxn.get_method("motif_pair"))->score(pose);
                    pose.snapshot(true);
                    accumulated_rejects = 0;
                } else {
                    pose.rollback(true);
                    accumulated_rejects += 1;
                }

            }
                

            // do fragment insertion again
            ires = frag_mover.pick_position(true);
            frag_mover.apply(pose,ires);


            cur_score = sfxn.score(pose) / protein_len;
            if (perturb_root) {
                cur_score_interchain = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1)) / protein_len;
            }
            delta = cur_score - best_score;
            pass = pass_metropolis(T, delta, runif(utils::rng()));

            if( pass ) {
            //    std::cout << "Accepted: " << cur_score << std::endl;
                best_score = cur_score;
                best_score_inter_chain = cur_score_interchain;
                best_motif_pair_dist = opts.motif_pair_weight * std::static_pointer_cast<scoring::MotifPairRelativePosScoreMethod>(sfxn.get_method("motif_pair"))->score(pose);
                pose.snapshot();
                accumulated_rejects = 0;
            } else {
                pose.rollback();
                accumulated_rejects += 1;
            }


            if(opts.early_stop && accumulated_rejects > maximum_allowed_rejects) {
                std::cout << "No acceptions after " << accumulated_rejects << " moves. Stop!" << std::endl;
                return false;
            }


            if(opts.verbose_level!=-1 && ii_inner%opts.verbose_level==0) {


                Real cur_sidechain_neighbors = sidechain_neighbors.compute(pose);

                std::cout << std::setprecision(3) << "Co-Optimization Stage  <===>  Outer: " << std::left << std::setw(3) << ii_outer
                                                  << "Inner: " << std::left << std::setw(7) << ii_inner 
                                                  << "Total Sc: " << std::left << std::setw(10) << best_score 
                                                  << "InterChain Sc: " << std::left << std::setw(10) << best_score_inter_chain
                                                  << "SN: " << std::left << std::setw(10) << cur_sidechain_neighbors
                                                  << "Motif Dist: " << std::setprecision(5) << std::left << std::setw(10) << best_motif_pair_dist;

                std::cout << std::endl;
            }

            
            // modify the weights

            if( ii_inner % M == 0 && ii_inner != 0 ) {
                T /= decrease_factor;
                perturb_ang_mag /= decrease_factor;
                perturb_ang_mag = perturb_ang_mag>5?perturb_ang_mag:5;
                perturb_trans_mag /= decrease_factor;
                perturb_trans_mag = perturb_trans_mag>0.5?perturb_trans_mag:0.5;
                loop_weight /= decrease_factor;
                loop_weight = loop_weight>=1?loop_weight:1;
                frag_mover.bias_loop_sampling_by_dssp(dssp, loop_weight, 1);
            }
            
            // after each move, check satisfication and return if good
            bool pass_requirements(true);
            if(best_score > opts.total_score_cutoff) pass_requirements = false;
            if(pass_requirements && perturb_root && best_score_inter_chain > opts.inter_chain_score_cutoff) pass_requirements=false;
            if(pass_requirements && best_motif_pair_dist > opts.motif_pair_dist_cutoff) pass_requirements = false;
            if(pass_requirements && sidechain_neighbors.compute(pose) < opts.sidechain_neighbor_cutoff) pass_requirements=false;

            if(pass_requirements) {
                std::cout << "Solution found! Eearly stop." << std::endl;

                Real temp_total_sc = sfxn.score(pose) / pose.size();
                Real temp_sn = sidechain_neighbors.compute(pose);
                Real temp_motif_dist = opts.motif_pair_weight * std::static_pointer_cast<scoring::MotifPairRelativePosScoreMethod>(sfxn.get_method("motif_pair"))->score(pose);

                std::cout <<std::setprecision(3)  << "Success <===> Total_Sc: " << std::left << std::setw(10) << temp_total_sc
                                                  << "SCNeighbors: " << std::left << std::setw(10) << temp_sn
                                                  << "MotifDist: " << std::left << std::setw(10) << temp_motif_dist
                                                  <<  std::endl;

                return true;
            }

        } // inner loop
    } // outer loop


    {
        Real temp_total_sc = sfxn.score(pose) / pose.size();
        Real temp_sn = sidechain_neighbors.compute(pose);
        Real temp_motif_dist = opts.motif_pair_weight * std::static_pointer_cast<scoring::MotifPairRelativePosScoreMethod>(sfxn.get_method("motif_pair"))->score(pose);

        std::cout <<std::setprecision(3)  << "Failed <===> Total_Sc: " << std::left << std::setw(10) << temp_total_sc
                                          << "SCNeighbors: " << std::left << std::setw(10) << temp_sn
                                          << "MotifDist: " << std::left << std::setw(10) << temp_motif_dist
                                          <<  std::endl;
    }

    std::cout << "No solution found! Too bad." << std::endl;
    return false;
}





bool mcmc_refine(scene::Pose & pose, scoring::ScoreFunction & sfxn, sampling::FragMover & frag_mover, Options const & opts)
{
    const bool perturb_root = opts.symmetry != "C1";


    Size protein_len = pose.size();
    std::string dssp = opts.fake_ss;
    if(dssp == ""){
        dssp = pose.dssp();
    }
    Size frag_len = frag_mover.frag_length();

    assert(frag_mover.protein_length() == protein_len);

    
    pose.update_coordinates();
    pose.clear_energies();


    // global parameters
    Size multi_cool_cycles = 1.0;
    Real best_score = 99999.0;
    Real best_score_inter_chain = 99999.0;

    // initialize the scores
    best_score = sfxn.score(pose) / protein_len;
    if(perturb_root)
        best_score_inter_chain = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1)) / protein_len;

    pose.snapshot();
    
    Size total_steps = protein_len * 5;
    Real decrease_factor = 1.6;
    Real increase_factor = 1.1;
    Size M = protein_len;


    for( Size ii_outer=0; ii_outer<multi_cool_cycles; ++ii_outer) {

        Real loop_weight(10.0);

        Real perturb_ang_mag(15);
        Real perturb_trans_mag(5);

        // initialize
        frag_mover.bias_loop_sampling_by_dssp(dssp, loop_weight, 1);

        for(Size ii_inner=0; ii_inner<total_steps; ++ii_inner) {
            //std::cout << "Step: " << ii << std::endl;
            Size ires = frag_mover.pick_position(true);

            // do fragment insertion
            frag_mover.apply(pose,ires);


            Real cur_score = sfxn.score(pose) / protein_len;
            Real cur_score_interchain = 0;
            if (perturb_root) {
                cur_score_interchain = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1)) / protein_len;
            }

            if( cur_score < best_score ) {
            //    std::cout << "Accepted: " << cur_score << std::endl;
                best_score = cur_score;
                best_score_inter_chain = cur_score_interchain;
                pose.snapshot();
            } else {
                pose.rollback();
            }

            // perturb root
            if( perturb_root ) {

                pose.perturb_root(perturb_ang_mag, perturb_trans_mag);
                cur_score_interchain =(std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1) ) / protein_len;
           
                if( cur_score_interchain < best_score_inter_chain ) {
                    //std::cout << "Accepted: " << delta << std::endl;
                    best_score_inter_chain = cur_score_interchain;
                    best_score = sfxn.score(pose) / protein_len;
                    pose.snapshot(true);
                } else {
                    pose.rollback(true);
                }

            }


            if( ii_inner % M == 0 && ii_inner != 0 ) {
                perturb_ang_mag /= decrease_factor;
                perturb_trans_mag /= decrease_factor;
                loop_weight /= decrease_factor;
                frag_mover.bias_loop_sampling_by_dssp(dssp, loop_weight, 1);
            }

        } // inner loop
    } // outer loop

    // this function always return true, and the score always goes down.
    pose.clear_energies();
    return true;


}

bool mcmc_extension(scene::Pose & pose, scoring::ScoreFunction & sfxn, sampling::FragMover & frag_mover, metric::SidechainNeighbors & sidechain_neighbors, Options const & opts)
{
    std::uniform_real_distribution<Real> runif(0,1);

    Size protein_len = pose.size();
    std::string dssp = opts.fake_ss;
    if(dssp == ""){
        dssp = pose.dssp();
    }
    
    Size frag_len = frag_mover.frag_length();

    assert(frag_mover.protein_length() == protein_len);

    pose.update_coordinates();

    // take a snapshoot right before sampling
    //pose.root_snapshot();
    // score the pose first before taking the snapshot
    // sfxn.score(pose);
    sfxn.score(pose);
    
    pose.snapshot();
    //pose.snapshot(1);

    // global parameters
    // Size multi_cool_cycles = 3.0;
    Size multi_cool_cycles = opts.mcmc_outer_cycles;
    Real best_score = 99999.0;
    Real best_score_inter_chain = 99999.0;
    Real best_score_rif =99999.0;
    
    // Size total_steps = protein_len * 10;
    Size total_steps = opts.mcmc_inner_cycles==-1 ? protein_len * 10 : opts.mcmc_inner_cycles;
    Real decrease_factor = 1.4;
    Real increase_factor = 1.1;
    // Size M = protein_len;
    Size M = total_steps / 10;

    const Size maximum_allowed_rejects = 20 * M;
    Size accumulated_rejects = 0;

    Size num_chains = pose.num_chains();

    Real T = 0.5;

    for( Size ii_outer=1; ii_outer<=multi_cool_cycles; ++ii_outer) {

        T *= 2;
        Real loop_weight(10.0);

        // initialize
        frag_mover.bias_loop_sampling_by_dssp(dssp, loop_weight, 1);

        for(Size ii_inner=0; ii_inner<total_steps; ++ii_inner) {
            //std::cout << "Step: " << ii << std::endl;
            Size ires = frag_mover.pick_position(true);
            frag_mover.apply(pose,ires);


            Real cur_score = sfxn.score(pose) / protein_len;
            Real delta = cur_score - best_score;
            bool pass = pass_metropolis(T, delta, runif(utils::rng()));

            if( pass ) {
            //    std::cout << "Accepted: " << cur_score << std::endl;
                best_score = cur_score;
                if (num_chains>1) {
                    best_score_inter_chain = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1)) / protein_len;
                }
                if(opts.rif_table != "") {
                    best_score_rif = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score(pose);
                }
                pose.snapshot();
                accumulated_rejects = 0;
            } else {
                pose.rollback();
                accumulated_rejects += 1;
            }

            if(opts.early_stop && accumulated_rejects > maximum_allowed_rejects) {
                std::cout << "No acceptions after " << accumulated_rejects << " moves. Stop!" << std::endl;
                return false;
            }

            if(opts.verbose_level!=-1 && ii_inner%opts.verbose_level==0) {


                Real tmp_sidechain_neighbors = sidechain_neighbors.compute(pose);
                Real tmp_score = sfxn.score(pose) / protein_len;
                Real tmp_rpx_score = std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->score(pose) / protein_len;
                Real tmp_clash_score = std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->score(pose) / protein_len;
                std::cout << std::setprecision(3) << " Outer: " << ii_outer 
                                                  << " Inner: " << ii_inner 
                                                  << " Total Sc: " << std::left << std::setw(10) << tmp_score 
                                                  << " Rpx: " << std::left << std::setw(10) << tmp_rpx_score
                                                  << " Clash: " << std::left << std::setw(10) << tmp_clash_score 
                                                  << " SN: " << std::left << std::setw(10) << tmp_sidechain_neighbors;
                if(opts.rif_table != "") {
                    Real tmp_rif_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score(pose);
                    std::cout <<std::setprecision(3) << "; Rif: " << std::setw(10) << tmp_rif_score;
                }
                if(opts.context_pdb != "") {
                    Real tmp_context_clash_score = opts.context_clash_weight * std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash"))->score(pose);
                    std::cout <<std::setprecision(3) << "; Context Clash: " << std::setw(10) << tmp_context_clash_score;
                }

                std::cout << std::endl;
            }


            if( ii_inner % M == 0 && ii_inner != 0 ) {
                T /= decrease_factor;
                loop_weight /= decrease_factor;
                loop_weight = loop_weight>=1?loop_weight:1;
                frag_mover.bias_loop_sampling_by_dssp(dssp, loop_weight, 1);
            }
            
            // after each move, check satisfication and return if good
            bool pass_requirements(true);
            if(best_score > opts.total_score_cutoff) pass_requirements = false;
            if(pass_requirements && num_chains>1 && best_score_inter_chain > opts.inter_chain_score_cutoff) pass_requirements=false;
            if(pass_requirements && opts.rif_table != "" && best_score_rif > opts.rif_score_cutoff) pass_requirements=false;
            if(pass_requirements && sidechain_neighbors.compute(pose) < opts.sidechain_neighbor_cutoff) pass_requirements=false;
            
            if(pass_requirements) {
                Real tmp_score = sfxn.score(pose) / protein_len;
                Real tmp_sidechain_neighbors = sidechain_neighbors.compute(pose);
                std::cout << std::setprecision(4) 
                          << "Success <===> " 
                          << "Total Sc: " << std::left << std::setw(10) << tmp_score 
                          << "SCN: " << std::left << std::setw(10) << tmp_sidechain_neighbors;
                if(opts.rif_table != "") {
                    Real tmp_rif_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score(pose);
                    std::cout <<std::setprecision(4) << "Rif: " << std::setw(10) << tmp_rif_score;
                }
                if(opts.context_pdb != "") {
                    Real tmp_context_clash_score = opts.context_clash_weight * std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash"))->score(pose);
                    std::cout <<std::setprecision(4) << "Context Clash: " << std::setw(10) << tmp_context_clash_score;
                }
              
                std::cout <<  std::endl;
                std::cout << "Solution found! Eearly stop." << std::endl;
                return true;
            }

        } // inner loop
    } // outer loop

    Real tmp_score = sfxn.score(pose) / protein_len;
    Real tmp_sidechain_neighbors = sidechain_neighbors.compute(pose);
    std::cout <<std::setprecision(4) 
              << "Failed <===> " 
              << "Total Sc: " << std::left << std::setw(10) << tmp_score 
              << "SCN: " << std::left << std::setw(10) << tmp_sidechain_neighbors;
    if(opts.rif_table != "") {
        Real tmp_rif_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score(pose);
        std::cout <<std::setprecision(4) << "Rif: " << std::setw(10) << tmp_rif_score;
    }
    if(opts.context_pdb != "") {
        Real tmp_context_clash_score = opts.context_clash_weight * std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash"))->score(pose);
        std::cout <<std::setprecision(4) << "Context Clash: " << std::setw(10) << tmp_context_clash_score;
    }
              
    std::cout <<  std::endl;

    std::cout << "No solution found! Too bad." << std::endl;
    return false;
}

bool mcmc_jump(scene::Pose & pose, scoring::ScoreFunction & sfxn, Options const & opts, EigenXform const & extra_ligand_xform)
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
    Size total_steps = opts.mcmc_inner_cycles==-1 ? protein_len * 10 : opts.mcmc_inner_cycles;
    Size M = total_steps / 10;


    Real best_score = 9999999.0;  
    Real best_rpx_score = 99999999.0;
    Real best_ligand1_score = 99999999.0;
    Real best_ligand2_score = 99999999.0;

    Real decrease_factor = 1.4;
    Real T = 0.125;

    for( Size ii_outer=0; ii_outer<multi_cool_cycles; ++ii_outer) {

        // local parameters for MCMC
        T *= 4.0; // multicycle cool good ???? perturb the pose out of local minima?? I don't know.

        Real perturb_ang_mag(25);
        Real perturb_trans_mag(8);

        for(Size ii_inner=0; ii_inner<total_steps; ++ii_inner) {

            pose.perturb_root(perturb_ang_mag, perturb_trans_mag);

            Real cur_score(0.0), cur_rpx_score(0.0), cur_ligand1_score(0.0), cur_ligand2_score(0.0);
            cur_rpx_score = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1) ) / protein_len;
            cur_ligand1_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score(pose)+std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash"))->score(pose);
            cur_score = cur_rpx_score + cur_ligand1_score/protein_len;

            if(opts.extra_ligand != "") {
                cur_ligand2_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif2"))->score(pose)+std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash2"))->score(pose);
                cur_score += cur_ligand2_score/protein_len;
            }
            // Real z_cen     = std::fabs(pose_center(pose)[2]);

            //
            Real delta = cur_score - best_score;

            bool pass = (cur_rpx_score<0) & (cur_ligand1_score<0);
            // TODO: ligand2 score ??????????????????????????????????????????????????

            if( pass ) {
             pass = pass_metropolis(T, delta, runif(utils::rng()));
            }

            if( pass ) {
            //    std::cout << "Accepted: " << cur_score << std::endl;
                best_score = cur_score;
                best_rpx_score = cur_rpx_score;
                best_ligand1_score = cur_ligand1_score;
                if(opts.extra_ligand != "") {
                    best_ligand2_score = cur_ligand2_score;
                }
                pose.snapshot(true);
            } else {
                pose.rollback(true);
            }

            bool need_reinitialize(false);
            Real cur_score_interchain(0), ligand1_score(0), ligand2_score(0);

            // std::cout << "LOG: " << best_score << " " << T << std::endl;
            // std::cout << "Inter Chain Score: " << best_rpx_score << "          " << (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1) ) / protein_len << std::endl;
            // std::cout << "Ligand      Score: " << best_ligand1_score << "          " << opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score(pose)+std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash"))->score(pose) << std::endl;

            // check inter chain interaction
            if( !need_reinitialize && best_rpx_score>opts.inter_chain_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor) {
                need_reinitialize = true;
            }
            // check the first ligand
            if( !need_reinitialize && best_ligand1_score>opts.rif_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor) {
                need_reinitialize = true;
            }
            // check the second ligand if needed
            if( !need_reinitialize && opts.extra_ligand != "" && best_ligand2_score>opts.rif_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor) {
                need_reinitialize = true;
            }
    
            if(need_reinitialize) {

                // std::cout << "LOG: " << best_score << " " << best_rpx_score << " " << best_ligand1_score << std::endl;

                bool perturb_root_pass = false;
                Size num_try(100);
                while(--num_try){
                    pose.random_root(true, true, true, true);

                    // check again
                    need_reinitialize = false;
                    if( !need_reinitialize ) {
                        cur_score_interchain =(std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1) ) / protein_len;
                        if (cur_score_interchain>=opts.inter_chain_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor) {
                            need_reinitialize = true;
                        }
                    }

                    // check the first ligand
                    if( !need_reinitialize ) {
                        ligand1_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score(pose)+std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash"))->score(pose);
                        if (ligand1_score>=opts.rif_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor) {
                            need_reinitialize = true;
                        }
                    }

                    // check the second ligand if needed
                    if( !need_reinitialize && opts.extra_ligand != "") {
                        ligand2_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif2"))->score(pose)+std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash2"))->score(pose);
                        if (ligand2_score>=opts.rif_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor) {
                            need_reinitialize = true;
                        }
                    }
                    if(!need_reinitialize){
                        best_score = sfxn.score(pose) / protein_len;
                        best_rpx_score = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1) ) / protein_len;
                        best_ligand1_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score(pose)+std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash"))->score(pose);
                        if(opts.extra_ligand != "") {
                            best_ligand2_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif2"))->score(pose)+std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash2"))->score(pose);
                        }
                        pose.snapshot(true);
                        perturb_root_pass = true;
                        break;
                    }  
                }
                if( !perturb_root_pass ) {
                    pose.rollback(true);
                }
            } /* else {

                pose.perturb_root(perturb_ang_mag, perturb_trans_mag);

                cur_score = sfxn.score(pose) / protein_len;

                delta = cur_score - best_score;

                pass = pass_metropolis(T, delta, runif(utils::rng()));

                if( pass ) {
                //    std::cout << "Accepted: " << cur_score << std::endl;
                    best_score = cur_score;
                    best_rpx_score = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1) ) / protein_len;
                    best_ligand1_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score(pose)+std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash"))->score(pose);
                    if(opts.extra_ligand != "") {
                        best_ligand2_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif2"))->score(pose)+std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash2"))->score(pose);
                    }
                    pose.snapshot(true);
                } else {
                    pose.rollback(true);
                }

            } */


            if( ii_inner % M == 0 && ii_inner != 0 ) {
                T /= decrease_factor;
                perturb_ang_mag /= decrease_factor;
                perturb_ang_mag = perturb_ang_mag>5?perturb_ang_mag:5;
                perturb_trans_mag /= decrease_factor;
                perturb_trans_mag = perturb_trans_mag>0.5?perturb_trans_mag:0.5;
            }

            if(opts.verbose_level!=-1 && ii_inner%opts.verbose_level==0) {

                Real tmp_ligand_neighbors(0.0);
                if(opts.ligand_name == "AMA") {
                    tmp_ligand_neighbors = utils::AMA_neighbors(pose);
                }  else if(opts.ligand_name == "HQA") {
                    tmp_ligand_neighbors = utils::HQA_neighbors(pose);
                } else {
                    std::cout << "Holy shit!!! I am died." << std::endl;
                    exit(0);
                }

                 std::cout <<std::setprecision(3) 
                          << "Outer: " << std::left << std::setw(4) << ii_outer 
                          << "Inner: " << std::left << std::setw(7) << ii_inner 
                          << "Total Score: " << std::left << std::setw(10) << best_score 
                          << "Inter Chain: " << std::left << std::setw(10) << best_rpx_score
                          << "SideChainNeighbors: " << std::left << std::setw(10) << utils::sidechain_neighbors(pose, false, true) 
                          << "Lig1 Sc: " << std::left << std::setw(10) << best_ligand1_score
                          << "Lig1 Neibhbor: " << std::left << std::setw(10) << tmp_ligand_neighbors;
                if(opts.extra_ligand != "") {
                    Real tmp_ligand_neighbors2(0.0);
                    if(opts.ligand_name == "AMA") {
                        tmp_ligand_neighbors2 = utils::AMA_neighbors(pose, extra_ligand_xform);
                    }  else if(opts.ligand_name == "HQA") {
                        tmp_ligand_neighbors2 = utils::HQA_neighbors(pose, extra_ligand_xform);
                    } else {
                        std::cout << "Holy shit!!! I am died." << std::endl;
                        exit(0);
                    }
                    std::cout << "Lig2 Sc: " << std::left << std::setw(10) << best_ligand2_score
                              << "Lig2 Neibhbor: " << std::left << std::setw(10) << tmp_ligand_neighbors2;
                }
                std::cout << std::endl;
            }


            // after each move, check satisfication and return if good
            bool pass_requirements(true);
            if(pass_requirements && best_rpx_score>opts.inter_chain_score_cutoff) {
                pass_requirements=false;
            }
            if(pass_requirements && best_ligand1_score>opts.rif_score_cutoff) {
                pass_requirements=false;
            }
            if(pass_requirements && opts.extra_ligand != "" && best_ligand2_score>opts.rif_score_cutoff) {
                pass_requirements=false;
            }
            if(pass_requirements && opts.ligand_neighbors_cutoff != -1) {
                Real tmp_ligand_neighbors(0.0);
                if(opts.ligand_name == "AMA") {
                    tmp_ligand_neighbors = utils::AMA_neighbors(pose);
                }  else if(opts.ligand_name == "HQA") {
                    tmp_ligand_neighbors = utils::HQA_neighbors(pose);
                } else {
                    std::cout << "Holy shit!!! I am died." << std::endl;
                    exit(0);
                }
                if(tmp_ligand_neighbors<opts.ligand_neighbors_cutoff) {
                    pass_requirements=false;
                }
            } 
            if(pass_requirements && opts.ligand2_neighbors_cutoff != -1) {
                Real tmp_ligand_neighbors2(0.0);
                if(opts.ligand_name == "AMA") {
                    tmp_ligand_neighbors2 = utils::AMA_neighbors(pose, extra_ligand_xform);
                }  else if(opts.ligand_name == "HQA") {
                    tmp_ligand_neighbors2 = utils::HQA_neighbors(pose, extra_ligand_xform);
                } else {
                    std::cout << "Holy shit!!! I am died." << std::endl;
                    exit(0);
                }
                if(tmp_ligand_neighbors2<opts.ligand2_neighbors_cutoff) {
                    pass_requirements = false;
                }
            } 
            if(pass_requirements && utils::sidechain_neighbors(pose,false,true) < opts.sidechain_neighbor_cutoff) {
                pass_requirements=false;
            }

            if(pass_requirements) {
                std::cout << "Solution found! Eearly stop." << std::endl;
                return true;
            }
        } // inner loop
    } // outer loop
    std::cout << "No solution found! Too bad." << std::endl;
    return false;
}



bool mcmc_repeat_peptide(scene::Pose & pose, scoring::ScoreFunction & sfxn, sampling::FragMover & frag_mover, Options const & opts)
{
    const bool perturb_root = opts.symmetry != "C1";
    std::shared_ptr<scoring::RpxScoreTargetMethod> peptide_rpx_sfxn=std::static_pointer_cast<scoring::RpxScoreTargetMethod>(sfxn.get_method("rpx_peptide"));
    std::shared_ptr<scoring::RepeatScoreMethod> repeat_paramter_sxfn=std::static_pointer_cast<scoring::RepeatScoreMethod>(sfxn.get_method("repeat"));
    std::vector<std::string> exclude_score_from_backbone;
    exclude_score_from_backbone.push_back("repeat");
    exclude_score_from_backbone.push_back("rpx_peptide");
    std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->only_count_two_repeat(true);
    std::uniform_real_distribution<Real> runif(0,1);
    Size peptide_len = 12;
    scene::Pose peptide_pose(peptide_len, 1, std::string("C1"));
    {
        // reset peptide root stub to ideal postion
        peptide_pose.set_dssp(std::string(peptide_len,'H'));
        peptide_pose.set_sequence(std::string(peptide_len,'A'));
        peptide_pose.set_root_index(1);
        peptide_pose.reset_coords();
        Eigen::MatrixXf peptide_point(peptide_len,3);
        for (Size i=1;i<=peptide_len;i++)peptide_point.row(i-1)=peptide_pose.xyz(i,ATOM_CA);
        Eigen::MatrixXf peptide_point_centered = peptide_point.rowwise() - peptide_point.colwise().mean();
        // std::cout<<peptide_point_centered<<std::endl;
        Eigen::JacobiSVD<Eigen::MatrixXf> svd(peptide_point_centered, Eigen::ComputeThinV);
        // std::cout<<svd.matrixV()<<std::endl;
        Eigen::Vector3f axis = svd.matrixV().leftCols(1).col(0);
        // std::cout<<svd.matrixV().leftCols(1)<<std::endl;
        Eigen::Vector3f center = peptide_point.colwise().mean();
        EigenXform correct_root;
        correct_root.translation() = center;
        correct_root.matrix().col(2) = axis;
        correct_root.matrix().col(1) = (correct_root.matrix().col(2).cross(Vec(0,0,1))).normalized();
        correct_root.matrix().col(0) = (correct_root.matrix().col(1).cross(correct_root.matrix().col(2))).normalized(); 
        // std::cout<<"test"<<std::endl;
        peptide_pose.conformation().split_root_from_residue(correct_root);
        peptide_pose.update_coordinates();
        // peptide_pose.dump_pdb("peptide.pdb");
        pose.set_target_pose(peptide_pose);
    }
    // std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->only_count_two_repeat();
    Size protein_len = pose.size();
    std::string dssp = opts.fake_ss;
    if(dssp == ""){
        dssp = pose.dssp();
    }
    Size frag_len = frag_mover.frag_length();

    assert(frag_mover.protein_length() == protein_len);

    Size root_index = utils::random_int(-8,8) + Size(protein_len/2);
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
    Size multi_cool_cycles = 3.0;
    Real best_score = 99999.0;
    Real best_score_inter_chain = 99999.0;
    
    Size total_steps = protein_len * 10;
    Real decrease_factor = 1.4;
    Real increase_factor = 1.1;
    Size M = protein_len;

    // for C/D-symmetry
    Real C_symmetry_score = -9e9;
    Real D_symmetry_score = -9e9;

    // inter chain pair scores
    // get the symmetry
    std::string symmetry = pose.symmetry();
    Size num_chains = pose.num_chains();
    Size num_chains_C_symmetry = symmetry.at(0)=='C'?num_chains:num_chains/2;
    bool pass_requirements(true);
    for( Size ii_outer=1; ii_outer<=multi_cool_cycles; ++ii_outer) {

        // local parameters for MCMC
        Real T = 1.5 / ii_outer; // multicycle cool good ???? perturb the pose out of local minima?? I don't know.
        Real loop_weight(10.0);

        Real perturb_ang_mag(30);
        Real perturb_trans_mag(10);

        // initialize
        frag_mover.bias_loop_sampling_by_dssp(dssp, loop_weight, 1);

        for(Size ii_inner=0; ii_inner<total_steps; ++ii_inner) {
            //std::cout << "Step: " << ii << std::endl;
            Size ires = frag_mover.pick_position(true);

            frag_mover.apply(pose,ires);
            Vec axis,axis_center;
            Real tmp_real;
            utils::get_repeat_parameters_from_stubs(pose,tmp_real,tmp_real,tmp_real,axis,axis_center);
            EigenXform new_root = EigenXform(pose.get_target_pose()->conformation().root());
            new_root.translation() = axis_center;
            new_root.matrix().col(2) = axis;
            new_root.matrix().col(1) = (new_root.matrix().col(2).cross(Vec(0,0,1))).normalized();
            // when i use the reverse cross value than above, the peptide becomes D-Chirality, there's some stange things happen
            new_root.matrix().col(0) = (new_root.matrix().col(1).cross(new_root.matrix().col(2))).normalized();
            pose.get_target_pose()->conformation().set_root(new_root);
            pose.get_target_pose()->conformation().update_coordinates();
            Real peptide_rpx = peptide_rpx_sfxn->score(pose);
            Real repeat_paramerter_score = repeat_paramter_sxfn->score(pose);
            Real backbone_score = sfxn.score(pose,exclude_score_from_backbone) / protein_len;
            Real cur_score = backbone_score + peptide_rpx*2/(200.0*pose.num_repeats());
            // Real cur_score = backbone_score;
            Real cur_C_symmetry_score =0;
            Real cur_D_symmetry_score =0;
            Real cur_score_interchain = 0;
            Real delta = cur_score - best_score;
            bool pass = pass_metropolis(T, delta, runif(utils::rng()));

            if( pass ) {
            //    std::cout << "Accepted: " << cur_score << std::endl;
                best_score = cur_score;
                best_score_inter_chain = cur_score_interchain;
                pose.snapshot();
                pose.get_target_pose()->conformation().snapshot();
                if(opts.dump_trajactory)pose.append_trajctory();
            } else {
                pose.rollback();
                pose.get_target_pose()->conformation().rollback();
            }

            if( ii_inner % M == 0 && ii_inner != 0 ) {
                T /= decrease_factor;
                perturb_ang_mag /= decrease_factor;
                perturb_trans_mag /= decrease_factor;
                loop_weight /= decrease_factor;
                frag_mover.bias_loop_sampling_by_dssp(dssp, loop_weight, 1);
            }

            // after each move, check satisfication and return if good
        
            // if(backbone_score > opts.total_score_cutoff || peptide_rpx>opts.peptide_score_cutoff ) pass_requirements = false;
            // if(pass_requirements && utils::sidechain_neighbors(pose) < opts.sidechain_neighbor_cutoff) pass_requirements=false;
            // if(pass_requirements) {
            //     Vec axis,axis_center;
            //     Real tmp_real;
            //     // std::cout<<(pose.stub(1).inverse() * pose.stub(1+pose.size()/pose.num_repeats())).matrix()<<std::endl;
            //     // std::cout<<(pose.stub(1+pose.size()/pose.num_repeats())* pose.stub(1).inverse()).translation().norm()<<std::endl;
            //     // std::cout<<(pose.conformation().center_vec(1,pose.size()/pose.num_repeats())-pose.conformation().center_vec(1+pose.size()/pose.num_repeats(),pose.size()/pose.num_repeats()*2)).norm()<<std::endl;
            //     // std::cout<<pose.stub(1+pose.size()/pose.num_repeats()).matrix()<<pose.stub(1).matrix()<<std::endl;
            //     utils::get_repeat_parameter(pose,tmp_real,tmp_real,tmp_real,axis,axis_center);
            //     char buf[128];
            //     Size anum(1), rnum(1);
            //     std::ofstream out("/home/chentong/wefold/build/test_axis.pdb");
            //     Vec xyz = axis_center;
            //     // snprintf(buf,128,"%s%5i %4s %3s %c%4i    %8.3f%8.3f%8.3f%6.2f%6.2f %11s\n",
			// 	// 			"HETATM",
			// 	// 			anum++,
			// 	// 			"BURR",
			// 	// 			"BUR",
			// 	// 			'B',
			// 	// 			rnum++,
			// 	// 			xyz[0],xyz[1],xyz[2],
			// 	// 			1.0,
			// 	// 			1.0,
			// 	// 			"B"
			// 	// 		);

            //     // out << buf;
            //     // xyz = pose.conformation().center_vec();
            //     // snprintf(buf,128,"%s%5i %4s %3s %c%4i    %8.3f%8.3f%8.3f%6.2f%6.2f %11s\n",
			// 	// 			"HETATM",
			// 	// 			anum++,
			// 	// 			"BURR",
			// 	// 			"BUR",
			// 	// 			'B',
			// 	// 			rnum++,
			// 	// 			xyz[0],xyz[1],xyz[2],
			// 	// 			1.0,
			// 	// 			1.0,
			// 	// 			"B"
			// 	// 		);
            //     // out << buf;
            //     std::cout << "Solution found! Eearly stop." << std::endl;
            //     // std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->only_count_two_repeat(false);
            //     return true;
            // }

        } // inner loop
    } // outer loop
    std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->only_count_two_repeat(false);
    Real peptide_rpx = peptide_rpx_sfxn->score(pose);
    Real backbone_score = sfxn.score(pose,exclude_score_from_backbone) / protein_len;
    Real cur_score = backbone_score + peptide_rpx/300;
    Real clash =std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash_inter"))->score(pose);

    Real radius,omega,rise;
    Vec tmp_vec;
    utils::get_repeat_parameters_from_coords(pose,rise,radius,omega,tmp_vec,tmp_vec);
    if(omega<opts.repeat_min_omega || omega>opts.repeat_max_omega || radius<opts.repeat_min_radius || radius>opts.repeat_max_radius|| rise<opts.repeat_min_rise || rise>opts.repeat_max_radius) pass_requirements=false;
    if(clash>opts.clash_cutoff)pass_requirements=false;
    if(backbone_score > opts.total_score_cutoff || peptide_rpx>opts.peptide_score_cutoff ) pass_requirements = false;
    if(pass_requirements && utils::sidechain_neighbors(pose) < opts.sidechain_neighbor_cutoff) pass_requirements=false;
    if(pass_requirements)return true;
    std::cout<<"peptide rpx :"<<peptide_rpx<<" backbone rpx :"<<backbone_score<<" sc_neighbour: "<< utils::sidechain_neighbors(pose) <<" interclash: "<<clash<<std::endl;
    std::cout << "No solution found! Too bad." << std::endl;
    // std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->only_count_two_repeat(false);
    return false;
}

bool mcmc_fold_on_target(scene::Pose & pose, scoring::ScoreFunction & sfxn, sampling::FragMover & frag_mover, Options const & opts)
{
    std::uniform_real_distribution<Real> runif(0,1);
    Size protein_len = pose.size();
    Size frag_len = frag_mover.frag_length();

    assert(frag_mover.protein_length() == protein_len);

    std::shared_ptr<scoring::RpxScoreMethod> sfxn_rpx = std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"));
    std::shared_ptr<scoring::ClashScoreMethod> sfxn_clash_intra =  std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash_intra"));
    std::shared_ptr<scoring::VoxelClashScoreMethod> sfxn_target_clash =  std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("clash_inter"));
    std::shared_ptr<scoring::Rpx1SideScoreMethod> sfxn_rpx1side =  std::static_pointer_cast<scoring::Rpx1SideScoreMethod>(sfxn.get_method("rpx1side"));
    assert(frag_mover.protein_length() == protein_len);

    sfxn.score(pose);
    pose.snapshot();

    Size multi_cool_cycles = opts.mcmc_outer_cycles;
    Size total_steps = opts.mcmc_inner_cycles==-1 ? protein_len * 10 : opts.mcmc_inner_cycles;
    Size M = total_steps / 10;


    Real best_score = 99999999.9;
    Real best_intra_score = 99999999.9;
    Real best_inter_score = 99999999.9;

    Size try_counter = 0;


    Real decrease_factor = 1.4;
    
    // inter chain pair scores
    // get the symmetry

    // protein folding stage
    Real T = 1.5;
    for( Size ii_outer=1; ii_outer<=multi_cool_cycles; ++ii_outer) {
        try_counter = 0;
        // local parameters for MCMC
        T *= 2.0;
        Real loop_weight(10.0);

        bool is_protein_paritially_folded = false;

        // initialize
        frag_mover.bias_loop_sampling_by_dssp(pose.dssp(), loop_weight, 1);

        for(Size ii_inner=0; ii_inner<total_steps; ++ii_inner) {
            
            // do fragment insertion
            Size ires = frag_mover.pick_position(true);
            frag_mover.apply(pose,ires);

            Real cur_intra_score = ( sfxn_rpx->score(pose) + sfxn_clash_intra->score(pose) ) / protein_len;
            Real cur_inter_score = ( opts.context_clash_weight * sfxn_target_clash->score(pose) ) / protein_len;
            Real cur_score = cur_inter_score + cur_intra_score;

            Real delta = cur_score - best_score;
            bool pass = pass_metropolis(T, delta, runif(utils::rng()));

            if( pass ) {
                try_counter=0;
            //    std::cout << "Accepted: " << cur_score << std::endl;
                best_score = cur_score;
                best_inter_score = cur_inter_score;
                best_intra_score = cur_intra_score;
                pose.snapshot();
            } else {
                try_counter++;
                if (opts.max_try != -1 && try_counter>opts.max_try){
                    std::cout<<"Too many tries, return false"<<std::endl;
                    return false;
                }
                pose.rollback();
            }

            if(opts.verbose_level!=-1 && ii_inner%opts.verbose_level==0) {

                Real cur_rpx = sfxn_rpx->score(pose)/protein_len;
                Real cur_clash = sfxn_clash_intra->score(pose)/protein_len;
                Real cur_target_clash = opts.context_clash_weight * sfxn_target_clash->score(pose) / protein_len;
                Real cur_rpx1side = opts.rpx1side_weight * sfxn_rpx1side->score(pose) / protein_len;
                Real cur_sidechain_neighbors = utils::sidechain_neighbors(pose);

                std::cout << std::setprecision(3) << "Folding <==> "
                                                  << " Outer: " << std::left << std::setw(8) << ii_outer 
                                                  << " Inner: " << std::left << std::setw(10) << ii_inner 
                                                  << " Total Sc: " << std::left << std::setw(10) << best_score 
                                                  << " RPX Sc: " << std::left << std::setw(10) << cur_rpx 
                                                  << " Clash: " << std::left << std::setw(10) << cur_clash 
                                                  << " RPX 1-side: " << std::left << std::setw(10) << cur_rpx1side
                                                  << " Target Clash: " << std::left << std::setw(10) << cur_target_clash
                                                  << " Sidechain Neighbor: " << std::left << std::setw(10) << cur_sidechain_neighbors
                                                  <<  std::endl;
            }

            if( ii_inner % M == 0 && ii_inner != 0 ) {
                T /= decrease_factor;
                loop_weight /= decrease_factor;
                loop_weight = loop_weight>=1?loop_weight:1;
                frag_mover.bias_loop_sampling_by_dssp(pose.dssp(), loop_weight, 1);
            }

            bool pass_requirements = true;
            if(pass_requirements && best_intra_score > opts.intra_score_cutoff) pass_requirements = false;
            if(pass_requirements && best_inter_score > 10.0) pass_requirements = false;
            if(pass_requirements && utils::sidechain_neighbors(pose) < 0.8 * opts.sidechain_neighbor_cutoff) pass_requirements=false;

            if(pass_requirements) 
            {
                is_protein_paritially_folded = true;    
                break;
            }
        } // inner loop
        if( is_protein_paritially_folded ) {
            std::cout << "The protein is kindof folded, and jump to the next stage!" << std::endl;
            std::cout << std::setprecision(3) << "Current Scores: " 
                                                  << "; Total score: " << best_score 
                                                  << "; Intra Score: " << best_intra_score 
                                                  << "; Inter Score: " << best_inter_score
                                                  <<  std::endl;
            break;
        }
    } // outer loop

    best_intra_score = ( sfxn_rpx->score(pose) + sfxn_clash_intra->score(pose) ) / protein_len;
    best_inter_score = ( opts.context_clash_weight * sfxn_target_clash->score(pose) 
            + opts.rpx1side_weight * sfxn_rpx1side->score(pose) ) / protein_len;
    best_score = best_inter_score + best_intra_score;
    // protein binding stage
    T = 0.1;
    for( Size ii_outer=1; ii_outer<=multi_cool_cycles; ++ii_outer) {
        try_counter=0;
        // local parameters for MCMC
        T *= 2.0;
        Real loop_weight(10.0);

        // initialize
        frag_mover.bias_loop_sampling_by_dssp(pose.dssp(), loop_weight, 1);

        for(Size ii_inner=0; ii_inner<total_steps; ++ii_inner) {
            
            // do fragment insertion
            Size ires = frag_mover.pick_position(true);
            frag_mover.apply(pose,ires);

            Real cur_intra_score = ( sfxn_rpx->score(pose) + sfxn_clash_intra->score(pose) ) / protein_len;
            Real cur_inter_score = ( opts.context_clash_weight * sfxn_target_clash->score(pose) 
            + opts.rpx1side_weight * sfxn_rpx1side->score(pose) ) / protein_len;
            Real cur_score = cur_inter_score + cur_intra_score;

            Real delta = cur_score - best_score;
            bool pass = pass_metropolis(T, delta, runif(utils::rng()));

            if( pass ) {
                try_counter = 0;
            //    std::cout << "Accepted: " << cur_score << std::endl;
                best_score = cur_score;
                best_inter_score = cur_inter_score;
                best_intra_score = cur_intra_score;
                pose.snapshot();
            } else {
                try_counter++;
                if (opts.max_try != -1 && try_counter>opts.max_try){
                    std::cout<<"Too many tries, return false"<<std::endl;
                    return false;
                }
                pose.rollback();
            }

            if(opts.verbose_level!=-1 && ii_inner%opts.verbose_level==0) {

                Real cur_rpx = std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->score(pose)/protein_len;
                Real cur_clash = std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash_intra"))->score(pose)/protein_len;
                Real cur_target_clash = opts.context_clash_weight * sfxn_target_clash->score(pose) / protein_len;
                Real cur_rpx1side = opts.rpx1side_weight * sfxn_rpx1side->score(pose) / protein_len;
                Real cur_sidechain_neighbors = utils::sidechain_neighbors(pose);

                std::cout << std::setprecision(3) << "Binding <==> "
                                                  << " Outer: " << std::left << std::setw(8) << ii_outer 
                                                  << " Inner: " << std::left << std::setw(10) << ii_inner 
                                                  << " Total Sc: " << std::left << std::setw(10) << best_score 
                                                  << " RPX Sc: " << std::left << std::setw(10) << cur_rpx 
                                                  << " Clash: " << std::left << std::setw(10) << cur_clash 
                                                  << " RPX 1-side: " << std::left << std::setw(10) << cur_rpx1side
                                                  << " Target Clash: " << std::left << std::setw(10) << cur_target_clash
                                                  << " Sidechain Neighbor: " << std::left << std::setw(10) << cur_sidechain_neighbors
                                                  <<  std::endl;
            }

            if( ii_inner % M == 0 && ii_inner != 0 ) {
                T /= decrease_factor;
                loop_weight /= decrease_factor;
                loop_weight = loop_weight>=1?loop_weight:1;
                frag_mover.bias_loop_sampling_by_dssp(pose.dssp(), loop_weight, 1);
            }

            bool pass_requirements(true);
            if(best_score > opts.total_score_cutoff) pass_requirements = false;
            if(pass_requirements && best_intra_score > opts.intra_score_cutoff) pass_requirements = false;
            if(pass_requirements && best_inter_score > opts.inter_score_cutoff) pass_requirements = false;
            if(pass_requirements && utils::sidechain_neighbors(pose) < opts.sidechain_neighbor_cutoff) pass_requirements=false;

            if(pass_requirements)return true;
        } // inner loop
    } // outer loop


    std::cout << "No solution found! Too bad." << std::endl;
    // std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->only_count_two_repeat(false);
    return false;
}


bool mcmc_fold_on_target_chentong(scene::Pose & pose, scoring::ScoreFunction & sfxn, sampling::FragMover & frag_mover, Options const & opts)
{
    std::vector<std::string> exclude_score_from_backbone;
    std::uniform_real_distribution<Real> runif(0,1);
    Size protein_len = pose.size();
    Size frag_len = frag_mover.frag_length();
    exclude_score_from_backbone.push_back("clash_inter");
    exclude_score_from_backbone.push_back("rpx1side");
    std::shared_ptr<scoring::VoxelClashScoreMethod> sfxn_target_clash =  std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("clash_inter"));
    std::shared_ptr<scoring::Rpx1SideScoreMethod> sfxn_rpx1side =  std::static_pointer_cast<scoring::Rpx1SideScoreMethod>(sfxn.get_method("rpx1side"));
    assert(frag_mover.protein_length() == protein_len);

    sfxn.score(pose);
    pose.snapshot();

    Size multi_cool_cycles = 2.0;
    Real best_score_backbone = 99999.0;
    Real best_score_target = 99999.0;
    Size total_steps = protein_len * 10;
    Real decrease_factor = 1.4;
    Real decrease_factor_target = 1.2;
    Size M = protein_len;

    // inter chain pair scores
    // get the symmetry
    std::string symmetry = pose.symmetry();
    bool pass_requirements(true);
    for( Size ii_outer=1; ii_outer<=multi_cool_cycles; ++ii_outer) {

        // local parameters for MCMC
        Real T = 1 /(pow(pow(decrease_factor,4),ii_outer-1)+ii_outer); // multicycle cool good ???? perturb the pose out of local minima?? I don't know.
        Real T_target = 5/(pow(pow(decrease_factor_target,8),ii_outer-1)+ii_outer);
        Real loop_weight(10.0);

        Real perturb_ang_mag(30);
        Real perturb_trans_mag(10);

        // initialize
        frag_mover.bias_loop_sampling_by_dssp(std::string(pose.dssp()), loop_weight, 1);

        for(Size ii_inner=0; ii_inner<total_steps; ++ii_inner) {
            //std::cout << "Step: " << ii << std::endl;
            Size ires = frag_mover.pick_position(true);

            frag_mover.apply(pose,ires);
            Real backbone_score = sfxn.score(pose,exclude_score_from_backbone) / protein_len;
            Real target_clash_score = sfxn_target_clash->hresl_score(pose);
            Real target_score =0;
            if(target_clash_score>30)target_score=9999;
            else target_score =  target_clash_score*10.0/protein_len + (ii_outer==1?0:(sfxn_rpx1side->score(pose)/ protein_len));
            // Real cur_score = backbone_score;
            Real delta_backbone = backbone_score - best_score_backbone;
            Real delta_target = target_score - best_score_target;
            bool pass_backbone = pass_metropolis(T, delta_backbone, runif(utils::rng()));
            bool pass_target = pass_metropolis(T_target, delta_target, runif(utils::rng()));
            if( pass_backbone && pass_target ) {
            //    std::cout << "Accepted: " << cur_score << std::endl;
                best_score_backbone = backbone_score;
                best_score_target = target_score;
                pose.snapshot();
                pose.get_target_pose()->conformation().snapshot();
                if(opts.dump_trajactory)pose.append_trajctory();
            } else {
                pose.rollback();
                pose.get_target_pose()->conformation().rollback();
            }

            if( ii_inner % M == 0 && ii_inner != 0 ) {
                T_target/=decrease_factor_target;
                T /= decrease_factor;
                perturb_ang_mag /= decrease_factor;
                perturb_trans_mag /= decrease_factor;
                loop_weight /= decrease_factor;
                frag_mover.bias_loop_sampling_by_dssp(pose.dssp(), loop_weight, 1);
            }
        } // inner loop
    } // outer loop
    Real backbone_score = sfxn.score(pose,exclude_score_from_backbone) / protein_len;
    Real target_score = sfxn_target_clash->hresl_score(pose) + sfxn_rpx1side->score(pose);
    Real clash =sfxn_target_clash->hresl_score(pose);

    
    if(clash>opts.clash_cutoff)pass_requirements=false;
    if(backbone_score > opts.backbone_score_cutoff || target_score>opts.rpx1side_cutoff ) pass_requirements = false;
    if(utils::sidechain_neighbors(pose) < opts.sidechain_neighbor_cutoff) pass_requirements=false;
    if(pass_requirements)return true;
    std::cout<<"backbone score :"<<backbone_score<<" target score :"<<target_score<<" sc_neighbour: "<< utils::sidechain_neighbors(pose) <<" interclash: "<<clash<<std::endl;
    std::cout << "No solution found! Too bad." << std::endl;
    // std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->only_count_two_repeat(false);
    return false;
}



}

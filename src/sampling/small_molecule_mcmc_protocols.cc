#include "sampling/mcmc_protocols.hh"

#include "utils/utils.hh"
#include "basic/types.hh"
#include "scene/Pose.hh"
#include "sampling/FragMover.hh"
#include "scoring/ScoreFunction.hh"
#include "utils/random_util.hh"
#include "utils/ligand_util.hh"
#include "metric/LigandNeighbors.hh"
#include "metric/SidechainNeighbors.hh"

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


// small molecule binder sampler
bool small_molecule_binder_mcmc(scene::Pose & pose,
                                scoring::ScoreFunction & sfxn,
                                sampling::FragMover & frag_mover,
                                metric::SidechainNeighbors & sidechain_neighbors,
                                metric::LigandNeighbors & ligand_neighbors,
                                Options const & opts)
{
    std::uniform_real_distribution<Real> runif(0,1);

    Size protein_len = pose.size();
    std::string dssp = opts.fake_ss;
    if(dssp == ""){
        dssp = pose.dssp();
    }
    Size frag_len = frag_mover.frag_length();

    assert(frag_mover.protein_length() == protein_len);

    Size root_index = utils::random_int(-35,35) + Size(protein_len/2);
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
    Real best_rif_score = 99999.0;
    
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

    std::vector<Size> anchor_seqpos;

    Real T = 0.5;
    if(opts.temperature != -1.0) {
        T = opts.temperature/2.0;
    }

    Size mcmc_failed_attempt_counter(0);
    Size mcmc_failed_attempt_exit_counter(0);
    const Size max_failed_attempts_before_change_root(Size(total_steps/10));
    const Size max_failed_attempts_before_quit(total_steps*3);

    for( Size ii_outer=0; ii_outer<multi_cool_cycles; ++ii_outer) {

        // local parameters for MCMC
        // Real T = 3.0 / (ii_outer+1) / (ii_outer+1); // multicycle cool good ???? perturb the pose out of local minima?? I don't know.
        T *= 2.0;
        Real loop_weight(10.0);

        Real change_root_prob(0.3);

        Real perturb_ang_mag(30);
        Real perturb_trans_mag(5);

        // initialize
        frag_mover.bias_loop_sampling_by_dssp(dssp, loop_weight, 1);

        if( opts.ramp_context_clash_weight ) {
            sfxn.get_method("target_context_clash")->set_weight(opts.context_clash_weight * 0.1);

            // not the first iter, clear the energy graph and recompute the scores
            if( ii_outer != 0 ) {
                //
                pose.clear_energies();
                best_score = sfxn.score(pose) / protein_len;
                if( num_chains>1 ) {
                    best_score_inter_chain = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1)) / protein_len;
                }
                best_rif_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score_with_anchor_seqpos(pose, anchor_seqpos);;
                pose.snapshot();
            }
        }

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
            Real rif_score = 0.0;
            Real delta = cur_score - best_score;
            bool pass = pass_metropolis(T, delta, runif(utils::rng()));

            if( pass ) {
            //    std::cout << "Accepted: " << cur_score << std::endl;
                best_score = cur_score;
                if( num_chains>1 ) {
                    best_score_inter_chain = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1)) / protein_len;
                }
                best_rif_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score_with_anchor_seqpos(pose, anchor_seqpos);;
                pose.snapshot();

                mcmc_failed_attempt_counter = 0;
                mcmc_failed_attempt_exit_counter = 0;
            } else {
                pose.rollback();
                mcmc_failed_attempt_counter += 1;
                mcmc_failed_attempt_exit_counter += 1;
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

            // early stop
            // stop
            // stop
            // stop
            // stop
            if( opts.early_stop && mcmc_failed_attempt_exit_counter > max_failed_attempts_before_quit ) {
                std::cout << "Got trapped in local minima, STOP!!!" << std::endl;
                return false;
            }

            if ( mcmc_failed_attempt_counter <= max_failed_attempts_before_change_root ) {
                rif_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score(pose);
            } else {
                rif_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score_with_anchor_seqpos(pose, anchor_seqpos);

                if( anchor_seqpos.size() > 1 && runif(utils::rng()) < change_root_prob) {

                    Size num_spans(1);
                    std::vector<Size> span_labels(anchor_seqpos.size(),0);
                    for(Size ii_seqpos=1; ii_seqpos<anchor_seqpos.size(); ++ii_seqpos) {
                        if(anchor_seqpos[ii_seqpos] - anchor_seqpos[ii_seqpos-1]>5) {
                            num_spans += 1;
                            span_labels[ii_seqpos] = span_labels[ii_seqpos-1] + 1;
                        } else {
                            span_labels[ii_seqpos] = span_labels[ii_seqpos-1];
                        }
                    }

                    if( num_spans == 1 ) {
                        pose.set_root_index(anchor_seqpos[utils::random_int(0, anchor_seqpos.size()-1)]);
                    } else {
                        // exclude the span containing the current root index
                        Size exclude_span(-1);
                        for(Size ii_seqpos=0; ii_seqpos<anchor_seqpos.size(); ++ii_seqpos) {
                            if( anchor_seqpos[ii_seqpos] == pose.root_index() ) exclude_span = span_labels[ii_seqpos];
                        }


                        std::vector<Size> span_counts;
                        std::vector<std::pair<Size, Size > > span_range;
                        Size ii_s(0), cur_span_label(0);
                        for(Size ii_temp=1; ii_temp<span_labels.size(); ++ii_temp) {
                            //
                            if(span_labels[ii_temp] != cur_span_label) {
                                span_range.push_back(std::pair<Size, Size>(ii_s,ii_temp-1));
                                if(cur_span_label == exclude_span) {
                                    span_counts.push_back(-99);
                                } else {
                                    //
                                    span_counts.push_back(ii_temp-ii_s);
                                }
                                ii_s = ii_temp;
                                cur_span_label = span_labels[ii_temp];
                            }
                        }
                        if(cur_span_label == exclude_span) {
                            span_counts.push_back(-99);
                        } else {
                            //
                            span_counts.push_back(span_labels.size()-ii_s);
                        }
                        span_range.push_back(std::pair<Size, Size>(ii_s, span_labels.size()-1));

                        Size max_span_size(-1), max_span_index(-1);
                        for(Size ii_temp=0; ii_temp<span_counts.size(); ++ii_temp) {
                            if ( span_counts[ii_temp] > max_span_size ) {
                                max_span_index = ii_temp;
                                max_span_size = span_counts[ii_temp];
                            } else if ( span_counts[ii_temp] == max_span_size && runif(utils::rng()) > 0.5 ) {
                                //
                                max_span_index = ii_temp;
                            } else {
                                //
                            }
                        }

                        pose.set_root_index(anchor_seqpos[utils::random_int(span_range[max_span_index].first, span_range[max_span_index].second)]);

                        /*
                        if( true ) {
                            std::cout << "The current root index is " << pose.root_index() << std::endl;
                            std::cout << "The anchoring seq position: ";
                            for(Size ii_temp : anchor_seqpos ) std::cout << ii_temp << " ";
                            std::cout << std::endl;

                            for(auto s : span_range) {
                                std::cout << s.first << "    " << s.second << std::endl;
                            }

                            std::cout << "The anchoring span: ";
                            for(Size ii_temp : span_labels ) std::cout << ii_temp << " ";
                            std::cout << std::endl;

                            std::stringstream dump_name;
                            dump_name << "Iter_" << ii_inner << ".pdb";
                            pose.dump_pdb(dump_name.str(),false,opts.gzip);
                        }
                        */
                        
                    }
                    

                    pose.snapshot(true);
                    mcmc_failed_attempt_counter = 0;

                }
            }



            if( (num_chains > 1 && cur_score_interchain>=opts.inter_chain_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor) 
                || (num_chains_C_symmetry>1 && cur_C_symmetry_score>=opts.C_symmetry_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor) 
                || (symmetry.at(0) == 'D' && cur_D_symmetry_score>=opts.D_symmetry_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor)
                || rif_score >= opts.rif_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor ) {

                Size num_try(100);
                bool perturb_root_pass = false;
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
                    rif_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score_with_anchor_seqpos(pose, anchor_seqpos);


                    if( (num_chains == 1 || cur_score_interchain<opts.inter_chain_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor) 
                        && (num_chains_C_symmetry==1 || cur_C_symmetry_score<opts.C_symmetry_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor) 
                        && ( symmetry.at(0) != 'D' || cur_D_symmetry_score<opts.D_symmetry_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor)
                        && rif_score < opts.rif_score_cutoff*opts.sampling_stage_cutoff_tolerance_factor ) {

                        pose.set_root_index(anchor_seqpos[utils::random_int(0, anchor_seqpos.size()-1)]);
                        best_score = sfxn.score(pose) / protein_len;
                        best_score_inter_chain = cur_score_interchain;
                        best_rif_score = rif_score;
                        perturb_root_pass = true;
                        pose.snapshot(true);
                        mcmc_failed_attempt_counter = 0;
                        mcmc_failed_attempt_exit_counter = 0;
                        break;
                    }

                }
                if( !perturb_root_pass ) {
                    pose.rollback(true);
                }
            }


            // perturb root to better wrap the small molecule
            pose.perturb_root(perturb_ang_mag, perturb_trans_mag);
            cur_score = sfxn.score(pose) / protein_len;
            delta = cur_score - best_score;
            pass = pass_metropolis(T, delta, runif(utils::rng()));
            if( pass ) {
            //    std::cout << "Accepted: " << cur_score << std::endl;
                best_score = cur_score;
                if( num_chains>1 ) {
                    best_score_inter_chain = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1)) / protein_len;
                }
                best_rif_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score_with_anchor_seqpos(pose, anchor_seqpos);;
                pose.snapshot(true);
                mcmc_failed_attempt_counter = 0;
                mcmc_failed_attempt_exit_counter = 0;
            } else {
                pose.rollback(true);
                mcmc_failed_attempt_counter += 1;
                mcmc_failed_attempt_exit_counter += 1;
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
                best_rif_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score_with_anchor_seqpos(pose, anchor_seqpos);;
                pose.snapshot();
                mcmc_failed_attempt_counter = 0;
                mcmc_failed_attempt_exit_counter = 0;
            } else {
                pose.rollback();
                mcmc_failed_attempt_counter += 1;
                mcmc_failed_attempt_exit_counter += 1;
            }                

            if(opts.verbose_level!=-1 && ii_inner%opts.verbose_level==0) {

                Real tmp_score = sfxn.score(pose) / protein_len;

                Real tmp_rpx = std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->score(pose)/pose.size();
                Real tmp_clash = std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->score(pose)/pose.size();

                Real tmp_sidechain_neighbors = sidechain_neighbors.compute(pose);
                Real tmp_ligand_neighbors = ligand_neighbors.is_ligand_loaded()?ligand_neighbors.compute(pose):0.0;

                Real tmp_rif_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score(pose);
                Real tmp_target_clash = opts.context_clash_weight * std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash"))->score(pose);

                std::cout <<std::setprecision(4) 
                                                  << "Accepted <===> Outer: " << std::left << std::setw(10) << ii_outer 
                                                  << "Inner: " << std::left << std::setw(10) << ii_inner 
                                                  << "Total Sc: " << std::left << std::setw(10) << tmp_score 
                                                  << "RPX: " << std::left << std::setw(10) << tmp_rpx 
                                                  << "Clash: " << std::left << std::setw(10) << tmp_clash 
                                                  << "Rif: " << std::left << std::setw(10) << tmp_rif_score
                                                  << "Target Clash: " << std::left << std::setw(10) << tmp_target_clash
                                                  << "SCN: " << std::left << std::setw(10) << tmp_sidechain_neighbors
                                                  << "LN: " << std::left << std::setw(10) << tmp_ligand_neighbors
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
                loop_weight /= decrease_factor;
                loop_weight = loop_weight>=1?loop_weight:1;
                frag_mover.bias_loop_sampling_by_dssp(dssp, loop_weight, 1);

                if( opts.ramp_context_clash_weight ) {
                    sfxn.get_method("target_context_clash")->set_weight(opts.context_clash_weight * 0.1 * (ii_inner / M + 1) );

                    // not the first iter, clear the energy graph and recompute the scores
                    if( ii_outer != 0 ) {
                        //
                        pose.clear_energies();
                        best_score = sfxn.score(pose) / protein_len;
                        if( num_chains>1 ) {
                            best_score_inter_chain = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2,-1)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2,-1)) / protein_len;
                        }
                        best_rif_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score_with_anchor_seqpos(pose, anchor_seqpos);;
                        pose.snapshot();
                    }
                }
            }

            // Real sc_neighbor = utils::sidechain_neighbors(pose);
            // Real target_clash_no_weight = std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash"))->score(pose);

            // Real rpx_score = std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->score(pose);

            // after each move, check satisfication and return if good
            bool pass_requirements(true);
            // std::cout << std::setprecision(5);
            // if(best_score > opts.total_score_cutoff) { pass_requirements = false; if (ii_inner % 200 ==0) std::cout << ii_outer << " <---> " << ii_inner << "     TOTAL FAILED:   " << best_score << "   <===>   " << rpx_score << "   <===>   " << sc_neighbor << "   <===>   " << best_rif_score << "   <===>   " << target_clash_no_weight << "    <===>   " << anchor_seqpos.size() << std::endl;}
            // if(pass_requirements && num_chains>1 && best_score_inter_chain > opts.inter_chain_score_cutoff) { pass_requirements=false; }
            // if(pass_requirements && utils::sidechain_neighbors(pose) < opts.sidechain_neighbor_cutoff) { pass_requirements=false; if (ii_inner % 200 ==0) std::cout << ii_outer << " <--->" << ii_inner << "     SIDECHAIN FAILED:   " << best_score << "   <===>   " << sc_neighbor << "   <===>   " << best_rif_score << "   <===>   " << target_clash_no_weight << "    <===>   " << anchor_seqpos.size() << std::endl;}
            // if(pass_requirements && best_rif_score > opts.rif_score_cutoff) { pass_requirements=false; if (ii_inner % 200 ==0)  std::cout << ii_outer << " <---> " << ii_inner << "     RIF FAILED:   " << best_score << "   <===>   " << sc_neighbor << "   <===>   " << best_rif_score << "   <===>   " << target_clash_no_weight << "    <===>   " << anchor_seqpos.size() << std::endl; }

            if(best_score > opts.total_score_cutoff) { pass_requirements = false; }
            if(pass_requirements && num_chains>1 && best_score_inter_chain > opts.inter_chain_score_cutoff) { pass_requirements=false; }
            if(pass_requirements && best_rif_score > opts.rif_score_cutoff) { pass_requirements=false; }
            if(pass_requirements && sidechain_neighbors.compute(pose) < opts.sidechain_neighbor_cutoff) { pass_requirements=false; }
            if(pass_requirements && opts.ligand_neighbors_cutoff != -1 && ligand_neighbors.compute(pose) < opts.ligand_neighbors_cutoff) { pass_requirements=false; }
            if(pass_requirements) {
                Real tmp_score = sfxn.score(pose) / protein_len;

                Real tmp_rpx = std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->score(pose)/pose.size();
                Real tmp_clash = std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->score(pose)/pose.size();

                Real tmp_sidechain_neighbors = sidechain_neighbors.compute(pose);
                Real tmp_ligand_neighbors = ligand_neighbors.is_ligand_loaded()?ligand_neighbors.compute(pose):0.0;

                Real tmp_rif_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score(pose);

                Real tmp_target_clash = opts.context_clash_weight * std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash"))->score(pose);

                std::cout <<std::setprecision(4) 
                                                  << "Success <===> "
                                                  << "Total Sc: " << std::left << std::setw(10) << tmp_score 
                                                  << "RPX: " << std::left << std::setw(10) << tmp_rpx 
                                                  << "Clash: " << std::left << std::setw(10) << tmp_clash 
                                                  << "Rif: " << std::left << std::setw(10) << tmp_rif_score
                                                  << "Target Clash: " << std::left << std::setw(10) << tmp_target_clash
                                                  << "SCN: " << std::left << std::setw(10) << tmp_sidechain_neighbors
                                                  << "LN: " << std::left << std::setw(10) << tmp_ligand_neighbors
                                                  <<  std::endl;
                std::cout << "Solution found! Eearly stop." << std::endl;
                return true;
            }

        } // inner loop
    } // outer loop

    Real tmp_score = sfxn.score(pose) / protein_len;

    Real tmp_rpx = std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->score(pose)/pose.size();
    Real tmp_clash = std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->score(pose)/pose.size();

    Real tmp_sidechain_neighbors = sidechain_neighbors.compute(pose);
    Real tmp_ligand_neighbors = ligand_neighbors.is_ligand_loaded()?ligand_neighbors.compute(pose):0.0;

    Real tmp_rif_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score(pose);

    Real tmp_target_clash = opts.context_clash_weight * std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash"))->score(pose);

    std::cout <<std::setprecision(4) 
                                      << "Failed <===> " 
                                      << "Total Sc: " << std::left << std::setw(10) << tmp_score 
                                      << "RPX: " << std::left << std::setw(10) << tmp_rpx 
                                      << "Clash: " << std::left << std::setw(10) << tmp_clash 
                                      << "Rif: " << std::left << std::setw(10) << tmp_rif_score
                                      << "Target Clash: " << std::left << std::setw(10) << tmp_target_clash
                                      << "SCN: " << std::left << std::setw(10) << tmp_sidechain_neighbors
                                      << "LN: " << std::left << std::setw(10) << tmp_ligand_neighbors
                                      <<  std::endl;

    std::cout << "No solution found! Too bad." << std::endl;
    return false;
}


// small molecule binder sampler
bool ncaa_protein_mcmc(scene::Pose & pose, scoring::ScoreFunction & sfxn, sampling::FragMover & frag_mover, metric::SidechainNeighbors & sidechain_neighbors,
                                metric::LigandNeighbors & ligand_neighbors, Options const & opts)
{
    std::uniform_real_distribution<Real> runif(0,1);

    Size protein_len = pose.size();
    std::string dssp = opts.fake_ss;
    if(dssp == ""){
        dssp = pose.dssp();
    }
    std::string sequence = pose.sequence();

    // only for AzoF ??
    Size cutting_point1(-1);
    Size cutting_point2(-2);
    if(opts.motif_neighbors_cutoff != -1) {

        Size ncaa_index = sequence.find('X');
        if (ncaa_index < pose.size()/2) {
            
            Size iwalker(0);
            while (true) {
                iwalker++;
                if(dssp[iwalker] != 'H') {
                    cutting_point1 = iwalker;
                    break;
                }
            }
            while (true) {
                iwalker++;
                if(dssp[iwalker] != 'L') {
                    cutting_point2 = iwalker;
                    break;
                }
            }
        } else {
            Size iwalker(pose.size()-1);
            while (true) {
                iwalker--;
                if(dssp[iwalker] != 'H') {
                    cutting_point2 = iwalker;
                    break;
                }
            }
            while (true) {
                iwalker--;
                if(dssp[iwalker] != 'L') {
                    cutting_point1 = iwalker;
                    break;
                }
            }
        }
    }

    Size frag_len = frag_mover.frag_length();

    assert(frag_mover.protein_length() == protein_len);

    pose.update_coordinates();

    std::cout << dssp << std::endl;

    // take a snapshoot right before sampling
    //pose.root_snapshot();
    // score the pose first before taking the snapshot
    // sfxn.score(pose);
    sfxn.score(pose);
    
    pose.snapshot();
    //pose.snapshot(1);

    // global parameters
    Size multi_cool_cycles = opts.mcmc_outer_cycles;
    Real best_score = 99999.0;
    Real best_rif_score = 99999.0;
    
    Size total_steps = opts.mcmc_inner_cycles==-1 ? protein_len * 10 : opts.mcmc_inner_cycles;
    Real decrease_factor = 1.4;
    Real increase_factor = 1.1;
    Size M = total_steps / 10;

    const Size maximum_allowed_rejects = 25 * M;
    Size accumulated_rejects = 0;


    Real T = 0.5;

    for( Size ii_outer=0; ii_outer<multi_cool_cycles; ++ii_outer) {

        // local parameters for MCMC
        // Real T = 3.0 / (ii_outer+1) / (ii_outer+1); // multicycle cool good ???? perturb the pose out of local minima?? I don't know.
        T *= 2.0;
        Real loop_weight(10.0);

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
            Real delta = cur_score - best_score;
            bool pass = pass_metropolis(T, delta, runif(utils::rng()));

            if( pass ) {
            //    std::cout << "Accepted: " << cur_score << std::endl;
                best_score = cur_score;
                best_rif_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score(pose);
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

                    cur_score = sfxn.score(pose) / protein_len;

                    Real cur_rpx = std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->score(pose)/pose.size();
                    Real cur_clash = std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->score(pose)/pose.size();

                    Real cur_sidechain_neighbors = sidechain_neighbors.compute(pose);

                    Real cur_rif_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score(pose);

                    Real cur_target_clash = opts.context_clash_weight * std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash"))->score(pose);

                    Real cur_ligand_neighbors = ligand_neighbors.compute(pose);

                    std::cout <<std::setprecision(3) 
                                                      << "Outer: " << std::left << std::setw(10) << ii_outer 
                                                      << "Inner: " << std::left << std::setw(10) << ii_inner 
                                                      << "Total score: " << std::left << std::setw(10) << cur_score 
                                                      << "RPX score: " << std::left << std::setw(10) << cur_rpx 
                                                      << "Clash: " << std::left << std::setw(10) << cur_clash 
                                                      << "Rif: " << std::left << std::setw(10) << cur_rif_score
                                                      << "Target Clash: " << std::left << std::setw(10) << cur_target_clash
                                                      << "SN: " << std::left << std::setw(10) << cur_sidechain_neighbors
                                                      << "LN: " << std::left << std::setw(10) << cur_ligand_neighbors;
                    if(opts.motif_neighbors_cutoff != -1) {
                        Real cur_motif_neighbors = sidechain_neighbors.motif_neighbors(pose, cutting_point1, cutting_point2);
                        std::cout << std::setprecision(3) << "MN: " << std::left << std::setw(10) << cur_motif_neighbors;
                    }
                    std::cout <<  std::endl;

                    /*

                    if(cur_clash > 20) {
                        scoring::TwoBodyEnergyOP clash_energy = std::static_pointer_cast<scoring::TwoBodyEnergy>(pose.energy_manager().get_energy_tables(CLASH)[0]);

                        for(Size itmp=1; itmp<=pose.size(); ++itmp) {
                            for(Size jtmp=itmp+1; jtmp<=pose.size(); ++jtmp) {
                                Real pair_clash = clash_energy->get_pair_score(itmp, jtmp);
                                if(pair_clash > 0) {
                                    std::cout << itmp << "    ----    " << jtmp << "         " << pair_clash << std::endl;
                                }
                            }
                        }
                        pose.dump_pdb("test.pdb.gz", false, true, false, false);
                        exit(0);
                        // pose.dump_pdb(name.str(),false,opts.gzip);
                    }

                    
                    scoring::OneBodyEnergyOP voxel_clash_energy = std::static_pointer_cast<scoring::OneBodyEnergy>(pose.energy_manager().get_energy_tables(TARGET_CONTEXT_CLASH)[0]);

                    for(Size itmp=1; itmp<=pose.size(); ++itmp) {
                        Real voxel_clash_ires = voxel_clash_energy->get_score_per_residue(itmp);
                        if(voxel_clash_ires != 0) {
                            std::cout << pose.sequence() << std::endl;
                            std::cout << itmp << "    ----    " << voxel_clash_energy->get_score_per_residue(itmp) << std::endl;
                        }
                    }
                    */
                    
            }


            if( ii_inner % M == 0 && ii_inner != 0 ) {
                T /= decrease_factor;
                //change_root_prob /= decrease_factor;
                loop_weight /= decrease_factor;
                loop_weight = loop_weight>=1?loop_weight:1;
                frag_mover.bias_loop_sampling_by_dssp(dssp, loop_weight, 1);
            }

            // Real sc_neighbor = utils::sidechain_neighbors(pose);
            // Real target_clash_no_weight = std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash"))->score(pose);

            // Real rpx_score = std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->score(pose);

            // after each move, check satisfication and return if good
            bool pass_requirements(true);
            // std::cout << std::setprecision(5);
            // if(best_score > opts.total_score_cutoff) { pass_requirements = false; if (ii_inner % 200 ==0) std::cout << ii_outer << " <---> " << ii_inner << "     TOTAL FAILED:   " << best_score << "   <===>   " << rpx_score << "   <===>   " << sc_neighbor << "   <===>   " << best_rif_score << "   <===>   " << target_clash_no_weight << "    <===>   " << anchor_seqpos.size() << std::endl;}
            // if(pass_requirements && num_chains>1 && best_score_inter_chain > opts.inter_chain_score_cutoff) { pass_requirements=false; }
            // if(pass_requirements && utils::sidechain_neighbors(pose) < opts.sidechain_neighbor_cutoff) { pass_requirements=false; if (ii_inner % 200 ==0) std::cout << ii_outer << " <--->" << ii_inner << "     SIDECHAIN FAILED:   " << best_score << "   <===>   " << sc_neighbor << "   <===>   " << best_rif_score << "   <===>   " << target_clash_no_weight << "    <===>   " << anchor_seqpos.size() << std::endl;}
            // if(pass_requirements && best_rif_score > opts.rif_score_cutoff) { pass_requirements=false; if (ii_inner % 200 ==0)  std::cout << ii_outer << " <---> " << ii_inner << "     RIF FAILED:   " << best_score << "   <===>   " << sc_neighbor << "   <===>   " << best_rif_score << "   <===>   " << target_clash_no_weight << "    <===>   " << anchor_seqpos.size() << std::endl; }

            if(best_score > opts.total_score_cutoff) { pass_requirements = false; }
            if(pass_requirements && best_rif_score > opts.rif_score_cutoff) { pass_requirements=false; }
            if(pass_requirements && opts.ligand_neighbors_cutoff != -1 && ligand_neighbors.compute(pose) < opts.ligand_neighbors_cutoff) { pass_requirements=false; }
            if(pass_requirements && sidechain_neighbors.compute(pose) < opts.sidechain_neighbor_cutoff) { pass_requirements=false; }
            
            // only for AzoF
            if(pass_requirements && opts.motif_neighbors_cutoff != -1) {
                if( sidechain_neighbors.motif_neighbors(pose, cutting_point1, cutting_point2) < opts.motif_neighbors_cutoff ) {
                    pass_requirements = false;
                }
            }


            if(pass_requirements) {
                std::cout << "Solution found! Eearly stop." << std::endl;

                Real tmp_score = sfxn.score(pose) / protein_len;

                Real tmp_rpx = std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->score(pose)/pose.size();
                Real tmp_clash = std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->score(pose)/pose.size();

                Real tmp_sidechain_neighbors = sidechain_neighbors.compute(pose);

                Real tmp_rif_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score(pose);

                Real tmp_target_clash = opts.context_clash_weight * std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash"))->score(pose);

                Real tmp_ligand_neighbors = ligand_neighbors.compute(pose);

                std::cout <<std::setprecision(3)  << "Success: Total_Sc: " << std::left << std::setw(10) << tmp_score 
                                                  << "RPX: " << std::left << std::setw(10) << tmp_rpx 
                                                  << "Clash: " << std::left << std::setw(10) << tmp_clash 
                                                  << "Rif: " << std::left << std::setw(10) << tmp_rif_score
                                                  << "Target_Clash: " << std::left << std::setw(10) << tmp_target_clash
                                                  << "SN: " << std::left << std::setw(10) << tmp_sidechain_neighbors
                                                  << "LN: " << std::left << std::setw(10) << tmp_ligand_neighbors;
                if(opts.motif_neighbors_cutoff != -1) {
                    Real tmp_motif_neighbors = sidechain_neighbors.motif_neighbors(pose, cutting_point1, cutting_point2);
                    std::cout << std::setprecision(3) << "MN: " << std::left << std::setw(10) << tmp_motif_neighbors;
                }
                std::cout <<  std::endl;

                return true;
            }

        } // inner loop
    } // outer loop

    std::cout << "No solution found! Too bad." << std::endl;

    Real tmp_score = sfxn.score(pose) / protein_len;

    Real tmp_rpx = std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->score(pose)/pose.size();
    Real tmp_clash = std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->score(pose)/pose.size();

    Real tmp_sidechain_neighbors = sidechain_neighbors.compute(pose);

    Real tmp_rif_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score(pose);

    Real tmp_target_clash = opts.context_clash_weight * std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash"))->score(pose);

    Real tmp_ligand_neighbors = ligand_neighbors.compute(pose);

    std::cout <<std::setprecision(3)  << "Failed: Total_Sc: " << std::left << std::setw(10) << tmp_score 
                                      << "RPX: " << std::left << std::setw(10) << tmp_rpx 
                                      << "Clash: " << std::left << std::setw(10) << tmp_clash 
                                      << "Rif: " << std::left << std::setw(10) << tmp_rif_score
                                      << "Target_Clash: " << std::left << std::setw(10) << tmp_target_clash
                                      << "SN: " << std::left << std::setw(10) << tmp_sidechain_neighbors
                                      << "LN: " << std::left << std::setw(10) << tmp_ligand_neighbors;
    if(opts.motif_neighbors_cutoff != -1) {
        Real tmp_motif_neighbors = sidechain_neighbors.motif_neighbors(pose, cutting_point1, cutting_point2);
        std::cout << std::setprecision(3) << "MN: " << std::left << std::setw(10) << tmp_motif_neighbors;
    }
    std::cout <<  std::endl;

    return false;
}


// small molecule binder sampler
bool GFP_protein_mcmc(scene::Pose & pose, scoring::ScoreFunction & sfxn, sampling::FragMover & frag_mover, metric::SidechainNeighbors & sidechain_neighbors,
                                metric::LigandNeighbors & ligand_neighbors, Options const & opts)
{
    std::uniform_real_distribution<Real> runif(0,1);

    Size protein_len = pose.size();
    std::string dssp = pose.dssp();
    std::string sequence = pose.sequence();


    Size frag_len = frag_mover.frag_length();

    assert(frag_mover.protein_length() == protein_len);

    pose.update_coordinates();

    std::cout << dssp << std::endl;

    // take a snapshoot right before sampling
    //pose.root_snapshot();
    // score the pose first before taking the snapshot
    // sfxn.score(pose);
    sfxn.score(pose);
    
    pose.snapshot();
    //pose.snapshot(1);

    // global parameters
    Size multi_cool_cycles = opts.mcmc_outer_cycles;
    Real best_score = 99999.0;
    Real best_rif_score = 99999.0;
    Real best_privileged_motif_score = 99999.0;
    
    Size total_steps = opts.mcmc_inner_cycles==-1 ? protein_len * 10 : opts.mcmc_inner_cycles;
    Real decrease_factor = 1.4;
    Real increase_factor = 1.1;
    Size M = total_steps / 10;

    const Size maximum_allowed_rejects = 30 * M;
    Size accumulated_rejects = 0;


    Real T = 0.5;

    for( Size ii_outer=0; ii_outer<multi_cool_cycles; ++ii_outer) {

        // local parameters for MCMC
        // Real T = 3.0 / (ii_outer+1) / (ii_outer+1); // multicycle cool good ???? perturb the pose out of local minima?? I don't know.
        T *= 2.0;
        Real loop_weight(10.0);

        // initialize
        frag_mover.bias_loop_sampling_by_dssp(dssp, loop_weight, 1);

        for(Size ii_inner=0; ii_inner<total_steps; ++ii_inner) {
            //std::cout << "Step: " << ii << std::endl;
            Size ires;
            bool is_good=false;
            Size num_try(1000);
            while(--num_try) {
                ires = frag_mover.pick_position(true);
                if(dssp.substr(ires-1, 7) != "LLLLLLL"){
                    is_good=true;
                    break;
                }
            }
            if(!is_good) continue;

            frag_mover.apply(pose,ires);


            Real cur_score = sfxn.score(pose) / protein_len;
            Real delta = cur_score - best_score;
            bool pass = pass_metropolis(T, delta, runif(utils::rng()));

            if( pass ) {
            //    std::cout << "Accepted: " << cur_score << std::endl;
                best_score = cur_score;
                best_rif_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score(pose);
                best_privileged_motif_score = opts.privileged_motif_weight * std::static_pointer_cast<scoring::PrivilegedMotifScoreMethod>(sfxn.get_method("privileged_motif"))->score(pose);
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

                    cur_score = sfxn.score(pose) / protein_len;

                    Real cur_rpx = std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->score(pose)/pose.size();
                    Real cur_clash = std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->score(pose)/pose.size();

                    Real cur_sidechain_neighbors = sidechain_neighbors.compute(pose);

                    Real cur_rif_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score(pose);

                    Real cur_privileged_motif_score = opts.privileged_motif_weight * std::static_pointer_cast<scoring::PrivilegedMotifScoreMethod>(sfxn.get_method("privileged_motif"))->score(pose);

                    Real cur_target_clash = opts.context_clash_weight * std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash"))->score(pose);

                    Real cur_ligand_neighbors = ligand_neighbors.compute(pose);

                    std::cout <<std::setprecision(3) 
                                                      << "Outer: " << std::left << std::setw(10) << ii_outer 
                                                      << "Inner: " << std::left << std::setw(10) << ii_inner 
                                                      << "TotalSc: " << std::left << std::setw(10) << cur_score 
                                                      << "RPX: " << std::left << std::setw(10) << cur_rpx 
                                                      << "Clash: " << std::left << std::setw(10) << cur_clash 
                                                      << "Rif: " << std::left << std::setw(10) << cur_rif_score
                                                      << "PvMotif: " << std::left << std::setw(10) << cur_privileged_motif_score
                                                      << "TargetClash: " << std::left << std::setw(10) << cur_target_clash
                                                      << "SN: " << std::left << std::setw(10) << cur_sidechain_neighbors
                                                      << "LN: " << std::left << std::setw(10) << cur_ligand_neighbors;
                    std::cout <<  std::endl;

                
                    
            }


            if( ii_inner % M == 0 && ii_inner != 0 ) {
                T /= decrease_factor;
                //change_root_prob /= decrease_factor;
                loop_weight /= decrease_factor;
                loop_weight = loop_weight>=1?loop_weight:1;
                frag_mover.bias_loop_sampling_by_dssp(dssp, loop_weight, 1);
            }

            // Real sc_neighbor = utils::sidechain_neighbors(pose);
            // Real target_clash_no_weight = std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash"))->score(pose);

            // Real rpx_score = std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->score(pose);

            // after each move, check satisfication and return if good
            bool pass_requirements(true);
            // std::cout << std::setprecision(5);
            // if(best_score > opts.total_score_cutoff) { pass_requirements = false; if (ii_inner % 200 ==0) std::cout << ii_outer << " <---> " << ii_inner << "     TOTAL FAILED:   " << best_score << "   <===>   " << rpx_score << "   <===>   " << sc_neighbor << "   <===>   " << best_rif_score << "   <===>   " << target_clash_no_weight << "    <===>   " << anchor_seqpos.size() << std::endl;}
            // if(pass_requirements && num_chains>1 && best_score_inter_chain > opts.inter_chain_score_cutoff) { pass_requirements=false; }
            // if(pass_requirements && utils::sidechain_neighbors(pose) < opts.sidechain_neighbor_cutoff) { pass_requirements=false; if (ii_inner % 200 ==0) std::cout << ii_outer << " <--->" << ii_inner << "     SIDECHAIN FAILED:   " << best_score << "   <===>   " << sc_neighbor << "   <===>   " << best_rif_score << "   <===>   " << target_clash_no_weight << "    <===>   " << anchor_seqpos.size() << std::endl;}
            // if(pass_requirements && best_rif_score > opts.rif_score_cutoff) { pass_requirements=false; if (ii_inner % 200 ==0)  std::cout << ii_outer << " <---> " << ii_inner << "     RIF FAILED:   " << best_score << "   <===>   " << sc_neighbor << "   <===>   " << best_rif_score << "   <===>   " << target_clash_no_weight << "    <===>   " << anchor_seqpos.size() << std::endl; }

            if(best_score > opts.total_score_cutoff) { pass_requirements = false; }
            if(pass_requirements && best_rif_score > opts.rif_score_cutoff) { pass_requirements=false; }
            if(pass_requirements && best_privileged_motif_score > opts.privileged_motif_score_cutoff) {pass_requirements=false; }
            if(pass_requirements && sidechain_neighbors.compute(pose) < opts.sidechain_neighbor_cutoff) { pass_requirements=false; }
            if(pass_requirements && opts.ligand_neighbors_cutoff != -1 && ligand_neighbors.compute(pose) < opts.ligand_neighbors_cutoff) { pass_requirements=false; }


            if(pass_requirements) {
                std::cout << "Solution found! Eearly stop." << std::endl;

                Real tmp_score = sfxn.score(pose) / protein_len;

                Real tmp_rpx = std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->score(pose)/pose.size();
                Real tmp_clash = std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->score(pose)/pose.size();

                Real tmp_sidechain_neighbors = sidechain_neighbors.compute(pose);

                Real tmp_rif_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score(pose);

                Real tmp_privileged_motif_score = opts.privileged_motif_weight * std::static_pointer_cast<scoring::PrivilegedMotifScoreMethod>(sfxn.get_method("privileged_motif"))->score(pose);

                Real tmp_target_clash = opts.context_clash_weight * std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash"))->score(pose);

                Real tmp_ligand_neighbors = ligand_neighbors.compute(pose);

                std::cout <<std::setprecision(3)  << "Success: Total_Sc: " << std::left << std::setw(10) << tmp_score 
                                                  << "RPX: " << std::left << std::setw(10) << tmp_rpx 
                                                  << "Clash: " << std::left << std::setw(10) << tmp_clash 
                                                  << "Rif: " << std::left << std::setw(10) << tmp_rif_score
                                                  << "PvMotif: " << std::left << std::setw(10) << tmp_privileged_motif_score
                                                  << "Target_Clash: " << std::left << std::setw(10) << tmp_target_clash
                                                  << "SN: " << std::left << std::setw(10) << tmp_sidechain_neighbors
                                                  << "LN: " << std::left << std::setw(10) << tmp_ligand_neighbors;
                std::cout <<  std::endl;

                return true;
            }

        } // inner loop
    } // outer loop

    std::cout << "No solution found! Too bad." << std::endl;

    Real tmp_score = sfxn.score(pose) / protein_len;

    Real tmp_rpx = std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->score(pose)/pose.size();
    Real tmp_clash = std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->score(pose)/pose.size();

    Real tmp_sidechain_neighbors = sidechain_neighbors.compute(pose);

    Real tmp_rif_score = opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score(pose);

    Real tmp_privileged_motif_score = opts.privileged_motif_weight * std::static_pointer_cast<scoring::PrivilegedMotifScoreMethod>(sfxn.get_method("privileged_motif"))->score(pose);

    Real tmp_target_clash = opts.context_clash_weight * std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash"))->score(pose);

    Real tmp_ligand_neighbors = ligand_neighbors.compute(pose);

    std::cout <<std::setprecision(3)  << "Failed: Total_Sc: " << std::left << std::setw(10) << tmp_score 
                                      << "RPX: " << std::left << std::setw(10) << tmp_rpx 
                                      << "Clash: " << std::left << std::setw(10) << tmp_clash 
                                      << "Rif: " << std::left << std::setw(10) << tmp_rif_score
                                      << "PvMotif: " << std::left << std::setw(10) << tmp_privileged_motif_score
                                      << "Target_Clash: " << std::left << std::setw(10) << tmp_target_clash
                                      << "SN: " << std::left << std::setw(10) << tmp_sidechain_neighbors
                                      << "LN: " << std::left << std::setw(10) << tmp_ligand_neighbors;
    std::cout <<  std::endl;

    return false;
}


// small molecule binder sampler
bool ncaa_protein_symm_mcmc(scene::Pose & pose, scoring::ScoreFunction & sfxn, sampling::FragMover & frag_mover, metric::SidechainNeighbors & sidechain_neighbors,
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
    Real best_score = 99999.0;
    Real best_score_inter_chain = 99999.0;
    
    Size total_steps = opts.mcmc_inner_cycles==-1 ? protein_len * 10 : opts.mcmc_inner_cycles;
    Real decrease_factor = 1.4;
    Real increase_factor = 1.1;
    Size M = total_steps / 10;

    const Size maximum_allowed_rejects = M * 10;
    Size accumulated_rejects = 0;


    Real T = 0.5;

    for( Size ii_outer=0; ii_outer<multi_cool_cycles; ++ii_outer) {

        // local parameters for MCMC
        // Real T = 3.0 / (ii_outer+1) / (ii_outer+1); // multicycle cool good ???? perturb the pose out of local minima?? I don't know.
        T *= 2.0;
        Real loop_weight(10.0);

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
            Real delta = cur_score - best_score;
            bool pass = pass_metropolis(T, delta, runif(utils::rng()));

            if( pass ) {
            //    std::cout << "Accepted: " << cur_score << std::endl;
                best_score = cur_score;
                best_score_inter_chain = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2) ) / pose.size();
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

                    cur_score = sfxn.score(pose) / protein_len;

                    Real cur_rpx = std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->score(pose)/pose.size();
                    Real cur_clash = std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->score(pose)/pose.size();

                    Real cur_sidechain_neighbors = sidechain_neighbors.compute(pose);

                    Real cur_target_clash = opts.context_clash_weight * std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash"))->score(pose);

                    Real cur_score_interchain = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2) ) / pose.size();

                    Real cur_ligand_neighbors = ligand_neighbors.compute(pose);

                    std::cout <<std::setprecision(3) 
                                                      << "Accepted <===> Outer: " << std::left << std::setw(10) << ii_outer 
                                                      << "Inner: " << std::left << std::setw(10) << ii_inner 
                                                      << "Total_Sc: " << std::left << std::setw(10) << cur_score 
                                                      << "RPX: " << std::left << std::setw(10) << cur_rpx 
                                                      << "Clash: " << std::left << std::setw(10) << cur_clash 
                                                      << "Score_InterChain: " << std::left << std::setw(10) << cur_score_interchain
                                                      << "Target_Clash: " << std::left << std::setw(10) << cur_target_clash
                                                      << "SCN: " << std::left << std::setw(10) << cur_sidechain_neighbors
                                                      << "LN: " << std::left << std::setw(10) << cur_ligand_neighbors
                                                      <<  std::endl;

                    /*

                    if(cur_clash > 20) {
                        scoring::TwoBodyEnergyOP clash_energy = std::static_pointer_cast<scoring::TwoBodyEnergy>(pose.energy_manager().get_energy_tables(CLASH)[0]);

                        for(Size itmp=1; itmp<=pose.size(); ++itmp) {
                            for(Size jtmp=itmp+1; jtmp<=pose.size(); ++jtmp) {
                                Real pair_clash = clash_energy->get_pair_score(itmp, jtmp);
                                if(pair_clash > 0) {
                                    std::cout << itmp << "    ----    " << jtmp << "         " << pair_clash << std::endl;
                                }
                            }
                        }
                        pose.dump_pdb("test.pdb.gz", false, true, false, false);
                        exit(0);
                        // pose.dump_pdb(name.str(),false,opts.gzip);
                    }

                    
                    scoring::OneBodyEnergyOP voxel_clash_energy = std::static_pointer_cast<scoring::OneBodyEnergy>(pose.energy_manager().get_energy_tables(TARGET_CONTEXT_CLASH)[0]);

                    for(Size itmp=1; itmp<=pose.size(); ++itmp) {
                        Real voxel_clash_ires = voxel_clash_energy->get_score_per_residue(itmp);
                        if(voxel_clash_ires != 0) {
                            std::cout << pose.sequence() << std::endl;
                            std::cout << itmp << "    ----    " << voxel_clash_energy->get_score_per_residue(itmp) << std::endl;
                        }
                    }
                    */
                    
            }


            if( ii_inner % M == 0 && ii_inner != 0 ) {
                T /= decrease_factor;
                //change_root_prob /= decrease_factor;
                loop_weight /= decrease_factor;
                loop_weight = loop_weight>=1?loop_weight:1;
                frag_mover.bias_loop_sampling_by_dssp(dssp, loop_weight, 1);
            }

            // Real sc_neighbor = utils::sidechain_neighbors(pose);
            // Real target_clash_no_weight = std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash"))->score(pose);

            // Real rpx_score = std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->score(pose);

            // after each move, check satisfication and return if good
            bool pass_requirements(true);
            // std::cout << std::setprecision(5);
            // if(best_score > opts.total_score_cutoff) { pass_requirements = false; if (ii_inner % 200 ==0) std::cout << ii_outer << " <---> " << ii_inner << "     TOTAL FAILED:   " << best_score << "   <===>   " << rpx_score << "   <===>   " << sc_neighbor << "   <===>   " << best_rif_score << "   <===>   " << target_clash_no_weight << "    <===>   " << anchor_seqpos.size() << std::endl;}
            // if(pass_requirements && num_chains>1 && best_score_inter_chain > opts.inter_chain_score_cutoff) { pass_requirements=false; }
            // if(pass_requirements && utils::sidechain_neighbors(pose) < opts.sidechain_neighbor_cutoff) { pass_requirements=false; if (ii_inner % 200 ==0) std::cout << ii_outer << " <--->" << ii_inner << "     SIDECHAIN FAILED:   " << best_score << "   <===>   " << sc_neighbor << "   <===>   " << best_rif_score << "   <===>   " << target_clash_no_weight << "    <===>   " << anchor_seqpos.size() << std::endl;}
            // if(pass_requirements && best_rif_score > opts.rif_score_cutoff) { pass_requirements=false; if (ii_inner % 200 ==0)  std::cout << ii_outer << " <---> " << ii_inner << "     RIF FAILED:   " << best_score << "   <===>   " << sc_neighbor << "   <===>   " << best_rif_score << "   <===>   " << target_clash_no_weight << "    <===>   " << anchor_seqpos.size() << std::endl; }

            if(best_score > opts.total_score_cutoff) { pass_requirements = false; }
            if(pass_requirements && best_score_inter_chain > opts.inter_chain_score_cutoff) { pass_requirements = false; }
            if(pass_requirements && sidechain_neighbors.compute(pose) < opts.sidechain_neighbor_cutoff) { pass_requirements=false; }
            if(pass_requirements && opts.ligand_neighbors_cutoff != -1 && ligand_neighbors.compute(pose) < opts.ligand_neighbors_cutoff) { pass_requirements=false; }
            
            if(pass_requirements) {
                std::cout << "Solution found! Eearly stop." << std::endl;

                Real tmp_score = sfxn.score(pose) / protein_len;

                Real tmp_rpx = std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->score(pose)/pose.size();
                Real tmp_clash = std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->score(pose)/pose.size();

                Real tmp_sidechain_neighbors = sidechain_neighbors.compute(pose);

                Real tmp_target_clash = opts.context_clash_weight * std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash"))->score(pose);

                Real tmp_score_interchain = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2) ) / pose.size();

                Real tmp_ligand_neighbors = ligand_neighbors.compute(pose);

                std::cout <<std::setprecision(3)  << "Success <===> Total_Sc: " << std::left << std::setw(10) << tmp_score 
                                                  << "RPX: " << std::left << std::setw(10) << tmp_rpx 
                                                  << "Clash: " << std::left << std::setw(10) << tmp_clash 
                                                  << "Score_InterChain: " << std::left << std::setw(10) << tmp_score_interchain
                                                  << "NCAA_Clash: " << std::left << std::setw(10) << tmp_target_clash
                                                  << "SCNeighbors: " << std::left << std::setw(10) << tmp_sidechain_neighbors
                                                  << "LigandNeighbors: " << std::left << std::setw(10) << tmp_ligand_neighbors
                                                  <<  std::endl;

                return true;
            }

        } // inner loop
    } // outer loop

    std::cout << "No solution found! Too bad." << std::endl;

    Real tmp_score = sfxn.score(pose) / protein_len;

    Real tmp_rpx = std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->score(pose)/pose.size();
    Real tmp_clash = std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->score(pose)/pose.size();

    Real tmp_sidechain_neighbors = sidechain_neighbors.compute(pose);

    Real tmp_target_clash = opts.context_clash_weight * std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash"))->score(pose);

    Real tmp_score_interchain = (std::static_pointer_cast<scoring::RpxScoreMethod>(sfxn.get_method("rpx"))->base_score(pose,2)+std::static_pointer_cast<scoring::ClashScoreMethod>(sfxn.get_method("clash"))->base_score(pose,2) ) / pose.size();

    Real tmp_ligand_neighbors = ligand_neighbors.compute(pose);

    std::cout <<std::setprecision(3)  << "Failed <===> Total_Sc: " << std::left << std::setw(10) << tmp_score 
                                      << "RPX: " << std::left << std::setw(10) << tmp_rpx 
                                      << "Clash: " << std::left << std::setw(10) << tmp_clash 
                                      << "Score_InterChain: " << std::left << std::setw(10) << tmp_score_interchain
                                      << "NCAA_Clash: " << std::left << std::setw(10) << tmp_target_clash
                                      << "SCNeighbors: " << std::left << std::setw(10) << tmp_sidechain_neighbors
                                      << "LigandNeighbors: " << std::left << std::setw(10) << tmp_ligand_neighbors
                                      <<  std::endl;

    return false;
}


}

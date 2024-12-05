#ifndef INCLUDED_apps_args_hh
#define INCLUDED_apps_args_hh

#include "basic/types.hh"
#include "utils/string_util.hh"

#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>

using basic::Size;
using basic::Real;

struct Options {

    public:

        Size len;
        Size num_repeats;
        std::string symmetry;
        Size num_frags;
        std::string ss;
        std::string fake_ss;
        std::string seq;
        std::string frag_path;
        std::string dssp_file;
        std::string rpx_db;
        Real rpx_cart_resl;
        Real rpx_ang_resl;
        std::string rpx_db_refine;
        Real rpx_cart_resl_refine;
        Real rpx_ang_resl_refine;
        Real rpx_inter_chain_weight;
        Real temperature;

        std::string seeding_pdb;
        Real redundancy_rmsd_cutoff;

        bool early_stop;

        Real CB_swelling_factor;

        std::string designable_residues;

        // mcmc sampling
        Size mcmc_outer_cycles;
        Size mcmc_inner_cycles;
        Size max_try;

        // disulfide
        std::string disulfide_hash_table;
        Real disulfide_hash_cart_resl;
        Real disulfide_hash_ang_resl;
        Real disulfide_weight;
        Size disulfide_requirement_num;
        Size disulfide_min_sequence_separation;

        Real gold_burial_cutoff;


        // heterodimer
        Size chain1_len;
        std::string chain1_pdb;
        Size chain2_len;
        std::string chain2_pdb;

        // motif pair
        std::string motif1_pdb;
        std::string motif1_ss;
        std::string motif2_pdb;
        std::string motif2_ss;
        Real motif_pair_weight;
        Real motif_pair_dist_cutoff;


        // fiber
        Real fiber_rotation_angle_expectation;
        Real fiber_raise_dist_expectation;

        // privileged pair
        std::string privileged_pair_hash_table;
        Real privileged_pair_hash_cart_resl;
        Real privileged_pair_hash_ang_resl;
        Real privileged_pair_weight;
        Real privileged_pair_score_cutoff;
        Size privileged_pair_min_sequence_separation;

        Real privileged_pair_burial_cutoff;

        // GFP
        Size GFP_upper_helix_max_len;
        Size GFP_upper_helix_min_len;
        Size GFP_lower_loop_max_len;
        Size GFP_lower_loop_min_len;
        Size GFP_lower_helix_max_len;
        Size GFP_lower_helix_min_len;

        // metal
        std::string metal_hash_table;
        Real metal_hash_cart_resl;
        Real metal_hash_ang_resl;
        Real metal_coordination_score_weight;
        Size metal_coordination_number;
        Real metal_radius;
        std::string metal_pdb;
        Real metal_coordination_score_cutoff;

        // interface metals for yannan's project
        Size interface_metal_chain1_len;
        Size interface_metal_chain2_len;
        std::string interface_metal_scaffold_pdb;
        std::string interface_metal_chain1_pdb;
        std::string interface_metal_chain2_pdb;
        std::string interface_metal_config1;
        std::string interface_metal_config2;
        std::string interface_metal_type;
        Size num_interface_metals;
        Real interface_metal_score_cutoff;
        Real interface_metal_distance_optimization_weight;
        Real interface_metal_radius;
        std::string exclude_interface_metal_names;

        // SMIF
        std::string SMIF_scaffold_pdb;
        std::string SMIF_config1;
        std::string SMIF_config2;
        Real small_molecule_radius;
        Real num_small_molecules;
        Real SMIF_score_cutoff;
        Real SMIF_optimization_weight;
        std::string exclude_small_molecule_names;

        // vaccine scaffolding
        Size motif_len;
        Size context_len;

        // pdb extension
        std::string N_extension_ss_str;
        std::string C_extension_ss_str;
        std::string contact_residues;

        bool use_ss;
        bool hallucinate;
        std::string output_prefix;
        Size nstruct;
        Size max_outputs;

        // rif
        std::string rif_table;
        Real rif_weight;
        Real rif_cart_resl;
        Real rif_ang_resl;
        std::string context_pdb;
        Real context_clash_weight;
        Real context_clash_radius;
        bool ramp_context_clash_weight;

        Real rif_score_cutoff;


        // privileged motifs
        std::vector<std::string> privileged_motif_hash_tables;
        Real privileged_motif_hash_cart_resl;
        Real privileged_motif_hash_ang_resl;
        Real privileged_motif_weight;
        Real privileged_motif_score_cutoff;

        std::string privileged_interface_motif_hash_table;
        Real privileged_interface_motif_hash_cart_resl;
        Real privileged_interface_motif_hash_ang_resl;
        Real privileged_interface_motif_weight;
        Real privileged_interface_motif_score_cutoff;

        std::string ligand_name;
        std::string ligand_pdb;

        std::string extra_ligand;

        bool perturb_ligand;

        Real ligand_neighbors_cutoff;
        Real ligand2_neighbors_cutoff;

        Real motif_neighbors_cutoff;


        std::string output_dir;
        Real backbone_score_cutoff;
        Real sidechain_neighbor_cutoff;
        Real total_score_cutoff;
        Real intra_chain_score_cutoff;
        Real inter_chain_score_cutoff;
        Real C_symmetry_score_cutoff;
        Real D_symmetry_score_cutoff;
        Real sampling_stage_cutoff_tolerance_factor;
        Real clash_cutoff;

        // option for hallucinate
        Size num_helix;
        Size min_helix_len;
        Size max_helix_len;
        Size min_loop_len;
        Size max_loop_len;
        Size max_len;
        Size segment_num;
        // option for hallucinate:repeat
        Real repeat_min_radius,repeat_max_radius;
        Real repeat_min_rise,repeat_max_rise;
        Real repeat_min_omega,repeat_max_omega;
	    Real repeat_weight;
	    Real repeat_score_cutoff;

	// peptide
        Real peptide_score_cutoff;
        bool peptide;
        Size peptide_len;
        Size scaffold_num_repeats;
        Real peptide_diversify_angle;
        Real peptide_diversify_cart;
        Real peptide_diversify_step;
        Real peptide_diversify_cutoff;
        Size peptide_output_num;
        Real extra_clash_distance;
        // option for load motif
        std::string motif_pdb;
        std::string motif_ss;
        bool Nter_extension;
        bool Cter_extension;
        Size motif_pos;


        // options that mimic rosetta
        std::vector<std::string> input_pdbs;
        std::string score_file;


        // conformational switch
        Size wiggling_pos;
        Real delta_rpx_score;
        Real delta_sidechain_neighbors;
        Real minimum_rmsd;

        // NCAA
        std::string ncaa;
        Size ncaa_pos;


        //option for dock
        std::string prerpx1side_path;
        std::string rpx1side_path;
        std::string rpx1side_selected_res; // 1,2,3
        Real rpx1side_weight;
        Real intra_score_cutoff;
        Real inter_score_cutoff;

        std::string scaffold_pdb;
        std::string target_pdb;
        std::string grid_pdb;
        std::vector<Size> select_res;
        Real up_biased_factor;
        Real down_biased_factor;
        std::vector<Size> up_biased_res;
        std::vector<Size> down_biased_res;
        Size output_num;
        Real cart_sample_resl;
        Real angle_sample_resl;
        Real rpx1side_cutoff;
        std::vector<std::pair<Size,std::string>> center_atoms;
        Real sample_iteration;
        std::vector<Real> cart_refine_resl;
        std::vector<Real> angle_refine_resl;
        bool dump_init_pos;
        //option for output
        bool dump_trajactory;
        bool gzip;
        bool dump_polyG;
        Size verbose_level;
        // constructor
        Options();
        void parse_args(int argc, char *argv[]);
        void parse_flags(std::vector<std::string> args);
};


#endif

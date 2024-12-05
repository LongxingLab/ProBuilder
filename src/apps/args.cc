#include "apps/args.hh"
#include "basic/assert.hh"

#include <string>
#include <fstream>

Options::Options() :
    len(3),
    num_repeats(1),
    symmetry("C1"),
    num_frags(50),
    ss(""),
    fake_ss(""),
    seq(""),
    frag_path(""),
    dssp_file(""),
    rpx_db(""),
    rpx_cart_resl(2.0),
    rpx_ang_resl(26),
    rpx_db_refine(""),
    rpx_cart_resl_refine(1.0),
    rpx_ang_resl_refine(16),
    rpx_inter_chain_weight(1.0),
    temperature(-1.0),

    seeding_pdb(""),
    redundancy_rmsd_cutoff(2.0),

    early_stop(false),

    CB_swelling_factor(1.0),

    designable_residues(""),

    mcmc_outer_cycles(3),
    mcmc_inner_cycles(-1),
    max_try(-1),

    disulfide_hash_table(""),
    disulfide_hash_cart_resl(1.0),
    disulfide_hash_ang_resl(16.0),
    disulfide_weight(5.0),
    disulfide_requirement_num(0),
    disulfide_min_sequence_separation(8),
    gold_burial_cutoff(9999.99),

    // heterodimer
    chain1_len(-1),
    chain1_pdb(""),
    chain2_len(-1),
    chain2_pdb(""),

    // motif pair
    motif1_pdb(""),
    motif1_ss(""),
    motif2_pdb(""),
    motif2_ss(""),
    motif_pair_weight(1.0),
    motif_pair_dist_cutoff(999999),

    // fiber
    fiber_rotation_angle_expectation(-1.0),
    fiber_raise_dist_expectation(-1.0),

    privileged_pair_hash_table(""),
    privileged_pair_hash_cart_resl(1.0),
    privileged_pair_hash_ang_resl(16.0),
    privileged_pair_weight(5.0),
    privileged_pair_score_cutoff(0.0),
    privileged_pair_min_sequence_separation(8),
    privileged_pair_burial_cutoff(0.0),

    // GFP
    GFP_upper_helix_max_len(7),
    GFP_upper_helix_min_len(1),
    GFP_lower_loop_max_len(5),
    GFP_lower_loop_min_len(0),
    GFP_lower_helix_max_len(7),
    GFP_lower_helix_min_len(4),

    metal_hash_table(""),
    metal_hash_cart_resl(0.5),
    metal_hash_ang_resl(15.0),
    metal_coordination_score_weight(1.0),
    metal_coordination_number(4),
    metal_radius(1.4),
    metal_pdb(""),
    metal_coordination_score_cutoff(0.0),

    // partially satisfied zinc
    interface_metal_chain1_len(1),
    interface_metal_chain2_len(1),
    interface_metal_scaffold_pdb(""),
    interface_metal_chain1_pdb(""),
    interface_metal_chain2_pdb(""),
    interface_metal_config1(""),
    interface_metal_config2(""),
    interface_metal_type(""),
    num_interface_metals(1),
    interface_metal_score_cutoff(1.5),
    interface_metal_distance_optimization_weight(1.0),
    interface_metal_radius(4.0),
    exclude_interface_metal_names(""),

    // smif
    SMIF_scaffold_pdb(""),
    SMIF_config1(""),
    SMIF_config2(""),
    small_molecule_radius(3.0),
    num_small_molecules(1),
    SMIF_score_cutoff(1.5),
    SMIF_optimization_weight(1.0),
    exclude_small_molecule_names(""),

    motif_len(-1),
    motif_ss(""),
    Nter_extension(false),
    Cter_extension(false),
    context_len(-1),

    C_extension_ss_str(""),
    N_extension_ss_str(""),
    contact_residues(""),

    use_ss(false),
    hallucinate(false),
    gzip(false),
    dump_polyG(false),
    output_prefix(""),
    nstruct(1),
    max_outputs(9e9),
    motif_pos(-1),
    rif_table(""),
    rif_weight(1.0),
    rif_cart_resl(1.0),
    rif_ang_resl(16.0),
    context_pdb(""),
    context_clash_weight(1.0),
    context_clash_radius(0.9),
    ramp_context_clash_weight(false),
    rif_score_cutoff(0.0),
    ligand_name(""),
    ligand_pdb(""),
    extra_ligand(""),
    perturb_ligand(false),
    ligand_neighbors_cutoff(-1.0),
    ligand2_neighbors_cutoff(-1.0),
    motif_neighbors_cutoff(-1.0),
    motif_pdb(""),
    output_dir(""),
    backbone_score_cutoff(0.0),
    sidechain_neighbor_cutoff(0.0),
    total_score_cutoff(0.0),
    intra_chain_score_cutoff(0.0),
    inter_chain_score_cutoff(0.0),
    C_symmetry_score_cutoff(0.0),
    D_symmetry_score_cutoff(0.0),
    clash_cutoff(0.1),
    sampling_stage_cutoff_tolerance_factor(1.0),


    privileged_motif_hash_cart_resl(1.0),
    privileged_motif_hash_ang_resl(16.0),
    privileged_motif_weight(1.0),
    privileged_motif_score_cutoff(0.0),

    privileged_interface_motif_hash_cart_resl(1.0),
    privileged_interface_motif_hash_ang_resl(16.0),
    privileged_interface_motif_weight(1.0),
    privileged_interface_motif_score_cutoff(0.0),

    // options for rpx scoring
    score_file("scores.sc"),

    // conformational switch
    wiggling_pos(1),
    delta_rpx_score(0.0),
    delta_sidechain_neighbors(0.0),
    minimum_rmsd(0.0),

    // ncaa
    ncaa(""),
    ncaa_pos(-1),

    // option for hallucinate
    num_helix(5),
    min_helix_len(5),
    max_helix_len(9999),
    min_loop_len(2),
    max_loop_len(6),
    max_len(1000),
    segment_num(-1),

    // repeat
    repeat_min_radius(0),
    repeat_max_radius(10000),
    repeat_min_rise(0),
    repeat_max_rise(10000),
    repeat_min_omega(0),
    repeat_max_omega(1000),
    repeat_weight(1.0),
    repeat_score_cutoff(9e9),


    extra_clash_distance(15),
    peptide(false),
    peptide_diversify_angle(16.0),
    peptide_diversify_cart(1.0),
    peptide_diversify_step(1),
    peptide_diversify_cutoff(0.5),
    peptide_output_num(100),

    // fold on traget
    rpx1side_path(""),
    rpx1side_selected_res(""),
    rpx1side_weight(1.0),
    intra_score_cutoff(0.0),
    inter_score_cutoff(0.0),

    //option for dock
    up_biased_factor(2.0),
    down_biased_factor(0.5),
    rpx1side_cutoff(0),
    //output
    dump_trajactory(false),
    verbose_level(-1)
{}

void Options::parse_args(int argc, char *argv[])
{

    std::vector<std::string> flags_temp(argv+1, argv+argc);
    std::vector<std::string> flags;

    // in case there is a flag file
    for(Size ii=0; ii<flags_temp.size(); ++ii) {
        if(flags_temp[ii][0] == '@') {
            // read in the file
            std::string flag_fname = flags_temp[ii].substr(1);
            std::ifstream f(flag_fname);
            if(f.fail()) {
                std::cout << "Failed to open the flag file: " << flag_fname << std::endl;
                exit(0);
            }
            for(std::string line; std::getline(f,line);) {
                utils::trim(line);
                if(line=="") {
                    continue;
                }

                if(line[0]=='#') {
                    continue;
                } else if(line[0]=='-') {
                    std::istringstream iss(line);
                    std::string token;
                    while (iss >> token) {
                        if (!token.empty()) {
                            flags.push_back(token);
                        }
                    }
                } else {
                    std::cout << "Unknown flag in flag file " << flag_fname << "  ===>  " << line << std::endl;
                    exit(0);
                }

            }
            f.close();

        } else {
            flags.push_back(flags_temp[ii]);
        }
    }

    parse_flags(flags);

    if(motif_ss == "" && motif_len != -1) {
        for(Size ii=0; ii<motif_len; ++ii) motif_ss += "H";
    }

    if (ss != "" && seq == "" && !hallucinate ) {
        // auto fill seq L => G, H/E => V
        // assert(ss.length() == len);
        for( Size ii=0; ii<len; ++ii ) {
            if(ss.at(ii) == 'L')
                seq += "G";
            else
                seq += "V";
        }
    }
}

void Options::parse_flags(std::vector<std::string> args) {

    for(Size ii=0; ii<args.size(); ++ii)
    {
        std::string val = args[ii];
        if(val[0]=='-')
        {
            if(val == "-len") { len = stoi(args[ii+1]); ++ii; continue;}
            if(val == "-num_repeats") { num_repeats = stoi(args[ii+1]); ++ii; continue;}
            if(val == "-symmetry") { symmetry = args[ii+1]; ++ii; continue;}
            if(val == "-num_frags") { num_frags = stoi(args[ii+1]); ++ii; continue;}
            if(val == "-ss") { ss = args[ii+1]; ++ii; continue;}
            if(val == "-fake_ss") { fake_ss = args[ii+1]; ++ii; continue;}
            if(val == "-seq") {seq = args[ii+1]; ++ii; continue;}
            if(val == "-frag_path") {frag_path = args[ii+1]; ++ii; continue;}
            if(val == "-dssp_file") {dssp_file = args[ii+1]; ++ii; continue;}
            if(val == "-rpx_db") {rpx_db = args[ii+1]; ++ii; continue;}
            if(val == "-rpx_cart_resl") {rpx_cart_resl = stof(args[ii+1]); ++ii; continue;}
            if(val == "-rpx_ang_resl") {rpx_ang_resl = stof(args[ii+1]); ++ii; continue;}
            if(val == "-rpx_db_refine") {rpx_db_refine = args[ii+1]; ++ii; continue;}
            if(val == "-rpx_cart_resl_refine") {rpx_cart_resl_refine = stof(args[ii+1]); ++ii; continue;}
            if(val == "-rpx_ang_resl_refine") {rpx_ang_resl_refine = stof(args[ii+1]); ++ii; continue;}
            if(val == "-rpx_inter_chain_weight") {rpx_inter_chain_weight = stof(args[ii+1]); ++ii; continue;}
            if(val == "-temperature") {temperature = stof(args[ii+1]); ++ii; continue;}

            if(val == "-seeding_pdb") { seeding_pdb = args[ii+1]; ++ii; continue;}
            if(val == "-redundancy_rmsd_cutoff") { redundancy_rmsd_cutoff = stof(args[ii+1]); ++ii; continue;}

            if(val == "-early_stop") {early_stop=true; continue;};

            if(val == "-CB_swelling_factor") {CB_swelling_factor = stof(args[ii+1]); ++ii; continue;}

            if(val == "-designable_residues") {designable_residues = args[ii+1]; ++ii; continue;}

            // sampling
            if(val == "-mcmc_outer_cycles") { mcmc_outer_cycles = stoi(args[ii+1]); ++ii; continue; }
            if(val == "-mcmc_inner_cycles") { mcmc_inner_cycles = stoi(args[ii+1]); ++ii; continue;}
            if(val == "-max_try") { max_try= stoi(args[ii+1]); ++ii; continue;}

            // disulfide bond
            if(val == "-disulfide_hash_table") {disulfide_hash_table = args[ii+1]; ++ii; continue;}
            if(val == "-disulfide_hash_cart_resl") {disulfide_hash_cart_resl = stof(args[ii+1]); ++ii; continue;}
            if(val == "-disulfide_hash_ang_resl") {disulfide_hash_ang_resl = stof(args[ii+1]); ++ii; continue;}
            if(val == "-disulfide_weight") {disulfide_weight = stof(args[ii+1]); ++ii; continue;}
            if(val == "-disulfide_requirement_num") {disulfide_requirement_num = stoi(args[ii+1]); ++ii; continue;}
            if(val == "-disulfide_min_sequence_separation") {disulfide_min_sequence_separation = stoi(args[ii+1]); ++ii; continue;}

            if(val == "-gold_burial_cutoff") {gold_burial_cutoff = stof(args[ii+1]); ++ii; continue;}

            // heterodimer
            if(val == "-chain1_len") {chain1_len = stoi(args[ii+1]); ++ii; continue;}
            if(val == "-chain1_pdb") {chain1_pdb = args[ii+1]; ++ii; continue;}
            if(val == "-chain2_len") {chain2_len = stoi(args[ii+1]); ++ii; continue;}
            if(val == "-chain2_pdb") {chain2_pdb = args[ii+1]; ++ii; continue;}

            // motif pair
            if(val == "-motif1_pdb") {motif1_pdb = args[ii+1]; ++ii; continue;}
            if(val == "-motif1_ss") {motif1_ss = args[ii+1]; ++ii; continue;}
            if(val == "-motif2_pdb") {motif2_pdb = args[ii+1]; ++ii; continue;}
            if(val == "-motif2_ss") {motif2_ss = args[ii+1]; ++ii; continue;}
            if(val == "-motif_pair_weight") {motif_pair_weight = stof(args[ii+1]); ++ii; continue;}
            if(val == "-motif_pair_dist_cutoff") {motif_pair_dist_cutoff = stof(args[ii+1]); ++ii; continue;}

            if(val == "-fiber_rotation_angle_expectation") {fiber_rotation_angle_expectation = stof(args[ii+1]); ++ii; continue;}
            if(val == "-fiber_raise_dist_expectation") {fiber_raise_dist_expectation = stof(args[ii+1]); ++ii; continue;}

            // conformational switch
            if(val == "-wiggling_pos") {wiggling_pos = stoi(args[ii+1]); ++ii; continue;}
            if(val == "-delta_rpx_score") {delta_rpx_score = stof(args[ii+1]); ++ii; continue;}
            if(val == "-delta_sidechain_neighbors") {delta_sidechain_neighbors = stof(args[ii+1]); ++ii; continue;}
            if(val == "-minimum_rmsd") {minimum_rmsd = stof(args[ii+1]); ++ii; continue;}

            // ncaa
            if(val == "-ncaa_pos") {ncaa_pos = stoi(args[ii+1]); ++ii; continue;}
            if(val == "-ncaa") {ncaa = args[ii+1]; ++ii; continue;}


            // privileged motif score
            if(val == "-privileged_motif_hash_tables") {
                ++ii;
                std::string temp_string(args[ii]);
                while (temp_string.at(0) != '-') {
                    privileged_motif_hash_tables.push_back(args[ii]);
                    ++ii;
                    if( ii<args.size() )
                        temp_string = args[ii];
                    else
                        break;
                }
                --ii;
                continue;
            }
            if(val == "-privileged_motif_hash_cart_resl") {privileged_motif_hash_cart_resl = stof(args[ii+1]); ++ii; continue;}
            if(val == "-privileged_motif_hash_ang_resl") {privileged_motif_hash_ang_resl = stof(args[ii+1]); ++ii; continue;}
            if(val == "-privileged_motif_weight") {privileged_motif_weight = stof(args[ii+1]); ++ii; continue;}
            if(val == "-privileged_motif_score_cutoff") {privileged_motif_score_cutoff = stof(args[ii+1]); ++ii; continue;}

            if(val == "-privileged_interface_motif_hash_table") {privileged_interface_motif_hash_table = args[ii+1]; ++ii; continue;}
            if(val == "-privileged_interface_motif_hash_cart_resl") {privileged_interface_motif_hash_cart_resl = stof(args[ii+1]); ++ii; continue;}
            if(val == "-privileged_interface_motif_hash_ang_resl") {privileged_interface_motif_hash_ang_resl = stof(args[ii+1]); ++ii; continue;}
            if(val == "-privileged_interface_motif_weight") {privileged_interface_motif_weight = stof(args[ii+1]); ++ii; continue;}
            if(val == "-privileged_interface_motif_score_cutoff") {privileged_interface_motif_score_cutoff = stof(args[ii+1]); ++ii; continue;}

            // privileged pair score cutoff
            if(val == "-privileged_pair_hash_table") {privileged_pair_hash_table = args[ii+1]; ++ii; continue;}
            if(val == "-privileged_pair_hash_cart_resl") {privileged_pair_hash_cart_resl = stof(args[ii+1]); ++ii; continue;}
            if(val == "-privileged_pair_hash_ang_resl") {privileged_pair_hash_ang_resl = stof(args[ii+1]); ++ii; continue;}
            if(val == "-privileged_pair_weight") {privileged_pair_weight = stof(args[ii+1]); ++ii; continue;}
            if(val == "-privileged_pair_score_cutoff") {privileged_pair_score_cutoff = stof(args[ii+1]); ++ii; continue;}
            if(val == "-privileged_pair_min_sequence_separation") {privileged_pair_min_sequence_separation = stoi(args[ii+1]); ++ii; continue;}
            if(val == "-privileged_pair_burial_cutoff") {privileged_pair_burial_cutoff = stof(args[ii+1]); ++ii; continue;}

            // metal
            if(val == "-metal_hash_table") {metal_hash_table = args[ii+1]; ++ii; continue;}
            if(val == "-metal_hash_cart_resl") {metal_hash_cart_resl = stof(args[ii+1]); ++ii; continue;}
            if(val == "-metal_hash_ang_resl") {metal_hash_ang_resl = stof(args[ii+1]); ++ii; continue;}
            if(val == "-metal_coordination_score_weight") {metal_coordination_score_weight = stof(args[ii+1]); ++ii; continue;}
            if(val == "-metal_coordination_number") {metal_coordination_number = stoi(args[ii+1]); ++ii; continue;}
            if(val == "-metal_radius") {metal_radius = stof(args[ii+1]); ++ii; continue;}
            if(val == "-metal_pdb") {metal_pdb = args[ii+1]; ++ii; continue;}
            if(val == "-metal_coordination_score_cutoff") {metal_coordination_score_cutoff = stof(args[ii+1]); ++ii; continue;}

            // GFP
            if(val == "-GFP_upper_helix_max_len") {GFP_upper_helix_max_len = stoi(args[ii+1]); ++ii; continue;}
            if(val == "-GFP_upper_helix_min_len") {GFP_upper_helix_min_len = stoi(args[ii+1]); ++ii; continue;}
            if(val == "-GFP_lower_loop_max_len")  {GFP_lower_loop_max_len = stoi(args[ii+1]); ++ii; continue;}
            if(val == "-GFP_lower_loop_min_len") {GFP_lower_loop_min_len = stoi(args[ii+1]); ++ii; continue;}
            if(val == "-GFP_lower_helix_max_len") {GFP_lower_helix_max_len = stoi(args[ii+1]); ++ii; continue;}
            if(val == "-GFP_lower_helix_min_len") {GFP_lower_helix_min_len = stoi(args[ii+1]); ++ii; continue;}

      
            // interface metal for yannan
            if(val == "-interface_metal_chain1_len") {interface_metal_chain1_len = stoi(args[ii+1]); ++ii; continue;}
            if(val == "-interface_metal_chain2_len") {interface_metal_chain2_len = stoi(args[ii+1]); ++ii; continue;}
            if(val == "-interface_metal_scaffold_pdb") {interface_metal_scaffold_pdb = args[ii+1]; ++ii; continue;}
            if(val == "-interface_metal_chain1_pdb") {interface_metal_chain1_pdb = args[ii+1]; ++ii; continue;}
            if(val == "-interface_metal_chain2_pdb") {interface_metal_chain2_pdb = args[ii+1]; ++ii; continue;}
            if(val == "-interface_metal_config1") {interface_metal_config1 = args[ii+1]; ++ii; continue;}
            if(val == "-interface_metal_config2") {interface_metal_config2 = args[ii+1]; ++ii; continue;}
            if(val == "-interface_metal_type") {interface_metal_type = args[ii+1]; ++ii; continue;}
            if(val == "-num_interface_metals") {num_interface_metals = stoi(args[ii+1]); ++ii; continue;}
            if(val == "-interface_metal_score_cutoff") {interface_metal_score_cutoff = stof(args[ii+1]); ++ii; continue;}
            if(val == "-interface_metal_distance_optimization_weight") {interface_metal_distance_optimization_weight = stof(args[ii+1]); ++ii; continue;}
            if(val == "-interface_metal_radius") {interface_metal_radius = stof(args[ii+1]); ++ii; continue;}
            if(val == "-exclude_interface_metal_names") {exclude_interface_metal_names = args[ii+1]; ++ii; continue;}

            // smif
            if(val == "-SMIF_scaffold_pdb") {SMIF_scaffold_pdb = args[ii+1]; ++ii; continue;}
            if(val == "-SMIF_config1") {SMIF_config1 = args[ii+1]; ++ii; continue;}
            if(val == "-SMIF_config2") {SMIF_config2 = args[ii+1]; ++ii; continue;}
            if(val == "-num_small_molecules") {num_small_molecules = stoi(args[ii+1]); ++ii; continue;}
            if(val == "-SMIF_score_cutoff") {SMIF_score_cutoff = stof(args[ii+1]); ++ii; continue;}
            if(val == "-SMIF_optimization_weight") {SMIF_optimization_weight = stof(args[ii+1]); ++ii; continue;}
            if(val == "-small_molecule_radius") {small_molecule_radius = stof(args[ii+1]); ++ii; continue;}
            if(val == "-exclude_small_molecule_names") {exclude_small_molecule_names = args[ii+1]; ++ii; continue;}



            if(val == "-N_extension_ss_str") {N_extension_ss_str = args[ii+1]; ++ii; continue;}
            if(val == "-C_extension_ss_str") {C_extension_ss_str = args[ii+1]; ++ii; continue;}
            if(val == "-contact_residues") {contact_residues = args[ii+1]; ++ii; continue;}

            if(val == "-motif_len") {motif_len = stoi(args[ii+1]); ++ii; continue;}
            if(val == "-motif_ss") {motif_ss = args[ii+1]; ++ii; continue;}
            if(val == "-Nter_extension") {Nter_extension = true; continue;}
            if(val == "-Cter_extension") {Cter_extension = true; continue;}
            if(val == "-context_len") {context_len = stoi(args[ii+1]); ++ii; continue;}

            if(val == "-use_ss") {use_ss = true; continue;}
            if(val == "-hallucinate") {hallucinate= true; continue;}
            if(val == "-output_prefix") {output_prefix = args[ii+1]; ++ii; continue;}
            if(val == "-nstruct") {nstruct = stoi(args[ii+1]); ++ii; continue;}
            if(val == "-max_outputs") {max_outputs = stoi(args[ii+1]); ++ii; continue;}
            if(val == "-rif_table") {rif_table = args[ii+1]; ++ii; continue;}
            if(val == "-rif_weight") {rif_weight = stof(args[ii+1]); ++ii; continue;}
            if(val == "-rif_cart_resl") {rif_cart_resl = stof(args[ii+1]); ++ii; continue;}
            if(val == "-rif_ang_resl") {rif_ang_resl = stof(args[ii+1]); ++ii; continue;}
            if(val == "-ligand_name") {ligand_name = args[ii+1]; ++ii; continue;}
            if(val == "-ligand_pdb") {ligand_pdb = args[ii+1]; ++ii; continue;}
            if(val == "-extra_ligand") {extra_ligand = args[ii+1]; ++ii; continue;}
            if(val == "-perturb_ligand") {perturb_ligand = true; continue;}
            if(val == "-ligand_neighbors_cutoff") {ligand_neighbors_cutoff = stof(args[ii+1]); ++ii; continue;}
            if(val == "-ligand2_neighbors_cutoff") {ligand2_neighbors_cutoff = stof(args[ii+1]); ++ii; continue;}
            if(val == "-motif_neighbors_cutoff") {motif_neighbors_cutoff = stof(args[ii+1]); ++ii; continue;}
            if(val == "-context_pdb") {context_pdb = args[ii+1]; ++ii; continue;}
            if(val == "-context_clash_weight") {context_clash_weight = stof(args[ii+1]); ++ii; continue;}
            if(val == "-context_clash_radius") {context_clash_radius = stof(args[ii+1]); ++ii; continue;}
            if(val == "-ramp_context_clash_weight") {ramp_context_clash_weight = true; continue;}
            if(val == "-rif_score_cutoff") {rif_score_cutoff = stof(args[ii+1]); ++ii; continue;}
            if(val == "-motif_pdb") {motif_pdb = args[ii+1]; ++ii; continue;}
            if(val == "-motif_pos") {motif_pos = stoi(args[ii+1]); ++ii; continue;}
            // if(val == "-motif_ss") {motif_ss = args[ii+1]; ++ii; continue;}
            if(val == "-clash_cutoff") {clash_cutoff = stof(args[ii+1]); ++ii; continue;}
            if(val == "-output_dir") {output_dir = args[ii+1]; ++ii; continue;}
            if(val == "-sidechain_neighbor_cutoff") {sidechain_neighbor_cutoff = stof(args[ii+1]); ++ii; continue;}
            if(val == "-total_score_cutoff") {total_score_cutoff = stof(args[ii+1]); ++ii; continue;}
            if(val == "-backbone_score_cutoff") {backbone_score_cutoff = stof(args[ii+1]); ++ii; continue;}
            if(val == "-intra_chain_score_cutoff") {intra_chain_score_cutoff = stof(args[ii+1]); ++ii; continue;}
            if(val == "-inter_chain_score_cutoff") {inter_chain_score_cutoff = stof(args[ii+1]); ++ii; continue;}
            if(val == "-C_symmetry_score_cutoff") {C_symmetry_score_cutoff = stof(args[ii+1]); ++ii; continue;}
            if(val == "-D_symmetry_score_cutoff") {D_symmetry_score_cutoff = stof(args[ii+1]); ++ii; continue;}
            if(val == "-sampling_stage_cutoff_tolerance_factor") {sampling_stage_cutoff_tolerance_factor = stof(args[ii+1]); ++ii; continue;}

            // options for scoring
            if(val == "-s") {
                ++ii;
                std::string temp_string(args[ii]);
                while (temp_string.at(0) != '-') {
                    input_pdbs.push_back(args[ii]);
                    ++ii;
                    if( ii<args.size() )
                        temp_string = args[ii];
                    else
                        break;
                }
                --ii;
                continue;
            }
            if(val == "-l") {
                std::ifstream f(args[ii+1]);
                std::string temp_fname;
                while(f >> temp_fname) input_pdbs.push_back(temp_fname);
                ++ii;
                continue;
            }
            if(val == "-score_file") {score_file = args[ii+1]; ++ii; continue;}

            // option for hallucinate
            if(val == "-num_helix") {num_helix = stoi(args[ii+1]); ++ii; continue;};
            if(val == "-min_helix_len") {min_helix_len = stoi(args[ii+1]); ++ii; continue;};
            if(val == "-max_helix_len") {max_helix_len = stoi(args[ii+1]); ++ii; continue;};
            if(val == "-min_loop_len") {min_loop_len = stoi(args[ii+1]); ++ii; continue;};
            if(val == "-max_loop_len") {max_loop_len = stoi(args[ii+1]); ++ii; continue;};
            if(val == "-max_len") {max_len = stoi(args[ii+1]); ++ii; continue;};
            if(val == "-segment_num") {segment_num = stoi(args[ii+1]); ++ii; continue;}
            // option for hallucinate:repeat
            if(val == "-repeat_min_radius") {repeat_min_radius = stof(args[ii+1]); ++ii; continue;};
            if(val == "-repeat_max_radius") {repeat_max_radius = stof(args[ii+1]); ++ii; continue;};
            if(val == "-repeat_min_rise") {repeat_min_rise = stof(args[ii+1]); ++ii; continue;};
            if(val == "-repeat_max_rise") {repeat_max_rise = stof(args[ii+1]); ++ii; continue;};
            if(val == "-repeat_min_omega") {repeat_min_omega = stof(args[ii+1]); ++ii; continue;};
            if(val == "-repeat_max_omega") {repeat_max_omega = stof(args[ii+1]); ++ii; continue;};
            if(val == "-repeat_weight") {repeat_weight = stof(args[ii+1]); ++ii; continue;};
            if(val == "-repeat_score_cutoff") {repeat_score_cutoff = stof(args[ii+1]); ++ii; continue;};
            // option for hallucinate:repeat peptide
            if(val == "-peptide") {peptide=true; continue;};
            if(val == "-scaffold_num_repeats") {scaffold_num_repeats = stoi(args[ii+1]); ++ii; continue;};
            if(val == "-peptide_output_num") {peptide_output_num = stoi(args[ii+1]); ++ii; continue;};
            if(val == "-extra_clash_distance") {extra_clash_distance = stof(args[ii+1]); ++ii; continue;};
            if(val == "-peptide_score_cutoff") {peptide_score_cutoff = stof(args[ii+1]); ++ii; continue;};
            if(val == "-peptide_diversify_angle") {peptide_diversify_angle = stof(args[ii+1]); ++ii; continue;};
            if(val == "-peptide_diversify_cart") {peptide_diversify_cart = stof(args[ii+1]); ++ii; continue;};
            if(val == "-peptide_diversify_step") {peptide_diversify_step = stoi(args[ii+1]); ++ii; continue;};
            if(val == "-peptide_diversify_cutoff") {peptide_diversify_cutoff = stof(args[ii+1]); ++ii; continue;};
            if(val == "-peptide_len") {peptide_len = stoi(args[ii+1]); ++ii; continue;};
            //option for dock
            if(val == "-prerpx1side_path") {prerpx1side_path = args[ii+1]; ++ii; continue;};
            if(val == "-rpx1side_weight") {rpx1side_weight = stof(args[ii+1]); ++ii; continue;}
            if(val == "-intra_score_cutoff") {intra_score_cutoff = stof(args[ii+1]); ++ii; continue;}
            if(val == "-inter_score_cutoff") {inter_score_cutoff = stof(args[ii+1]); ++ii; continue;}
            if(val == "-rpx1side_path") {rpx1side_path = args[ii+1]; ++ii; continue;};
            if(val == "-rpx1side_selected_res") {rpx1side_selected_res = args[ii+1]; ++ii; continue;};
            if(val == "-rpx1side_cutoff") {rpx1side_cutoff = stof(args[ii+1]); ++ii; continue;};
            if(val == "-scaffold_pdb") {scaffold_pdb = args[ii+1]; ++ii; continue;};
            if(val == "-target_pdb") {target_pdb = args[ii+1]; ++ii; continue;};
            if(val == "-grid_pdb") {grid_pdb = args[ii+1]; ++ii; continue;};
            if(val == "-up_biased_factor") {up_biased_factor = stof(args[ii+1]); ++ii; continue;};
            if(val == "-down_biased_factor") {down_biased_factor = stof(args[ii+1]); ++ii; continue;};
            if(val == "-select_res") {
                ii++;
                while(ii<args.size()&&args[ii].at(0)!='-'){
                    select_res.push_back(std::stoi(args[ii])); ++ii;
                }
                ii--;
                continue;
            };
            if(val == "-up_biased_res") {
                ii++;
                while(ii<args.size()&&args[ii].at(0)!='-'){
                    up_biased_res.push_back(std::stoi(args[ii])); ++ii;
                }
                ii--;
                continue;
            };
            if(val == "-down_biased_res") {
                ii++;
                while(ii<args.size()&&args[ii].at(0)!='-'){
                    down_biased_res.push_back(std::stoi(args[ii])); ++ii;
                }
                ii--;
                continue;
            };
            if(val == "-cart_sample_resl") {cart_sample_resl = stof(args[ii+1]); ++ii; continue;};
            if(val == "-angle_sample_resl") {angle_sample_resl = stof(args[ii+1]); ++ii; continue;};
            if(val == "-output_num") {output_num = stoi(args[ii+1]); ++ii; continue;};
            if(val == "-gzip") {gzip=true; continue;};
            if(val == "-dump_polyG") {dump_polyG=true; continue;};
            if(val == "-verbose_level") {verbose_level=stoi(args[ii+1]); ++ii; continue;};
            if(val == "-dump_trajactory") {dump_trajactory=true; continue;};
            std::cout << "Unknown argument:   " << val << std::endl;
            exit(0);
        }
    }
}

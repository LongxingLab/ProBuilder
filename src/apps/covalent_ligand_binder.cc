#include "basic/types.hh"
#include "scene/Residue.hh"
#include "scene/Conformation.hh"
#include "scene/Pose.hh"
#include "basic/macros.hh"
#include "utils/hash_util.hh"
#include "sampling/FragMover.hh"
#include "sampling/SSDatabase.hh"
#include "sampling/mcmc_protocols.hh"
#include "scoring/ScoreFunction.hh"
#include "utils/utils.hh"
#include "basic/assert.hh"
#include "utils/random_util.hh"
#include "basic/assert.hh"
#include "utils/dssp_util.hh"

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

// for monomers only

int main(int argc, char *argv[])
{
    Options opts;
    opts.parse_args(argc, argv);

    std::cout << "Protein length: " << opts.len << std::endl
              << "Num Repeats: "    << opts.num_repeats << std::endl
              << "Symmetry: "       << opts.symmetry << std::endl;


    // repeats == 1
    assert(opts.num_repeats == 1);

    assert(opts.num_repeats == 1 && opts.symmetry == "C1");
    std::cout << "Prepare the pose ...\n";
    Pose pose(opts.len, opts.num_repeats, opts.symmetry);
    if(!opts.hallucinate) {
        pose.set_dssp(opts.ss);
        pose.set_sequence(opts.seq);
    }

    SSDatabase ss_db(opts.len);
    if(opts.dssp_file != "") {
        ss_db.load_database(opts.dssp_file);
    }

    assert(opts.ncaa != "");
    assert(opts.ncaa_pos == -1 || (opts.ncaa_pos>=1 && opts.ncaa_pos<=opts.len) );

    std::cout << "Prepare the score function ...\n";
    std::cout << "Use secondary informtion for rpx scoring: " << opts.use_ss << std::endl;
    ScoreFunction sfxn;
    sfxn.regist_method("clash",ScoreMethodOP(new ClashScoreMethod()));
    std::static_pointer_cast<ClashScoreMethod>(sfxn.get_method("clash"))->set_CB_swelling_factor(opts.CB_swelling_factor);
    sfxn.regist_method("rpx",ScoreMethodOP(new RpxScoreMethod(opts.rpx_db, opts.rpx_cart_resl, opts.rpx_ang_resl, opts.use_ss, false)));

    //
    // load the rif table and context pdb into voxel array
    assert(opts.rif_table != "" && opts.context_pdb != "");

    Size num_ligands(1);
    std::cout << "Loading rif table " << opts.rif_table << std::endl;
    sfxn.regist_method("target_rif",ScoreMethodOP(new RifScoreMethod(opts.rif_table, opts.rif_cart_resl, opts.rif_ang_resl)));
    sfxn.get_method("target_rif")->set_weight(opts.rif_weight);
    sfxn.get_method("target_rif")->set_score_type(TARGET_RIF);
    

    std::cout << "Loading ncaa side chain pdb " << opts.context_pdb << " for clash check" << std::endl;
    sfxn.regist_method("target_context_clash",ScoreMethodOP(new VoxelClashScoreMethod(opts.context_pdb,"ALL",opts.context_clash_radius,0.25F,false,"MOLECULE")));
    sfxn.get_method("target_context_clash")->set_score_type(TARGET_CONTEXT_CLASH);
    sfxn.get_method("target_context_clash")->set_weight(opts.context_clash_weight);
    // std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash"))->visualize_voxel_grid("clash_grid.pdb");
    


    std::cout << "Load the fragment_library: " << opts.frag_path << std::endl;
    FragMapMover frag_mover(opts.len, 7);
    frag_mover.load_fragments(opts.frag_path);


    metric::LigandNeighbors ligand_neighbors;
    if( opts.ligand_pdb != "" ) {
        ligand_neighbors.load_ligand(opts.ligand_pdb);
    }
    metric::SidechainNeighbors sidechain_neighbors;
    
    auto t1 = Clock::now();

    std::string hallucinate_dssp_prefix = "";
    Size insert_pos;

    for(Size itry = 1; itry <= opts.nstruct; ++itry) {
        
        if(opts.hallucinate && (itry-1) % 20 == 0){

            // prepare a new pos
            pose = Pose(opts.len, opts.num_repeats, opts.symmetry);
            pose.energy_manager().add_energy_onebody(TARGET_RIF);
            pose.energy_manager().add_energy_onebody(TARGET_CONTEXT_CLASH);


            std::string sequence = "";
            std::string dssp, dssp_prefix;
            // std::string dssp = random_helix_dssp(opts.min_helix_len,opts.max_helix_len,opts.min_loop_len,opts.max_loop_len,opts.max_len,opts.segment_num);
            bool dssp_pass;

            if(opts.motif_ss == "") {
                dssp_pass = random_bundle_dssp(opts.len, opts.num_helix, opts.min_helix_len, opts.max_helix_len, opts.min_loop_len, opts.max_loop_len, dssp_prefix, dssp);
                if(opts.ncaa_pos == -1) {
                    insert_pos = utils::random_int(1, opts.len);
                } else{
                    insert_pos = opts.ncaa_pos;
                }
            } else if(opts.motif_ss == "AZOF_TER_HELIX") {
                dssp_pass = random_bundle_dssp_AzoF(opts.len, opts.num_helix, opts.min_helix_len, opts.max_helix_len, opts.min_loop_len, opts.max_loop_len, dssp_prefix, dssp, insert_pos);
            } else if (opts.motif_ss == "FROM_DB") {
                dssp_pass = true;
                std::cout << "Dssp config from database: " << sequence << std::endl;
                ss_db.fetch_dssp_config(dssp_prefix, dssp, insert_pos);
            }
            else {
                dssp_pass = random_motif_bundle_dssp(opts.len,opts.motif_ss,false,opts.num_helix, opts.min_helix_len, opts.max_helix_len, opts.min_loop_len, opts.max_loop_len, dssp_prefix, dssp,insert_pos);
            }
            
            
            std::cout << "dssp : " << dssp  << std::endl;
            for( Size ii=0; ii<dssp.length(); ++ii ) {
                if(dssp.at(ii) == 'L')
                    sequence += "G";
                else
                    sequence += "V";
            }

            // for FDI insert_pos += 6;
            // mutate that residue to X
            sequence[insert_pos-1] = 'X';

            std::cout << "sequence: " << sequence << std::endl;
            if (dssp_pass) {
                pose.set_dssp(dssp);
                pose.set_sequence(sequence);
                // load the ncaa    
                pose.load_pdb(opts.ncaa, insert_pos, true, true, true, true);

                hallucinate_dssp_prefix = dssp_prefix;

            } else {
                itry--;
                continue;
            }
        }

        pose.reset_coords();
        
        std::cout << "mcmc ... " << itry << std::endl;
        // not perturbing the jumps between chains
        bool success = sampling::ncaa_protein_mcmc(pose, sfxn, frag_mover, sidechain_neighbors, ligand_neighbors, opts);

        if( !success ) continue;

        Real scs = sidechain_neighbors.compute(pose);
        Real total_sc = sfxn.score(pose) / pose.size();

        std::stringstream name;

        if(opts.output_dir != "")
            name << opts.output_dir << "/";
        // prefix, itry, time
        if(opts.output_prefix != "") name << opts.output_prefix << "_";
        if( hallucinate_dssp_prefix != "" ) name << hallucinate_dssp_prefix << "_";

        // the ncaa pos
        name << "_NCAA_" << insert_pos;
        
        name << "_" << std::chrono::high_resolution_clock::now().time_since_epoch().count();
        name << "_SN_" << std::setprecision(5) << scs;
        if(opts.ligand_pdb != "") {
            Real ln = ligand_neighbors.compute(pose);
            name << "_LN_" << ln;
        }
        name << "_TotalSc_" << total_sc;

        name << "_RPX_" << std::static_pointer_cast<RpxScoreMethod>(sfxn.get_method("rpx"))->score(pose)/pose.size();


        name << "_Clash_" << std::static_pointer_cast<ClashScoreMethod>(sfxn.get_method("clash"))->score(pose)/pose.size();

        name << "_RIF_" << opts.rif_weight * std::static_pointer_cast<scoring::RifScoreMethod>(sfxn.get_method("target_rif"))->score(pose);

        name << "_NCAAClash_" << + opts.context_clash_weight * std::static_pointer_cast<scoring::VoxelClashScoreMethod>(sfxn.get_method("target_context_clash"))->score(pose);

        name << ".pdb"<<(opts.gzip?".gz":"");

        std::cout << "Dump model: " << name.str() << std::endl;

        std::string original_seq = "";
        if(opts.dump_polyG) {
            original_seq = pose.sequence();
            std::string polyG_str(pose.size(), 'G');
            polyG_str[insert_pos-1] = 'X';
            pose.set_sequence(polyG_str);
        }

        pose.dump_pdb(name.str(),false,opts.gzip);

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

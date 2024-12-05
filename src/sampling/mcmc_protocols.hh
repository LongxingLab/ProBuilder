#ifndef INCLUDED_sampling_mcmc_protocols_hh
#define INCLUDED_sampling_mcmc_protocols_hh

#include "apps/args.hh"
#include "basic/types.hh"
#include "scene/Pose.hh"
#include "scene/SmallMoleculeInteractionField.hh"
#include "scene/InterfaceMetal.hh"
#include "scene/InterfaceZinc.hh"
#include "scene/InterfaceCalcium.hh"
#include "sampling/FragMover.hh"
#include "scoring/ScoreFunction.hh"
#include "utils/random_util.hh"
#include "metric/LigandNeighbors.hh"
#include "metric/SidechainNeighbors.hh"

#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>

namespace sampling {

using namespace basic;

inline
bool
pass_metropolis(
    Real const & temperature,
    Real const & deltaE,
    Real const & random_uniform
)
{
    if ( deltaE < 0 ) {
        return true;
    } else { //evaluate prob of substitution
        Real lnprob = deltaE / temperature;
        if ( lnprob < 10.0 ) {
            Real probability = std::exp(-lnprob);
            if ( probability > random_uniform ) return true;
        }
    }
    return false;
}


bool mcmc(scene::Pose & pose, scoring::ScoreFunction & sfxn, sampling::FragMover & frag_mover, Options const & opts);
bool mcmc_motif_pair(scene::Pose & pose, scoring::ScoreFunction & sfxn, sampling::FragMover & frag_mover, metric::SidechainNeighbors & sidechain_neighbors, Options const & opts);
bool mcmc_refine(scene::Pose & pose, scoring::ScoreFunction & sfxn, sampling::FragMover & frag_mover, Options const & opts);

bool mcmc_extension(scene::Pose & pose, scoring::ScoreFunction & sfxn, sampling::FragMover & frag_mover, metric::SidechainNeighbors & sidechain_neighbors, Options const & opts);
// vaccine
bool mcmc_vaccine(scene::Pose & pose, scoring::ScoreFunction & sfxn, sampling::FragMover & frag_mover, std::string & fake_dssp, Options const & opts);

// for ligand binding?
bool mcmc_jump(scene::Pose & pose, scoring::ScoreFunction & sfxn, Options const & opts, EigenXform const & extra_ligand_xform=EigenXform::Identity());

// kinda like rpxdock, but using mcmc for sampling
bool mcmc_pH_assembly(scene::Pose & pose, scoring::ScoreFunction & sfxn, Options const & opts);

// light
bool mcmc_light_heterodimer(scene::Pose & pose, scoring::ScoreFunction & sfxn, metric::SidechainNeighbors & sidechain_neighbors, Options const & opts);
bool mcmc_light_assembly(scene::Pose & pose, scoring::ScoreFunction & sfxn, metric::SidechainNeighbors & sidechain_neighbors, Options const & opts);

bool interface_metal_heterodimer_mcmc(scene::Pose & pose, scene::InterfaceMetalOP & interface_metal_p, scoring::ScoreFunction & sfxn, Options const & opts);
bool interface_metal_homooligomer_mcmc(scene::Pose & pose, scene::InterfaceMetalOP & interface_metal_p, scoring::ScoreFunction & sfxn, metric::SidechainNeighbors & sidechain_neighbors, Options const & opts);
bool interface_metal_fiber_mcmc(scene::Pose & pose, scene::InterfaceMetalOP & interface_metal_p1, scene::InterfaceMetalOP & interface_metal_p2, Size & metal1_chain_idx, Size & metal2_chain_idx, scoring::ScoreFunction & sfxn, metric::SidechainNeighbors & sidechain_neighbors, Options const & opts);


bool SMIF_homooligomer_mcmc(scene::Pose & pose, scene::SmallMoleculeInteractionField & smif, scoring::ScoreFunction & sfxn, metric::SidechainNeighbors & sidechain_neighbors, Options const & opts);


bool metal_coordination_mcmc(scene::Pose & pose, scoring::ScoreFunction & sfxn, sampling::FragMover & frag_mover, metric::SidechainNeighbors & sidechain_neighbors,
                                metric::LigandNeighbors & ligand_neighbors, Options const & opts);

bool small_molecule_binder_mcmc(scene::Pose & pose, scoring::ScoreFunction & sfxn, sampling::FragMover & frag_mover, metric::SidechainNeighbors & sidechain_neighbors,
                                metric::LigandNeighbors & ligand_neighbors, Options const & opts);
bool ncaa_protein_mcmc(scene::Pose & pose, scoring::ScoreFunction & sfxn, sampling::FragMover & frag_mover, metric::SidechainNeighbors & sidechain_neighbors,
                                metric::LigandNeighbors & ligand_neighbors, Options const & opts);
bool GFP_protein_mcmc(scene::Pose & pose, scoring::ScoreFunction & sfxn, sampling::FragMover & frag_mover, metric::SidechainNeighbors & sidechain_neighbors,
                                metric::LigandNeighbors & ligand_neighbors, Options const & opts);
bool ncaa_protein_symm_mcmc(scene::Pose & pose, scoring::ScoreFunction & sfxn, sampling::FragMover & frag_mover, metric::SidechainNeighbors & sidechain_neighbors,
                                metric::LigandNeighbors & ligand_neighbors, Options const & opts);

bool mcmc_repeat_peptide(scene::Pose & pose, scoring::ScoreFunction & sfxn, sampling::FragMover & frag_mover, Options const & opts);
bool mcmc_fold_on_target(scene::Pose & pose, scoring::ScoreFunction & sfxn, sampling::FragMover & frag_mover, Options const & opts );
}

#endif

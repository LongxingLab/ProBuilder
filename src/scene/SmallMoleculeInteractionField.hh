#ifndef INCLUDED_scene_SmallMoleculeInteractionField_hh
#define INCLUDED_scene_SmallMoleculeInteractionField_hh

#include "basic/types.hh"

#include <set>
#include <vector>
#include <memory>

namespace scene {

using namespace basic;

struct SmallMolecule_distance {
    Real distance;
    Size idx;
    Size jdx;
};


class SmallMoleculeInteractionField
{
    public:
        SmallMoleculeInteractionField();
        ~SmallMoleculeInteractionField();

        // small molecule positions in json format
        void load_SMIF_configs(std::string json_fname1, std::string json_fname2);
        void exclude_small_molecules(std::string const & small_molecule_names);

        void update_SMIF_position(EigenXform const & ref_xform_chain1, EigenXform const & ref_xform_chain2);
        void update_chain1_SMIF_position(EigenXform const & ref_xform);
        void update_chain2_SMIF_position(EigenXform const & ref_xform);

        Real small_molecule_radius();
        void set_small_molecule_radius(Real radius);
        Real ligand_neighbors_cutoff();
        void set_ligand_neighbors_cutoff(Real cutoff);

        Real topN_pair_distance(Size topN); // the sum distance of top N metal matches
        std::vector<std::string> report_topN_pairs(Size topN);

        void rollback(bool chain1_only=false);
        void snapshot(bool chain1_only=false);

    protected:

        void parse_json(std::string const & json_fname, std::vector<std::string> & cluster_names, std::vector<EigenXform> & xforms, Size & num_metals, std::vector<bool> & intra_chain_cluster_compatibility, std::vector<Real> & ligand_neighbors);
        Real xform_magnitude(EigenXform const & x);
        std::vector<SmallMolecule_distance> collect_topN_pairs(Size topN);

        //
        void load_chain1_SMIF_configs(std::string json_fname);
        void load_chain2_SMIF_configs(std::string json_fname);
        void compute_pairwise_distance();

    protected:
        std::vector<EigenXform> _relative_xforms_chain1;
        std::vector<EigenXform> _current_xforms_chain1;
        std::vector<EigenXform> _current_xforms_chain1_snapshot;


        std::vector<EigenXform> _relative_xforms_chain2;
        std::vector<EigenXform> _current_xforms_chain2;        
        std::vector<EigenXform> _current_xforms_chain2_snapshot;

        std::vector<std::string> _cluster_names_chain1;
        std::vector<std::string> _cluster_names_chain2;

        std::vector<bool> _intra_chain_cluster_compatibility_chain1;
        std::vector<bool> _intra_chain_cluster_compatibility_chain2;

        std::vector<bool> _inter_chain_cluster_match;

        // ligand neighbors
        std::vector<Real> _ligand_neighbors_chain1;
        std::vector<Real> _ligand_neighbors_chain2;

        // exclude metal clusters with name ...
        std::set<std::string> _excluded_small_molecules;

        Real _ligand_neighbors_cutoff;
        Size _num_small_molecules_chain1;
        Size _num_small_molecules_chain2;
        Real _small_molecule_radius;

        std::vector<SmallMolecule_distance> _distances;
        std::vector<SmallMolecule_distance> _distances_snapshot;

        bool _distance_updated;
        bool _distance_updated_snapshot;
};

typedef std::shared_ptr<SmallMoleculeInteractionField> SmallMoleculeInteractionFieldOP;


}

#endif

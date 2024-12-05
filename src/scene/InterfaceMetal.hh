#ifndef INCLUDED_scene_InterfaceMetal_hh
#define INCLUDED_scene_InterfaceMetal_hh

#include "basic/types.hh"

#include <set>
#include <vector>
#include <memory>

namespace scene {

using namespace basic;

struct Metal_distance {
    Real distance;
    Size idx;
    Size jdx;
};


class InterfaceMetal
{
    public:
        InterfaceMetal();
        ~InterfaceMetal();

        // metal positions in json format
        void load_metal_configs(std::string json_fname1, std::string json_fname2);
        void exclude_metals(std::string const & metal_names);

        void update_metal_position(EigenXform const & ref_xform_chain1, EigenXform const & ref_xform_chain2);
        void update_chain1_metal_position(EigenXform const & ref_xform);
        void update_chain2_metal_position(EigenXform const & ref_xform);

        Real metal_radius();
        void set_metal_radius(Real radius);

        Real topN_pair_distance(Size topN); // the sum distance of top N metal matches
        std::vector<std::string> report_topN_pairs(Size topN);

        void rollback(bool chain1_only=false);
        void snapshot(bool chain1_only=false);

    protected:

        void parse_json(std::string const & json_fname, std::vector<std::string> & cluster_names, std::vector<EigenXform> & xforms, Size & num_metals, std::vector<bool> & intra_chain_cluster_compatibility, std::vector<Size> & cluster_tokens, std::vector<std::set<Size> > & matching_clusters);
        Real xform_magnitude(EigenXform const & x);
        std::vector<Metal_distance> collect_topN_pairs(Size topN);

        // every new metal needs to reimplement these three functions
        virtual void load_chain1_metal_configs(std::string json_fname) = 0;
        virtual void load_chain2_metal_configs(std::string json_fname) = 0;
        virtual void compute_pairwise_distance() = 0;

    protected:
        std::vector<EigenXform> _relative_xforms_chain1;
        std::vector<EigenXform> _current_xforms_chain1;
        std::vector<EigenXform> _current_xforms_chain1_snapshot;


        std::vector<EigenXform> _relative_xforms_chain2;
        std::vector<EigenXform> _current_xforms_chain2;        
        std::vector<EigenXform> _relative_xforms_chain2_flipped;
        std::vector<EigenXform> _current_xforms_chain2_flipped;
        std::vector<EigenXform> _current_xforms_chain2_snapshot;
        std::vector<EigenXform> _current_xforms_chain2_flipped_snapshot;

        std::vector<std::string> _cluster_names_chain1;
        std::vector<std::string> _cluster_names_chain2;

        std::vector<bool> _intra_chain_cluster_compatibility_chain1;
        std::vector<bool> _intra_chain_cluster_compatibility_chain2;

        // exclude metal clusters with name ...
        std::set<std::string> _excluded_metals;

        Size _num_metals_chain1;
        Size _num_metals_chain2;
        Real _metal_radius;

        // temporary args
        std::vector<Size> _cluster_tokens_chain1;
        std::vector<Size> _cluster_tokens_chain2;
        std::vector<std::set<Size> > _matching_clusters_chain1;
        std::vector<std::set<Size> > _matching_clusters_chain2;

        std::vector<bool> _inter_chain_cluster_match;

        std::vector<Metal_distance> _distances;
        std::vector<Metal_distance> _distances_snapshot;

        bool _distance_updated;
        bool _distance_updated_snapshot;
};

typedef std::shared_ptr<InterfaceMetal> InterfaceMetalOP;


}

#endif

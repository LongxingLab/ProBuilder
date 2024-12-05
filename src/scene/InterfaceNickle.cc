#include <algorithm>
#include <numeric>
#include <fstream>
#include <cmath>

#include "utils/utils.hh"
#include "basic/types.hh"
#include "scene/InterfaceNickle.hh"
#include "utils/random_util.hh"

// external/json
#include <json/json.hpp>

namespace scene {

using namespace basic;

InterfaceNickle::InterfaceNickle() : InterfaceMetal() {}
InterfaceNickle::~InterfaceNickle() {}

void InterfaceNickle::load_chain1_metal_configs(std::string json_fname)
{
    // parse the first json file
    parse_json(json_fname, _cluster_names_chain1, _relative_xforms_chain1, _num_metals_chain1, _intra_chain_cluster_compatibility_chain1, _cluster_tokens_chain1, _matching_clusters_chain1);
    _current_xforms_chain1.resize(_num_metals_chain1);
    _current_xforms_chain1_snapshot.resize(_num_metals_chain1);
}

void InterfaceNickle::load_chain2_metal_configs(std::string json_fname)
{
    // Actually the chain2 xform is not used
    // not need to load the _relative_xforms_chain2
    // for compatibility, I will just leave it here.
    parse_json(json_fname, _cluster_names_chain2, _relative_xforms_chain2, _num_metals_chain2, _intra_chain_cluster_compatibility_chain2, _cluster_tokens_chain2, _matching_clusters_chain2);
    _current_xforms_chain2.resize(_num_metals_chain2);
    _current_xforms_chain2_snapshot.resize(_num_metals_chain2);
}

void InterfaceNickle::compute_pairwise_distance()
{
    if(!_distance_updated) {
        // TODO: can I directly store the inverse transform? Probabely

        for(Size idx=0; idx<_num_metals_chain1; ++idx) {
            for(Size jdx=0; jdx<_num_metals_chain2; ++jdx) {
                // only fill in diagonal
                // other positions the dist is 9e9
                if(idx==jdx && _inter_chain_cluster_match[idx*_num_metals_chain2+jdx]) {

                    Real dist(9e9);

                    Vec y_axis = _current_xforms_chain1[idx].rotation().col(1);
                    Vec z_axis = _current_xforms_chain1[idx].rotation().col(2);
                    Vec metal  = _current_xforms_chain1[idx].translation();

                    // distance between two lines
                    Vec axis1 = (z_axis + y_axis).normalized();
                    Vec axis2 = (z_axis - y_axis).normalized();

                    // angle
                    Real cos_theta1 = std::fabs(axis1[2]);
                    Real err_rot1 = std::sqrt( std::max( 0.0, 1.0 - cos_theta1*cos_theta1 ) ) * _metal_radius;

                    Real cos_theta2 = std::fabs(axis2[2]);
                    Real err_rot2 = std::sqrt( std::max( 0.0, 1.0 - cos_theta2*cos_theta2 ) ) * _metal_radius;

                    Real err_rot = std::min(err_rot1, err_rot2);

                    Real err_trans2 = metal[0]*metal[0] + metal[1]*metal[1];

                    dist = std::sqrt( err_trans2 + err_rot*err_rot );

                    _distances[idx*_num_metals_chain2+jdx].distance = dist;
                } else {

                    _distances[idx*_num_metals_chain2+jdx].distance = 9e9;
                }
                _distances[idx*_num_metals_chain2+jdx].idx = idx;
                _distances[idx*_num_metals_chain2+jdx].jdx = jdx;
            }
        }

        _distance_updated = true;
    }
}

}

#include <algorithm>
#include <numeric>
#include <fstream>
#include <cmath>

#include "utils/utils.hh"
#include "basic/types.hh"
#include "scene/InterfaceCalcium.hh"
#include "utils/random_util.hh"

// external/json
#include <json/json.hpp>

namespace scene {

using namespace basic;

InterfaceCalcium::InterfaceCalcium() : InterfaceMetal() {}
InterfaceCalcium::~InterfaceCalcium() {}

void InterfaceCalcium::load_chain1_metal_configs(std::string json_fname)
{
    // parse the first json file
    parse_json(json_fname, _cluster_names_chain1, _relative_xforms_chain1, _num_metals_chain1, _intra_chain_cluster_compatibility_chain1, _cluster_tokens_chain1, _matching_clusters_chain1);
    _current_xforms_chain1.resize(_num_metals_chain1);
    _current_xforms_chain1_snapshot.resize(_num_metals_chain1);
}

void InterfaceCalcium::load_chain2_metal_configs(std::string json_fname)
{
    // parse chain2 metal config file
    std::vector<EigenXform> relative_xforms_chain2_origin;
    parse_json(json_fname, _cluster_names_chain2, relative_xforms_chain2_origin, _num_metals_chain2, _intra_chain_cluster_compatibility_chain2, _cluster_tokens_chain2, _matching_clusters_chain2);
    // the relative_xforms of chain2
    for(Size idx=0; idx<_num_metals_chain2; ++idx) {
        //
        EigenXform const & tempx = relative_xforms_chain2_origin[idx];
        Vec x_axis = tempx.rotation().col(0);
        Vec z_axis = tempx.rotation().col(2);
        EigenXform tempy, tempz;
        Mat m = AngleAxis( 3.14159265359, z_axis ) * tempx.rotation();
        for(Size i = 0; i < 9; ++i) tempy.data()[i] = m.data()[i];

        tempy.translation() = tempx.translation();

        _relative_xforms_chain2.push_back(tempy);
    }
    _current_xforms_chain2.resize(_num_metals_chain2);
    _current_xforms_chain2_snapshot.resize(_num_metals_chain2);
}

void InterfaceCalcium::compute_pairwise_distance()
{
    if(!_distance_updated) {
        // TODO: can I directly store the inverse transform? Probabely

        for(Size idx=0; idx<_num_metals_chain1; ++idx) {
            EigenXform y = _current_xforms_chain1[idx].inverse(Eigen::Isometry);
            for(Size jdx=0; jdx<_num_metals_chain2; ++jdx) {
                if(_inter_chain_cluster_match[idx*_num_metals_chain2+jdx]) {

                    Real dist(9e9);

                    Vec t = _current_xforms_chain1[idx].translation() - _current_xforms_chain2[jdx].translation();
                    Real d2 = t.squaredNorm();

                    // TODO: best cutoff??
                    if( d2 > 64 ) {
                        dist = std::sqrt(d2) + _metal_radius;
                    } else {
                        dist = xform_magnitude(y*_current_xforms_chain2[jdx]);
                    }

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

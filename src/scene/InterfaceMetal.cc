#include <algorithm>
#include <numeric>
#include <fstream>
#include <cmath>

#include "utils/utils.hh"
#include "basic/types.hh"
#include "scene/InterfaceMetal.hh"
#include "utils/random_util.hh"

// external/json
#include <json/json.hpp>

namespace scene {

using namespace basic;

InterfaceMetal::InterfaceMetal() : 
    _num_metals_chain1(0),
    _num_metals_chain2(0),
    _metal_radius(4.0),
    _distance_updated(false),
    _distance_updated_snapshot(false)
{}

InterfaceMetal::~InterfaceMetal() {}

// example of one entry
/*
"Zn_Cluster0028": {
        "avg_dist": 0.0035325423814356327,
        "max_dist": 0.007065084762871265,
        "center_xform": "0.5899181483454792_0.7125986009982159_-0.37973650352320204_0.6464643676878301_-0.1350218608472628_0.7509014039161976_0.4838185605809426_-0.6884564844263783_-0.5403214501467222_5.057025581907076_-11.311456195554017_-1.736",
        "center_match": "new_UM_104_H30D26_17_bcov_4helix_11072_ZNX_HD_1.pdb",
        "descriptions": [
            "new_UM_104_H30D26_17_bcov_4helix_11072_ZNX_HD_1.pdb",
            "new_UM_104_H30D26_9_bcov_4helix_11072_ZNX_HD_1.pdb"
        ],
        "residue_index": [
            26,
            30
        ],
        "cluster_token": 1,
        "matching_clusters": [
            2
        ]
*/
// The following information must be included in the json file
// cluster_name : {
//      center_xform: the transform of this metal
//      residue_index: the residue indexes of for this metal, to be used to check whether two clusters on the same chain are compatible or not
//      cluster_token: the identifier of the cluster, to be used to control pairing => HH--HD, HD--HD, HH--HH
//      matching_cluster: to be used to control pairing along with cluster_token
//    }


// The convention of the xform
// atom1, Zn, atom2
// vec1 = atom1 - Zn
// vec2 = atom2 - Zn
// x = vec1 + vec2; x.normalize()
// z = vec1.cross(vec2)
// y = z.cross(x)

// convert the matrix to str
// def xform_to_str(xform):
//    return '_'.join([str(ii) for ii in xform[:3,:].flatten(order='F')])

void InterfaceMetal::parse_json(std::string const & json_fname, std::vector<std::string> & cluster_names, std::vector<EigenXform> & xforms, Size & num_metals, std::vector<bool> & intra_chain_cluster_compatibility, std::vector<Size> & cluster_tokens, std::vector<std::set<Size> > & matching_clusters)
{
    std::ifstream jfile( json_fname );
    if ( ! jfile ) {
        std::cout << "Error: failed to open file " << json_fname << " " << std::endl;
        exit (0);
    }
    std::stringstream buffer;
    buffer << jfile.rdbuf();
    nlohmann::json  j = nlohmann::json::parse( buffer.str() );
    jfile.close();

    std::vector<std::set<Size> > cluster_residue_index;

    cluster_tokens.clear();
    matching_clusters.clear();

    for(auto metal_cluster : j.items()){
        std::string xform_str = metal_cluster.value()["center_xform"];
        std::string cluster_name = metal_cluster.key();

        // check whether the cluster has been excluded or not
        // 
        if( _excluded_metals.find(cluster_name) != _excluded_metals.end() ) {
            std::cout << "       ====> excludeing metal cluster " << cluster_name << " while parsing the interface metal config file." << std::endl;
            continue;
        }

        std::vector<std::string> splited_str;
        std::string delimiters = "_";
        std::string::size_type lastPos = xform_str.find_first_not_of(delimiters,0);
        std::string::size_type pos = xform_str.find_first_of(delimiters, lastPos);
        while( std::string::npos != pos || std::string::npos != lastPos ) {
            splited_str.push_back(xform_str.substr(lastPos, pos-lastPos));
            lastPos = xform_str.find_first_not_of(delimiters, pos);
            pos = xform_str.find_first_of(delimiters, lastPos);
        }

        if(splited_str.size() != 12) {
            std::cout << "Bad Xform for interface metal cluster xform!!!!!!!!!!!!" << std::endl;
            exit(0);
        }

        EigenXform tempx(EigenXform::Identity());
        for(Size idx=0; idx<12; ++idx) {
            tempx.data()[idx] = std::stof(splited_str[idx]);
        }

        xforms.push_back(tempx);       
        cluster_names.push_back(cluster_name);
        ++num_metals;

        std::set<Size> res_idx;
        for(Size idx : metal_cluster.value()["residue_index"]) {
            res_idx.insert(idx);
        }
        cluster_residue_index.push_back(res_idx);

        cluster_tokens.push_back(metal_cluster.value()["cluster_token"]);
        std::set<Size> m_tokens;
        for(Size i_token : metal_cluster.value()["matching_clusters"]) {
            m_tokens.insert(i_token);
        }
        matching_clusters.push_back(m_tokens);

    }

    // two metal clusters on the same chain compatible??? 
    intra_chain_cluster_compatibility.resize(num_metals*num_metals);
    for(Size idx=0; idx<num_metals; ++idx) {
        for(Size jdx=0; jdx<num_metals; ++jdx) {
            // check 
            bool flag(true);
            std::set<Size> & s1 = cluster_residue_index[idx];
            std::set<Size> & s2 = cluster_residue_index[jdx];
            std::set<Size> r;
            std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
                 std::inserter(r, r.begin()));
            if( r.size() != 0 ) flag = false;
            intra_chain_cluster_compatibility[idx*num_metals+jdx] = flag;
        }
    }
}

void InterfaceMetal::load_metal_configs(std::string json_fname1, std::string json_fname2)
{

    load_chain1_metal_configs(json_fname1);
    load_chain2_metal_configs(json_fname2);

    // do cluster match???
    _inter_chain_cluster_match.resize(_num_metals_chain1*_num_metals_chain2);

    // initialize distances
    _distances.resize(_num_metals_chain1*_num_metals_chain2);
    _distances_snapshot.resize(_num_metals_chain1*_num_metals_chain2);

    for(Size idx=0; idx < _num_metals_chain1; ++idx) {
        for(Size jdx=0; jdx < _num_metals_chain2; ++jdx) {
            _distances[idx*_num_metals_chain2+jdx].idx = idx;
            _distances[idx*_num_metals_chain2+jdx].jdx = jdx;
            _distances[idx*_num_metals_chain2+jdx].distance = 9e9;

            _distances_snapshot[idx*_num_metals_chain2+jdx].idx = idx;
            _distances_snapshot[idx*_num_metals_chain2+jdx].jdx = jdx;
            _distances_snapshot[idx*_num_metals_chain2+jdx].distance = 9e9;


            if( _matching_clusters_chain2[jdx].find( _cluster_tokens_chain1[idx] ) != _matching_clusters_chain2[jdx].end()
                && _matching_clusters_chain1[idx].find( _cluster_tokens_chain2[jdx] ) != _matching_clusters_chain1[jdx].end() ) {
                _inter_chain_cluster_match[idx*_num_metals_chain2+jdx] = true;
            } else {

                _inter_chain_cluster_match[idx*_num_metals_chain2+jdx] = false;
            }
        }
    }

}

void InterfaceMetal::exclude_metals(std::string const & metal_names)
{
    std::string delimiters = ":";
    std::string::size_type lastPos = metal_names.find_first_not_of(delimiters,0);
    std::string::size_type pos = metal_names.find_first_of(delimiters, lastPos);
    while( std::string::npos != pos || std::string::npos != lastPos ) {
        _excluded_metals.insert(metal_names.substr(lastPos, pos-lastPos));
        lastPos = metal_names.find_first_not_of(delimiters, pos);
        pos = metal_names.find_first_of(delimiters, lastPos);
    }

}

void InterfaceMetal::update_metal_position(EigenXform const & ref_xform_chain1, EigenXform const & ref_xform_chain2)
{
    update_chain1_metal_position(ref_xform_chain1);
    update_chain2_metal_position(ref_xform_chain2);
    _distance_updated = false;
}
void InterfaceMetal::update_chain1_metal_position(EigenXform const & ref_xform)
{
    for(Size idx=0; idx<_num_metals_chain1; ++idx) {
        _current_xforms_chain1[idx] = ref_xform * _relative_xforms_chain1[idx];
    }
    _distance_updated = false;
}
void InterfaceMetal::update_chain2_metal_position(EigenXform const & ref_xform)
{
    for(Size idx=0; idx<_num_metals_chain2; ++idx) {
        _current_xforms_chain2[idx] = ref_xform * _relative_xforms_chain2[idx];
    }

    // might be zero
    for(Size idx=0; idx<_current_xforms_chain2_flipped.size(); ++idx) {
        _current_xforms_chain2_flipped[idx] = ref_xform * _relative_xforms_chain2_flipped[idx];
    }
    _distance_updated = false;
}


// convert the relative xform into a real number
Real InterfaceMetal::xform_magnitude(
    EigenXform const & x
){
    Real err_trans2 = x.translation().squaredNorm();
    Real cos_theta = (x.rotation().trace()-1.0)/2.0;

    Real err_rot = std::sqrt( std::max( 0.0, 1.0 - cos_theta*cos_theta ) ) * _metal_radius;
    if( cos_theta < 0 ) err_rot = _metal_radius;
    Real err = std::sqrt( err_trans2 + err_rot*err_rot );
    return err;
}


struct customLess{
    bool operator()(Metal_distance a, Metal_distance b) const { return a.distance < b.distance; }
};

std::vector<Metal_distance> InterfaceMetal::collect_topN_pairs(Size topN)
{
    compute_pairwise_distance();

    std::vector<Metal_distance> results;

    if(topN == 1) {
        std::vector<Metal_distance>::iterator result = std::min_element(_distances.begin(), _distances.end(), customLess());
        results.push_back(*result);
        return results;
    }


    // std::nth_element(_distances.begin(), _distances.begin()+200, _distances.end(), customLess());
    std::partial_sort(_distances.begin(), _distances.begin()+100, _distances.end(), customLess());

    // std::nth_element(_distances.begin(), _distances.begin()+topN, _distances.end(), customLess());
    // std::sort(_distances.begin(), _distances.end(), customLess());

    Size counter(0);
    for(auto p : _distances) {
        bool is_conflict(false);
        for(auto q : results) {
            if (p.idx == q.idx || p.jdx == q.jdx) {
                is_conflict = true;
                break;
            }
            if(_intra_chain_cluster_compatibility_chain1[p.idx*_num_metals_chain1+q.idx] == false || 
               _intra_chain_cluster_compatibility_chain2[p.jdx*_num_metals_chain2+q.jdx] == false) {
                is_conflict = true;
                break;
            }
        }
        if( !is_conflict ) {
            results.push_back(p);
            counter++;
        }

        if(counter>=topN) break;
    }

    return results;
}

Real InterfaceMetal::topN_pair_distance(Size topN)
{
    std::vector<Metal_distance> pairs = collect_topN_pairs(topN);

    if(pairs.size() < topN) return 9e9; 

    Real val(0.0);
    for(auto p : pairs)
        val += p.distance;

    return val;
}

std::vector<std::string> InterfaceMetal::report_topN_pairs(Size topN)
{

    std::vector<Metal_distance> pairs = collect_topN_pairs(topN); 
    std::vector<std::string> str_results;
    
    // safe guard
    if(pairs.size() < topN) {
        str_results.push_back("HOLY_SHIT");
        return str_results;
    }

    for(auto p : pairs)
        str_results.push_back(_cluster_names_chain1[p.idx] + "_" + _cluster_names_chain2[p.jdx]);
    return str_results;
}

void InterfaceMetal::rollback(bool chain1_only)
{
    //
    std::copy(_current_xforms_chain1_snapshot.begin(), _current_xforms_chain1_snapshot.end(), _current_xforms_chain1.begin());
    if(!chain1_only) {
        std::copy(_current_xforms_chain2_snapshot.begin(), _current_xforms_chain2_snapshot.end(), _current_xforms_chain2.begin());
        if(_current_xforms_chain2_flipped.size() != 0) {
            std::copy(_current_xforms_chain2_flipped_snapshot.begin(), _current_xforms_chain2_flipped_snapshot.end(), _current_xforms_chain2_flipped.begin());
        }
    }
    std::copy(_distances_snapshot.begin(), _distances_snapshot.end(), _distances.begin());

    _distance_updated = _distance_updated_snapshot;
}

void InterfaceMetal::snapshot(bool chain1_only)
{
    //
    std::copy(_current_xforms_chain1.begin(), _current_xforms_chain1.end(), _current_xforms_chain1_snapshot.begin());
    if(!chain1_only) {
        std::copy(_current_xforms_chain2.begin(), _current_xforms_chain2.end(), _current_xforms_chain2_snapshot.begin());
        if(_current_xforms_chain2_flipped.size() != 0) {
            std::copy(_current_xforms_chain2_flipped.begin(), _current_xforms_chain2_flipped.end(), _current_xforms_chain2_flipped_snapshot.begin());
        }
    }
    std::copy(_distances.begin(), _distances.end(), _distances_snapshot.begin());

    _distance_updated_snapshot = _distance_updated;
}


Real InterfaceMetal::metal_radius()
{
    return _metal_radius;
}
void InterfaceMetal::set_metal_radius(Real radius)
{
    _metal_radius = radius;
}


}

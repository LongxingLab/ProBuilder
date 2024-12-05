#include <algorithm>
#include <numeric>
#include <fstream>
#include <cmath>

#include "utils/utils.hh"
#include "basic/types.hh"
#include "scene/SmallMoleculeInteractionField.hh"
#include "utils/random_util.hh"

// external/json
#include <json/json.hpp>

namespace scene {

using namespace basic;

SmallMoleculeInteractionField::SmallMoleculeInteractionField() : 
    _ligand_neighbors_cutoff(6.0),
    _num_small_molecules_chain1(0),
    _num_small_molecules_chain2(0),
    _small_molecule_radius(4.0),
    _distance_updated(false),
    _distance_updated_snapshot(false)
{}

SmallMoleculeInteractionField::~SmallMoleculeInteractionField() {}

// example of one entry
/*
"SMcluster_0001": {
    "avg_dist": 0.18374095857143402,
    "max_dist": 0.36748191714286804,
    "center_xform": "0.11107458996824077_0.12306858094794722_-0.9861625423057021_-0.46768406510004135_-0.8690841731877916_-0.1611344629990514_-0.8768888474116248_0.47911045104229455_-0.03897595396089208_14.193551535988862_6.883736983610457_24.593633119625867",
    "center_match": "WU65aa3H05095_000000419",
    "descriptions": [
        "WU65aa3H05095_000000419",
        "WU65aa3H05095_000000587"
    ]
}
*/


// convert the matrix to str
// def xform_to_str(xform):
//    return '_'.join([str(ii) for ii in xform[:3,:].flatten(order='F')])

void SmallMoleculeInteractionField::parse_json(std::string const & json_fname, std::vector<std::string> & cluster_names, std::vector<EigenXform> & xforms, Size & num_small_molecules, std::vector<bool> & intra_chain_cluster_compatibility, std::vector<Real> & ligand_neighbors)
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

    for(auto small_molecule_cluster : j.items()){
        std::string xform_str = small_molecule_cluster.value()["center_xform"];
        std::string cluster_name = small_molecule_cluster.key();

        // check whether the cluster has been excluded or not
        // 
        if( _excluded_small_molecules.find(cluster_name) != _excluded_small_molecules.end() ) {
            std::cout << "       ====> excludeing small molecule cluster " << cluster_name << " while parsing the SMIF config file." << std::endl;
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

        ligand_neighbors.push_back(small_molecule_cluster.value()["ligand_neighbors"]);

        xforms.push_back(tempx);       
        cluster_names.push_back(cluster_name);
        ++num_small_molecules;
    }

    // two small molecule clusters on the same chain compatible??? 
    intra_chain_cluster_compatibility.resize(num_small_molecules*num_small_molecules);
    for(Size idx=0; idx<num_small_molecules; ++idx) {
        for(Size jdx=0; jdx<num_small_molecules; ++jdx) {
            // check 
            bool flag(true);
            Vec t = xforms[idx].translation() - xforms[jdx].translation();            
            if( t.squaredNorm() <= _small_molecule_radius * _small_molecule_radius ) flag = false;         
            intra_chain_cluster_compatibility[idx*num_small_molecules+jdx] = flag;
        }
    }
}

void SmallMoleculeInteractionField::load_SMIF_configs(std::string json_fname1, std::string json_fname2)
{

    load_chain1_SMIF_configs(json_fname1);
    std::cout << "Number of small molecule clusters on chain1: " << _num_small_molecules_chain1 << std::endl;
    load_chain2_SMIF_configs(json_fname2);
    std::cout << "Number of small molecule clusters on chain2: " << _num_small_molecules_chain2 << std::endl;

    // initialize distances
    _distances.resize(_num_small_molecules_chain1*_num_small_molecules_chain2);
    _distances_snapshot.resize(_num_small_molecules_chain1*_num_small_molecules_chain2);

    for(Size idx=0; idx < _num_small_molecules_chain1; ++idx) {
        for(Size jdx=0; jdx < _num_small_molecules_chain2; ++jdx) {
            _distances[idx*_num_small_molecules_chain2+jdx].idx = idx;
            _distances[idx*_num_small_molecules_chain2+jdx].jdx = jdx;
            _distances[idx*_num_small_molecules_chain2+jdx].distance = 9e9;

            _distances_snapshot[idx*_num_small_molecules_chain2+jdx].idx = idx;
            _distances_snapshot[idx*_num_small_molecules_chain2+jdx].jdx = jdx;
            _distances_snapshot[idx*_num_small_molecules_chain2+jdx].distance = 9e9;
        }
    }


    Size match_pair_counter(0);
    Size not_match_pair_counter(0);
    _inter_chain_cluster_match.resize(_num_small_molecules_chain1*_num_small_molecules_chain2);
    for(Size idx=0; idx < _num_small_molecules_chain1; ++idx) {
        for(Size jdx=0; jdx < _num_small_molecules_chain2; ++jdx) {
            if( _ligand_neighbors_chain1[idx]+_ligand_neighbors_chain2[jdx] >= _ligand_neighbors_cutoff ) {
                _inter_chain_cluster_match[idx*_num_small_molecules_chain2+jdx] = true;
                match_pair_counter++;
            } else {
                _inter_chain_cluster_match[idx*_num_small_molecules_chain2+jdx] = false;
                not_match_pair_counter++;
            }
        }
    }
    std::cout << "Pairs match: " << match_pair_counter << "     ;    Pairs not match: " << not_match_pair_counter << std::endl;

}

void SmallMoleculeInteractionField::exclude_small_molecules(std::string const & small_molecule_names)
{
    std::string delimiters = ":";
    std::string::size_type lastPos = small_molecule_names.find_first_not_of(delimiters,0);
    std::string::size_type pos = small_molecule_names.find_first_of(delimiters, lastPos);
    while( std::string::npos != pos || std::string::npos != lastPos ) {
        _excluded_small_molecules.insert(small_molecule_names.substr(lastPos, pos-lastPos));
        lastPos = small_molecule_names.find_first_not_of(delimiters, pos);
        pos = small_molecule_names.find_first_of(delimiters, lastPos);
    }
}

void SmallMoleculeInteractionField::update_SMIF_position(EigenXform const & ref_xform_chain1, EigenXform const & ref_xform_chain2)
{
    update_chain1_SMIF_position(ref_xform_chain1);
    update_chain2_SMIF_position(ref_xform_chain2);
    _distance_updated = false;
}
void SmallMoleculeInteractionField::update_chain1_SMIF_position(EigenXform const & ref_xform)
{
    for(Size idx=0; idx<_num_small_molecules_chain1; ++idx) {
        _current_xforms_chain1[idx] = ref_xform * _relative_xforms_chain1[idx];
    }
    _distance_updated = false;
}
void SmallMoleculeInteractionField::update_chain2_SMIF_position(EigenXform const & ref_xform)
{
    for(Size idx=0; idx<_num_small_molecules_chain2; ++idx) {
        _current_xforms_chain2[idx] = ref_xform * _relative_xforms_chain2[idx];
    }
    _distance_updated = false;
}


// convert the relative xform into a real number
Real SmallMoleculeInteractionField::xform_magnitude(
    EigenXform const & x
){
    Real err_trans2 = x.translation().squaredNorm();
    Real cos_theta = (x.rotation().trace()-1.0)/2.0;

    Real err_rot = std::sqrt( std::max( 0.0, 1.0 - cos_theta*cos_theta ) ) * _small_molecule_radius;
    if( cos_theta < 0 ) err_rot = _small_molecule_radius;
    Real err = std::sqrt( err_trans2 + err_rot*err_rot );
    return err;
}


struct customLess{
    bool operator()(SmallMolecule_distance a, SmallMolecule_distance b) const { return a.distance < b.distance; }
};

std::vector<SmallMolecule_distance> SmallMoleculeInteractionField::collect_topN_pairs(Size topN)
{
    compute_pairwise_distance();
    std::vector<SmallMolecule_distance> results;

    if(topN == 1) {
        std::vector<SmallMolecule_distance>::iterator result = std::min_element(_distances.begin(), _distances.end(), customLess());
        results.push_back(*result);
        return results;
    }


    // std::nth_element(_distances.begin(), _distances.begin()+200, _distances.end(), customLess());
    std::partial_sort(_distances.begin(), _distances.begin()+100, _distances.end(), customLess());

    Size counter(0);
    for(auto p : _distances) {
        bool is_conflict(false);
        for(auto q : results) {
            if (p.idx == q.idx || p.jdx == q.jdx) {
                is_conflict = true;
                break;
            }
            if(_intra_chain_cluster_compatibility_chain1[p.idx*_num_small_molecules_chain1+q.idx] == false || 
               _intra_chain_cluster_compatibility_chain2[p.jdx*_num_small_molecules_chain2+q.jdx] == false) {
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

Real SmallMoleculeInteractionField::topN_pair_distance(Size topN)
{
    std::vector<SmallMolecule_distance> pairs = collect_topN_pairs(topN);

    if(pairs.size() < topN) return 9e9; 

    Real val(0.0);
    for(auto p : pairs)
        val += p.distance;

    return val;
}

std::vector<std::string> SmallMoleculeInteractionField::report_topN_pairs(Size topN)
{

    std::vector<SmallMolecule_distance> pairs = collect_topN_pairs(topN); 
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

void SmallMoleculeInteractionField::rollback(bool chain1_only)
{
    //
    std::copy(_current_xforms_chain1_snapshot.begin(), _current_xforms_chain1_snapshot.end(), _current_xforms_chain1.begin());
    if(!chain1_only) {
        std::copy(_current_xforms_chain2_snapshot.begin(), _current_xforms_chain2_snapshot.end(), _current_xforms_chain2.begin());
    }
    std::copy(_distances_snapshot.begin(), _distances_snapshot.end(), _distances.begin());

    _distance_updated = _distance_updated_snapshot;
}

void SmallMoleculeInteractionField::snapshot(bool chain1_only)
{
    //
    std::copy(_current_xforms_chain1.begin(), _current_xforms_chain1.end(), _current_xforms_chain1_snapshot.begin());
    if(!chain1_only) {
        std::copy(_current_xforms_chain2.begin(), _current_xforms_chain2.end(), _current_xforms_chain2_snapshot.begin());
    }
    std::copy(_distances.begin(), _distances.end(), _distances_snapshot.begin());

    _distance_updated_snapshot = _distance_updated;
}


Real SmallMoleculeInteractionField::small_molecule_radius()
{
    return _small_molecule_radius;
}
void SmallMoleculeInteractionField::set_small_molecule_radius(Real radius)
{
    _small_molecule_radius = radius;
}

Real SmallMoleculeInteractionField::ligand_neighbors_cutoff()
{
    return _ligand_neighbors_cutoff;
}
void SmallMoleculeInteractionField::set_ligand_neighbors_cutoff(Real cutoff)
{
    _ligand_neighbors_cutoff = cutoff;
}

void SmallMoleculeInteractionField::load_chain1_SMIF_configs(std::string json_fname)
{
    // parse the first json file
    parse_json(json_fname, _cluster_names_chain1, _relative_xforms_chain1, _num_small_molecules_chain1, _intra_chain_cluster_compatibility_chain1, _ligand_neighbors_chain1);
    _current_xforms_chain1.resize(_num_small_molecules_chain1);
    _current_xforms_chain1_snapshot.resize(_num_small_molecules_chain1);
}

void SmallMoleculeInteractionField::load_chain2_SMIF_configs(std::string json_fname)
{
    parse_json(json_fname, _cluster_names_chain2, _relative_xforms_chain2, _num_small_molecules_chain2, _intra_chain_cluster_compatibility_chain2, _ligand_neighbors_chain2);
    // the relative_xforms of chain2
    _current_xforms_chain2.resize(_num_small_molecules_chain2);
    _current_xforms_chain2_snapshot.resize(_num_small_molecules_chain2);
}

void SmallMoleculeInteractionField::compute_pairwise_distance()
{
    if(!_distance_updated) {
        // TODO: can I directly store the inverse transform? Probabely
        for(Size idx=0; idx<_num_small_molecules_chain1; ++idx) {
            EigenXform y = _current_xforms_chain1[idx].inverse(Eigen::Isometry);
            for(Size jdx=0; jdx<_num_small_molecules_chain2; ++jdx) {

                if(_inter_chain_cluster_match[idx*_num_small_molecules_chain2+jdx]) {
                    Real dist(9e9);

                    Vec t = _current_xforms_chain1[idx].translation() - _current_xforms_chain2[jdx].translation();
                    Real d2 = t.squaredNorm();

                    // TODO: best cutoff??
                    if( d2 > 16 ) {
                        dist = std::sqrt(d2) + _small_molecule_radius;
                    } else {
                        dist = xform_magnitude(y*_current_xforms_chain2[jdx]);
                    }

                    _distances[idx*_num_small_molecules_chain2+jdx].distance = dist;
                } else {
                    _distances[idx*_num_small_molecules_chain2+jdx].distance = 9e9;
                }
                _distances[idx*_num_small_molecules_chain2+jdx].idx = idx;
                _distances[idx*_num_small_molecules_chain2+jdx].jdx = jdx;
            }
        }
        _distance_updated = true;
    }
}


}

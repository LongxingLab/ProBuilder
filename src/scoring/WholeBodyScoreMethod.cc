#include "scoring/WholeBodyScoreMethod.hh"
#include "scoring/Energy.hh"
#include "utils/hash_util.hh"
#include "utils/utils.hh"
#include "utils/math_util.hh"

#include <gzip/gzstream.hh>
#include <cmath>
#include <string>
#include <fstream>
#include <Eigen/Geometry>
#include <parallel_hashmap/phmap_dump.h>
#include <map>
namespace scoring {

using namespace basic;
using namespace scene;

Real  RepeatScoreMethod::score(scene::Pose & pose) const{
    Real radius_score(0),rise_score(0),omega_score(0);
    Real radius,rise,omega;
    Vec tmp_vec;
    utils::get_repeat_parameters_from_stubs(pose,rise,radius,omega,tmp_vec,tmp_vec);
    omega = omega / 3.14159 * 180;
    radius_score +=  radius<_radius_control.first?pow(_radius_control.first-radius,2):(radius>_radius_control.second?pow(radius-_radius_control.second,2):0);
    rise_score +=  rise<_rise_control.first?pow(_rise_control.first-rise,2):(rise>_rise_control.second?pow(rise-_rise_control.second,2):0);
    omega_score +=  omega<_omega_control.first?pow(_omega_control.first-omega,2):(omega>_omega_control.second?pow(omega-_omega_control.second,2):0);
    
    return radius_score*_radius_multiplier+rise_score*_rise_multiplier+omega_score*_omega_multiplier;
    };
    

MotifPairRelativePosScoreMethod::MotifPairRelativePosScoreMethod(std::string motif1_pdb, std::string motif2_pdb):
    _motif1_pdb(motif1_pdb),
    _motif2_pdb(motif2_pdb),
    _motif1_chain_idx(1),
    _motif2_chain_idx(2),
    _motif1_insert_pos(-1),
    _motif2_insert_pos(-1),
    _motif1_representative_res_index(-1),
    _motif2_representative_res_index(-1),
    _motif_radius(10.0), // should this be correlated with the size of the motifs??
    _inv_relative_xform(EigenXform::Identity())
{
    //
    initialize();
}

MotifPairRelativePosScoreMethod::~MotifPairRelativePosScoreMethod() {}

void MotifPairRelativePosScoreMethod::initialize()
{
    EigenXform motif1_res_xform, motif2_res_xform;

    load_motif(_motif1_pdb, _motif1_representative_res_index, motif1_res_xform);
    load_motif(_motif2_pdb, _motif2_representative_res_index, motif2_res_xform);

    _inv_relative_xform = (motif1_res_xform.inverse(Eigen::Isometry) * motif2_res_xform).inverse(Eigen::Isometry);
}

void MotifPairRelativePosScoreMethod::load_motif(std::string fname, Size & rep_res_index, EigenXform & xform)
{
    std::vector<std::string> raw_pdb_lines;

    // load the pdb file into pdblines
    if   ( fname.substr(fname.size()-4,4)==".pdb" )         {  } //
    else if ( fname.substr(fname.size()-7,7)==".pdb.gz" )   { std::cout << "Does not support gz format, gunzip your motif file!" << std::endl; exit(0);}
    else                                                    { std::cout << "Unknown Format!" << std::endl; exit(0); }

    std::ifstream pdb_file(fname, std::ios_base::in);
    if(pdb_file.fail()) {
        std::cout << "Failed to load motif pdb: " << fname << std::endl;
        exit(0);
    }

    Size motif_len(0);
    std::vector<std::string> pdb_lines;
    for(std::string line;std::getline(pdb_file,line);) {
        if( line.substr(0,4) == "ATOM" && (line.substr(13,2)=="N " || line.substr(13,3)=="CA " || line.substr(13,2)=="C ") )
        {
            pdb_lines.push_back(line);
            if( line.substr(13,3)=="CA " ){
                ++motif_len;
            }
        }
    }
    pdb_file.close();
    // num * [N, CA, C]
    assert(pdb_lines.size() == 3*motif_len);

    // mid residue
    rep_res_index = Size((motif_len-1)/2)+1;
    Vec N(std::stof(pdb_lines[(rep_res_index-1)*3+0].substr(30,8)),std::stof(pdb_lines[(rep_res_index-1)*3+0].substr(38,8)),std::stof(pdb_lines[(rep_res_index-1)*3+0].substr(46,8)));
    Vec CA(std::stof(pdb_lines[(rep_res_index-1)*3+1].substr(30,8)),std::stof(pdb_lines[(rep_res_index-1)*3+1].substr(38,8)),std::stof(pdb_lines[(rep_res_index-1)*3+1].substr(46,8)));
    Vec C(std::stof(pdb_lines[(rep_res_index-1)*3+2].substr(30,8)),std::stof(pdb_lines[(rep_res_index-1)*3+2].substr(38,8)),std::stof(pdb_lines[(rep_res_index-1)*3+2].substr(46,8)));

    xform = utils::xform_from_3points(N,CA,C);
}

void MotifPairRelativePosScoreMethod::set_motif_insert_pos(Size motif1_insert_pos, Size motif2_insert_pos)
{
    _motif1_insert_pos = motif1_insert_pos;
    _motif2_insert_pos = motif2_insert_pos;
}

void MotifPairRelativePosScoreMethod::set_motif_chain(Size motif1_chain_idx, Size motif2_chain_idx) {
    _motif1_chain_idx = motif1_chain_idx;
    _motif1_chain_idx = motif2_chain_idx;
}

Real MotifPairRelativePosScoreMethod::score(scene::Pose & pose) const {

    EigenXform stub1 = pose.stub(_motif1_insert_pos+_motif1_representative_res_index-1, _motif1_chain_idx); // all 1-index, so minus 1
    EigenXform stub2 = pose.stub(_motif2_insert_pos+_motif2_representative_res_index-1, _motif2_chain_idx); // all 1-index, so minus 1

    EigenXform cur_relative_xform = stub1.inverse(Eigen::Isometry) * stub2;

    Real dist = utils::xform_magnitude(_inv_relative_xform * cur_relative_xform, _motif_radius);

    if( dist <= 5 ) {
        return dist;
    } else {
        return std::log(dist-4)+5;
    }
}

// end
} // namespace scoring




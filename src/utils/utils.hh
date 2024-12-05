#ifndef INCLUDED_utils_utils_hh
#define INCLUDED_utils_utils_hh

#include "apps/args.hh"
#include "basic/types.hh"
#include "scene/Pose.hh"
#include "sampling/FragMover.hh"
#include "utils/random_util.hh"
#include "dssp/dssp.h"
#include "dssp/structure.h"

#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>

namespace utils {

using namespace basic;

std::string random_helix_dssp(Size min_helix_len,Size max_helix_len,Size min_loop_len,Size max_loop_len,Size max_len=65,Size segment_num=-1);
bool check_break(scene::Pose const & pose);

// very naive method
bool random_bundle_dssp(Size length, Size num_helix, Size min_helix_len, Size max_helix_len, Size min_loop_len, Size max_loop_len, std::string & dssp_prefix, std::string & dssp);
bool random_bundle_dssp_AzoF(Size length, Size num_helix, Size min_helix_len, Size max_helix_len, Size min_loop_len, Size max_loop_len, std::string & dssp_prefix, std::string & dssp, Size & insert_pos);
bool random_bundle_dssp_HQA(Size length, Size num_helix, Size min_helix_len, Size max_helix_len, Size min_loop_len, Size max_loop_len, std::string & dssp_prefix, std::string & dssp, Size & insert_pos);
bool random_motif_bundle_dssp(Size length,std::string motif_ss, bool side_require,Size num_helix, Size min_helix_len, Size max_helix_len, Size min_loop_len, Size max_loop_len, std::string & dssp_prefix, std::string & dssp,Size & insert_pos);
bool random_motif_pair_bundle_dssp(Size length,std::string motif1_ss, std::string motif2_ss, Size num_helix, Size min_helix_len, Size max_helix_len, Size min_loop_len, Size max_loop_len, std::string & dssp_prefix, std::string & dssp, std::string & fake_dssp, Size & insert_pos1, Size & insert_pos2);
bool random_motif_bundle_dssp_GFP(Size length,std::string motif_ss,Size num_helix, Size min_helix_len, Size max_helix_len, Size min_loop_len, Size max_loop_len, 
                                    Size GFP_upper_helix_max_len,
                                    Size GFP_upper_helix_min_len,
                                    Size GFP_lower_loop_max_len,
                                    Size GFP_lower_loop_min_len,
                                    Size GFP_lower_helix_max_len,
                                    Size GFP_lower_helix_min_len,
                                    std::string & dssp_prefix, std::string & dssp,Size & insert_pos);
bool random_motif_bundle_Nter_extension_dssp(Size length,std::string motif_ss,Size num_helix, Size min_helix_len, Size max_helix_len, Size min_loop_len, Size max_loop_len, std::string & dssp_prefix, std::string & dssp);
bool random_repeat_bundle_dssp(Size length, Size num_repeats, Size num_helix, Size min_helix_len, Size max_helix_len, Size min_loop_len, Size max_loop_len, std::string & dssp_prefix, std::string & dssp);

Eigen::Array<Real, Eigen::Dynamic, Eigen::Dynamic> per_res_sidechain_neighbors(scene::Pose const & pose, bool including_intra_chain=true, bool including_inter_chain=true);
Real sidechain_neighbors(scene::Pose const & pose, bool including_intra_chain=true, bool including_inter_chain=true);
Real ligand_neighbors(scene::Pose const & pose, Eigen::Matrix<Real, Eigen::Dynamic, 3> const & ligand);

Real get_dihedral(Vec const & a, Vec const & b, Vec const & c, Vec const & d);
Real get_angle(Vec const & a, Vec const & b, Vec const & c);
EigenXform xform_from_3points(Vec const & a, Vec const & b, Vec const & c);
void get_repeat_parameters_from_coords(scene::Pose & pose, Real & rise_out, Real & radius_out, Real & omega_out,Vec & axis_out,Vec & axis_center);
void get_repeat_parameters_from_stubs(scene::Pose & pose,Real & rise_out, Real & radius_out, Real & omega_out,Vec & axis,Vec & axis_center,bool debug=false);
void print_xform(EigenXform const & x);
void print_vec(Vec const & v);

Real rmsd_no_super(const scene::Pose & pose1, const scene::Pose & pose2, Size start_res=1, Size end_res=-1, bool CA_only=false);

std::string get_dssp_from_pose(scene::Pose pose ,Size len=-1,bool reduece_ss=true);

}

#endif

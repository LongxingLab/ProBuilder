#ifndef INCLUDED_utils_dssp_util_hh
#define INCLUDED_utils_dssp_util_hh

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

std::string get_dssp_from_pose(scene::Pose pose ,Size len=-1,bool reduece_ss=true);
std::string normalize_bundle_dssp(std::string dssp_str);

}

#endif

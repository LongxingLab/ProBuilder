#ifndef INCLUDED_utils_ligand_util_hh
#define INCLUDED_utils_ligand_util_hh

#include "basic/types.hh"
#include "scene/Pose.hh"

namespace utils {

using namespace basic;

EigenXform random_rotate_pNPA(std::string & xform_str);
Real pNPA_neighbors(scene::Pose const & pose, EigenXform const & ligand_xform);

EigenXform random_place_AMA(std::string & xform_str);
Real AMA_neighbors(scene::Pose const & pose, EigenXform const & ligand_xform=EigenXform::Identity());

Real HQA_neighbors(scene::Pose const & pose, EigenXform const & ligand_xform=EigenXform::Identity());

EigenXform random_rotate_DTG(std::string & xform_str);

}

#endif

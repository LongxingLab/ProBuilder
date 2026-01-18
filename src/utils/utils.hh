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

bool check_break(scene::Pose const & pose);

Eigen::Array<Real, Eigen::Dynamic, Eigen::Dynamic> per_res_sidechain_neighbors(scene::Pose const & pose, bool including_intra_chain=true, bool including_inter_chain=true);
Real sidechain_neighbors(scene::Pose const & pose, bool including_intra_chain=true, bool including_inter_chain=true);
Real ligand_neighbors(scene::Pose const & pose, Eigen::Matrix<Real, Eigen::Dynamic, 3> const & ligand);

void get_repeat_parameters_from_coords(scene::Pose & pose, Real & rise_out, Real & radius_out, Real & omega_out,Vec & axis_out,Vec & axis_center);
void get_repeat_parameters_from_stubs(scene::Pose & pose,Real & rise_out, Real & radius_out, Real & omega_out,Vec & axis,Vec & axis_center,bool debug=false);
void print_xform(EigenXform const & x);
void print_vec(Vec const & v);

Real rmsd_no_super(const scene::Pose & pose1, const scene::Pose & pose2, Size start_res=1, Size end_res=-1, bool CA_only=false);

}

#endif

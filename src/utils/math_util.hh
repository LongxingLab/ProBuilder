#ifndef INCLUDED_utils_math_util_hh
#define INCLUDED_utils_math_util_hh

#include "basic/types.hh"

#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <random>

namespace utils {

using namespace basic;

Real get_dihedral(Vec const & a, Vec const & b, Vec const & c, Vec const & d);
Real get_angle(Vec const & a, Vec const & b, Vec const & c);
EigenXform xform_from_3points(Vec const & a, Vec const & b, Vec const & c);

Real xform_magnitude(EigenXform const & x, Real lever);

}

#endif

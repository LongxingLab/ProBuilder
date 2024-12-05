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

Real xform_magnitude(EigenXform const & x, Real rg);

}

#endif

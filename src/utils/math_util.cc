#include "basic/types.hh"
#include "utils/math_util.hh"

#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>

namespace utils {

using namespace basic;

Real xform_magnitude(
    EigenXform const & x, Real rg
){
    Real err_trans2 = x.translation().squaredNorm();
    Real cos_theta = (x.rotation().trace()-1.0)/2.0;

    Real err_rot = std::sqrt( std::max( 0.0, 1.0 - cos_theta*cos_theta ) ) * rg;
    if( cos_theta < 0 ) err_rot = rg;
    Real err = std::sqrt( err_trans2 + err_rot*err_rot );
    return err;
}


}

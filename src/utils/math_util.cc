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

Real
get_dihedral(Vec const & a, Vec const & b, Vec const & c, Vec const & d)
{
    Vec b0 = -1.0*(b-a);
    Vec b1 = (c - b);
    Vec b2 = d - c;

    b1.normalize();

    Vec v = b0 - b0.dot(b1) * b1;
    Vec w = b2 - b2.dot(b1) * b1;

    return std::atan2(b1.cross(v).dot(w), v.dot(w));
}

Real
get_angle(Vec const & a, Vec const & b, Vec const & c)
{
    Vec v = a - b;
    v.normalize();
    Vec w = c - b;
    w.normalize();

    return std::acos(v.dot(w));
}


EigenXform
xform_from_3points(Vec const & a, Vec const & b, Vec const & c)
{
    EigenXform X;
    X.translation() = b;

    Vec e1 = (c - b).normalized();
    Vec e3 = e1.cross(a-b).normalized();
    Vec e2 = e3.cross(e1).normalized();
    X.matrix().col(0) = e1;
    X.matrix().col(1) = e2;
    X.matrix().col(2) = e3;

    return X;

}

//Real xform_magnitude(
//    EigenXform const & x, Real rg
//){
//    Real err_trans2 = x.translation().squaredNorm();
//    Real cos_theta = (x.rotation().trace()-1.0)/2.0;
//
//    Real err_rot = std::sqrt( std::max( 0.0, 1.0 - cos_theta*cos_theta ) ) * rg;
//    if( cos_theta < 0 ) err_rot = rg;
//    Real err = std::sqrt( err_trans2 + err_rot*err_rot );
//    return err;
//}

// convert the relative xform into a real number
//Real xform_magnitude(
//    EigenXform const & x,
//    Real lever
//){
//    Real err_trans2 = x.translation().squaredNorm();
//    Real cos_theta = (x.rotation().trace()-1.0)/2.0;
//
//    Real err_rot = std::sqrt( std::max( 0.0, 1.0 - cos_theta*cos_theta ) ) * lever;
//    if( cos_theta < 0 ) err_rot = lever;
//    Real err = std::sqrt( err_trans2 + err_rot*err_rot );
//    return err;
//}


// convert the relative xform into a real number
Real xform_magnitude(
    EigenXform const & x,
    Real lever
){
    Real err_trans2 = x.translation().squaredNorm();
    Real cos_theta = (x.rotation().trace()-1.0)/2.0;
    // make sure in the range of -1.0 to 1.0
    cos_theta = std::max(Real(-1.0), std::min(cos_theta, Real(1.0)));
    Real theta = std::acos(cos_theta);
    Real err_rot = theta * lever;
    Real err = std::sqrt( err_trans2 + err_rot*err_rot );
    return err;
}

}

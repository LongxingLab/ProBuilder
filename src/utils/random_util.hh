#ifndef INCLUDED_utils_random_util_hh
#define INCLUDED_utils_random_util_hh

#include "basic/types.hh"

#include <Eigen/Geometry>

#include <random>

namespace utils {

using namespace basic;

std::mt19937 & rng();

Size random_int(Size lower=0, Size upper=999);
Real random_real(Real lower=0, Real upper=1);

void
rand_xform(
	EigenXform & x,
    bool rand_ori, bool rand_x, bool rand_y, bool rand_z,
	Real const & x_cart_lb=1.5, Real const & x_cart_ub=12.0,
    Real const & y_cart_lb=1.5, Real const & y_cart_ub=12.0,
    Real const & z_cart_lb=0.0, Real const & z_cart_ub=40.0
);

EigenXform rand_roll(Real const & angle_mag, Real const & trans_mag);

Real random_ang();

}

#endif

#include "utils/random_util.hh"

#include <random>
#include <ctime>
#include <chrono>

//static std::mt19937 _mt_rand(time(0));

namespace utils {

using namespace basic;

static std::mt19937 _mt_rand(std::chrono::high_resolution_clock::now().time_since_epoch().count());
// static std::mt19937 _mt_rand(200);

std::mt19937 & rng() {
    return _mt_rand;
}

Size random_int(Size lower, Size upper)
{
    std::uniform_int_distribution<Size> rand_int(lower, upper);
    return rand_int(rng());
}

Real random_real(Real lower, Real upper)
{
    std::uniform_real_distribution<Real> runif(lower,upper);
    return runif(rng());
}

void
rand_xform(
    EigenXform & x,
    bool rand_ori, bool rand_x, bool rand_y, bool rand_z,
    Real const & x_cart_lb, Real const & x_cart_ub,
    Real const & y_cart_lb, Real const & y_cart_ub,
    Real const & z_cart_lb, Real const & z_cart_ub
){
	std::uniform_real_distribution<Real> runif;

    if (rand_ori) {
        std::normal_distribution<Real> rnorm;
        Eigen::Quaterniond qrand( rnorm(rng()), rnorm(rng()), rnorm(rng()), rnorm(rng()) );
        qrand.normalize();
        Eigen::Matrix3d m = qrand.matrix();
        for(Size i = 0; i < 9; ++i) x.data()[i] = m.data()[i];
    }

    Real x_sign = runif(rng()) > 0.5? 1 : -1;
    Real y_sign = runif(rng()) > 0.5? 1 : -1;
    Real z_sign = runif(rng()) > 0.5? 1 : -1;

    if(rand_x) {
        x(0,3) = x_sign * ( (x_cart_ub-x_cart_lb)*runif(rng())+x_cart_lb );        
    }
    if(rand_y) {
        x(1,3) = y_sign * ( (y_cart_ub-y_cart_lb)*runif(rng())+y_cart_lb );        
    }
    if(rand_z) {
        x(2,3) = z_sign * ( (z_cart_ub-z_cart_lb)*runif(rng())+z_cart_lb );        
    }
}

EigenXform
rand_roll(
        Real const & ang_mag,  // in degree
        Real const & trans_mag
){

    std::normal_distribution<Real> rnorm;
	std::uniform_real_distribution<Real> runif;

	Real ang = (1.0 - runif(rng())*runif(rng())) * ang_mag / 180.0 * 3.1415926;
	Vec axis( rnorm(rng()), rnorm(rng()), rnorm(rng()) );
	axis.normalize();
    // assert(axis.norm()>0.01);
	AngleAxis aa( ang, axis );
	EigenXform x( aa );

	Vec delta( 9e9, 9e9, 9e9 );
	while( delta.squaredNorm() > trans_mag*trans_mag ){
		delta = 2.0 * trans_mag * Vec( runif(rng())-0.5, runif(rng())-0.5, runif(rng())-0.5 );
	}

	x.translation() = delta;

    return x;
}

Real random_ang() {
    std::uniform_real_distribution<Real> runif(-3.14,3.14);
    return runif(rng());
}

}

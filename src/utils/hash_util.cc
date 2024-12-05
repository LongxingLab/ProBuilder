#include "utils/hash_util.hh"
#include "basic/macros.hh"

#include <cmath>
#include <iostream>

// @deprecated
// to make it compatible with ss hash
/*
uint64_t xform_hash64(EigenXform const & transform, Real cw, Real aw, Real dist_cutoff)
{
    Vec e = euler_angles(transform);
    Vec t = transform.translation();
    if(fabs(t.x())>dist_cutoff || fabs(t.y())>dist_cutoff || fabs(t.z())>dist_cutoff)
        return (uint64_t)0;
	uint64_t x = (uint64_t)fabs((t.x()+cw/2.0)/cw);
	uint64_t y = (uint64_t)fabs((t.y()+cw/2.0)/cw);
	uint64_t z = (uint64_t)fabs((t.z()+cw/2.0)/cw);
	uint64_t a = (uint64_t)(fmod(e.x()+aw/2.0,360.0)/aw);
	uint64_t b = (uint64_t)(fmod(e.y()+aw/2.0,360.0)/aw);
	uint64_t c = (uint64_t)(fmod(e.z()+aw/2.0,360.0)/aw);

	//assert(a < 1024);
	//assert(b < 1024);
	//assert(c < 1024);
	uint64_t k = x ^ (y<<10) ^ (z<<20) ^ (a<<30) ^ (b<<40) ^ (c<<50);
	k ^= (((uint64_t)(t.x()<0.0))<<60) ^ (((uint64_t)(t.y()<0.0))<<61) ^ (((uint64_t)(t.z()<0.0))<<62);
	return k;
}
*/

namespace utils {

using namespace basic;

uint64_t xform_hash64(EigenXform const & transform, Real cw, Real aw)
{
    Vec e = euler_angles(transform);
    Vec t = transform.translation();
	uint64_t x = (uint64_t)fabs((t.x()+cw/2.0)/cw);
	uint64_t y = (uint64_t)fabs((t.y()+cw/2.0)/cw);
	uint64_t z = (uint64_t)fabs((t.z()+cw/2.0)/cw);
	uint64_t a = (uint64_t)(fmod(e.x()+aw/2.0,360.0)/aw);
	uint64_t b = (uint64_t)(fmod(e.y()+aw/2.0,360.0)/aw);
	uint64_t c = (uint64_t)(fmod(e.z()+aw/2.0,360.0)/aw);

	//assert(a < 1024);
	//assert(b < 1024);
	//assert(c < 1024);
	uint64_t k = x ^ (y<<8) ^ (z<<16) ^ (a<<24) ^ (b<<34) ^ (c<<44);
	k ^= (((uint64_t)(t.x()<0.0))<<60) ^ (((uint64_t)(t.y()<0.0))<<61) ^ (((uint64_t)(t.z()<0.0))<<62);
	return k;
}

uint64_t xform_ss_hash64(EigenXform const & transform, char ss1, char ss2, Real cw, Real aw)
{
    uint64_t ss = ss1=='H'?(ss2=='H'?0:(ss2=='E'?1:2)):
                 (ss1=='E'?(ss2=='H'?3:(ss2=='E'?4:5)):(ss2=='H'?6:(ss2=='E'?7:8)));

    Vec e = euler_angles(transform);
    Vec t = transform.translation();
	uint64_t x = (uint64_t)fabs((t.x()+cw/2.0)/cw);
	uint64_t y = (uint64_t)fabs((t.y()+cw/2.0)/cw);
	uint64_t z = (uint64_t)fabs((t.z()+cw/2.0)/cw);
	uint64_t a = (uint64_t)(fmod(e.x()+aw/2.0,360.0)/aw);
	uint64_t b = (uint64_t)(fmod(e.y()+aw/2.0,360.0)/aw);
	uint64_t c = (uint64_t)(fmod(e.z()+aw/2.0,360.0)/aw);

	//assert(a < 1024);
	//assert(b < 1024);
	//assert(c < 1024);
    // 8 bits for distance
    // 10 bits for angle
	uint64_t k = x ^ (y<<8) ^ (z<<16) ^ (a<<24) ^ (b<<34) ^ (c<<44);
    k ^= (ss<<54);
	k ^= (((uint64_t)(t.x()<0.0))<<60) ^ (((uint64_t)(t.y()<0.0))<<61) ^ (((uint64_t)(t.z()<0.0))<<62);
	return k;
}

Vec euler_angles(EigenXform const & transform)
{
    Vec e;
    const Real FLOAT_PRECISION = 1e-5;
    if( transform(2,2) >= 1-FLOAT_PRECISION) {
        e(0) = std::acos(transform(0,0));
        e(1) = 0.0;
        e(2) = 0.0;
        return e;
    }
    if( transform(2,2) <= -1+FLOAT_PRECISION) {
        e(0) = std::acos(transform(0,0));
        e(1) = 0.0;
        e(2) = M_PI;
        return e;
    }
    Real pos_sin_theta = std::sqrt(1-transform(2,2)*transform(2,2));
    e(2) = std::asin(pos_sin_theta);
    if(transform(2,2)<0){
        e(2) = M_PI-e(2);
    }
    e(0) = std::atan2(transform(2,0),-transform(2,1));
    e(1) = std::atan2(transform(0,2),transform(1,2));

    e(0) += e(0)<0? 2*M_PI : 0.0;
    e(1) += e(1)<0? 2*M_PI : 0.0;
    
    return e/M_PI*180;
}

}

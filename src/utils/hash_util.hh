#ifndef INCLUDED_utils_hash_util_hh
#define INCLUDED_utils_hash_util_hh

#include "basic/types.hh"

namespace utils {

using namespace basic;

uint64_t xform_hash64(EigenXform const & x, Real cw=2.0, Real aw=22.5);
uint64_t xform_ss_hash64(EigenXform const & x, char ss1, char ss2, Real cw=2.0, Real aw=22.5);
Vec euler_angles(EigenXform const & x);

}

#endif

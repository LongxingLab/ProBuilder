#include "scene/pose_util.hh"

#include <cmath>
#include <iostream>

namespace scene {

using namespace basic;

Real pose_rmsd_no_super(scene::Pose & pose1, scene::Pose & pose2)
{
    Real v(0.0);

    for(Size ires=1; ires<=pose1.size(); ++ires) {
        Vec t = pose1.conformation().xyz(ires, ATOM_CA) - pose2.conformation().xyz(ires, ATOM_CA);
        v += t.squaredNorm();
    }

    return std::sqrt(v/pose1.size());
}

}

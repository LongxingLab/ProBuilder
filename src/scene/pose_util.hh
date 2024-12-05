#ifndef INCLUDED_scene_pose_util_hh
#define INCLUDED_scene_pose_util_hh

#include "basic/types.hh"
#include "scene/Pose.hh"

namespace scene {

using namespace basic;

Real pose_rmsd_no_super(scene::Pose & pose1, scene::Pose & pose2);

}

#endif

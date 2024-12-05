#ifndef INCLUDED_scoring_ClashScoreMethod_hh
#define INCLUDED_scoring_ClashScoreMethod_hh

#include "basic/types.hh"
#include "scene/Pose.hh"
#include "scene/Conformation.hh"
#include "scene/Residue.hh"
#include "scoring/Energy.hh"
#include "scoring/ScoreMethod.hh"
#include "basic/VoxelArray.hh"
#include "utils/utils.hh"
#include "utils/string_util.hh"
#include <parallel_hashmap/phmap_dump.h>
#include <iostream>
#include <set>

namespace scoring {
using namespace basic;
using namespace scene;


class ClashScoreMethod:public BaseScoreMethod{
    public:

    ClashScoreMethod(Size sequence_separation=5, Real clash_multiplier=5.0,bool scaffold_intra_clash=true,bool scaffold_target_inter_clash=false,Real interface_extra_punish=5.0);
    virtual Real score(scene::Pose & pose) const override{return base_score(pose);};
    Real base_score(scene::Pose & pose,Size start_chain=1,Size end_chain=-1) const;
    void set_CB_swelling_factor(Real v);
    protected:

    Real _R1R2_pow[5][5];
    Real _CB_swelling_factor;
    Size _min_sequence_separation;
    Real _squared_dist_cutoff;
    Real _clash_multiplier;
    Real _interface_extra_punish;
    bool _scaffold_intra_clash,_scaffold_target_inter_clash;

};

class VoxelClashScoreMethod:public BaseScoreMethod{
    public:
    
    VoxelClashScoreMethod(std::string fname, std::string mode="ALL", Real ball_radius=1.5f, Real grid_resolution = 0.25f,bool no_cache=false,std::string pdb_mode="PROTEIN"):
    BaseScoreMethod(),_no_cache(no_cache),grid_pose(get_grid_pose(fname,pdb_mode)),_clash_multiplier(30),_is_relative_xform_set(false)
    {
        vander_waal_radius = {{ATOM_N,1.54},{ATOM_C,1.70},{ATOM_O,1.52},{ATOM_S,1.8},{ATOM_CA,1.70},{ATOM_CB,1.70},{ATOM_P,1.9}};
        initialize_voxel_grid(mode,ball_radius,grid_resolution);
    };
    VoxelClashScoreMethod(scene::Pose in_pose, std::string mode="ALL", Real ball_radius=1.5f, Real grid_resolution = 0.25f,bool no_cache=false,std::string pdb_mode="PROTEIN"):
    BaseScoreMethod(),_no_cache(no_cache),grid_pose(in_pose),_clash_multiplier(30),_is_relative_xform_set(false)
    {
        vander_waal_radius = {{ATOM_N,1.54},{ATOM_C,1.70},{ATOM_O,1.52},{ATOM_S,1.8},{ATOM_CA,1.70},{ATOM_CB,1.70},{ATOM_P,1.9}};
        initialize_voxel_grid(mode,ball_radius,grid_resolution);
    };
    Pose get_grid_pose(std::string fname, std::string mode);
    void initialize_voxel_grid(std::string mode, Real ball_radius, Real grid_resolution);
    virtual Real score(scene::Pose & pose) const override{return score(pose,99999);};
    bool clash(Vec pos){return _voxel_array->at(pos);}
    Real score(std::vector<EigenXform> & stubs,Size cutoff) const;
    Real score(scene::Pose & pose,Size cutoff) const;
    Real hresl_score(scene::Pose & pose);
    void visualize_voxel_grid( std::string fname );
    void set_relative_pos(EigenXform const & xform);
    protected:
    std::string _pdb_mode;
    std::unordered_map<ATOM_TYPE,Real> vander_waal_radius;
    VoxelArrayOP _voxel_array;
    VoxelIndicateArrayOP _indicate_array;
    VoxelArrayOP _indicate_voxel_array;
    bool _no_cache;
    Real _clash_multiplier;
    EigenXform _relative_xform;
    bool _is_relative_xform_set;
    scene::Pose grid_pose;


};

}

#endif

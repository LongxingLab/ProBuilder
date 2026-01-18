#ifndef INCLUDED_scoring_WholeBodyScoreMethod_hh
#define INCLUDED_scoring_WholeBodyScoreMethod_hh

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

class RepeatScoreMethod:public BaseScoreMethod{
    public:

    RepeatScoreMethod(std::pair<Real,Real> radius_control = std::pair<Real,Real>(0,10000),
                       std::pair<Real,Real> rise_control = std::pair<Real,Real>(0,10000),
                       std::pair<Real,Real> omega_control = std::pair<Real,Real>(0,10000)):
        _radius_control(radius_control),
	_rise_control(rise_control),
	_omega_control(omega_control),
	_radius_multiplier(20.0),
	_rise_multiplier(20.0),
	_omega_multiplier(2.0){};
    virtual Real score(scene::Pose & pose) const override;
    protected:

    std::pair<Real,Real> _radius_control;
    std::pair<Real,Real> _rise_control;
    std::pair<Real,Real> _omega_control;
    // multiplier
    Real _radius_multiplier;
    Real _rise_multiplier;
    Real _omega_multiplier;
};

class MotifPairRelativePosScoreMethod:public BaseScoreMethod {

public:
    MotifPairRelativePosScoreMethod(std::string motif1_pdb, std::string motif2_pdb);
    ~MotifPairRelativePosScoreMethod();

    void initialize();
    void load_motif(std::string pdb, Size & rep_res_index, EigenXform & xform);
    void set_motif_insert_pos(Size motif1_insert_pos, Size motif2_insert_pos);
    void set_motif_chain(Size motif1_chain_idx, Size motif2_chain_idx);

    virtual Real score(scene::Pose & pose) const override;

protected:
    std::string _motif1_pdb;
    std::string _motif2_pdb;

    Size _motif1_chain_idx;
    Size _motif2_chain_idx;

    Size _motif1_insert_pos;
    Size _motif2_insert_pos;

    Size _motif1_representative_res_index; // the residue used to calculate the relative position
    Size _motif2_representative_res_index; // the residue index of moitf2 used to calculate the relative position

    Real _motif_radius;
    EigenXform _inv_relative_xform;

};


}

#endif

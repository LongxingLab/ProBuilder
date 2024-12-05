#ifndef INCLUDED_scoring_ScoreMethod_hh
#define INCLUDED_scoring_ScoreMethod_hh

#include "basic/types.hh"
#include "scene/Pose.hh"
#include "scene/Conformation.hh"
#include "scene/Residue.hh"
#include "scoring/Energy.hh"
#include "basic/VoxelArray.hh"
#include "utils/utils.hh"
#include "utils/string_util.hh"
#include <parallel_hashmap/phmap_dump.h>
#include <iostream>
#include <set>

namespace scoring {
using namespace basic;
using namespace scene;

class BaseScoreMethod{
    public:

    BaseScoreMethod():_weight(1.0),_energy_type(ENERGY_TYPE_UNDEFINED) {};
    virtual Real score(scene::Pose & pose) const=0;
    Real weighted_score(scene::Pose & pose) const{
        if(_weight == 0 ) {
            return 0.0;
        } else {
            return score(pose)*_weight;    
        }
    };
    void set_weight(Real weight){_weight=weight;};
    void set_energy_type(ENERGY_TYPE energy_type){_energy_type=energy_type;};
    void set_score_type(SCORE_TYPE score_type){_score_type = score_type;};
    void set_designable_res(const std::string & designable_res_str);
    void set_designable_res(std::vector<Size> designable_res);
    protected:

    Real _weight;
    ENERGY_TYPE _energy_type;
    SCORE_TYPE _score_type;

    std::set<Size> _designable_res;
};

typedef std::shared_ptr<BaseScoreMethod> ScoreMethodOP;

}

#endif

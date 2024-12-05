#ifndef INCLUDED_scoring_ScoreFunction_hh
#define INCLUDED_scoring_ScoreFunction_hh

#include "basic/types.hh"
#include "scene/Pose.hh"
#include "scene/Conformation.hh"
#include "scene/Residue.hh"
#include "basic/VoxelArray.hh"
#include "utils/utils.hh"
#include "utils/string_util.hh"

#include "scoring/Energy.hh"
#include "scoring/ScoreMethod.hh"
#include "scoring/ClashScoreMethod.hh"
#include "scoring/PrivilegedMotifScoreMethod.hh"
#include "scoring/PrivilegedPairScoreMethod.hh"
#include "scoring/RifScoreMethod.hh"
#include "scoring/RpxScoreMethod.hh"
#include "scoring/WholeBodyScoreMethod.hh"

#include <parallel_hashmap/phmap_dump.h>
#include <iostream>
#include <set>


namespace scoring {
using namespace basic;
using namespace scene;

class ScoreFunction{
    public:
        void regist_method(std::string method_name,ScoreMethodOP methodOP);
        
        void delete_method(std::string method_name);

        ScoreMethodOP get_method(std::string method_name);

        std::vector<ScoreMethodOP> get_all_methods () const;

        Real score(scene::Pose & pose);

        Real score(scene::Pose & pose,std::vector<std::string> exclude_methods);


    private:
        std::unordered_map<std::string,ScoreMethodOP> _method_map;
};

}

#endif

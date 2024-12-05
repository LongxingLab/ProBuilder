#ifndef INCLUDED_scoring_RifScoreMethod_hh
#define INCLUDED_scoring_RifScoreMethod_hh

#include "basic/types.hh"
#include "scoring/ScoreMethod.hh"
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

class RifScoreMethod:public BaseScoreMethod{
    public:

    RifScoreMethod(std::string rif_table_path,Real cart_resl=1.0,Real angle_resl=16.0);
    virtual Real score(scene::Pose & pose) const override;
    Real score_with_anchor_seqpos(scene::Pose & pose, std::vector<Size> & anchor_seqpos, bool return_seqpos=true) const;
    void unserialize_phmap(std::string fname);
    void set_relative_pos(EigenXform const & xform);
    protected:

    std::string _rif_table_path;
    hash_table<uint64_t, Real> _rif_table;
    wefoldHasher _hasher; 
    EigenXform _relative_xform;
    bool _is_relative_xform_set;


};


}

#endif

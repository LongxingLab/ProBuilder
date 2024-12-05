#ifndef INCLUDED_scoring_PrivilegedPairScoreMethod_hh
#define INCLUDED_scoring_PrivilegedPairScoreMethod_hh

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


class DisulfideScoreMethod:public BaseScoreMethod{
    public:

    DisulfideScoreMethod(std::string disulfide_hash_table_path, Size disulfide_min_sequence_separation=8, Real cart_resl=1.0,Real angle_resl=16.0);
    ~DisulfideScoreMethod();
    virtual Real score(scene::Pose & pose) const override{return base_score(pose);};
    Real base_score(scene::Pose & pose,Size start_chain=1,Size end_chain=-1) const;
    bool check_disulfide_exist(EigenXform const & stub1, EigenXform const & stub2) const;
    void unserialize_phmap(std::string const & fname);

    void set_scnb_burial_cutoff(Real cutoff) {_scnb_burial_cutoff = cutoff;};

    protected:

    std::string _disulfide_hash_table_path;
    hash_set_table<uint64_t, int32_t> _hash_table;
    wefoldHasher _hasher; 
    Size _min_sequence_separation;
    Real _squared_dist_cutoff;
    Real _scnb_burial_cutoff;

};

class PrivilegedPairScoreMethod:public BaseScoreMethod{
    public:

    PrivilegedPairScoreMethod(std::string hash_table_path, Size min_sequence_separation=8, Real cart_resl=1.0,Real angle_resl=16.0);
    ~PrivilegedPairScoreMethod();
    virtual Real score(scene::Pose & pose) const override{return base_score(pose);};
    Real base_score(scene::Pose & pose,Size start_chain=1,Size end_chain=-1) const;
    bool check_pair_exist(EigenXform const & stub1, EigenXform const & stub2) const;
    void unserialize_phmap(std::string const & fname);

    void set_scnb_burial_cutoff(Real cutoff) {_scnb_burial_cutoff = cutoff;};

    protected:

    std::string _hash_table_path;
    hash_set_table<uint64_t, uint8_t> _hash_table;
    wefoldHasher _hasher; 
    Size _min_sequence_separation;
    Real _squared_dist_cutoff;
    Real _scnb_burial_cutoff;
    bool _ignore_twin_pair;

};


}

#endif

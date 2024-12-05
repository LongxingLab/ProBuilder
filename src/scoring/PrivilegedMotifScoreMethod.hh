#ifndef INCLUDED_scoring_PrivilegedMotifScoreMethod_hh
#define INCLUDED_scoring_PrivilegedMotifScoreMethod_hh

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

class PrivilegedMotifScoreMethod:public BaseScoreMethod{

public:
    //
    PrivilegedMotifScoreMethod(std::vector<std::string> motif_hash_table_paths, Real cart_resl=1.0, Real angle_resl=16.0);
    ~PrivilegedMotifScoreMethod();

    virtual Real score(scene::Pose & pose) const override;

    void unserialize_phmaps(std::vector<std::string> const & fnames);
    void set_relative_pos(EigenXform const & xform);
protected:
    //
    std::vector<std::string> _motif_hash_table_paths;
    Size _num_hash_tables;
    std::vector<hash_set_table<uint64_t, int32_t> > _hash_tables;
    wefoldHasher _hasher;
    EigenXform _relative_xform;
    bool _is_relative_xform_set;
};


// code for NCAA at the interface, so I only consider the interaction between chain A and chain B
// TODO: more chains
// at current stage, one NCAA per interface
// if you want more, let me know
// longxing, 2023-09-15
class PrivilegedInterfaceMotifScoreMethod:public BaseScoreMethod {

public:
    //
    PrivilegedInterfaceMotifScoreMethod(std::string motif_hash_table_path, Real cart_resl=1.0, Real angle_resl=16.0);
    ~PrivilegedInterfaceMotifScoreMethod();

    virtual Real score(scene::Pose & pose) const override;

    void unserialize_phmap(std::string const & fname);

protected:
    //
    std::string _motif_hash_table_path;
    hash_set_table<uint64_t, int32_t> _hash_table;
    wefoldHasher _hasher;
};

// metal coordination score is a zero (whole) body energy

class MetalCoordinationScoreMethod:public BaseScoreMethod{
    public:

    MetalCoordinationScoreMethod(std::string metal_hash_table_path, Size coordination_number, Real cart_resl=0.5,Real angle_resl=15.0);
    ~MetalCoordinationScoreMethod();
    virtual Real score(scene::Pose & pose) const override;
    Real score_with_tokens(scene::Pose & pose, std::vector<Size> & coordination_res, Size & token) const;
    void unserialize_phmap(std::string const & fname);
    Size coordination_number() const {return _coordination_number;};

    protected:

    std::string _metal_hash_table_path;
    Size _coordination_number;
    hash_table<uint64_t, int32_t> _hash_table;
    wefoldHasher _hasher;

};


}

#endif

#ifndef INCLUDED_scoring_RpxScoreMethod_hh
#define INCLUDED_scoring_RpxScoreMethod_hh

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


class RpxScoreMethod:public BaseScoreMethod{
    public:

    RpxScoreMethod(std::string rpx_table_path,Real cart_resl=1.0,Real angle_resl=16.0,bool use_ss=true,bool include_target=true);
    ~RpxScoreMethod();
    virtual Real score(scene::Pose & pose) const override{return base_score(pose);};
    Real base_score(scene::Pose & pose,Size start_chain=1,Size end_chain=-1) const;
    Real score_rpx_pair(EigenXform const & stub1, EigenXform const & stub2, char ss1, char ss2) const;
    void unserialize_phmap(std::string const & fname);
    void only_count_two_repeat(bool option=true){_only_count_two_repeat=option;}
    void set_inter_chain_weight(Real weight) {_inter_chain_weight = weight;}
    protected:

    std::string _rpx_table_path;
    hash_table<uint64_t, Real> _rpx_table;
    wefoldHasher _hasher; 
    Size _min_sequence_separation;
    Real _squared_dist_cutoff;
    bool _use_ss;
    bool _only_count_two_repeat = false;
    bool _include_target;
    Real _inter_chain_weight;

};

class RpxScoreInterfaceMethod:public RpxScoreMethod{
    public:
        //
        RpxScoreInterfaceMethod(Size first_chain_len, std::string rpx_table_path,Real cart_resl=1.0,Real angle_resl=16.0,bool use_ss=true):RpxScoreMethod(rpx_table_path,cart_resl,angle_resl,use_ss), _first_chain_len(first_chain_len){};

        virtual Real score(scene::Pose & pose) const;
    private:
        //
        Size _first_chain_len;

};

// this class is for rpx betwwen target pose and origin pose
class RpxScoreTargetMethod:public RpxScoreMethod{
    public:

    RpxScoreTargetMethod(std::string rpx_table_path,Real cart_resl=1.0,Real angle_resl=16.0,bool use_ss=true):RpxScoreMethod(rpx_table_path,cart_resl,angle_resl,use_ss){};

    virtual Real score(scene::Pose & pose) const;

    private:
    // use this for some score that has no energy cache
    
};



class Rpx1SideScoreMethod:public BaseScoreMethod{
    public:

    Rpx1SideScoreMethod(std::string rpx_table_path,Real cart_resl=3.5,Real angle_resl=38.0,bool use_ss=false,bool no_chache=false);
    virtual Real score(scene::Pose & pose) const override;
    std::vector<Real> score(std::vector<EigenXform> & stubs,scene::Pose & target_pose) const;
    Real score_rpx_pair(EigenXform const & stub1, EigenXform const & stub2, char aa, char ss1='H', char ss2='H') const;
    void set_selected_target_res(std::vector<Size> & select_res){_select_target_res = select_res;}
    void set_target_residue_weight(std::vector<Real> & weight){_res_weight = weight;}
    void unserialize_phmap(std::string fname);

    protected:
    std::string _rpx_table_path;
    hash_table<uint64_t, Real> _rpx_table;
    wefoldHasher _hasher;
    std::vector<Size> _select_target_res;
    std::vector<Real> _res_weight;
    Real _squared_dist_cutoff;
    bool _use_ss;
    bool _no_cache;

};


class PreRpx1SideScoreMethod:public BaseScoreMethod{
    public:

    PreRpx1SideScoreMethod(std::string rpx1side_table_path,Real cart_resl=0.5,Real angle_resl=16.0,bool no_cache=false);
    virtual Real score(scene::Pose & pose) const override;

    Real score(std::vector<EigenXform> & stubs) const;
    void unserialize_phmap(std::string fname);

    protected:

    std::string _rpx1side_table_path;
    hash_table<uint64_t, Real> _rpx1side_table;
    wefoldHasher _hasher; 
    bool _no_cache;
};



}

#endif

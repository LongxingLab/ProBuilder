#include "scoring/RpxScoreMethod.hh"
#include "scoring/Energy.hh"
#include "utils/hash_util.hh"

#include <cmath>
#include <string>
#include <fstream>
#include <Eigen/Geometry>
#include <parallel_hashmap/phmap_dump.h>
#include <map>
namespace scoring {

using namespace basic;
using namespace scene;

// rpx score

RpxScoreMethod::RpxScoreMethod(std::string rpx_table_path,Real cart_resl,Real angle_resl,bool use_ss,bool include_target):
    BaseScoreMethod(),
    _rpx_table_path(rpx_table_path),
    _hasher(cart_resl,angle_resl,512),
    _use_ss(use_ss),
    _min_sequence_separation(5),
    _squared_dist_cutoff(10.0*10.0),
    _include_target(include_target),
    _inter_chain_weight(1.0)
{
    set_score_type(RPX);
    unserialize_phmap(_rpx_table_path);
}

RpxScoreMethod::~RpxScoreMethod() {}

Real RpxScoreMethod::base_score(scene::Pose &pose,Size start_chain,Size end_chain) const {
    Size full_len = pose.size(),len = pose.size();
    std::string sequence = pose.sequence();
    if(pose.num_repeats()!=1 && _only_count_two_repeat && pose.num_chains()==1)len = 2*len/pose.num_repeats();
    EnergyManager & energy_manager = pose.energy_manager();
    std::vector<EnergyOP> energyOPs = energy_manager.get_energy_tables(_score_type);
    Real total_rpx = 0;
    std::string const dssp = pose.dssp();
    for(auto energy_base:energyOPs){
        TwoBodyEnergyOP energy_rpx = std::static_pointer_cast<TwoBodyEnergy>(energy_base);
        if(energy_rpx == nullptr){
            std::cout<<"energy op from base to twobody fail!! check energy initialize is right"<<std::endl;
            exit(0);
        }
        if(energy_rpx->jchain()<start_chain) continue;
        if(end_chain!=-1&&energy_rpx->jchain()>end_chain)continue;
        if(_energy_type != ENERGY_TYPE_UNDEFINED && _energy_type != energy_rpx->energy_type() ) continue;
    // std::string sequence = pose.sequence();
        for(Size ires=1; ires<=len; ++ires) {
            char ss1 = dssp.at(ires-1);
            EigenXform const & stub1 = pose.stub(ires,energy_rpx->ichain());

            for(Size jres=(energy_rpx->energy_type()==TWO_BODY_INTRA_CHAIN?ires+_min_sequence_separation:1); jres<=len; ++jres) {

                if(sequence[ires-1] == 'X' || sequence[jres-1] == 'X') {
                    energy_rpx->set_pair_score(ires, jres, 0.0);
                    continue;
                }

                // move this to the outer loop???
                // to be safe, each time reset the energy graph to zero
                if(_designable_res.size() != 0 && (_designable_res.count(ires) == 0 || _designable_res.count(jres) == 0) ) {
                    energy_rpx->set_pair_score(ires, jres, 0.0);
                    continue;
                }

                if( !energy_rpx->is_changed_pair(ires, jres) ) {
                    continue;
                }
                EigenXform const & stub2 = pose.stub(jres,energy_rpx->jchain());

                Vec t = stub1.translation() - stub2.translation();

                /*
                if(std::isinf(t.squaredNorm())) {
                    // pose.dump_pdb("inf.pdb",false,false,false);
                    pose.dump_pdb("inf.pdb",false,false,false,false);
                    std::vector<EigenXform> jumps = pose.jumps();
                    utils::print_xform(jumps[0]);
                    std::cout << "Bad!!!" << std::endl;
                    exit(0);
                }
                */
                
                // std::cout << t.squaredNorm() << std::endl;

                Real pair_dist = t.squaredNorm();

                if(pair_dist > _squared_dist_cutoff || std::isinf(pair_dist) || std::isnan(pair_dist) ) {
                    energy_rpx->set_pair_score(ires, jres, 0.0);
                } else {
                    char ss2 = dssp.at(jres-1);
                    
                    Real this_pair_rpx_score = score_rpx_pair(stub1, stub2, ss1, ss2);

                    if(energy_rpx->ichain() != energy_rpx->jchain()) {
                        this_pair_rpx_score *= _inter_chain_weight;
                    }

                    energy_rpx->set_pair_score(ires, jres, this_pair_rpx_score);
                }
            }

        }
        total_rpx+= energy_rpx->weighted_score();
    }
    // the code below do not compat with muti chain now
    if(pose.get_target_pose() != nullptr && _include_target){
        Real inter_score=0;
        scene::PoseOP target_poseOP = pose.get_target_pose();
        std::string target_dssp = target_poseOP->dssp();
        Size target_len = target_poseOP->size();
        EigenXform transform = pose.conformation().root().inverse()*target_poseOP->conformation().root();
        for (Size ires =1;ires<=full_len;ires++){
            char ss1 = dssp.at(ires-1);
            EigenXform const & stub1 = pose.stub(ires);

            for(Size jres=1; jres<=target_len; ++jres) {

                EigenXform const & stub2 = transform*target_poseOP->stub(jres);

                Vec t = stub1.translation() - stub2.translation();

                if(t.squaredNorm() > _squared_dist_cutoff) {
                    continue;
                } else {
                    char ss2 = target_dssp.at(jres-1);
                    
                    Real this_pair_rpx_score = score_rpx_pair(stub1, stub2, ss1, ss2);
                    inter_score+= this_pair_rpx_score;
                    
                }
            }
        }
        total_rpx+= inter_score;
    }
    return total_rpx;
}

Real RpxScoreInterfaceMethod::score(scene::Pose &pose) const  {
    Size full_len = pose.size();
    Real total_rpx = 0;
    std::string const dssp = pose.dssp();

    for(Size ires=1; ires<=_first_chain_len; ++ires) {
        char ss1 = dssp.at(ires-1);
        EigenXform const & stub1 = pose.stub(ires);
        for(Size jres=_first_chain_len+1; jres<=full_len; ++jres) {
            EigenXform const & stub2 = pose.stub(jres);

            Vec t = stub1.translation() - stub2.translation();

            if(t.squaredNorm() > _squared_dist_cutoff) {
                continue;
            } else {
                char ss2 = dssp.at(jres-1);
                
                total_rpx += score_rpx_pair(stub1, stub2, ss1, ss2);
            }
        } // second chain loop
    } // first chain loop

    return total_rpx;

}

Real RpxScoreTargetMethod::score(scene::Pose &pose) const  {
    Size full_len = pose.size();
    Real total_rpx = 0;
    std::string const dssp = pose.dssp();
    // the code below do not compat with muti chain now
    if(pose.get_target_pose() != nullptr){
        Real inter_score=0;
        scene::PoseOP target_poseOP = pose.get_target_pose();
        std::string target_dssp = target_poseOP->dssp();
        Size target_len = target_poseOP->size();

        for (Size ires =1;ires<=full_len;ires++){
            char ss1 = dssp.at(ires-1);
            EigenXform const & stub1 = pose.stub(ires);

            for(Size jres=1; jres<=target_len; ++jres) {

                EigenXform const & stub2 = target_poseOP->stub(jres);

                Vec t = stub1.translation() - stub2.translation();

                if(t.squaredNorm() > _squared_dist_cutoff) {
                    continue;
                } else {
                    char ss2 = target_dssp.at(jres-1);
                    
                    Real this_pair_rpx_score = score_rpx_pair(stub1, stub2, ss1, ss2);
                    inter_score+= this_pair_rpx_score;
                    
                }
            }
        }
        total_rpx+= inter_score;
    }

    return total_rpx;
}

Real RpxScoreMethod::score_rpx_pair(EigenXform const & stub1, EigenXform const & stub2, char ss1, char ss2) const{
    uint64_t hash1(0), hash2(0);
    Real this_pair_score(0.0);
    Size counts(0);

    if (_use_ss) {
        // hash1 = xform_ss_hash64(stub1.inverse(Eigen::Isometry)*stub2, ss1, ss2);
        // hash2 = xform_ss_hash64(stub2.inverse(Eigen::Isometry)*stub1, ss2, ss1);
        hash1 = _hasher.get_key(stub1.inverse(Eigen::Isometry)*stub2)^(u_int64_t(ss1=='H'?1:(ss1=='L'?2:3))<<60)^((u_int64_t(ss2=='H'?1:(ss2=='L'?2:3)))<<62);
        hash2 = _hasher.get_key(stub2.inverse(Eigen::Isometry)*stub1)^(u_int64_t(ss2=='H'?1:(ss2=='L'?2:3))<<60)^(u_int64_t((ss1=='H'?1:(ss1=='L'?2:3)))<<62);
    } else {
        hash1 = _hasher.get_key(stub1.inverse(Eigen::Isometry)*stub2);
        hash2 = _hasher.get_key(stub2.inverse(Eigen::Isometry)*stub1);
    }

    // too far apart
    if( hash1 == 0 || hash2 == 0) {

    } else {

        hash_table<uint64_t,Real>::const_iterator got1 = _rpx_table.find(hash1);
        hash_table<uint64_t,Real>::const_iterator got2 = _rpx_table.find(hash2);
        if(got1 != _rpx_table.end()){
            counts += 1;
            this_pair_score += got1->second;
        }
        if(got2 != _rpx_table.end()){
            counts += 1;
            this_pair_score += got2->second;
        }
    }

    if(counts==0) {
        return 0.0;
    } else {
        return this_pair_score / counts;
    }
}

void RpxScoreMethod::unserialize_phmap(std::string const & fname){
    #ifdef USE_PHMAP
        phmap::BinaryInputArchive ar_in(fname.c_str());
        _rpx_table.clear();
        _rpx_table.load(ar_in);
    #else
        std::cout << "Not implemented yet!!! You have to modify the rifdock code to dump rif table.!" << std::endl;
        exit(-1);
    #endif
    assert(_rpx_table.size()!=0);
}


Rpx1SideScoreMethod::Rpx1SideScoreMethod(std::string rpx_table_path,Real cart_resl,Real angle_resl,bool use_ss,bool no_chache):
    _rpx_table_path(rpx_table_path),
    _hasher(wefoldHasher(cart_resl,angle_resl,512)),
    _use_ss(use_ss),
    _no_cache(no_chache),
    _squared_dist_cutoff(15*15)
{
    set_score_type(RPX1SIDE);
    unserialize_phmap(_rpx_table_path);

};


Real Rpx1SideScoreMethod::score_rpx_pair(EigenXform const & scaffold_stub, EigenXform const & target_stub, char aa_type,char ss1, char ss2) const{
    uint64_t  hash(0);
    Real this_pair_score(0.0);
    
    if (_use_ss) {
        // hash1 = xform_ss_hash64(stub1.inverse(Eigen::Isometry)*stub2, ss1, ss2);
        // hash2 = xform_ss_hash64(stub2.inverse(Eigen::Isometry)*stub1, ss2, ss1);
        std::cout<<"It is not recommended to use ss here; Not implemented yet!!"<<std::endl;
        exit(-1);
    } else {
        hash = _hasher.get_key(target_stub.inverse(Eigen::Isometry)*scaffold_stub) | (uint64_t(aa_type)<<50);
    }

    hash_table<uint64_t,Real>::const_iterator got = _rpx_table.find(hash);
    if(got!=_rpx_table.end())this_pair_score = got->second;

    return this_pair_score;
}

Real Rpx1SideScoreMethod::score(scene::Pose &pose) const {
    scene::PoseOP target_poseOP = pose.get_target_pose();
    Size len = pose.size();
    OneBodyEnergyOP energy_rpx1side;
    EnergyManager & energy_manager = pose.energy_manager();
    std::vector<EnergyOP> energyOPs = energy_manager.get_energy_tables(_score_type);

    assert(_no_cache || (energyOPs.size()==1 && energyOPs[0]->energy_type()==ONE_BODY) );
    assert( _res_weight.size()==0 || target_poseOP->size()==_res_weight.size() );
    
    if(!_no_cache)energy_rpx1side= std::static_pointer_cast<OneBodyEnergy>(energyOPs[0]);
    
    // iterate through each scaffold residue and each target residue
    Real total_rpx = 0;
    for(Size ires=1; ires<=len; ++ires) {

        // unchanged residue
        if( !_no_cache && !energy_rpx1side->is_changed_residue(ires) ) {
            continue;
        }

        Real res_score=0;
        EigenXform const & stub1 = pose.stub(ires,1);

        if(_select_target_res.size() == 0) {
            for(Size jres = 1; jres <= target_poseOP->size(); ++jres) {
                Vec  res1_cb_xyz = pose.xyz(ires, ATOM_CB);
                Vec  res2_ca_xyz = target_poseOP->xyz(jres, ATOM_CA);
                Vec t = res1_cb_xyz-res2_ca_xyz;
                if(t.squaredNorm() > _squared_dist_cutoff) {
                    continue;
                }
                EigenXform const & stub2 = target_poseOP->stub(jres);
                char aa2 = target_poseOP->sequence().at(jres-1);
                res_score += score_rpx_pair(stub1, stub2, aa2)*(_res_weight.size()==0?1:_res_weight[jres-1]);
            }
        } else {
            for(Size jres:_select_target_res) {
                Vec  res1_cb_xyz = pose.xyz(ires, ATOM_CB);
                Vec  res2_ca_xyz = target_poseOP->xyz(jres, ATOM_CA);
                Vec t = res1_cb_xyz-res2_ca_xyz;
                if(t.squaredNorm() > _squared_dist_cutoff) {
                    continue;
                }
                EigenXform const & stub2 = target_poseOP->stub(jres);
                char aa2 = target_poseOP->sequence().at(jres-1);
                res_score += score_rpx_pair(stub1, stub2, aa2)*(_res_weight.size()==0?1:_res_weight[jres-1]);
            }
        }
        if(!_no_cache)energy_rpx1side->set_score_per_residue(ires,res_score);
        total_rpx +=res_score;
    }
    if(!_no_cache)return energy_rpx1side->weighted_score();
    return total_rpx;
}

std::vector<Real> Rpx1SideScoreMethod::score(std::vector<EigenXform> & stubs, scene::Pose & target_pose) const {
    std::vector<Real> total_rpx(stubs.size(),0);
    assert(target_pose.size()==_res_weight.size());
    

    for(Size ires=1; ires<=stubs.size(); ++ires) {
        Real res_score=0;
        EigenXform const & stub1 = stubs[ires-1];

        if(_select_target_res.size() == 0) {
            for(Size jres = 1; jres <= target_pose.size(); ++jres) {
                EigenXform const & stub2 = target_pose.stub(jres);
                char aa2 = target_pose.sequence().at(jres-1);
                res_score += score_rpx_pair(stub1, stub2, aa2)*(_res_weight.size()==0?1:_res_weight[jres-1]);
            }
        } else {
            for(Size jres:_select_target_res) {
                EigenXform const & stub2 = target_pose.stub(jres);
                char aa2 = target_pose.sequence().at(jres-1);
                res_score += score_rpx_pair(stub1, stub2, aa2)*(_res_weight.size()==0?1:_res_weight[jres-1]);
            }
        }
        total_rpx[ires-1] = res_score;
    }
    return total_rpx;
}


void Rpx1SideScoreMethod::unserialize_phmap(std::string fname){
    #ifdef USE_PHMAP
        phmap::BinaryInputArchive ar_in(fname.c_str());
        _rpx_table.clear();
        _rpx_table.load(ar_in);
    #else
        std::cout << "Not implemented yet!!! You have to modify the rifdock code to dump rif table.!" << std::endl;
        exit(-1);
    #endif
    assert(_rpx_table.size()!=0);
}

// pre rpx 1side


PreRpx1SideScoreMethod::PreRpx1SideScoreMethod(std::string rpx1side_table_path,Real cart_resl,Real angle_resl,bool no_cache):
    _rpx1side_table_path(rpx1side_table_path),
    _hasher(wefoldHasher(cart_resl,angle_resl,512)),
    _no_cache(no_cache)
{
    set_score_type(PREBUILD_RPX1SIDE);
    unserialize_phmap(_rpx1side_table_path);
};


Real PreRpx1SideScoreMethod::score(scene::Pose &pose)const{
    Size len = pose.size();

    OneBodyEnergyOP energy_rpx1side;
    EnergyManager & energy_manager = pose.energy_manager();
    std::vector<EnergyOP> energyOPs = energy_manager.get_energy_tables(_score_type);
    assert(_no_cache || (energyOPs.size()==1 && energyOPs[0]->energy_type()==ONE_BODY) );
     if(!_no_cache)energy_rpx1side= std::static_pointer_cast<OneBodyEnergy>(energyOPs[0]);

    Real total_rpx = 0;
    for(Size ires=1; ires<=len; ++ires) {

        // unchanged residue
        if( !_no_cache && !energy_rpx1side->is_changed_residue(ires) ) {
            continue;
        }

        Real res_score=0;
        EigenXform const & stub = pose.stub(ires,1);

        uint64_t key = _hasher.get_key(stub);
        hash_table<uint64_t,Real>::const_iterator got = _rpx1side_table.find(key);
        if(got != _rpx1side_table.end()){
            // std::cout<<got->second<<std::endl;
            total_rpx+=got->second;
            if(!_no_cache)energy_rpx1side->set_score_per_residue(ires,got->second);
        }
    }
    if(!_no_cache)return energy_rpx1side->weighted_score();
    return total_rpx;
}

Real PreRpx1SideScoreMethod::score(std::vector<EigenXform> & stubs)const{
    Real rpx_total = 0;
    for(EigenXform const & stub:stubs) {
        uint64_t key = _hasher.get_key(stub);
        hash_table<uint64_t,Real>::const_iterator got = _rpx1side_table.find(key);
        if(got != _rpx1side_table.end()){
            // std::cout<<got->second<<std::endl;
            rpx_total+=got->second;
        }
    }
    return rpx_total;
}

void PreRpx1SideScoreMethod::unserialize_phmap(std::string fname){
    #ifdef USE_PHMAP
        phmap::BinaryInputArchive ar_in(fname.c_str());
        _rpx1side_table.clear();
        _rpx1side_table.load(ar_in);
    #else
        std::cout << "Not implemented yet!!!" << std::endl;
        exit(-1);
    #endif
    assert(_rpx1side_table.size()!=0);
}


} // namespace scoring


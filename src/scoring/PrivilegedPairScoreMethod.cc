#include "scoring/PrivilegedPairScoreMethod.hh"
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


// disulfide score

DisulfideScoreMethod::DisulfideScoreMethod(std::string disulfide_hash_table_path, Size disulfide_min_sequence_separation, Real cart_resl,Real angle_resl):
    BaseScoreMethod(),
    _disulfide_hash_table_path(disulfide_hash_table_path),
    _min_sequence_separation(disulfide_min_sequence_separation),
    _hasher(cart_resl,angle_resl,512),
    _squared_dist_cutoff(8.5*8.5),
    _scnb_burial_cutoff(99999.99)
{
    set_score_type(DISULFIDE);
    unserialize_phmap(_disulfide_hash_table_path);
}

DisulfideScoreMethod::~DisulfideScoreMethod() {}

Real DisulfideScoreMethod::base_score(scene::Pose &pose,Size start_chain,Size end_chain) const {
    Size len = pose.size();
    EnergyManager & energy_manager = pose.energy_manager();
    std::vector<EnergyOP> energyOPs = energy_manager.get_energy_tables(_score_type);
    Real total_disulfide_score = 0;

    auto scnb_vals = utils::per_res_sidechain_neighbors(pose);

    for(auto energy_base:energyOPs){
        TwoBodyEnergyOP energy_disulfide = std::static_pointer_cast<TwoBodyEnergy>(energy_base);
        if(energy_disulfide == nullptr){
            std::cout<<"energy op from base to twobody fail!! check energy initialize is right"<<std::endl;
            exit(0);
        }
        if(energy_disulfide->jchain()<start_chain) continue;
        if(end_chain!=-1&&energy_disulfide->jchain()>end_chain)continue;
        if(_energy_type != ENERGY_TYPE_UNDEFINED && _energy_type != energy_disulfide->energy_type() ) continue;
    // std::string sequence = pose.sequence();

        std::set<Size> occupied;

        for(Size ires=1; ires<=len; ++ires) {
            EigenXform const & stub1 = pose.stub(ires,energy_disulfide->ichain());

            for(Size jres=(energy_disulfide->energy_type()==TWO_BODY_INTRA_CHAIN?ires+_min_sequence_separation:1); jres<=len; ++jres) {

                // move this to the outer loop???
                // to be safe, each time reset the energy graph to zero
                if(_designable_res.size() != 0 && (_designable_res.count(ires) == 0 || _designable_res.count(jres) == 0) ) {
                    energy_disulfide->set_pair_score(ires, jres, 0.0);
                    continue;
                }
                // clear the rest to zero
                // have to check this in the first place
                // must check this befor is_changed_pair
                if (occupied.count(ires) != 0 || occupied.count(jres) != 0 || scnb_vals(ires-1,0)>_scnb_burial_cutoff || scnb_vals(jres-1,0)>_scnb_burial_cutoff ) {
                    energy_disulfide->set_pair_score(ires, jres, 0.0);
                    continue;
                }

                if( !energy_disulfide->is_changed_pair(ires, jres) ) {
                    continue;
                }

                EigenXform const & stub2 = pose.stub(jres,energy_disulfide->jchain());

                Vec t = stub1.translation() - stub2.translation();

                if(t.squaredNorm() > _squared_dist_cutoff) {
                    energy_disulfide->set_pair_score(ires, jres, 0.0);
                } else {
                    
                    bool find_a_disulfide_bond = check_disulfide_exist(stub1, stub2);

                    if (find_a_disulfide_bond) {
                        energy_disulfide->set_pair_score(ires, jres, -1.0);

                        occupied.insert(ires);
                        occupied.insert(jres);

                    } else {
                        energy_disulfide->set_pair_score(ires, jres, 0.0);
                    }
                }
            }

        }
        total_disulfide_score+= energy_disulfide->weighted_score();
    }
    return total_disulfide_score;
}

bool DisulfideScoreMethod::check_disulfide_exist(EigenXform const & stub1, EigenXform const & stub2) const{
    uint64_t hash1(0), hash2(0);
    
    hash1 = _hasher.get_key(stub1.inverse(Eigen::Isometry)*stub2);
    hash2 = _hasher.get_key(stub2.inverse(Eigen::Isometry)*stub1);


    // too far apart
    if( hash1 == 0 || hash2 == 0) {
        return false;
    } else {
        auto got1 = _hash_table.find(hash1);
        auto got2 = _hash_table.find(hash2);
        if(got1 != _hash_table.end() || got2 != _hash_table.end()){
            return true;
        }
    }

    return false;
}

void DisulfideScoreMethod::unserialize_phmap(std::string const & fname){
    #ifdef USE_PHMAP
        phmap::BinaryInputArchive ar_in(fname.c_str());
        _hash_table.clear();

        size_t size;

        ar_in.load( & size );
        _hash_table.reserve( size );

        for ( size_t idx = 0; idx < size; idx++ ) {
            uint64_t k; phmap::flat_hash_set<int32_t> v;

            ar_in.load( & k );
            v.load( ar_in );

            _hash_table.insert_or_assign( std::move( k ), std::move( v ) );
        }
    #else
        std::cout << "Not implemented yet!!!" << std::endl;
        exit(-1);
    #endif
    assert(_hash_table.size()!=0);
}


// privileged pair score

PrivilegedPairScoreMethod::PrivilegedPairScoreMethod(std::string hash_table_path, Size min_sequence_separation, Real cart_resl,Real angle_resl):
    BaseScoreMethod(),
    _hash_table_path(hash_table_path),
    _min_sequence_separation(min_sequence_separation),
    _hasher(cart_resl,angle_resl,512),
    _squared_dist_cutoff(10.0*10.0),
    _scnb_burial_cutoff(99999.99),
    _ignore_twin_pair(true)
{
    set_score_type(PRIVILEGED_PAIR);
    unserialize_phmap(_hash_table_path);
}

PrivilegedPairScoreMethod::~PrivilegedPairScoreMethod() {}

Real PrivilegedPairScoreMethod::base_score(scene::Pose &pose,Size start_chain,Size end_chain) const {
    Size len = pose.size();
    EnergyManager & energy_manager = pose.energy_manager();
    std::vector<EnergyOP> energyOPs = energy_manager.get_energy_tables(_score_type);
    Real total_score = 0;

    // inter chain only
    auto scnb_vals = utils::per_res_sidechain_neighbors(pose, false, true);

    for(auto energy_base:energyOPs){
        TwoBodyEnergyOP energy_table = std::static_pointer_cast<TwoBodyEnergy>(energy_base);
        if(energy_table == nullptr){
            std::cout<<"energy op from base to twobody fail!! check energy initialize is right"<<std::endl;
            exit(0);
        }
        if(energy_table->jchain()<start_chain) continue;
        if(end_chain!=-1&&energy_table->jchain()>end_chain)continue;
        if(_energy_type != ENERGY_TYPE_UNDEFINED && _energy_type != energy_table->energy_type() ) continue;
    // std::string sequence = pose.sequence();

        std::set<Size> occupied;

        for(Size ires=1; ires<=len; ++ires) {
            EigenXform const & stub1 = pose.stub(ires,energy_table->ichain());

            for(Size jres=(energy_table->energy_type()==TWO_BODY_INTRA_CHAIN?ires+_min_sequence_separation:1); jres<=len; ++jres) {

                if(_ignore_twin_pair && ires == jres) {
                    energy_table->set_pair_score(ires, jres, 0.0);
                    continue;
                }

                // move this to the outer loop???
                // to be safe, each time reset the energy graph to zero
                if(_designable_res.size() != 0 && (_designable_res.count(ires) == 0 || _designable_res.count(jres) == 0) ) {
                    energy_table->set_pair_score(ires, jres, 0.0);
                    continue;
                }

                // clear the rest to zero
                // have to check this in the first place
                // must check this befor is_changed_pair
                if (occupied.count(ires) != 0 || occupied.count(jres) != 0 || scnb_vals(ires-1,0)<_scnb_burial_cutoff || scnb_vals(jres-1,0)<_scnb_burial_cutoff ) {
                    energy_table->set_pair_score(ires, jres, 0.0);
                    continue;
                }

                if( !energy_table->is_changed_pair(ires, jres) ) {
                    continue;
                }

                EigenXform const & stub2 = pose.stub(jres,energy_table->jchain());

                Vec t = stub1.translation() - stub2.translation();

                if(t.squaredNorm() > _squared_dist_cutoff) {
                    energy_table->set_pair_score(ires, jres, 0.0);
                } else {
                    
                    bool find_a_privileged_pair = check_pair_exist(stub1, stub2);

                    if (find_a_privileged_pair) {
                        energy_table->set_pair_score(ires, jres, -1.0);

                        occupied.insert(ires);
                        occupied.insert(jres);

                    } else {
                        energy_table->set_pair_score(ires, jres, 0.0);
                    }
                }
            }

        }
        total_score+= energy_table->weighted_score();
    }
    return total_score;
}

bool PrivilegedPairScoreMethod::check_pair_exist(EigenXform const & stub1, EigenXform const & stub2) const{
    uint64_t hash1(0), hash2(0);
    
    hash1 = _hasher.get_key(stub1.inverse(Eigen::Isometry)*stub2);
    hash2 = _hasher.get_key(stub2.inverse(Eigen::Isometry)*stub1);


    // too far apart
    if( hash1 == 0 || hash2 == 0) {
        return false;
    } else {
        auto got1 = _hash_table.find(hash1);
        auto got2 = _hash_table.find(hash2);
        if(got1 != _hash_table.end() || got2 != _hash_table.end()){
            return true;
        }
    }

    return false;
}

void PrivilegedPairScoreMethod::unserialize_phmap(std::string const & fname){
    #ifdef USE_PHMAP
        phmap::BinaryInputArchive ar_in(fname.c_str());
        _hash_table.clear();

        size_t size;

        ar_in.load( & size );
        _hash_table.reserve( size );

        for ( size_t idx = 0; idx < size; idx++ ) {
            uint64_t k; phmap::flat_hash_set<uint8_t> v;

            ar_in.load( & k );
            v.load( ar_in );

            _hash_table.insert_or_assign( std::move( k ), std::move( v ) );
        }
    #else
        std::cout << "Not implemented yet!!!" << std::endl;
        exit(-1);
    #endif
    assert(_hash_table.size()!=0);
}


} // namespace scoring


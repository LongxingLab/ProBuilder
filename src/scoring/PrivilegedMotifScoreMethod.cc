#include "scoring/PrivilegedMotifScoreMethod.hh"
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


// privileged Motif Score

PrivilegedMotifScoreMethod::PrivilegedMotifScoreMethod(std::vector<std::string> motif_hash_table_paths, Real cart_resl,Real angle_resl):
    BaseScoreMethod(),
    _motif_hash_table_paths(motif_hash_table_paths),
    _hasher(cart_resl,angle_resl,512)
{
    set_score_type(PRIVILEGED_MOTIF);
    _num_hash_tables = _motif_hash_table_paths.size();
    unserialize_phmaps(_motif_hash_table_paths);
    _is_relative_xform_set = false;

    if(_num_hash_tables>20) {
        std::cout << "What the fuck!!!" << std::endl;
        exit(-1);
    }  
}

PrivilegedMotifScoreMethod::~PrivilegedMotifScoreMethod() {}


Real PrivilegedMotifScoreMethod::score(scene::Pose & pose) const {
    Size len = pose.size();
    std::string sequence = pose.sequence();
    EnergyManager & energy_manager = pose.energy_manager();
    std::vector<EnergyOP> energyOPs = energy_manager.get_energy_tables(_score_type);

    assert(energyOPs.size()==1 and energyOPs[0]->energy_type()==ONE_BODY);
    OneBodyEnergyOP privileged_motif_energy = std::static_pointer_cast<OneBodyEnergy>(energyOPs[0]);

    for(Size ires=1; ires<=len; ++ires) {
        if (sequence[ires-1] == 'X') {
            privileged_motif_energy->set_score_per_residue(ires, 0.0);
            continue;
        }
        if (!privileged_motif_energy->is_changed_residue(ires))continue;

        EigenXform const & stub = pose.stub(ires);
        uint64_t key(0);
        if(!_is_relative_xform_set) {
            key = _hasher.get_key(stub);
        } else {
            key = _hasher.get_key(_relative_xform * stub);
        }

        Size token(0);
        for(Size i_table=0; i_table<_num_hash_tables; ++i_table) {
            auto got = _hash_tables[i_table].find(key);
            if (got != _hash_tables[i_table].end()) {
                //
                token = token | ( 1 << i_table );
            }
        }
        privileged_motif_energy->set_score_per_residue(ires, Real(token));
    }

    std::vector<Size> num_hits_per_motif(_num_hash_tables, 0);
    for(Size ires=1; ires<=len; ++ires) {
        Size this_token = privileged_motif_energy->get_score_per_residue(ires); // implicately convert Real to Size
        if( this_token != 0 ) {
            for(Size i_table=0; i_table<_num_hash_tables; ++i_table) {
                num_hits_per_motif[i_table] = num_hits_per_motif[i_table] + (this_token >> i_table) & 1;
            }
        }
    }

    Real score(0.0);
    Real cap = 0.99 / _num_hash_tables;
    for(Size i_table=0; i_table<_num_hash_tables; ++i_table) {
        if(num_hits_per_motif[i_table]>0) {
            score += 1.0 + (1.0/(1+std::exp(-(num_hits_per_motif[i_table]-1)))-0.5)*2*cap;
        }
    }

    return -1.0 * score;
}


void PrivilegedMotifScoreMethod::unserialize_phmaps(std::vector<std::string> const & fnames){
    #ifdef USE_PHMAP

        _hash_tables.resize(_num_hash_tables);

        for(Size i_table=0; i_table<_num_hash_tables; ++i_table) {

            std::cout << "Loading table:    " << fnames[i_table] << std::endl;
            phmap::BinaryInputArchive ar_in(fnames[i_table].c_str());
            _hash_tables[i_table].clear();

            size_t size;

            ar_in.load( & size );
            _hash_tables[i_table].reserve( size );

            for ( size_t idx = 0; idx < size; idx++ ) {
                uint64_t k; phmap::flat_hash_set<int32_t> v;

                ar_in.load( & k );
                v.load( ar_in );

                _hash_tables[i_table].insert_or_assign( std::move( k ), std::move( v ) );
            }
            assert(_hash_tables[i_table].size()!=0);
        }
    #else
        std::cout << "Not implemented yet!!!" << std::endl;
        exit(-1);
    #endif
}

void PrivilegedMotifScoreMethod::set_relative_pos(EigenXform const & xform)
{
    //
    _is_relative_xform_set = true;
    _relative_xform = xform.inverse();
}



/////////////////////////////////////////////////////////
// Privileged Interface Motif Score Method
////////////////////////////////////////////////////////

PrivilegedInterfaceMotifScoreMethod::PrivilegedInterfaceMotifScoreMethod(std::string motif_hash_table_path, Real cart_resl, Real angle_resl):
    BaseScoreMethod(),
    _motif_hash_table_path(motif_hash_table_path),
    _hasher(cart_resl,angle_resl,512)
{
    set_score_type(PRIVILEGED_INTERFACE_MOTIF);
    unserialize_phmap(_motif_hash_table_path); 
}

PrivilegedInterfaceMotifScoreMethod::~PrivilegedInterfaceMotifScoreMethod() {}


Real PrivilegedInterfaceMotifScoreMethod::score(scene::Pose & pose) const
{
    // use the stub of the first residue as the relative xform
    EigenXform relative_xform = pose.stub(1).inverse();

    Size num_hits(0);
    for(Size ires=1; ires<=pose.size(); ++ires) {

        if(_designable_res.size() != 0 && _designable_res.count(ires) == 0) {
            continue;
        }

        EigenXform const & stub = pose.stub(ires, 2);
        uint64_t key = _hasher.get_key(relative_xform * stub);

        if( _hash_table.find(key) != _hash_table.end() ) {
            ++num_hits;
        } 
    }

    if(num_hits>0) {
        return -1 - (1.0/(1+std::exp(-(num_hits-1)))-0.5)*2.0;
    } else {
        return 0.0;
    }

}


void PrivilegedInterfaceMotifScoreMethod::unserialize_phmap(std::string const & fname)
{
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


// Zinc binding score method

MetalCoordinationScoreMethod::MetalCoordinationScoreMethod(std::string metal_hash_table_path, Size coordination_number, Real cart_resl,Real angle_resl):
    BaseScoreMethod(),
    _metal_hash_table_path(metal_hash_table_path),
    _coordination_number(coordination_number),
    _hasher(cart_resl,angle_resl,512)
{
    set_score_type(METAL_COORDINATION);
    unserialize_phmap(_metal_hash_table_path);
}

MetalCoordinationScoreMethod::~MetalCoordinationScoreMethod() {}

Real MetalCoordinationScoreMethod::score(scene::Pose & pose) const {
    std::vector<Size> temp_coordination_res;
    Size temp_token;

    return score_with_tokens(pose, temp_coordination_res, temp_token);
}

Real MetalCoordinationScoreMethod::score_with_tokens(scene::Pose & pose, std::vector<Size> & coordination_res, Size & token) const {

    Size len = pose.size();
    EnergyManager & energy_manager = pose.energy_manager();
    std::vector<EnergyOP> energyOPs = energy_manager.get_energy_tables(_score_type);

    assert(energyOPs.size()==1 and energyOPs[0]->energy_type()==ONE_BODY);
    OneBodyEnergyOP metal_energy = std::static_pointer_cast<OneBodyEnergy>(energyOPs[0]);

    for(Size ires=1; ires<=len; ++ires) {
        if (!metal_energy->is_changed_residue(ires))continue;
        EigenXform const & stub = pose.stub(ires);
        uint64_t key = _hasher.get_key(stub);

        auto got = _hash_table.find(key);
        if(got != _hash_table.end()){
            metal_energy->set_score_per_residue(ires,got->second);
        } else {

            metal_energy->set_score_per_residue(ires,0.0);
        }
    }

    Size num_hits(0);
    Size num_satisfied(0);

    token = 0;
    coordination_res.clear();
    for(Size ires=1; ires<=len; ++ires) {
        Size this_token = metal_energy->get_score_per_residue(ires); // implicately convert Real to Size
        token |= this_token;
        if( this_token != 0 ) {
            coordination_res.push_back(ires);
            ++num_hits;
        }
    }

    for(Size ibit=1; ibit<=_coordination_number; ++ibit) {
        if( token & (1<<ibit) ) ++num_satisfied;
    }

    Size sat_score = std::min(num_hits, num_satisfied);
    Size over_sat_score = std::max(num_hits, num_satisfied) - sat_score;

    return -(sat_score + 0.2 * over_sat_score);
}


void MetalCoordinationScoreMethod::unserialize_phmap(std::string const & fname){
    #ifdef USE_PHMAP

        // key => xform hash
        // val => int32_t format, and each bit represents whether this
        // xform can satisfy the coordination site

        phmap::BinaryInputArchive ar_in(fname.c_str());
        _hash_table.clear();
        _hash_table.load(ar_in);
        assert(_hash_table.size()!=0);
    #else
        std::cout << "Not implemented yet!!!" << std::endl;
        exit(-1);
    #endif
    assert(_hash_table.size()!=0);
}


} // namespace scoring


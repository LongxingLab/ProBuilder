#include "scoring/RifScoreMethod.hh"
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


RifScoreMethod::RifScoreMethod(std::string rif_table_path,Real cart_resl,Real angle_resl):
    BaseScoreMethod(),
    _rif_table_path(rif_table_path),
    _hasher(cart_resl,angle_resl,512.0),
    _is_relative_xform_set(false)
    {
        set_score_type(TARGET_RIF);
        unserialize_phmap(_rif_table_path);
    };

Real RifScoreMethod::score(scene::Pose &pose) const {
    std::vector<Size> dummy_vec;
    return score_with_anchor_seqpos(pose, dummy_vec, false);
}

Real RifScoreMethod::score_with_anchor_seqpos(scene::Pose & pose, std::vector<Size> & anchor_seqpos, bool return_seqpos) const {
    Size len = pose.size();
    std::string sequence = pose.sequence();
    EnergyManager & energy_manager = pose.energy_manager();
    std::vector<EnergyOP> energyOPs = energy_manager.get_energy_tables(_score_type);

    assert(energyOPs.size()==1 and energyOPs[0]->energy_type()==ONE_BODY);
    OneBodyEnergyOP rif_energy = std::static_pointer_cast<OneBodyEnergy>(energyOPs[0]);

    for(Size ires=1; ires<=len; ++ires) {
        if (sequence[ires-1] == 'X') {
            rif_energy->set_score_per_residue(ires, 0.0);
            continue;
        }
        if(_designable_res.size() != 0 && _designable_res.count(ires) == 0 ) {
            rif_energy->set_score_per_residue(ires, 0.0);
            continue;
        }
        if (!rif_energy->is_changed_residue(ires))continue;
        EigenXform const & stub = pose.stub(ires);
        uint64_t key(0);
        if(!_is_relative_xform_set) {
            key = _hasher.get_key(stub);
        } else {
            key = _hasher.get_key(_relative_xform * stub);
        }
        hash_table<uint64_t,Real>::const_iterator got = _rif_table.find(key);
        if(got != _rif_table.end()){
            rif_energy->set_score_per_residue(ires,got->second);
        } else {

            rif_energy->set_score_per_residue(ires,0.0);
        }
    }

    if( return_seqpos ) {
        anchor_seqpos.clear();
        for(Size ires=1; ires<=len; ++ires) {
            if( rif_energy->get_score_per_residue(ires) < 0.0 ) {
                anchor_seqpos.push_back(ires);
            }
        }
    }

    return rif_energy->weighted_score();
}

void RifScoreMethod::unserialize_phmap(std::string fname){
    // #ifdef USE_PHMAP
    //     phmap::BinaryInputArchive ar_in(fname.c_str());
    //     _rif_table.clear();
    //     _rif_table.load(ar_in);
    // #else
    //     std::cout << "Not implemented yet!!! You have to modify the rifdock code to dump rif table.!" << std::endl;
    //     exit(-1);
    // #endif

    phmap::BinaryInputArchive ar_in(fname.c_str());
    _rif_table.clear();
    _rif_table.load(ar_in);
    assert(_rif_table.size()!=0);
}

void RifScoreMethod::set_relative_pos(EigenXform const & xform)
{
    //
    _is_relative_xform_set = true;
    _relative_xform = xform.inverse();
}


} // namespace scoring


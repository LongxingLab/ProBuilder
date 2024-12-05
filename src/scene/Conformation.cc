#include <algorithm>

#include "utils/utils.hh"
#include "basic/types.hh"
#include "scene/Conformation.hh"
#include "utils/random_util.hh"

namespace scene {

using namespace basic;

Conformation::Conformation() : _nres(0), _symmetry("C1"), _num_chains(1), _root_moved(true), _root_index(1), _root_index_snapshot(1),_root_local(EigenXform::Identity()) {}
Conformation::Conformation(Size nres, std::string symmetry):
    _nres(nres),
    _symmetry(symmetry),
    _root_moved(true),
    _root_index(1),
    _root_index_snapshot(1),
    _root_local(EigenXform::Identity())
{
    initialize();
}

Conformation::~Conformation(){}

void
Conformation::initialize()
{
    _root = EigenXform::Identity();
    // _sym_operator = EigenXform::Identity();

    _residues.clear();
    _residues_snapshot.clear();
    for( Size ires=1; ires<=_nres; ++ires )
    {
        _residues.emplace_back(ires);
        _residues_snapshot.emplace_back(ires);
    }

    // _residues_twin.resize(_nres);
    // _residues_twin_snapshot.resize(_nres);

    if(_symmetry.at(0) == 'C') {
        _num_chains = std::stoi(_symmetry.substr(1));

        _residues_twin.resize(_num_chains-1);
        _residues_twin_snapshot.resize(_num_chains-1);

        // calculate the symmetric operator and the weight of each chain
        Size MAX_EFFECTIVE_CHAINS = 3;
        for(Size ichain=1; ichain < _num_chains; ++ichain)
        {
            _residues_twin[ichain-1].resize(_nres);
            _residues_twin_snapshot[ichain-1].resize(_nres);

            EigenXform x(EigenXform::Identity());
            x.translation() = Vec(0.0, 0.0, 0.0);
            x.rotate( AngleAxis( 360.0/_num_chains/180.0*3.14159265359*ichain, Vec(0,0,1) ) );
            _jumps.push_back(x);

            Real weight(0.0);
            if( ichain <= Size(_num_chains/2) && ichain <= MAX_EFFECTIVE_CHAINS ) {
                weight = 1.0;
            }
            if( ichain * 2 == _num_chains && ichain <= MAX_EFFECTIVE_CHAINS ) {
                weight = 0.5;
            }
            _chain_weights.push_back(weight);
        }


        // weights for each chain
    } else if (_symmetry.at(0) == 'D') {
        Size sym_num = stoi(_symmetry.substr(1));
        _num_chains = sym_num * 2;

        _residues_twin.resize(_num_chains-1);
        _residues_twin_snapshot.resize(_num_chains-1);

        Size MAX_EFFECTIVE_CHAINS = 3;
        for(Size ichain=1; ichain < sym_num; ++ichain)
        {
            _residues_twin[ichain-1].resize(_nres);
            _residues_twin_snapshot[ichain-1].resize(_nres);

            EigenXform x(EigenXform::Identity());
            x.translation() = Vec(0.0, 0.0, 0.0);
            x.rotate( AngleAxis( 360.0/sym_num/180.0*3.14159265359*ichain, Vec(0,0,1) ) );
            _jumps.push_back(x);

            Real weight(0.0);
            if( ichain <= Size(sym_num/2) && ichain <= MAX_EFFECTIVE_CHAINS ) {
                weight = 1.0;
            }
            if( ichain * 2 == _num_chains && ichain <= MAX_EFFECTIVE_CHAINS ) {
                weight = 0.5;
            }
            _chain_weights.push_back(weight);
        }

        for(Size ichain=0; ichain<sym_num; ++ichain)
        {
            _residues_twin[ichain+sym_num-1].resize(_nres);
            _residues_twin_snapshot[ichain+sym_num-1].resize(_nres);

            EigenXform x(EigenXform::Identity());
            x.translation() = Vec(0.0, 0.0, 0.0);
            x.rotate( AngleAxis( 360.0/sym_num/180.0*3.14159265359*ichain, Vec(0,0,1) ) );
            x.rotate( AngleAxis( 3.14159265359, Vec(0,1,0) ) );
            _jumps.push_back(x);

            Real weight(0.5);
            // stupid
            //if( ichain  <= MAX_EFFECTIVE_CHAINS || sym_num-ichain <= MAX_EFFECTIVE_CHAINS ) {
            //    weight = 0.5;
            //}
            _chain_weights.push_back(weight);
        }
    } else if(_symmetry.at(0) == 'H') {

        // some definitions for helix
        _num_chains = std::stoi(_symmetry.substr(1))+1;

        if(_num_chains <= 3) {
            std::cout << "What the fuck???" << std::endl;
            exit(0);
        }

        const Real max_shift = 50.0 / (_num_chains - 1) + 10;
        const Real avg_angle = 360.0 / (_num_chains - 1);

        _residues_twin.resize(_num_chains-1);
        _residues_twin_snapshot.resize(_num_chains-1);

        Real shift = utils::random_real(5.0, max_shift);
        Real angle = utils::random_real(avg_angle-25, avg_angle+25);
        if(utils::random_real(-1.0,1.0)<0) {
            angle *= -1;
        }

        // calculate the symmetric operator and the weight of each chain
        for(Size ichain=1; ichain < _num_chains; ++ichain)
        {
            _residues_twin[ichain-1].resize(_nres);
            _residues_twin_snapshot[ichain-1].resize(_nres);

            EigenXform x(EigenXform::Identity());
            x.translation() = Vec(0.0, 0.0, shift * ichain);
            x.rotate( AngleAxis( ichain * angle / 180.0 * 3.14159265359, Vec(0,0,1) ) );
            _jumps.push_back(x);
            _jumps_snapshot.push_back(x);

            Real weight(0.0);
            _chain_weights.push_back(weight);
        }
        
        // n+1
        _chain_weights[0] = 1.0;
        // Hx+1
        _chain_weights[_num_chains-2] = 1.0;
        // Hx
        _chain_weights[_num_chains-3] = 1.0;
        // Hx-1
        _chain_weights[_num_chains-4] = 1.0;


    } else {
        // not implemented error
        std::cout << "Not Implemented Error!!!!!!!!" << std::endl;
        exit(0);
    }

    if(_num_chains>1) random_root();
}

void
Conformation::set_sequence( std::string const & sequence )
{
    for(Size ii = 1; ii <= _nres; ++ii)
        _residues[ii-1].set_name(sequence.at(ii-1));
}

void
Conformation::update_twin_residue(Size ires)
{
    for(Size ichain=0; ichain < _num_chains-1; ++ichain) 
    {
        if (_chain_weights[ichain] != 0) {
            _residues_twin[ichain][ires-1]._N = _jumps[ichain] * _residues[ires-1].xyz(ATOM_N);
            _residues_twin[ichain][ires-1]._CA = _jumps[ichain] * _residues[ires-1].stub(ATOM_CA);
            _residues_twin[ichain][ires-1]._C = _jumps[ichain] * _residues[ires-1].xyz(ATOM_C);
            _residues_twin[ichain][ires-1]._O = _jumps[ichain] * _residues[ires-1].xyz(ATOM_O);
            _residues_twin[ichain][ires-1]._CB = _jumps[ichain] * _residues[ires-1].xyz(ATOM_CB);
        }
    }
}

void
Conformation::update_coordinates()
{
    bool upstream_moved(false), root_res_updated(false);

    EigenXform forward_stub, reverse_stub;
    {
        // update the global xform of root
        root_res_updated = _root_moved || _residues[_root_index-1].is_moved();
        _residues[_root_index-1].update_global_xform_root(_root*_root_local, _root_moved, forward_stub, reverse_stub);

        if(_num_chains>1 && root_res_updated) {
            update_twin_residue(_root_index);
        }
    }

    {
        upstream_moved = false;
        // update the global xform of the residues before root
        for(Size ires=_root_index-1; ires>=1; --ires) {
            upstream_moved = upstream_moved || _residues[ires-1].is_moved() || root_res_updated;
            reverse_stub = _residues[ires-1].update_global_xform_backward(reverse_stub, upstream_moved);

            // update the multiple chain units
            if(_num_chains>1 && upstream_moved) {
                update_twin_residue(ires);
            }
        }
    }

    {
        // update the global xform of the residues after root
        upstream_moved = false;
        for(Size ires=_root_index+1; ires<=_nres; ++ires) {
            upstream_moved = upstream_moved || _residues[ires-1].is_moved() || root_res_updated;
            forward_stub = _residues[ires-1].update_global_xform_forward(forward_stub, upstream_moved);

            // update the multiple chain units
            if(_num_chains>1 && upstream_moved) {
                update_twin_residue(ires);
            }
        }
    }


    // reset
    _root_moved = false;
    
}

void
Conformation::random_root(bool perturb_ori, bool perturb_x, bool perturb_y, bool perturb_z)
{

    // best root ???
    // rand_root(_root, 3.0, 10.0);
    Real x_lb(2.5),x_ub(50),y_lb(2.5),y_ub(50),z_lb(0.0),z_ub(60.0);
    if(_symmetry.at(0) == 'C')
    {
        if(_symmetry == "C1") { x_ub = 60.0; y_ub = 60.0; }
        else                  { x_ub = _num_chains * 15.0; y_ub = _num_chains * 15.0; }
    } 
    else if(_symmetry.at(0) == 'D')
    {
        x_ub = _num_chains / 2 * 6.0;
        y_ub = _num_chains / 2 * 6.0;
    }
    else if(_symmetry.at(0) == 'H')
    {
        x_lb = 15.0;
        x_ub = _num_chains  * 7.0;
        y_lb = 0.0;
        y_ub = 0.0;
        z_lb = 0.0;
        z_ub = 0.0;
    }
    else
    {
        // not implemented error
        std::cout << "Not Implemented Error!!!!!!!!" << std::endl;
        exit(0);
    }

    // EigenXform _root, bool rand_ori, bool rand_x, bool rand_y, bool rand_z, x_cart_lb, x_cart_ub, y_cart_lb, y_cart_ub, z_cart_lb, z_cart_ub,
    utils::rand_xform(_root,perturb_ori, perturb_x, perturb_y, perturb_z, x_lb, x_ub, y_lb, y_ub, z_lb, z_ub);
    
    _root_moved = true;
    update_coordinates();
}

void
Conformation::perturb_root(Real ang_mag, Real trans_mag, Size mode)
{
    Vec v(0.0,0.0,0.0);
    switch (mode){
        case 1:
            v = _root.translation();
            break;
        case 2:
            v = stub(Size(_nres/2+0.5)).translation();
            break;
        case 3:
            for(Size ires=1; ires<=_nres; ++ires) {
                v += stub(ires).translation();
            }
            v /= _nres;
            break;
        case 4:
            v = Vec(0.0,0.0,0.0);
        default:
            std::cerr << "Unknown mode!!!!" << std::endl;
            std::exit(-1);
    }

    // different symmetry
    if(_symmetry.at(0) == 'C')
    {
        EigenXform x = utils::rand_roll(ang_mag, trans_mag);
        EigenXform y(EigenXform::Identity());
        y.translation() = v;
        set_root( y * x * y.inverse() * _root);
    } 
    else if(_symmetry.at(0) == 'D')
    {
        EigenXform x = utils::rand_roll(ang_mag, trans_mag);
        EigenXform y(EigenXform::Identity());
        y.translation() = v;
        set_root( y * x * y.inverse() * _root);
    }
    else if(_symmetry.at(0) == 'H')
    {
        EigenXform x = utils::rand_roll(ang_mag, 0.0);
        EigenXform y(EigenXform::Identity());
        y.translation() = v;
        EigenXform z(EigenXform::Identity());
        Real x_shift = utils::random_real(-trans_mag, trans_mag);
        z.translation() = Vec(x_shift, 0.0, 0.0);
        set_root( z * y * x * y.inverse() * _root);
    }
    else
    {
        // not implemented error
        std::cout << "Not Implemented Error!!!!!!!!" << std::endl;
        exit(0);
    }

    update_coordinates();
}

void
Conformation::random_jump(Real ang_mean, Real trans_mean)
{
    // different symmetry
    if(_symmetry.at(0) == 'C')
    {
        // not implemented error
        std::cout << "What the fuck???" << std::endl;
        exit(0);
    } 
    else if(_symmetry.at(0) == 'D')
    {
        // not implemented error
        std::cout << "What the fuck???" << std::endl;
        exit(0);
    }
    else if(_symmetry.at(0) == 'H')
    {
        const Real max_shift = 50.0 / (_num_chains - 1)+10.0;
        const Real avg_angle = 360.0 / (_num_chains - 1);

        Real shift = utils::random_real(5.0, max_shift);
        Real angle = utils::random_real(avg_angle-25, avg_angle+25);

        if( ang_mean != -1.0 ) {
            angle = utils::random_real(ang_mean-25, ang_mean+25);
        }
        if( trans_mean != -1.0 ) {
            shift = utils::random_real(trans_mean-10, trans_mean+10);
        }

        if(utils::random_real(-1.0,1.0)<0) {
            angle *= -1;
        }

        // calculate the symmetric operator and the weight of each chain
        for(Size ichain=1; ichain < _num_chains; ++ichain)
        {
            EigenXform x(EigenXform::Identity());
            x.translation() = Vec(0.0, 0.0, shift * ichain);
            x.rotate( AngleAxis( ichain * angle / 180.0 * 3.14159265359, Vec(0,0,1) ) );
            _jumps[ichain-1] = x;
        }
    }
    else
    {
        // not implemented error
        std::cout << "Not Implemented Error!!!!!!!!" << std::endl;
        exit(0);
    }

    // explicitely set the root moved flag
    _root_moved = true;
    update_coordinates();
}

void
Conformation::perturb_jump(Real ang_mag, Real trans_mag)
{
    // different symmetry
    if(_symmetry.at(0) == 'C')
    {
        // not implemented error
        std::cout << "What the fuck???" << std::endl;
        exit(0);
    } 
    else if(_symmetry.at(0) == 'D')
    {
        // not implemented error
        std::cout << "What the fuck???" << std::endl;
        exit(0);
    }
    else if(_symmetry.at(0) == 'H')
    {
        AngleAxis angle_axis(_jumps[0].rotation());
        Real rotation_angle = angle_axis.angle() / 3.14159265359 * 180;
        Real z_shift = _jumps[0].translation()[2];

        z_shift = utils::random_real(-trans_mag, trans_mag) + z_shift;
        rotation_angle = utils::random_real(-ang_mag, ang_mag) + rotation_angle;

        // calculate the symmetric operator and the weight of each chain
        for(Size ichain=1; ichain < _num_chains; ++ichain)
        {
            EigenXform x(EigenXform::Identity());
            x.translation() = Vec(0.0, 0.0, z_shift * ichain);
            x.rotate( AngleAxis( ichain * rotation_angle / 180.0 * 3.14159265359, Vec(0,0,1) ) );
            _jumps[ichain-1] = x;
        }
    }
    else
    {
        // not implemented error
        std::cout << "Not Implemented Error!!!!!!!!" << std::endl;
        exit(0);
    }

    // explicitely set the root moved flag
    _root_moved = true;
    update_coordinates();
}

void
Conformation::snapshot()
{
    _root_index_snapshot = _root_index;
    _root_snapshot = _root;
    // _sym_operator_snapshot = _sym_operator;

    std::copy(_residues.begin(), _residues.end(), _residues_snapshot.begin());
    for(Size ichain=0; ichain<_num_chains-1; ++ichain)
    {
        if( _chain_weights[ichain] != 0 ) {
            std::copy(_residues_twin[ichain].begin(), _residues_twin[ichain].end(), _residues_twin_snapshot[ichain].begin());
        }
    }

    if(_symmetry.at(0) == 'H') {
        std::copy(_jumps.begin(), _jumps.end(), _jumps_snapshot.begin());
    }
	// std::copy(_residues_twin.begin(), _residues_twin.end(), _residues_twin_snapshot.begin());
}
void
Conformation::rollback()
{
    _root_index = _root_index_snapshot;
    _root = _root_snapshot;
    // _sym_operator = _sym_operator_snapshot;

    std::copy(_residues_snapshot.begin(), _residues_snapshot.end(), _residues.begin());
    for(Size ichain=0; ichain<_num_chains-1; ++ichain)
    {
        if( _chain_weights[ichain] != 0 ) {
            std::copy(_residues_twin_snapshot[ichain].begin(), _residues_twin_snapshot[ichain].end(), _residues_twin[ichain].begin());
        }
    }
	    // std::copy(_residues_twin_snapshot.begin(), _residues_twin_snapshot.end(), _residues_twin.begin());    	
    if(_symmetry.at(0) == 'H') {
        std::copy(_jumps_snapshot.begin(), _jumps_snapshot.end(), _jumps.begin());
    }
}

void
Conformation::set_stub(Size ires, ATOM_TYPE atom_name, EigenXform const & stub)
{
    _residues[ires-1].set_stub(stub, atom_name);
}


// TODO:: stupid, can not reture a const reference
// unknown behavios will happen
// especially when returning a product of two xforms
// that will generate a temporary variable
// returning a const reference of a temporary variable will cause unknown behaviors
// FIX it!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
const EigenXform
Conformation::stub(Size ires, Size ii_unit) const {

    if(ii_unit == 1) {
        return _residues[ires-1].stub(ATOM_CA);
    } 

    if ( _num_chains > 1) {
    	if ( _chain_weights[ii_unit-2] != 0 ) {
    		return _residues_twin[ii_unit-2][ires-1]._CA; 
    	} else {
    		// EigenXform r(_residues_twin[ires-1]._CA);
    		// for(Size ii=0;ii<ii_unit-2;++ii) r = _sym_operator * r;
    		return _jumps[ii_unit-2] * _residues[ires-1].stub(ATOM_CA);
    	}
    } else {
    	std::cerr << "This is a single chain protein!!!!" << std::endl;
        std::exit(-1);
    }
}

// TODO:: stupid, can not reture a const reference
// unknown behavios will happen
// especially when returning a product of two xforms
// that will generate a temporary variable
// returning a const reference of a temporary variable will cause unknown behaviors
// FIX it!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
const EigenXform Conformation::stub(Size ires, Size ii_unit, ATOM_TYPE atom_name) const {
    if(ii_unit==1){
        return _residues[ires-1].stub(atom_name);
    } else {
        // EigenXform r(_residues[ires-1].stub(atom_name));
        // for(Size ii=0;ii<ii_unit-1;++ii) r = _sym_operator * r;
        // return r;
        return _jumps[ii_unit-2] * _residues[ires-1].stub(atom_name);
    }
}

Vec
Conformation::xyz(Size ires, ATOM_TYPE atom_name, Size ii_unit) const {

	if(ii_unit == 1) {
		return _residues[ires-1].xyz(atom_name);
    }

    if ( _num_chains > 1 && ii_unit>1 && ii_unit<=_num_chains) {
    	// if (ii_unit == 2) {
    	// 	switch (atom_name) {
        //         case ATOM_N:
        //             return _residues_twin[ires-1]._N;
        //         case ATOM_CA:
        //             return _residues_twin[ires-1]._CA.translation();
        //         case ATOM_C:
        //             return _residues_twin[ires-1]._C;
        //         case ATOM_O:
        //             return _residues_twin[ires-1]._O;
        //         case ATOM_CB:
        //             return _residues_twin[ires-1]._CB;
        //         default:
    	//             std::cerr << "Unknown atom type!!!!" << std::endl;
        //         	exit(-1);
        //     }
    	// } else {
    	// 	Vec r(_residues[ires-1].xyz(atom_name));
        //     for(Size ii=0;ii<ii_unit-1;++ii) r = _sym_operator * r;
    	// 	return r;
    	// }

        if ( _chain_weights[ii_unit-2] != 0 ) {
            switch (atom_name) {
                case ATOM_N:
                    return _residues_twin[ii_unit-2][ires-1]._N;
                case ATOM_CA:
                    return _residues_twin[ii_unit-2][ires-1]._CA.translation();
                case ATOM_C:
                    return _residues_twin[ii_unit-2][ires-1]._C;
                case ATOM_O:
                    return _residues_twin[ii_unit-2][ires-1]._O;
                case ATOM_CB:
                    return _residues_twin[ii_unit-2][ires-1]._CB;
                default:
                    std::cerr << "Unknown atom type!!!!" << std::endl;
                exit(-1);
            }
        } else {
            // EigenXform r(_residues_twin[ires-1]._CA);
            // for(Size ii=0;ii<ii_unit-2;++ii) r = _sym_operator * r;
            return _jumps[ii_unit-2] * _residues[ires-1].xyz(atom_name);
        }
    } else {
    	std::cerr << "chain number:" <<ii_unit <<"is illegal,chain number should be between 1 and "<<_num_chains<< std::endl;
        std::exit(-1);
    }
}


void Conformation::extend_conformation(Size Nter_extension, Size Cter_extension)
{
    // N terminus extension
    for(Size i_tmp=0; i_tmp<Nter_extension; ++i_tmp) {
        _residues.emplace(_residues.begin(), i_tmp+1);
        _residues_snapshot.emplace(_residues_snapshot.begin(), i_tmp+1);
    }
    _nres += Nter_extension;

    // C terminus extension
    for(Size i_tmp=0; i_tmp<Cter_extension; ++i_tmp) {
        _residues.emplace_back(_nres+i_tmp+1);
        _residues_snapshot.emplace_back(_nres+i_tmp+1);
    }
    _nres += Cter_extension;

    // reorder the residues
    for(Size i_tmp=1; i_tmp<=_nres; ++i_tmp) {
        _residues[i_tmp-1].set_seqpos(i_tmp);
        _residues_snapshot[i_tmp-1].set_seqpos(i_tmp);
    }

    // reset root
    _root_index += Nter_extension;
    _root_index_snapshot = _root_index;
    _root_moved = true;

    // twin residues
    for(Size ichain=1; ichain < _num_chains; ++ichain)
    {
        _residues_twin[ichain-1].resize(_nres);
        _residues_twin_snapshot[ichain-1].resize(_nres);
    }

    update_coordinates();

}


}

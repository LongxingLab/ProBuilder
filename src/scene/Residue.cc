#include "scene/Residue.hh"
#include "basic/macros.hh"

#include <cmath>

namespace scene {

using namespace basic;

Residue::Residue(int ires):
    _ires(ires),
    _name('V'),
    _natoms(5),
    _psudo_cb(true),
    _is_omega_movable(true),
    _is_omega_moved(true),
    _is_phi_movable(true),
    _is_phi_moved(true),
    _is_psi_movable(true),
    _is_psi_moved(true),
    _omega(3.14),
    _phi(-1.04),
    _psi(-0.87),
    _CX_global_xform(EigenXform::Identity()),
    _N_local_xform(EigenXform::Identity()),
    _N_global_xform(EigenXform::Identity()),
    _CA_local_xform(EigenXform::Identity()),
    _CA_global_xform(EigenXform::Identity()),
    _C_local_xform(EigenXform::Identity()),
    _C_global_xform(EigenXform::Identity()),
    _O_local_pos(IDL_O_X, IDL_O_Y, IDL_O_Z),
    _O_global_pos(0.0,0.0,0.0),
    _CB_local_pos(PSEUDO_CB_X, PSEUDO_CB_Y, PSEUDO_CB_Z),
    _CB_global_pos(0.0,0.0,0.0),
    _H_local_pos(IDL_H_X, IDL_H_Y, IDL_H_Z),
    _H_global_pos(0.0,0.0,0.0)
{
    initialize();
}

Residue::~Residue(){}

void
Residue::initialize()
{
    // N
    _N_local_xform(0,0) = -idl_cosW_N;
    _N_local_xform(0,1) = -idl_sinW_N;
    _N_local_xform(0,2) = 0.0;
    _N_local_xform(0,3) = idl_R_N;
    // CA
    _CA_local_xform(0,0) = -idl_cosW_CA;
    _CA_local_xform(0,1) = -idl_sinW_CA;
    _CA_local_xform(0,2) = 0.0;
    _CA_local_xform(0,3) = idl_R_CA;
    // C
    _C_local_xform(0,0) = -idl_cosW_C;
    _C_local_xform(0,1) = -idl_sinW_C;
    _C_local_xform(0,2) = 0.0;
    _C_local_xform(0,3) = idl_R_C;


}

void
Residue::update_local_xform(){

    // update_omega local xform
    if( _is_omega_moved ) {
        Real sin_omega = sin(_omega);
        Real cos_omega = cos(_omega);
        _N_local_xform(1,0) = idl_sinW_N*cos_omega;
        _N_local_xform(1,1) = -idl_cosW_N*cos_omega;
        _N_local_xform(1,2) = -sin_omega;
        _N_local_xform(2,0) = idl_sinW_N*sin_omega;
        _N_local_xform(2,1) = -idl_cosW_N*sin_omega;
        _N_local_xform(2,2) = cos_omega;
    }

    // update phi local xform
    if( _is_phi_moved ) {
        Real sin_phi = sin(_phi);
        Real cos_phi = cos(_phi);       
        _CA_local_xform(1,0) = idl_sinW_CA*cos_phi;
        _CA_local_xform(1,1) = -idl_cosW_CA*cos_phi;
        _CA_local_xform(1,2) = -sin_phi;
        _CA_local_xform(2,0) = idl_sinW_CA*sin_phi;
        _CA_local_xform(2,1) = -idl_cosW_CA*sin_phi;
        _CA_local_xform(2,2) = cos_phi;
    }

    // update psi local xform
    if( _is_psi_moved ) {
        Real sin_psi = sin(_psi);
        Real cos_psi = cos(_psi);
        _C_local_xform(1,0) = idl_sinW_C*cos_psi;
        _C_local_xform(1,1) = -idl_cosW_C*cos_psi;
        _C_local_xform(1,2) = -sin_psi;
        _C_local_xform(2,0) = idl_sinW_C*sin_psi;
        _C_local_xform(2,1) = -idl_cosW_C*sin_psi;
        _C_local_xform(2,2) = cos_psi;
    }
}

void 
Residue::update_local_xform_from_global_xform(bool update_local_N, bool update_local_CA, bool update_local_C)
{
    if(update_local_N){
        //
        _N_local_xform = _CX_global_xform.inverse(Eigen::Isometry) * _N_global_xform;
    }
    if(update_local_CA){
        //
        _CA_local_xform = _N_global_xform.inverse(Eigen::Isometry) * _CA_global_xform;
    }
    if(update_local_C){
        //
        _C_local_xform = _CA_global_xform.inverse(Eigen::Isometry) * _C_global_xform;
    }
}

void
Residue::update_peripherial_atoms(bool update_H, bool update_CB, bool update_O)
{
    // update the O/CB/H coordinates
    //if( update_H )
    //    _H_global_pos = _N_global_xform * _H_local_pos;
    if( update_CB )
        _CB_global_pos = _CA_global_xform * _CB_local_pos;
    if( update_O )
        _O_global_pos = _C_global_xform * _O_local_pos;
}

// update from CA
void
Residue::update_global_xform_root(EigenXform const & root, bool root_moved, EigenXform & forward_stub, EigenXform & reverse_stub){
    //
    if( is_moved() || root_moved ) {

        // update local xform
        update_local_xform();

        _CA_global_xform = root;
        _C_global_xform = _CA_global_xform * _C_local_xform;

        _N_global_xform = _CA_global_xform * _CA_local_xform.inverse(Eigen::Isometry);
        _CX_global_xform = _N_global_xform * _N_local_xform.inverse(Eigen::Isometry);    

        // update O/CB/H
        update_peripherial_atoms();

        // reset
        _is_omega_moved = false;
        _is_phi_moved = false;
        _is_psi_moved = false;
    }
    reverse_stub = _CX_global_xform;
    forward_stub = _C_global_xform;
}

void
Residue::update_global_by_ca_stub(){

    _C_global_xform = _CA_global_xform * _C_local_xform;

    _N_global_xform = _CA_global_xform * _CA_local_xform.inverse(Eigen::Isometry);
    _CX_global_xform = _N_global_xform * _N_local_xform.inverse(Eigen::Isometry);    

    // update O/CB/H
    update_peripherial_atoms();

}

// update from C
const EigenXform &
Residue::update_global_xform_backward(EigenXform const & upper_stub, bool upstream_moved){
    //

    if( is_moved() || upstream_moved ) {
        // update local xform
        update_local_xform();

        // 
        _C_global_xform = upper_stub;
        _CA_global_xform = _C_global_xform * _C_local_xform.inverse(Eigen::Isometry);
        _N_global_xform = _CA_global_xform * _CA_local_xform.inverse(Eigen::Isometry);
        _CX_global_xform =  _N_global_xform * _N_local_xform.inverse(Eigen::Isometry);

        // update O/CB/H
        update_peripherial_atoms();

        //reset
        _is_omega_moved = false;
        _is_phi_moved = false;
        _is_psi_moved = false;
    }
    return _CX_global_xform;
}

// update from N
const EigenXform &
Residue::update_global_xform_forward(EigenXform const & upper_stub, bool upstream_moved){
    if( is_moved() || upstream_moved) {
        // update local xform
        update_local_xform();

        //
        _CX_global_xform = upper_stub;
        _N_global_xform = _CX_global_xform * _N_local_xform;
        _CA_global_xform = _N_global_xform * _CA_local_xform;
        _C_global_xform = _CA_global_xform * _C_local_xform;

        // update O/CB/H
        update_peripherial_atoms();

        // reset
        _is_omega_moved = false;
        _is_phi_moved = false;
        _is_psi_moved = false;
    }
    return _C_global_xform;
}

}

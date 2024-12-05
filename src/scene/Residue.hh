#ifndef INCLUDED_scene_Residue_hh
#define INCLUDED_scene_Residue_hh

#include "basic/types.hh"
#include "basic/macros.hh"

#include <string>
#include <iostream>
#include <vector>

namespace scene {

using namespace basic;

typedef struct AtomInfo{
    std::string atom_name;
    std::string res_name;
    Vec local_postion;
    Vec global_postion;
    std::string element;
    ATOM_TYPE atom_type;
    Size ires;
};

class Residue
{
    public:
        Residue(Size ires);
        ~Residue();


    inline char name() const {
        return _name;
    }

    inline Size natoms(){
        return _natoms;
    }

    inline void set_name(char const & name) {
        _name = name;
        if (_name == 'G' ) {
            if(_psudo_cb)_natoms=4;
            _CB_local_pos(0,0) = 0.0;
            _CB_local_pos(1,0) = 0.0;
            _CB_local_pos(2,0) = 0.0;
        } else {
            if(_psudo_cb)_natoms=5;
            _CB_local_pos(0,0) = PSEUDO_CB_X;
            _CB_local_pos(1,0) = PSEUDO_CB_Y;
            _CB_local_pos(2,0) = PSEUDO_CB_Z;
        }
    }
    inline Size seqpos() const {
        return _ires;
    }
    inline void set_seqpos(Size ires) {
        _ires = ires;

        for(Size i_tmp=0; i_tmp<_atom_infos.size(); ++i_tmp) {
            _atom_infos[i_tmp].ires = _ires;
        }
    }
    // torsion angles
    inline Real omega() const {
        return _omega;
    }
    inline void set_omega(Real value) {
        if(_is_omega_movable){
            _is_omega_moved = true;
            _omega = value;
        }
    }
    inline Real phi() const {
        return _phi;
    }
    inline void set_phi(Real value) {
        if(_is_phi_movable){
            _is_phi_moved = true;
            _phi = value;
        }
    }
    inline Real psi() const {
        return _psi;
    }
    inline void set_psi(Real value) {
        if(_is_psi_movable) {
            _is_psi_moved = true;
            _psi = value;
        }
    }
    inline std::vector<AtomInfo> atom_info(){
        return _atom_infos;
    }
    inline void push_atom_info(AtomInfo atom_info){
        _psudo_cb = false;
        _atom_infos.push_back(atom_info);
        //  only keep N CA C
        _natoms = _atom_infos.size()+3;
    }
    inline void freeze_torsions( bool freeze_omega = true, bool freeze_phi = true, bool freeze_psi = true) {
        if( freeze_omega ){
            _is_omega_movable = false;
            _is_omega_moved = false;
        }
        if( freeze_phi ) {
            _is_phi_movable = false;
            _is_phi_moved = false;
        }
        if( freeze_psi ) {
            _is_psi_movable = false;
            _is_psi_moved = false;
        }
    }

    // could the status be changed from outside???
    inline bool is_moved() const {
        return _is_omega_moved || _is_phi_moved || _is_psi_moved;
    }
    
    void update_local_xform();
    void update_global_by_ca_stub();
    void update_local_xform_from_global_xform(bool update_local_N=true, bool update_local_CA=true, bool update_local_C=true);

    // three different ways to update the global xform
    void update_global_xform_root(EigenXform const & root, bool root_moved, EigenXform & forward_stub, EigenXform & reverse_stub);
    const EigenXform & update_global_xform_backward(EigenXform const & upper_stub, bool upstream_moved);
    const EigenXform & update_global_xform_forward(EigenXform const & upper_stub, bool upstream_moved);

    // O/CB/H
    void update_peripherial_atoms(bool update_H=true, bool update_CB=true, bool update_O=true);

    // what's the point of setting the global xform of CX ?????????
    // Should I allow this ??????????????
    inline void set_stub(EigenXform const & stub, ATOM_TYPE atom=ATOM_CA) {
        switch(atom) {
            case ATOM_CX:
                _CX_global_xform = stub;
                break;
            case ATOM_N:
                _N_global_xform = stub;
                update_peripherial_atoms(true, false, false);
                break;
            case ATOM_CA:
                _CA_global_xform = stub;
                update_peripherial_atoms(false, true, true);
                break;
            case ATOM_C:
                _C_global_xform = stub;
                update_peripherial_atoms(false, false, true);
                break;
            default:
                std::cerr << "Unknown atom type!!!!" << std::endl;
                std::exit(-1);
        }
    }

    // 
    inline const EigenXform stub(ATOM_TYPE atom=ATOM_CA) const {
        if ( atom == ATOM_CX ) {
            return _CX_global_xform;
        } else if (atom == ATOM_N) {
            return _N_global_xform;
        } else if ( atom == ATOM_CA ) {
            return _CA_global_xform;
        } else if ( atom == ATOM_C ) {
            return _C_global_xform;
        } else {
            std::cerr << "Unknown atom type!!!!" << std::endl;
            std::exit(-1);
        }
    }

    inline Vec xyz(ATOM_TYPE atom) const {
        if (atom == ATOM_N) {
            return _N_global_xform.translation();
        } else if ( atom == ATOM_CA ) {
            return _CA_global_xform.translation();
        } else if ( atom == ATOM_C ) {
            return _C_global_xform.translation();
        } else if ( atom == ATOM_O ) {
            return _O_global_pos;
        } else if ( atom == ATOM_CB ) {
            return _CB_global_pos;
        } else if ( atom == ATOM_H) {
            return _H_global_pos;
        } else {
            std::cerr << "Unknown atom type!!!!" << std::endl;
            std::exit(-1);
        }
    }
    inline Vec xyz(Size iatom) const {
        assert(0<iatom<=_natoms);
        if(iatom<=3){
            if(iatom==1)return _N_global_xform.translation();
            if(iatom==2)return _CA_global_xform.translation();
            if(iatom==3)return _C_global_xform.translation();
        }
        else{
            if(_psudo_cb){
                if(iatom==4)return _O_global_pos;
                else if(iatom==5)return _CB_global_pos;
                else {
                    std::cerr<<"out range atom for ideal residue"<<std::endl;
                    exit(1);
                }
            }else{
                return _CA_global_xform*_atom_infos[iatom-4].local_postion;
            }
        }
    }
    inline Vec xyz(std::string atom_name) const {
        if(atom_name=="N")return _N_global_xform.translation();
        else if(atom_name=="CA")return _CA_global_xform.translation();
        else if(atom_name=="C")return _C_global_xform.translation();
        else{
            for(auto & atom_info : _atom_infos){
                if(atom_info.atom_name==atom_name)return _CA_global_xform*atom_info.local_postion;
            }
            std::cerr<<"try to get xyz of a none exist atom "<<atom_name<<"on residue "<<_ires<<std::endl;
            exit(0);
        }
    }

    inline ATOM_TYPE atom_type(Size iatom) const {
        assert(0<iatom<=_natoms);
        if(iatom<=3){
            if(iatom==1)return ATOM_N;
            if(iatom==2)return ATOM_C;
            if(iatom==3)return ATOM_C;
        }
        else{
            if(_psudo_cb){
                if(iatom==4)return ATOM_O;
                else if(iatom==5)return ATOM_CB;
                else {
                    std::cerr<<"out range atom for ideal residue"<<std::endl;
                    exit(1);
                }
            }else{
                return _atom_infos[iatom-4].atom_type;
            }
        }
    }
    inline std::string atom_name(Size iatom) const {
        assert(0<iatom<=_natoms);
        if(iatom<=3){
            if(iatom==1)return "N";
            if(iatom==2)return "CA";
            if(iatom==3)return "C";
        }
        else{
            if(_psudo_cb){
                if(iatom==4)return "O";
                else if(iatom==5)return "CB";
                else {
                    std::cerr<<"out range atom for ideal residue"<<std::endl;
                    exit(1);
                }
            }else{
                return _atom_infos[iatom-4].atom_name;
            }
        }
    }
    private:
        void initialize();
        
    private:
        Size _ires;
        char _name;
        bool _is_omega_movable;
        bool _is_omega_moved;
        bool _is_phi_movable;
        bool _is_phi_moved;
        bool _is_psi_movable;
        bool _is_psi_moved;

        Real _omega;
        Real _phi;
        Real _psi;
        bool _psudo_cb;
        Size _natoms;
        EigenXform _CX_global_xform;
        EigenXform _N_local_xform;
        EigenXform _N_global_xform;
        EigenXform _CA_local_xform;
        EigenXform _CA_global_xform;
        EigenXform _C_local_xform;
        EigenXform _C_global_xform;
        Vec _O_local_pos;
        Vec _O_global_pos;
        Vec _CB_local_pos;
        Vec _CB_global_pos;
        Vec _H_local_pos;
        Vec _H_global_pos;
        std::vector<AtomInfo> _atom_infos;
};

}

#endif

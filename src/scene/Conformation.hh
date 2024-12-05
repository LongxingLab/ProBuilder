#ifndef INCLUDED_scene_Conformation_hh
#define INCLUDED_scene_Conformation_hh

#include "basic/types.hh"
#include "scene/Residue.hh"

#include <vector>
#include <memory>

namespace scene {

using namespace basic;

struct ResidueLite {
    Vec _N;
    EigenXform _CA;
    Vec _C;
    Vec _O;
    Vec _CB;

    ResidueLite():_N(0,0,0),_CA(EigenXform::Identity()),_C(0,0,0),_O(0,0,0),_CB(0,0,0){};
};

class Conformation
{
    public:
        Conformation();
        Conformation(Size nres, std::string symmetry);
        ~Conformation();

        inline Size size() const {
            return _nres;
        }
        inline void size(Size nres) {
            _nres = nres;
        }
        inline Size num_chains () const {
            return _num_chains;
        }

        void set_sequence( std::string const & sequence);

        // torsion angles
        inline Real omega(Size ires) const {
            return _residues[ires-1].omega();
        }
        inline void set_omega(Size ires, Real value) {
            _residues[ires-1].set_omega(value);
        }
        inline Real phi(Size ires) const {
            return _residues[ires-1].phi();
        }
        inline void set_phi(Size ires, Real value) {
            _residues[ires-1].set_phi(value);
        }
        inline Real psi(Size ires) const {
            return _residues[ires-1].psi();
        }
        inline void set_psi(Size ires, Real value) {
            _residues[ires-1].set_psi(value);
        }

        inline void freeze_torsions( Size ires, bool freeze_omega, bool freeze_phi, bool freeze_psi ) {
            _residues[ires-1].freeze_torsions(freeze_omega, freeze_phi, freeze_psi);
        }
        inline void update_local_xform_from_global_xform(Size ires, bool update_local_N, bool update_local_CA, bool update_local_C) {
            _residues[ires-1].update_local_xform_from_global_xform(update_local_N, update_local_CA, update_local_C);
        }

        inline const EigenXform & jump(Size ichain) const {
            std::cout << "Why you want to get jump?? Let me know!!!" << std::endl;
            exit(1);
            return _jumps[ichain-2];
        }

        inline const EigenXform & root () const {
            return _root;
        }

        inline void set_root_index(Size residue_index, bool move_to_current_root=false) {
            _root_index = residue_index;
            if( !move_to_current_root ) {
                _root = _residues[residue_index-1].stub(ATOM_CA);
            } else {
                _root_moved = true;
            }
        }
        inline Size root_index() const {
            return _root_index;
        }

        inline void set_root(EigenXform const & root) {
            _root = root;
            _root_moved = true;
        }
        inline Vec center_vec(Size first_residue=1,Size last_residue=-1){
            assert((first_residue>0 && first_residue<=_nres)&&(last_residue==-1 || (last_residue >0 && last_residue<=_nres) ));
            Vec v(0.0,0.0,0.0);
            if(last_residue==-1)last_residue=_nres;
            for(Size ires=first_residue; ires<=last_residue; ++ires) {
                v += stub(ires).translation();
            }
            // a bug wasted 3h to find
            // v /= _nres;
            v /= (last_residue-first_residue+1);
            return v;
        }
        inline Size center_residue(){
            Vec v = center_vec();
            Real min_dis=99999,best_index;
            for(Size ires=1; ires<=_nres; ++ires) {
                if((stub(ires).translation()-v).norm()<min_dis)min_dis=(stub(ires).translation()-v).norm(),best_index=ires;
            }
            return best_index;
        };
        inline Real max_radius(){
            Vec v = center_vec();
            Size max_dis=0;
            for(Size ires=1; ires<=_nres; ++ires) {
                if((stub(ires).translation()-v).norm()>max_dis)max_dis=(stub(ires).translation()-v).norm();
            }
            return max_dis;
        }
        inline Residue & residue(Size ires){
            assert(ires>0&&ires<=_nres);
            return _residues[ires-1];
        }
        inline std::vector<AtomInfo>  atom_info (Size ires){
            return _residues[ires-1].atom_info();
        }
        inline void push_atom_info (Size ires,AtomInfo atom_info){
            _residues[ires-1].push_atom_info(atom_info);
        }
        // this funtion split the initial root wich is the root reisude CA stub to somewhere else
        // add this only need to be call only once since _root_local means tranform from root to ca stub !!!!! 
        inline void split_root_from_residue(EigenXform newroot){
            _root_local = newroot.inverse() * _root;
            _root_moved = true;
            _root  = newroot;
        }
        void random_root(bool perturb_ori=true, bool perturb_x=true, bool perturb_y=true, bool perturb_z=true);
        void perturb_root(Real ang_mag=10.0, Real trans_mag=5.0, Size mode=1);
        void random_jump(Real ang_mean=-1.0, Real trans_mean=-1.0);
        void perturb_jump(Real ang_mag=10.0, Real trans_mag=5.0);
        std::vector<EigenXform> jumps() {
            return _jumps;
        }

        void initialize();
        void update_coordinates();

        const EigenXform stub(Size ires, Size ii_unit=1) const;
        const EigenXform stub(Size ires, Size ii_unit, ATOM_TYPE atom_name) const;
        void set_stub(Size ires, ATOM_TYPE atom_name, EigenXform const & stub);
        void snapshot();
        void rollback();

        Vec xyz(Size ires, ATOM_TYPE atom_name, Size ii_unit=1) const;

        inline std::vector<Real> const & chain_weights() const { return _chain_weights; }

        void update_twin_residue(Size ires);

        void extend_conformation(Size Nter_extension, Size Cter_extension);

    private:
        Size _nres;
        std::string _symmetry;
        Size _num_chains;
        bool _root_moved;
        Size _root_index;
        Size _root_index_snapshot;
        EigenXform _root;
        // this update is to split identity of root and ca stub of root residue to represent a peptide which root is in the axis
        EigenXform _root_local;
        EigenXform _root_snapshot;
        std::vector<EigenXform> _jumps;
        std::vector<EigenXform> _jumps_snapshot;
        std::vector<Real> _chain_weights;
        // EigenXform _sym_operator;
        // EigenXform _sym_operator_snapshot;
        std::vector<Residue> _residues;
        std::vector<Residue> _residues_snapshot;

        std::vector<std::vector<ResidueLite> > _residues_twin;
        std::vector<std::vector<ResidueLite> > _residues_twin_snapshot;

};

typedef std::shared_ptr<Conformation> ConformationOP;

}

#endif

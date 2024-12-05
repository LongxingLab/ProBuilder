#ifndef INCLUDED_scene_Pose_hh
#define INCLUDED_scene_Pose_hh

#include "basic/types.hh"
#include "scene/Residue.hh"
#include "scene/Conformation.hh"
#include "scoring/Energy.hh"

#include <gzip/gzstream.hh>
#include <memory>
#include <set>

namespace scene {

using namespace basic;

class Pose
{
    public:
        Pose(Size nres, Size num_repeats=1, std::string symmetry="C1");
        Pose(const Pose& pose,bool copy_enery=true,bool copy_trajctory=false);
        Pose(std::string pdb_name,char chain=' ' ,std::string symmetry="C1",Size repeat_num=1);
        ~Pose();
        void initialize();
        inline Size nres() const {
            return _nres;
        }
        inline std::string symmetry() const {
            return _symmetry;
        }
        inline Size num_repeats() const {
            return _num_repeats;
        }
        inline Size num_chains() const {
            return _num_chains;
        }
        inline Size size() const {
            return _conformation.size();
        }
        inline Size chains() const {
            return _conformation.num_chains();
        }
        // torsion angles
        inline Real omega(Size ires) const {
            return _conformation.omega(ires);
        }
        inline void set_omega(Size ires, Real value) {
            _conformation.set_omega(ires, value);
        }
        inline Real phi(Size ires) const {
            return _conformation.phi(ires);
        }
        inline void set_phi(Size ires, Real value) {
            _conformation.set_phi(ires, value);
        }
        inline Real psi(Size ires) const {
            return _conformation.psi(ires);
        }
        inline void set_psi(Size ires, Real value) {
            _conformation.set_psi(ires, value);
        }

        inline void freeze_torsions( Size ires, bool freeze_omega, bool freeze_phi, bool freeze_psi ) {
            _conformation.freeze_torsions(ires, freeze_omega, freeze_phi, freeze_psi);
        }

        void freeze_segment(Size start_pos, Size end_pos);

        // sequence and ss getters and setters
        inline std::string sequence() const {
            return _sequence;
        }
        inline void set_sequence(std::string const & seq, Size ires=1) {
            assert(seq.length() + ires - 1 <= size());
            _sequence.replace(ires-1, seq.length(), seq);
            _conformation.set_sequence(_sequence);
        }
        inline std::string dssp() const {
            return _ss;
        }
        inline void set_dssp(std::string const ss) {
            assert(ss.length() == _nres);
            _ss = ss;
            if(ss.size()==0){
                std::string seq;
                for(auto s:ss)seq.push_back(s=='H'?'V':(s=='L'?'G':'A'));
            }
        }
        inline std::vector<std::string> pdb_info() {return _pdb_info;};
        
        void set_stub(Size ires, ATOM_TYPE atom_name, EigenXform & stub){_conformation.set_stub(ires,atom_name,stub);}
        // inline const EigenXform & sym_operator() const {
        // 	return _conformation.sym_operator();
        //}
        inline std::vector<EigenXform> get_stubs(){
            std::vector<EigenXform> stubs;
            stubs.resize(_nres);
            for(Size i =1;i<=_nres;i++){
                stubs[i-1] = stub(i);
            }
            return stubs;
        }
        inline void update_global_by_ca_stub(){
            for(Size i =0;i<_nres;i++){
                _conformation.residue(i+1).update_global_by_ca_stub();
            }
        }

        inline void set_voxel_clash_free(Size ires, bool clash_free) {
            _voxel_clash_free_residues[ires-1] = clash_free;
        }
        inline std::vector<bool> voxel_clash_free_residues() const {
            return _voxel_clash_free_residues;
        }
        inline const EigenXform stub(Size ires, Size ii_chain=1) const {
            return _conformation.stub(ires, ii_chain);
        }
        inline const EigenXform stub(Size ires, Size ii_chain, ATOM_TYPE atom_name) const {
            return _conformation.stub(ires, ii_chain, atom_name);
        }
        inline Vec xyz(Size ires, ATOM_TYPE atom, Size ii_chain=1) const {
            return _conformation.xyz(ires, atom, ii_chain);
        }
        inline void clear_energies() {
            _energy_manager.marker_changed_segment_moved(1, _nres, _conformation.root_index());
        }
        inline void update_coordinates() {
            _conformation.update_coordinates();
        }

        inline Size root_index() const {
            return _conformation.root_index();
        }

        inline void set_root_index(Size residue_index, bool move_to_current_root=false) {
            _conformation.set_root_index(residue_index, move_to_current_root);

            // reset the energy graph only if the residue (residue_index) is physically
            // moved to the current root xform
            if(move_to_current_root) {
                _energy_manager.marker_changed_segment_moved(1, _nres, _conformation.root_index());
            }
        }

        void reset_coords();
        void insert_fragment(Size ires, std::vector<Real> const & torsions, bool update_coordinates=true);
        void perturb_root(Real ang_mag=10.0, Real trans_mag=5.0, Size mode=1);
        void random_root(bool perturb_ori=true, bool perturb_x=true, bool perturb_y=true, bool perturb_z=true);

        void random_jump(Real ang_mean=-1.0, Real trans_mean=-1.0);
        void perturb_jump(Real ang_mag=10.0, Real trans_mag=5.0);
        std::vector<EigenXform> jumps() {
            return _conformation.jumps();
        }

        Conformation & conformation() { return _conformation;}
        scoring::EnergyManager & energy_manager(){return _energy_manager;};

        // pose extension
        void extend_pose(Size Nter_extension, Size Cter_extension);


        inline void snapshot(bool root_only=false) {
            _conformation.snapshot();
            if(root_only) {
                _energy_manager.snapshot(ONE_BODY);
                _energy_manager.snapshot(TWO_BODY_INTER_CHAIN);
            } else {
                 _energy_manager.snapshot();
            }
        }

        inline void rollback(bool root_only=false) {
            _conformation.rollback();
            if(root_only) {
                _energy_manager.rollback(ONE_BODY);
                _energy_manager.rollback(TWO_BODY_INTER_CHAIN);
            } else {
                 _energy_manager.rollback();
            }
        }

        inline std::vector<Real> const & chain_weights() const { return _conformation.chain_weights(); }
        inline void set_target_pose(std::string fname,char chain){_target_poseOP=std::make_shared<Pose>(fname,chain);}
        inline void set_target_pose(Pose pose){_target_poseOP=std::make_shared<Pose>(pose);}
        inline void set_target_pose(std::shared_ptr<Pose> pose){_target_poseOP=pose;}
        inline std::shared_ptr<Pose> get_target_pose(){return _target_poseOP;}
        inline void append_trajctory(){_trajactory.push_back(Pose(*this,false,false));}
        inline void pop_trajctory(){_trajactory.pop_back();}
        inline void clear_trajctory(){_trajactory.clear();}
        inline std::vector<Pose> & get_trajctory(){return _trajactory;}
        // freeze torsions
        // freeze_cart == true && freeze_torsions == false: the loaded pdb is only used to initialize the coords of the pose
        // freeze_cart == true && freeze_torsions == true:  this is used to fix any load segments
        // freeze_cart == false && freeze_torsions == true: torsions are fixed but the loaded pdb is allowed to move freely
        // freeze_cart == false && freeze_torsions == false: the torsions from the input pdb file is used to inialize the torison angles
        void load_pdb(std::string const & fname, Size insert_pos=1, bool freeze_cart=true, bool freeze_torsions=true, bool update=true,bool voxel_clash_free = false);
        void dump_pdb(std::string const & fname , bool dump_target=false,bool dump_gzip=false,bool dump_trajactory=false,bool append=false, std::string const & extra_info="");
        std::stringstream dump_pdb( bool dump_target,bool dump_gzip=false,bool dump_trajactory=false,bool append=false);



    private:
        Size _nres;
        Size _num_repeats;
        Size _num_chains;
        bool _freeze_root;
        std::string _symmetry;
        std::string _sequence;
        std::string _ss;
        Conformation _conformation;
        scoring::EnergyManager _energy_manager;
        std::set<std::pair<Size, Size> > _fixed_segments;
        std::shared_ptr<Pose> _target_poseOP;
        std::vector<Pose> _trajactory;

        std::vector<std::string> _pdb_info;
        
        //this variable tell voxelclash scorefunction which residues ignore when calculate clash with target
        std::vector<bool> _voxel_clash_free_residues;

};

typedef std::shared_ptr<Pose> PoseOP;

}

#endif

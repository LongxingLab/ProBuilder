#ifndef INCLUDED_metric_LigandNeighbors_hh
#define INCLUDED_metric_LigandNeighbors_hh

#include "basic/types.hh"
#include "scene/Pose.hh"

namespace metric {

using namespace basic;

class LigandNeighbors
{
    public:
        LigandNeighbors();
        LigandNeighbors(const std::string & ligand_pdb);
        LigandNeighbors(Real dist_midpoint, Real dist_exponent, Real angle_cutoff, Real angle_shift_factor, Real angle_exponent);
        ~LigandNeighbors();

        // load ligand
        void load_ligand(const std::string & ligand_pdb);
        bool is_ligand_loaded() {return _is_ligand_loaded;}
        void set_use_angle_component(bool flag);
        void set_ligand_xform(const EigenXform & xform);
        void clear_ligand_xform();

        // compute
        Real compute(const scene::Pose & pose);

    private:
        Real _dist_midpoint;
        Real _dist_exponent;
        Real _angle_cutoff;
        Real _angle_shift_factor;
        Real _angle_exponent;
        bool _is_ligand_loaded;
        bool _use_angle_component;
        EigenXform _ligand_xform;
        Eigen::Matrix<Real, Eigen::Dynamic, 3> _ligand_atoms;
        Eigen::Matrix<Real, Eigen::Dynamic, 3> _working_coords;

};

typedef std::shared_ptr<LigandNeighbors> LigandNeighborsOP;

}

#endif

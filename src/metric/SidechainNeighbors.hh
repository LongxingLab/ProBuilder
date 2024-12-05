#ifndef INCLUDED_metric_SidechainNeighbors_hh
#define INCLUDED_metric_SidechainNeighbors_hh

#include "basic/types.hh"
#include "scene/Pose.hh"

namespace metric {

using namespace basic;

class SidechainNeighbors
{
    public:
        SidechainNeighbors();
        SidechainNeighbors(Real dist_midpoint, Real dist_exponent, Real angle_shift_factor, Real angle_exponent);
        ~SidechainNeighbors();

        // compute
        Real compute(const scene::Pose & pose);
        Real compute(const scene::Pose & pose, bool including_intra_chain, bool including_inter_chain);
        Eigen::Array<Real, Eigen::Dynamic, Eigen::Dynamic> per_res_sidechain_neighbors(const scene::Pose & pose, bool including_intra_chain, bool including_inter_chain);
        Real motif_neighbors(const scene::Pose & pose, Size cutting_point1, Size cutting_point2);
        Eigen::Array<Real, Eigen::Dynamic, Eigen::Dynamic> sidechain_neighbors_matrix(const scene::Pose & pose, bool including_intra_chain, bool including_inter_chain);
        void set_global_flag(bool including_intra_chain, bool including_inter_chain);
        void clear_global_flag();

    private:
        Real _dist_midpoint;
        Real _dist_exponent;
        Real _angle_shift_factor;
        Real _angle_exponent;
        bool _global_flag_set;
        bool _including_intra_chain;
        bool _including_inter_chain;
};

typedef std::shared_ptr<SidechainNeighbors> SidechainNeighborsOP;

}

#endif

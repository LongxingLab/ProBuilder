#include "metric/SidechainNeighbors.hh"

namespace metric {

using namespace basic;

// constructor
SidechainNeighbors::SidechainNeighbors() :
    _dist_midpoint(9.0),
    _dist_exponent(1.0),
    _angle_shift_factor(0.5),
    _angle_exponent(2.0),
    _global_flag_set(false),
    _including_intra_chain(true),
    _including_inter_chain(true)
{}

SidechainNeighbors::SidechainNeighbors(Real dist_midpoint, Real dist_exponent, Real angle_shift_factor, Real angle_exponent) :
    _dist_midpoint(dist_midpoint),
    _dist_exponent(dist_exponent),
    _angle_shift_factor(angle_shift_factor),
    _angle_exponent(angle_exponent)
{}

// destructor
SidechainNeighbors::~SidechainNeighbors() {}

void SidechainNeighbors::set_global_flag(bool including_intra_chain, bool including_inter_chain)
{
    _global_flag_set = true;
    _including_intra_chain = including_intra_chain;
    _including_inter_chain = including_inter_chain;
}

void SidechainNeighbors::clear_global_flag()
{
    _global_flag_set = false;
}

Real SidechainNeighbors::compute(const scene::Pose & pose)
{
    if(_global_flag_set) {
        return compute(pose, _including_intra_chain, _including_inter_chain);
    } else {
        return compute(pose, true, true);
    }
}

Real SidechainNeighbors::compute(const scene::Pose & pose, bool including_intra_chain, bool including_inter_chain)
{
    auto matrix = sidechain_neighbors_matrix(pose, including_intra_chain, including_inter_chain);
    return matrix.sum() / pose.size();
}

Eigen::Array<Real, Eigen::Dynamic, Eigen::Dynamic> SidechainNeighbors::per_res_sidechain_neighbors(const scene::Pose & pose, bool including_intra_chain, bool including_inter_chain)
{
    auto matrix = sidechain_neighbors_matrix(pose, including_intra_chain, including_inter_chain);
    return matrix.rowwise().sum();
}

Real SidechainNeighbors::motif_neighbors(const scene::Pose & pose, Size cutting_point1, Size cutting_point2) 
{
    auto matrix = sidechain_neighbors_matrix(pose, true, true);

    if (cutting_point1+cutting_point2<pose.size()) {
        //
        Real s = 0.0;
        for(Size idx=0; idx<cutting_point1; ++idx){
            for(Size jdx=cutting_point2; jdx<matrix.cols(); ++jdx){
                s += matrix(idx, jdx);
            }
        }
        return s / cutting_point1;
    } else {
        //
        Real s = 0.0;
        for(Size idx=cutting_point2; idx<matrix.rows(); ++idx){
            for(Size jdx=0; jdx<cutting_point1; ++jdx){
                s += matrix(idx, jdx);
            }
        }
        return s / (matrix.rows()-cutting_point2);
    }

}
        
Eigen::Array<Real, Eigen::Dynamic, Eigen::Dynamic> SidechainNeighbors::sidechain_neighbors_matrix(const scene::Pose & pose, bool including_intra_chain, bool including_inter_chain)
{
    Size nres = pose.size();
    Size nchains = pose.num_chains();
    Size total_res = nres * nchains;

    Eigen::Matrix<Real, Eigen::Dynamic, 3> CA(total_res, 3), 
                                           CB(total_res, 3), 
                                           CA_CB(total_res, 3);
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> dist_map(nres, total_res), angle_map(nres, total_res);

    Real total_score=0;

    const std::string sequence = pose.sequence(); 
    for(Size ichain=1; ichain<=nchains; ++ichain) {
        for(Size idx=1; idx<=nres; ++idx) {
            Size global_idx = (ichain-1) * nres + idx - 1;
            CA.row(global_idx) = pose.xyz(idx, ATOM_CA, ichain);

            // get CB, depends on the residue identity
            if( sequence.at(idx-1) == 'G' ) {
                Vec N = pose.xyz(idx, ATOM_N, ichain);
                Vec C = pose.xyz(idx, ATOM_C, ichain);
                Vec N_CA = Vec(CA.row(global_idx)) - N;
                Vec CA_C = C - Vec(CA.row(global_idx));
                CB.row(global_idx) = -0.58273431*N_CA.cross(CA_C) + 0.56802827*N_CA - 0.54067466*CA_C;
                CB.row(global_idx) += CA.row(global_idx);
            } else {
                CB.row(global_idx) = pose.xyz(idx, ATOM_CB, ichain);
            }
        }
    }

    CA_CB = CB - CA;
    CA_CB.rowwise().normalize();

    // intra_chain / inter_chain
    Size start_res(0), end_res(total_res);
    if (!including_intra_chain) start_res = nres;
    if(!including_inter_chain)  end_res   = nres;

    for(Size idx=0; idx<nres; ++idx){
        for(Size jdx=0; jdx<total_res; ++jdx){
            // necessory ????
            // not efficient, not smart
            angle_map(idx, jdx) = 0.0;
            dist_map(idx, jdx)  = 0.0;
        }
    }

    for(Size idx=0; idx<nres; ++idx)
    {
        for(Size jdx=start_res; jdx<end_res; ++jdx)
        {
            if(idx == jdx) continue;

            dist_map(idx, jdx) = 1.0/(1.0+std::exp(_dist_exponent*((CB.row(jdx)-CB.row(idx)).norm() - _dist_midpoint)));

            angle_map(idx, jdx) = (CA_CB.row(idx).dot((CB.row(jdx)-CB.row(idx)).normalized()) + _angle_shift_factor) / (1+_angle_shift_factor);
            if (angle_map(idx, jdx) > 0.0) {
                angle_map(idx, jdx) = std::pow( angle_map(idx, jdx), _angle_exponent);
            } else {
                angle_map(idx, jdx) = 0.0;
            }
            /*
            // This actually works pretty well
            // To make the code consistent with the python script       
            angle_map(idx, jdx) = CA_CB.row(idx).dot((CB.row(jdx)-CB.row(idx)).normalized());
            if(angle_map(idx,jdx) > 0.0) {
                angle_map(idx,jdx) = std::pow( (angle_map(idx,jdx) + angle_shift_factor) / (1+angle_shift_factor), 2.0 );
            }
            else {
                angle_map(idx, jdx) = 0.0;
            }
            */

        }
    }

    return angle_map.array() * dist_map.array();
    // return (angle_map.array() * dist_map.array()).sum() / nres;
}


}

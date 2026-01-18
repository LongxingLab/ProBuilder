#include "utils/utils.hh"
#include "basic/types.hh"
#include "scene/Pose.hh"
#include "sampling/FragMover.hh"
#include "utils/random_util.hh"

#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include <algorithm>

namespace utils {

using namespace basic;

bool check_break(scene::Pose const & pose)
{
    Size len = pose.size();

    bool flag = false;
    for(Size ires=1; ires<len; ++ires)
    {
        Vec atom_C = pose.xyz(ires, ATOM_C);
        Vec atom_N = pose.xyz(ires+1, ATOM_N);

        Real dist = (atom_C - atom_N).squaredNorm();
        
        if(dist > 2.5) {
            flag = true;
            break;
        }
    }
    return flag;
}

Real sidechain_neighbors(scene::Pose const & pose, bool including_intra_chain, bool including_inter_chain)
{
    auto per_res_scnb = per_res_sidechain_neighbors(pose, including_intra_chain, including_inter_chain);
    return per_res_scnb.sum() / pose.size();
}

Eigen::Array<Real, Eigen::Dynamic, Eigen::Dynamic>
per_res_sidechain_neighbors(scene::Pose const & pose, bool including_intra_chain, bool including_inter_chain)
{

    const Real dist_midpoint(9.0), dist_exponent(1.0), angle_shift_factor(0.5), angle_exponent(2.0);
    
    Size nres = pose.size();
    Size nchains = pose.num_chains();
    Size total_res = nres * nchains;

    Eigen::Matrix<Real, Eigen::Dynamic, 3> N(total_res, 3), 
                                           CA(total_res, 3),
                                           C(total_res, 3), 
                                           CB(total_res, 3), 
                                           N_CA(total_res, 3), 
                                           CA_C(total_res, 3), 
                                           CA_CB(total_res, 3);
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> dist_map(nres, total_res), angle_map(nres, total_res);

    Real total_score=0;

    for(Size ichain=1; ichain<=nchains; ++ichain) {
        for(Size idx=1; idx<=nres; ++idx) {
            Size global_idx = (ichain-1) * nres + idx - 1;
            CA.row(global_idx) = pose.xyz(idx, ATOM_CA, ichain);
            C.row(global_idx)  = pose.xyz(idx, ATOM_C, ichain);
            N.row(global_idx)  = pose.xyz(idx, ATOM_N, ichain);
        }
    }

    N_CA = CA - N;
    CA_C = C  - CA;
    for(Size idx=0; idx<total_res; ++idx)
        CA_CB.row(idx) = -0.58273431*N_CA.row(idx).cross(CA_C.row(idx)) + 0.56802827*N_CA.row(idx) - 0.54067466*CA_C.row(idx);
    CB = CA_CB + CA;
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

            dist_map(idx, jdx) = 1.0/(1.0+std::exp(dist_exponent*((CB.row(jdx)-CB.row(idx)).norm() - dist_midpoint)));

            angle_map(idx, jdx) = (CA_CB.row(idx).dot((CB.row(jdx)-CB.row(idx)).normalized()) + angle_shift_factor) / (1+angle_shift_factor);
            if (angle_map(idx, jdx) > 0.0) {
                angle_map(idx, jdx) = std::pow( angle_map(idx, jdx), angle_exponent);
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

    return (angle_map.array() * dist_map.array()).rowwise().sum();
    // return (angle_map.array() * dist_map.array()).sum() / nres;
}

Real ligand_neighbors(scene::Pose const & pose, Eigen::Matrix<Real, Eigen::Dynamic, 3> const & ligand)
{
    const Real dist_midpoint(5.0), dist_exponent(1.0), angle_shift_factor(0.5), angle_exponent(1.0);
    
    Size nres = pose.size();
    Size nchains = pose.num_chains();
    Size total_res = nres * nchains;

    Size ligand_natoms = ligand.rows();

    std::string const sequence = pose.sequence();

    Eigen::Matrix<Real, Eigen::Dynamic, 3> N(total_res, 3), 
                                           CA(total_res, 3),
                                           C(total_res, 3), 
                                           CB(total_res, 3),
                                           /*
                                           N_CA(total_res, 3), 
                                           CA_C(total_res, 3), */
                                           CA_CB(total_res, 3);
                                        
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> dist_map(total_res, ligand_natoms), angle_map(total_res, ligand_natoms);

    Real total_score=0;

    for(Size ichain=1; ichain<=nchains; ++ichain) {
        for(Size idx=1; idx<=nres; ++idx) {
            Size global_idx = (ichain-1) * nres + idx - 1;
            CA.row(global_idx) = pose.xyz(idx, ATOM_CA, ichain);
            C.row(global_idx)  = pose.xyz(idx, ATOM_C, ichain);
            N.row(global_idx)  = pose.xyz(idx, ATOM_N, ichain);

            if( sequence.at(idx-1) == 'G' ) {
                Vec N_CA = CA.row(global_idx) - N.row(global_idx);
                Vec CA_C = C.row(global_idx) - CA.row(global_idx);
                CB.row(global_idx) = -0.58273431*N_CA.cross(CA_C) + 0.56802827*N_CA - 0.54067466*CA_C;
                CB.row(global_idx) += CA.row(global_idx);
            } else {
                CB.row(global_idx) = pose.xyz(idx, ATOM_CB, ichain);
            }
        }
    }

    /*
    N_CA = CA - N;
    CA_C = C  - CA;
    for(Size idx=0; idx<total_res; ++idx)
        CA_CB.row(idx) = -0.58273431*N_CA.row(idx).cross(CA_C.row(idx)) + 0.56802827*N_CA.row(idx) - 0.54067466*CA_C.row(idx);
    CB = CA_CB + CA;
    */
    CA_CB = CB - CA;
    CA_CB.rowwise().normalize();


    // initialize to 0
    for(Size idx=0; idx<total_res; ++idx){
        for(Size jdx=0; jdx<ligand_natoms; ++jdx){
            // necessory ????
            // not efficient, not smart
            angle_map(idx, jdx) = 0.0;
            dist_map(idx, jdx)  = 0.0;
        }
    }

    // fill in the numbers
    for(Size idx=0; idx<total_res; ++idx)
    {
        for(Size jdx=0; jdx<ligand_natoms; ++jdx)
        {

            dist_map(idx, jdx) = 1.0/(1.0+std::exp(dist_exponent*((ligand.row(jdx)-CB.row(idx)).norm() - dist_midpoint)));

            angle_map(idx, jdx) = (CA_CB.row(idx).dot((ligand.row(jdx)-CB.row(idx)).normalized()) + angle_shift_factor) / (1+angle_shift_factor);
            if (angle_map(idx, jdx) > 0.0) {
                angle_map(idx, jdx) = std::pow( angle_map(idx, jdx), angle_exponent);
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


    return (angle_map.array() * dist_map.array()).sum() / ligand_natoms;
}


void print_xform(EigenXform const & x)
{
    std::cout << x(0,0) << " " << x(0,1) << " " << x(0,2) << " " << x(0,3) << std::endl
              << x(1,0) << " " << x(1,1) << " " << x(1,2) << " " << x(1,3) << std::endl
              << x(2,0) << " " << x(2,1) << " " << x(2,2) << " " << x(2,3) << std::endl;
}

void print_vec(Vec const & v)
{
    std::cout << v(0) << " " << v(1) << " " << v(2) << std::endl;
}


void get_repeat_parameters_from_stubs(scene::Pose & pose,Real & rise_out, Real & radius_out, Real & omega_out,Vec & axis_out,Vec & axis_center,bool debug )
{
    using namespace Eigen;
    assert(pose.num_repeats()>1);
    assert(pose.size()%pose.num_repeats() ==0);
	const Real FLOAT_PRECISION = 1e-6;

    Size residue_offset = pose.size()/pose.num_repeats();
	EigenXform sym_operator =pose.stub(1+residue_offset)*pose.stub(1).inverse() ;
	Eigen::AngleAxis<Real> angle_axis(sym_operator.rotation());

	// axis ang
	Vec axis = angle_axis.axis();
	Real ang = angle_axis.angle();

	Real translation_along_axis_of_rotation = axis.dot(sym_operator.translation());
	Vec cen(0.0,0.0,0.0);
	// these points lie on a circle and the plane is perpendicular with the axis 
    Vec c_A = pose.conformation().center_vec(1,residue_offset);
    Vec c_B = pose.conformation().center_vec(residue_offset+1,residue_offset*2);
    Vec p0 = c_A;
	Vec p1 = sym_operator * p0;
	Vec p2 = sym_operator * p1;
    if(debug){
        char buf[128];
        Size anum(1), rnum(1);
        std::ofstream out("/home/chentong/wefold/build/test_point.pdb");
        Vec xyz = p0;
        snprintf(buf,128,"%s%5i %4s %3s %c%4i    %8.3f%8.3f%8.3f%6.2f%6.2f %11s\n",
                    "HETATM",
                    anum++,
                    "BURR",
                    "BUR",
                    'B',
                    rnum++,
                    xyz[0],xyz[1],xyz[2],
                    1.0,
                    1.0,
                    "B"
                );

        out << buf;
        xyz = p1;
        snprintf(buf,128,"%s%5i %4s %3s %c%4i    %8.3f%8.3f%8.3f%6.2f%6.2f %11s\n",
                    "HETATM",
                    anum++,
                    "BURR",
                    "BUR",
                    'B',
                    rnum++,
                    xyz[0],xyz[1],xyz[2],
                    1.0,
                    1.0,
                    "B"
                );

        out << buf;
        xyz = p2;
        snprintf(buf,128,"%s%5i %4s %3s %c%4i    %8.3f%8.3f%8.3f%6.2f%6.2f %11s\n",
                    "HETATM",
                    anum++,
                    "BURR",
                    "BUR",
                    'B',
                    rnum++,
                    xyz[0],xyz[1],xyz[2],
                    1.0,
                    1.0,
                    "B"
                );

        out << buf;
    }
	p1 -= axis * (p1-p0).dot(axis);
	p2 -= axis * (p2-p0).dot(axis);
	Real d = (p1-p0).norm();

    Real l = d / (2.0*std::tan(ang/2.0));
    Vec tocen = (p1-p0).normalized().cross(axis)*l;
    if(tocen.dot(p2-p1)<0.0){
        tocen = -tocen;
    }
    cen = (p0+p1)/2.0 + tocen;
    Real radius = (p0-cen).norm();
    radius_out =radius;
    omega_out = ang;
    rise_out = translation_along_axis_of_rotation;
    axis_out = axis;
    
    axis_center = cen + translation_along_axis_of_rotation * axis*(float(pose.num_repeats())/2.0 -0.5);
    return;
}
void get_repeat_parameters_from_coords(scene::Pose & pose, Real & rise_out, Real & radius_out, Real & omega_out,Vec & axis_out,Vec & axis_center)  {
	using Eigen::MatrixXd;
	using namespace Eigen;
	using namespace std;
	Size seg_size = pose.size()/pose.num_repeats();
    Size startAtRepeat = pose.num_repeats()/2;
	Size startResOffset = seg_size*(startAtRepeat-1);
	MatrixXf A(3, seg_size);
	MatrixXf B(3, seg_size);
	for ( Size i = 1; i <= seg_size; i++ ) {
		Size offset =i+startResOffset;
		Vec coord = pose.xyz(offset,ATOM_CA);
		A(0,i-1) = coord[0];
		A(1,i-1) = coord[1];
		A(2,i-1) = coord[2];

	}
	for ( Size i = seg_size+1; i <= seg_size*2 ; i++ ) {
		Size offset =i+startResOffset;
		Vec coord = pose.xyz(offset,ATOM_CA);
		B(0,i-1-seg_size) = coord[0];
		B(1,i-1-seg_size) = coord[1];
		B(2,i-1-seg_size) = coord[2];
	}
	Vector3f c_A(A.row(0).mean(), A.row(1).mean(), A.row(2).mean());
	Vector3f c_B(B.row(0).mean(), B.row(1).mean(), B.row(2).mean());
	MatrixXf x_A(A.rows(),A.cols());
	MatrixXf x_B(B.rows(),B.cols());
	for ( int i=0; i<A.cols(); i++ ) {
		x_A.col(i)=A.col(i)-c_A;
		x_B.col(i)=B.col(i)-c_B;
	}
	Matrix3f cov= (x_B * x_A.transpose()) / x_A.cols();
	JacobiSVD<MatrixXf> svd(cov, ComputeFullU | ComputeFullV);

	Matrix3f Rt=svd.matrixU() * svd.matrixV().transpose();
	Matrix3f R;
	R<< 1,0,0, 0,1,0, 0,0,Rt.determinant();
	Matrix3f H= svd.matrixU() * R * svd.matrixV().transpose();

	Real acos_fix_tmp = (H.trace()-1)/2;
	if ( acos_fix_tmp <= -1.0 ) {
		acos_fix_tmp = -1.0;
	}
	if ( acos_fix_tmp >= 1.0 ) {
		acos_fix_tmp = 1.0;
	}
	Real omega = acos(acos_fix_tmp);
	Matrix3f I = Matrix3f::Identity();
	Matrix3f N = 0.5*(H+H.transpose()) - cos(omega)*I;


	Vector3f hN(0,0,0);//initializiation for the testing server
	Real scalar = 0;
	Real max_scalar = -10000000;
	for ( Size i = 0; i<=2; i++ ) {
		scalar = N.col(i).norm();
		if ( scalar > max_scalar ) {
			max_scalar = scalar;
			hN = N.col(i)/N.col(i).norm();
		}
	}

	Real sin_omega = (H(1,0)-H(0,1)) / (2*hN(2));
	if ( sin_omega < 0 ) hN = -1 * hN;

	Vector3f t = c_B - H*c_A;
	Real L = t.dot(hN) ;
	Real rise=abs(L);

	Matrix3f Ncross;
	Ncross << 0,-1*hN(2),hN(1), hN(2),0,-1*hN(0), -1*hN(1),hN(0),0 ;
	Matrix3f R0t= (1-cos(omega))*I - sin(omega)*Ncross;
	Vector3f R0 = R0t.inverse() * (t-L*hN);
	Vector3f pA= (c_A-R0)-(hN*(hN.dot(c_A-R0)));
	Vector3f pB= (c_B-R0)-(hN*(hN.dot(c_B-R0)));


	Real direction = L * hN.dot(pA.cross(pB));
	Real radius=pA.norm();
	rise_out = rise;
	radius_out = radius;
	omega_out = omega;
    axis_out = hN;
    axis_center = R0+(hN*(hN.dot(c_A-R0))) + hN*(float(pose.num_repeats())/2+0.5-(pose.num_repeats()/2));
}

Real rmsd_no_super(const scene::Pose & pose1, const scene::Pose & pose2, Size start_res, Size end_res, bool CA_only)
{

    if( end_res == -1 ) end_res = pose1.size();
    assert( pose1.size() == pose2.size() && start_res < end_res && end_res <= pose1.size() );

    std::vector<ATOM_TYPE> atom_types = {ATOM_CA, ATOM_N, ATOM_C};

    Real dev(0.0);
    Size counts(0);
    for(Size ires=start_res; ires<=end_res; ++ires) {
        for( Size iatom=0; iatom<3; ++iatom ) {

            Vec  iatom_xyz = pose1.xyz(ires, atom_types[iatom]);
            Vec  jatom_xyz = pose2.xyz(ires, atom_types[iatom]);
            dev += (iatom_xyz - jatom_xyz).squaredNorm();
            ++counts;

            if(CA_only) break;
        }
    }

    return std::sqrt(dev / counts);
}

}

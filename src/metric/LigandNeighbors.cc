#include "metric/LigandNeighbors.hh"
#include "basic/assert.hh"

#include <fstream>

namespace metric {

using namespace basic;

// constructor
LigandNeighbors::LigandNeighbors() :
    _dist_midpoint(7.0),
    _dist_exponent(0.7),
    _angle_cutoff(-0.5),
    _angle_shift_factor(2.0),
    _angle_exponent(2.0),
    _is_ligand_loaded(false),
    _use_angle_component(false)
{}

LigandNeighbors::LigandNeighbors(const std::string & ligand_pdb) :
    _dist_midpoint(7.0),
    _dist_exponent(0.7),
    _angle_cutoff(-0.5),
    _angle_shift_factor(2.0),
    _angle_exponent(2.0),
    _use_angle_component(false)
{
    load_ligand(ligand_pdb);
}

LigandNeighbors::LigandNeighbors(Real dist_midpoint, Real dist_exponent, Real angle_cutoff, Real angle_shift_factor, Real angle_exponent) :
    _dist_midpoint(dist_midpoint),
    _dist_exponent(dist_exponent),
    _angle_cutoff(angle_cutoff),
    _angle_shift_factor(angle_shift_factor),
    _angle_exponent(angle_exponent),
    _is_ligand_loaded(false),
    _use_angle_component(false)
{}

// destructor
LigandNeighbors::~LigandNeighbors() {}

void LigandNeighbors::load_ligand(const std::string & ligand_pdb)
{
    // parse the pdb and extract the coords of the atoms
    // all atoms excluding H
    ALWAYS_ASSERT_MSG(ligand_pdb.substr(ligand_pdb.size()-4,4) == ".pdb", "Only support pdb format!!!");
    
    std::ifstream f(ligand_pdb);
    std::vector<Vec> atoms_xyz;
    {
        for(std::string line; std::getline(f, line);)
        {
            if( (line.substr(0,4) == "ATOM" || line.substr(0,6) == "HETATM") && (line.substr(13,1) != "H") && (line.substr(77,1) != "H") ) {
                atoms_xyz.push_back(Vec(std::stof(line.substr(30,8)),std::stof(line.substr(38,8)),std::stof(line.substr(46,8))));
            }
        }
        _ligand_atoms.resize(atoms_xyz.size(), 3);
        for(Size iatom=0; iatom<atoms_xyz.size(); ++iatom) {
            _ligand_atoms.row(iatom) = atoms_xyz[iatom];
        }
    }
    f.close();

    _working_coords = _ligand_atoms;

    std::cout << "The ligand coordinates are loaded successfully!!!!" << std::endl;

    _is_ligand_loaded = true;
}

void LigandNeighbors::set_ligand_xform(const EigenXform & xform)
{
    _ligand_xform = xform;
    _working_coords = ( ( _ligand_xform.rotation() * ( _ligand_atoms.transpose() ) ).colwise() + _ligand_xform.translation() ).transpose();
}
void LigandNeighbors::clear_ligand_xform()
{
    _ligand_xform = EigenXform::Identity();
    _working_coords = _ligand_atoms;
}

void LigandNeighbors::set_use_angle_component(bool flag)
{
    _use_angle_component = flag;
}

Real LigandNeighbors::compute(const scene::Pose & pose)
{
    Size nres = pose.size();
    Size nchains = pose.num_chains();
    Size total_res = nres * nchains;

    Size ligand_natoms = _ligand_atoms.rows();

    Eigen::Matrix<Real, Eigen::Dynamic, 3> CA(total_res, 3),
                                           CB(total_res, 3),
                                           CA_CB(total_res, 3);
                                        
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> dist_map(total_res, ligand_natoms), angle_map(total_res, ligand_natoms);

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

            dist_map(idx, jdx) = 1.0/(1.0+std::exp(_dist_exponent*((_working_coords.row(jdx)-CB.row(idx)).norm() - _dist_midpoint)));

            if( _use_angle_component ) {
                Real cos_ang = CA_CB.row(idx).dot((_working_coords.row(jdx)-CA.row(idx)).normalized());
                if(cos_ang > _angle_cutoff) {
                    cos_ang = 1.0;
                } 
                angle_map(idx, jdx) = (cos_ang + _angle_shift_factor) / (1+_angle_shift_factor);
                if (angle_map(idx, jdx) > 0.0) {
                    angle_map(idx, jdx) = std::pow( angle_map(idx, jdx), _angle_exponent);
                } else {
                    angle_map(idx, jdx) = 0.0;
                }
            }

        }
    }

    if( _use_angle_component ) {
        return (angle_map.array() * dist_map.array()).sum() / ligand_natoms;
    } else {
        return dist_map.array().sum() / ligand_natoms;
    }
}


}

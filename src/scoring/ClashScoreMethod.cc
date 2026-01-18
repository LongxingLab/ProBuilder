#include "scoring/ClashScoreMethod.hh"
#include "scoring/Energy.hh"
#include "utils/hash_util.hh"
#include "utils/math_util.hh"

#include <cmath>
#include <string>
#include <fstream>
#include <Eigen/Geometry>
#include <parallel_hashmap/phmap_dump.h>
#include <map>
namespace scoring {

using namespace basic;
using namespace scene;

Pose
VoxelClashScoreMethod::get_grid_pose( std::string fname,std::string mode){
    assert( mode=="PROTEIN" || mode=="MOLECULE" );
    _pdb_mode = mode;
    if(mode=="PROTEIN")return Pose(fname,' ',"C1",1);
    else if(mode=="MOLECULE"){
        std::ifstream pdb_file(fname);
        std::vector<std::vector<AtomInfo>> all_atom_info;
        Size segment_len = 0;
        std::map<std::string,ATOM_TYPE> aa1toatom= {{"N",ATOM_N},{"C",ATOM_C},{"O",ATOM_O},{"S",ATOM_S},{"H",ATOM_H}};
        if(!pdb_file.is_open()){
            std::cout<<"Error: can not open file "<<fname<<std::endl;
            exit(1);
        }

        Size cur_residue_number=0;
        std::string cur_residue_origin_number;
        for(std::string line;std::getline(pdb_file,line);)
        {   
            if( line.substr(0,4) == "ATOM" ||  line.substr(0,6) == "HETATM" ){
                if(cur_residue_origin_number!=line.substr(22,4)){
                    cur_residue_origin_number=line.substr(22,4);
                    cur_residue_number++;
                    all_atom_info.push_back(std::vector<AtomInfo>());
                }
                AtomInfo atom_info;
                atom_info.atom_name = utils::trim(line.substr(12,4));
                atom_info.res_name = line.substr(17,3);
                atom_info.local_postion = Vec(std::stof(line.substr(30,8)),std::stof(line.substr(38,8)),std::stof(line.substr(46,8)));
                atom_info.element = utils::trim(line.substr(77,2));
                atom_info.atom_type = aa1toatom[line.substr(77,1)];
                all_atom_info.back().push_back(atom_info);
            }
        }
        pdb_file.close();
        Pose pose = Pose(all_atom_info.size());
        pose.conformation() = Conformation(all_atom_info.size(),"C1");
        for(Size ii=1; ii<=all_atom_info.size(); ++ii){
            EigenXform ca_stub = utils::xform_from_3points(Vec(1,0,0),Vec(0,1,0),Vec(0,0,1));
            pose.conformation().set_stub(ii, ATOM_CA, ca_stub );
            for(auto & atom_info : all_atom_info[ii-1]){
                atom_info.local_postion = ca_stub.inverse()*atom_info.local_postion;
                pose.conformation().push_atom_info(ii,atom_info);
            }
        }
        return pose;
    }
    else{
        std::cout<<"Error: voxel clash load mode is not supported"<<std::endl;
        exit(1);
    }
}
void
VoxelClashScoreMethod::initialize_voxel_grid(std::string mode, Real ball_radius, Real grid_resolution)
{
    assert( mode=="ALL" || mode=="BB_CB" || mode=="BB_ONLY" );
    
    Real ball_radius_(ball_radius);
    Real grid_resolution_(grid_resolution);

    // load pdb file
    std::vector<Vec> atoms_xyz;
    std::vector<ATOM_TYPE> atoms_type;
    std::vector<Size> atoms_ires,atoms_iatom;

    for(Size ires=1;ires<=grid_pose.size();ires++)
    {   
        for(Size iatom=1;iatom<=grid_pose.conformation().residue(ires).natoms();iatom++){
            std::string atom_name = grid_pose.conformation().residue(ires).atom_name(iatom);

            if( mode == "ALL" ) {
                // nothing
            }

            if( mode == "BB_CB" ) {
                if( atom_name!="O" && atom_name!="C" && atom_name!="N" && atom_name!="CA" && atom_name!="CB"){
                    continue;
                }
            }

            if( mode == "BB_ONLY" ) {
                if( atom_name!="O" && atom_name!="C" && atom_name!="N" && atom_name!="CA"){
                    continue;
                }
            }

            atoms_xyz.push_back(grid_pose.conformation().residue(ires).xyz(iatom));
            atoms_type.push_back(grid_pose.conformation().residue(ires).atom_type(iatom));
            atoms_ires.push_back(ires);
            atoms_iatom.push_back(iatom);
        }
    }

    // the upper bound and lower bound
    Vec ub{ -1e9, -1e9, -1e9};
    Vec lb{  1e9,  1e9,  1e9};
    //
    for( Vec & xyz : atoms_xyz ) {
        for ( int i = 0; i < 3; i++ ) {
            ub[ i ] = std::max<Real>( ub[ i ], xyz[ i ] );
            lb[ i ] = std::min<Real>( lb[ i ], xyz[ i ] );
        }        
    }

    // Pad the edges of the box
    for ( Size i = 0; i < 3; i++ ) {
        ub[ i ] += ball_radius_ * 2 + grid_resolution_ * 2; // should theoretically be ok without doing
        lb[ i ] -= ball_radius_ * 2 + grid_resolution_ * 2; // double resolution. But I don't trust it.
    }

    // Make the voxel array at the right resolution
    Vec cs{ grid_resolution_, grid_resolution_, grid_resolution_ };

    _indicate_array = std::make_shared< VoxelIndicateArray >( lb, ub, cs );
    _voxel_array = std::make_shared< VoxelArray >( lb, ub, cs );
    _indicate_voxel_array = std::make_shared< VoxelArray >( lb, ub, cs );
    _VoxelArray<Real,Real> indicate_dis_array = _VoxelArray<Real,Real>( lb, ub, cs );
    typedef typename VoxelArray::Bounds Bounds;

    Real double_radius_sq_ = ball_radius_ * ball_radius_ * 4;
    // fill in the grid
    for( Size i =0;i< atoms_xyz.size();i++ ) {
        // ball_radius = vander_waal_radius[atoms_type[i]];
        Size double_radius_steps = (Size) std::ceil( (ball_radius_*2)  / grid_resolution_ ) + 1;

        Real double_radius_sq = ball_radius_*ball_radius_ * 4;
        Real double_radius_sq_max = 1.8*1.8 * 4;
        Bounds centered_xyz = _voxel_array->indices_to_center( _voxel_array->floats_to_index( atoms_xyz[i] ) );

        Vec worker { 0, 0, 0};
        for ( Size ix = -double_radius_steps; ix <= double_radius_steps; ix++ ) {
            worker[0] = ix * grid_resolution_ + centered_xyz[0];

            for ( Size iy = -double_radius_steps; iy <= double_radius_steps; iy++ ) {
                worker[1] = iy * grid_resolution_ + centered_xyz[1];

                for ( Size iz = -double_radius_steps; iz <= double_radius_steps; iz++ ) {
                    worker[2] = iz * grid_resolution_ + centered_xyz[2];
                    // if ( xyz.distance_squared( worker ) > double_radius_sq ) continue;
                    if((atoms_xyz[i]-worker).squaredNorm()<double_radius_sq){
                        (*_voxel_array)[ worker ] = true;
                        // Real dis = (atoms_xyz[(*_indicate_array)[ worker ].first]-worker).norm()-vander_waal_radius[atoms_type[(*_indicate_array)[ worker ].first]];
                        // // temporary indicate for postion at upper array to fast caculate
                        // if((*_indicate_voxel_array)[worker]){
                        //     if(indicate_dis_array[worker]>dis){
                        //         indicate_dis_array[worker]=dis;
                        //         (*_indicate_array)[ worker ].first = i;
                        //     }
                        // }else{
                        //     (*_indicate_array)[ worker ]=std::pair<Size,Size>(i,0);
                        //     (*_indicate_voxel_array)[worker] = true;
                        //     indicate_dis_array[worker] = dis;
                        // }
                    }
                    if((atoms_xyz[i]-worker).squaredNorm()<double_radius_sq_max){
                        Real dis = pow((vander_waal_radius[atoms_type[(*_indicate_array)[ worker ].first]]+1.8),2)- (atoms_xyz[(*_indicate_array)[ worker ].first]-worker).squaredNorm();
                        // temporary indicate for postion at upper array to fast caculate
                        if((*_indicate_voxel_array)[worker]){
                            if(indicate_dis_array[worker]>dis){
                                indicate_dis_array[worker]=dis;
                                (*_indicate_array)[ worker ].first = i;
                            }
                        }else{
                            (*_indicate_array)[ worker ]=std::pair<Size,Size>(i,0);
                            (*_indicate_voxel_array)[worker] = true;
                            indicate_dis_array[worker] = dis;
                        }
                    }
                } // z
            } // y
        } // x
    }
    typedef std::array< Size, 3> Indices;
    for(Size xx=0; xx<(*_indicate_array).shape()[0]; ++xx) {
        for(Size yy=0; yy<(*_indicate_array).shape()[1]; ++yy) {
            for(Size zz=0; zz<(*_indicate_array).shape()[2]; ++zz) {
                Indices idx = Indices{ { xx, yy, zz } };
                if(0 != _indicate_array->operator()(idx).first) {
                _indicate_array->operator()(idx)=std::pair<Size,Size>(atoms_ires[_indicate_array->operator()(idx).first],atoms_iatom[_indicate_array->operator()(idx).first]);
                }
            }
        }
    }
    
}

void
VoxelClashScoreMethod::visualize_voxel_grid( std::string fname ) {
    assert(_voxel_array != nullptr);

    std::ofstream of(fname);
    _voxel_array->visualize(of);
    of.close();
}

// update here tommorrow
Real VoxelClashScoreMethod::score(scene::Pose &pose,Size cutoff) const {

    Size len = pose.size();
    std::string sequence = pose.sequence();
    std::vector<bool> voxel_clash_free_residues = pose.voxel_clash_free_residues();
    EnergyManager & energy_manager = pose.energy_manager();
    std::vector<EnergyOP> energyOPs = energy_manager.get_energy_tables(_score_type);
    assert(_no_cache||(energyOPs.size()==1 and energyOPs[0]->energy_type()==ONE_BODY));
    OneBodyEnergyOP voxel_clash_energy = _no_cache?nullptr:std::static_pointer_cast<OneBodyEnergy>(energyOPs[0]);
    Real total_clash=0;
    Size ichain=_no_cache?1:voxel_clash_energy->ichain();
    if(!_no_cache&&voxel_clash_energy == nullptr){
        std::cout<<"energy op from base to onebody fail!! check energy initialize is right"<<std::endl;
        exit(0);
    }
    for(Size ires=1; ires<=len; ++ires) {
        Real clash=0;

        // if the amino acid is 'X', then skip
        if(sequence[ires-1] == 'X') {
            voxel_clash_energy->set_score_per_residue(ires, 0.0);
            continue;
        }
        if(_designable_res.size() != 0 && _designable_res.count(ires) == 0 ) {
            voxel_clash_energy->set_score_per_residue(ires, 0.0);
            continue;
        }
        if(voxel_clash_free_residues[ires-1]){
            voxel_clash_energy->set_score_per_residue(ires, 0.0);
            continue;
        }

        if (!_no_cache&&!voxel_clash_energy->is_changed_residue(ires))continue;
        ATOM_TYPE atom_list[] = {ATOM_N,ATOM_CA,ATOM_C,ATOM_O,ATOM_CB};
        for(auto atom:atom_list){
            Vec xyz = pose.xyz(ires, atom,ichain);
            // relative xform set
            // ligand shifts
            if(_is_relative_xform_set) xyz = _relative_xform * xyz;
            if ( _voxel_array->at( xyz ) )clash+=1;
        }
        total_clash+=clash;

        if(!_no_cache)voxel_clash_energy->set_score_per_residue(ires,clash*_clash_multiplier);
    }
    if(_no_cache)return total_clash;
    else return voxel_clash_energy->weighted_score();
}

Real VoxelClashScoreMethod::hresl_score(scene::Pose &pose)  {
    std::unordered_map<ATOM_TYPE,Real> vander_waal_radius = {{ATOM_N,1.54},{ATOM_C,1.70},{ATOM_O,1.52},{ATOM_S,1.8},{ATOM_CA,1.70},{ATOM_CB,1.70}};
    Size len = pose.size();
    EnergyManager energy_manager = pose.energy_manager();
    std::vector<EnergyOP> energyOPs = energy_manager.get_energy_tables(_score_type);
    assert(_no_cache||(energyOPs.size()==1 and energyOPs[0]->energy_type()==ONE_BODY));
    OneBodyEnergyOP voxel_clash_energy = _no_cache?nullptr:std::static_pointer_cast<OneBodyEnergy>(energyOPs[0]);
    Real total_clash=0;
    Size ichain=_no_cache?1:voxel_clash_energy->ichain();
    if(!_no_cache&&voxel_clash_energy == nullptr){
        std::cout<<"energy op from base to onebody fail!! check energy initialize is right"<<std::endl;
        exit(0);
    }
    for(Size ires=1; ires<=len; ++ires) {
        Real clash=0;
        if (!_no_cache&&!voxel_clash_energy->is_changed_residue(ires))continue;

        // use vitual O and CB here , notice xyz(Size iatom) get real atoms
        ATOM_TYPE atom_list[] = {ATOM_CA,ATOM_C,ATOM_N,ATOM_O,ATOM_CB};
        for(auto atom:atom_list){
            Vec xyz = pose.xyz(ires, atom,ichain);
            if(_indicate_voxel_array->at( xyz ) && (_indicate_array->at( xyz ).first>0) ){
                std::pair<Size,Size> pair = _indicate_array->at( xyz );
                Vec grid_atom_vec =  grid_pose.conformation().residue(pair.first).xyz(pair.second);
                ATOM_TYPE grid_atom_type = grid_pose.conformation().residue(pair.first).atom_type(pair.second);
                Real according_dis= pow((vander_waal_radius[grid_atom_type]+vander_waal_radius[atom]),2);
                Real dist = (xyz-grid_atom_vec).squaredNorm();
                clash+= according_dis<dist?0:according_dis-dist;
            }
        }
        total_clash+=clash;
        if(!_no_cache)voxel_clash_energy->set_score_per_residue(ires,clash);
    }
    if(_no_cache)return total_clash;
    else return voxel_clash_energy->weighted_score();
}

Real VoxelClashScoreMethod::score(std::vector<EigenXform> & stubs,Size cutoff)const{
    Real clash_total = 0;
    for(auto stub:stubs) {
        if ( _voxel_array->at( stub.translation() ) )clash_total+=1;
        if(clash_total>=cutoff)return clash_total;
    }
    return clash_total;
}

void VoxelClashScoreMethod::set_relative_pos(EigenXform const & xform)
{
    //
    _is_relative_xform_set = true;
    _relative_xform = xform.inverse();
}



// clash score


ClashScoreMethod::ClashScoreMethod(Size sequence_separation, Real clash_multiplier,bool scaffold_intra_clash,bool scaffold_target_inter_clash,Real interface_extra_punish) :
    BaseScoreMethod(),
    _min_sequence_separation(sequence_separation),
    _clash_multiplier(clash_multiplier),
    _squared_dist_cutoff(36),
    _scaffold_intra_clash(scaffold_intra_clash),
    _scaffold_target_inter_clash(scaffold_target_inter_clash),
    _interface_extra_punish(interface_extra_punish),
    _CB_swelling_factor(1.0),
    _R1R2_pow{6.96,5.80,8.88,4.97,6.96,5.80,14.28,5.76,7.67,14.13,8.88,5.76,8.46,7.56,7.56,4.97,7.67,7.56,9.30,5.85,6.96,14.13,7.56,5.85,13.10}
    {
        set_score_type(CLASH);
    }

void ClashScoreMethod::set_CB_swelling_factor(Real v) {

    // the CB-CB got swelled twice
    for(Size i=0; i<5; i++) _R1R2_pow[i][4] = _R1R2_pow[i][4] / _CB_swelling_factor * v;
    for(Size j=0; j<5; j++) _R1R2_pow[4][j] = _R1R2_pow[4][j] / _CB_swelling_factor * v;
    _CB_swelling_factor = v;
}

Real ClashScoreMethod::base_score(scene::Pose & pose,Size start_chain,Size end_chain) const{

    std::vector<ATOM_TYPE> atom_types = {ATOM_N, ATOM_CA, ATOM_C, ATOM_O, ATOM_CB};

    Size len = pose.size();
    std::string const sequence = pose.sequence();
    EnergyManager & energy_manager = pose.energy_manager();
    Conformation & conformation = pose.conformation();
    std::vector<EnergyOP> energyOPs = energy_manager.get_energy_tables(_score_type);
    assert(energyOPs.size()>0 && (energyOPs[0]->energy_type()==TWO_BODY_INTER_CHAIN || energyOPs[0]->energy_type()==TWO_BODY_INTRA_CHAIN) );
    Real clash_total = 0;
    // handle situation of intra and inter clash
    for(auto energy_base:energyOPs){
        if(!_scaffold_intra_clash)break;
        TwoBodyEnergyOP energy_clash= std::static_pointer_cast<TwoBodyEnergy>(energy_base);
        if(energy_clash == nullptr){
            std::cout<<"energy op from base to twobody fail!! check energy initialize is right"<<std::endl;
            exit(0);
        }
        if(energy_clash->jchain()<start_chain) continue;
        if(end_chain!=-1&&energy_clash->jchain()>end_chain)continue;
        if(_energy_type != ENERGY_TYPE_UNDEFINED && _energy_type != energy_clash->energy_type())continue;
        for(Size ires=1; ires<=len; ++ires) {
            for(Size jres=(energy_clash->energy_type()==TWO_BODY_INTRA_CHAIN?ires+_min_sequence_separation:1); jres<=len; ++jres) {

                if( !energy_clash->is_changed_pair(ires, jres) ) {
                    continue;
                }
                Vec  res1_cb_xyz = pose.xyz(ires, ATOM_CB,energy_clash->ichain());
                Vec  res2_cb_xyz = pose.xyz(jres, ATOM_CB,energy_clash->jchain());
                Vec t = res1_cb_xyz-res2_cb_xyz;
                if(t.squaredNorm() > _squared_dist_cutoff) {
                    energy_clash->set_pair_score(ires, jres, 0.0);
                } else {         
                    

                    // values based on Rosetta vdw radii
                    //std::vector<Real> atom_vdw = {1.80,2.0,1.91,1.54,2.2}; // CB 2.4

                    // values from Chentong based on the de novo scaffold set
                    // 75% coverage
                    // N--N  N--CA   N--C    N--O    N--CB   CA--CA  CA--C   CA--O   CA--CB  C--C    C--O    C--CB   O--O    O--CB   CB--CB
                    // 2.64  2.41    2.98    2.23    2.64    3.78    2.40    2.77    3.76    2.91    2.75    2.75    3.05    2.42    3.62

                    // clash check is important for model quality, I recalculated the cutoff values based on Chentong's script
                    // chash check is also vital for the efficiency of the mcmc protocol
                    // soft clash is better??
                    // /storage/caolongxingLab/software/database/pdb/PISCES_181019/biounit/relaxed30.chains/*pdb

                    // All, also CB for Gly (this doesn't make sense, and also incompatible with wefold) !!!!!!!!!!!!!!!!!!! Bad !!!!!!!
                    // 75% coverage
                    // N--N  N--CA   N--C    N--O    N--CB   CA--CA  CA--C   CA--O   CA--CB  C--C    C--O    C--CB   O--O    O--CB   CB--CB
                    // 2.57  2.36    2.99    2.20    2.60    3.04    2.34    2.61    2.81    2.84    2.67    2.62    2.91    2.06    2.64

                    // All, Gly (CB == CA) !!!!!!!!!!!!!!!!!!!!!!! Bad!!!!!!!!
                    // 75% coverage
                    // N--N  N--CA   N--C    N--O    N--CB   CA--CA  CA--C   CA--O   CA--CB  C--C    C--O    C--CB   O--O    O--CB   CB--CB
                    // 2.57 2.35    2.99    2.20    2.38    3.04    2.34    2.61    3.52    2.84    2.67    2.36    2.92    2.66    3.48

                    // Others (no Gly)
                    // 75% coverage
                    // N--N  N--CA   N--C    N--O    N--CB   CA--CA  CA--C   CA--O   CA--CB  C--C    C--O    C--CB   O--O    O--CB   CB--CB
                    // 2.60  2.36    2.99    2.21    2.76    3.06    2.35    2.61    3.59    2.84    2.67    2.96    2.94    2.89    3.50

                    // GLY
                    // 75% coverage
                    // N--N  N--CA   N--C    N--O    N--CB   CA--CA  CA--C   CA--O   CA--CB  C--C    C--O    C--CB   O--O    O--CB   CB--CB
                    // 2.75  2.42    2.99    2.24    2.90     3.80    2.42    2.74    3.43    3.14    3.12    2.79    3.43    2.46    3.69

                    // Others v.s. Gly
                    // 75% coverage
                    // N(G)--N(O)   N(G)--CA(O) N(G)--C(O)  N(G)--O(O)  N(G)--CB(O) CA(G)--N(O) CA(G)--CA(O)    CA(G)--C(O) CA(G)--O(O) CA(G)--CB(O)    C(G)--N(O)  C(G)--CA(O) C(G)--C(O)  C(G)--O(O)  C(G)--CB(O) O(G)--N(O)  O(G)--CA(O) O(G)--C(O)  O(G)--O(O)  O(G)--CB(O)
                    //    2.62        2.39        2.99         2.22       2.97         2.38         3.74           2.38       2.66         3.68            2.99        2.38       2.90        2.80         3.08       2.22        2.66         2.81        3.07       3.09    

                    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //
                    //    manually tunned values used in wefold
                    //
                    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    // N--N  N--CA   N--C    N--O    N--CB   CA--CA  CA--C   CA--O   CA--CB  C--C    C--O    C--CB   O--O    O--CB   CB--CB
                    // 2.60  2.36    2.99    2.21    2.76    3.06    2.35    2.74    3.59    2.84    2.67    2.96    2.94    2.89    3.50
                    

                    // Real R1R2_pow[5][5] = {6.96,5.80,8.88,4.97,6.96,5.80,14.28,5.76,7.67,14.13,8.88,5.76,8.46,7.56,7.56,4.97,7.67,7.56,9.30,5.85,6.96,14.13,7.56,5.85,13.10};
                    // Real R1R2_pow[5][5] = {6.96,5.80,8.88,4.97,6.96,5.80,14.28,5.76,7.67,14.13,8.88,5.76,8.46,7.56,7.56,4.97,7.67,7.56,9.30,5.85,6.96,14.13,7.56,5.85,13.10};
                    Real clash=0;
                    for( Size iatom=0; iatom<5; ++iatom ) {
                        for( Size jatom=0; jatom<5; ++jatom ) {
                            Vec  iatom_xyz = pose.xyz(ires, atom_types[iatom], energy_clash->ichain());
                            Vec  jatom_xyz = pose.xyz(jres, atom_types[jatom], energy_clash->jchain());
                            Real dist = (iatom_xyz - jatom_xyz).squaredNorm();
                            Real s = _R1R2_pow[iatom][jatom] - dist;
                            if(s>0.0) clash += s;
                        }
                    }
                    energy_clash->set_pair_score(ires,jres,clash*10);
                }
            }
            // calculate clash between sidechain of motif residues and backbone of generate scaffold
            // Tooooooooooooooooooooooooo slowwwwwwwwwwwwwwwwwwwwwwww
            // TODO: bug!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            const bool sc_bb_clash = false;
            if(sc_bb_clash) {
                for(Size jres=ires+1; jres<=len; ++jres){
                    Real clash_backbone = energy_clash->get_pair_score(ires,jres);
                    Real clash = 0;
                    for( Size iatom=0; iatom<5; ++iatom ) {
                        for (AtomInfo const & j_atom_info : conformation.atom_info(jres)){
                            Vec iatom_xyz = pose.xyz(ires, atom_types[iatom], energy_clash->ichain());
                            Vec jatom_xyz = j_atom_info.global_postion;
                            Real dist = (iatom_xyz - jatom_xyz).squaredNorm();
                            Size j_atom_index = std::distance(atom_types.begin(),std::find(atom_types.begin(), atom_types.end(), j_atom_info.atom_type));
                            Real s = _R1R2_pow[iatom][j_atom_index] - dist;
                            if(s>0.0) clash += s;
                        }
                    }
                    energy_clash->set_pair_score(ires,jres,clash_backbone+clash*10);
                }
            }

        }
        clash_total+= energy_clash->weighted_score();
    }
    std::shared_ptr<scene::Pose> target_op = pose.get_target_pose();
    if(target_op!= nullptr && _scaffold_target_inter_clash){
        std::string const target_sequence = target_op->sequence();
        Size target_len = target_op->size();
        for(Size ires=1; ires<=len; ++ires) {
            for(Size jres=1; jres<=target_len; ++jres) {
                Vec  res1_cb_xyz = pose.xyz(ires, ATOM_CB);
                Vec  res2_cb_xyz = target_op->xyz(jres, ATOM_CB);
                Vec t = res1_cb_xyz-res2_cb_xyz;
                if(t.squaredNorm() > _squared_dist_cutoff+20) {
                    continue;
                } else {         
                    std::vector<ATOM_TYPE> atom_types = {ATOM_N, ATOM_CA, ATOM_C, ATOM_O, ATOM_CB};
                    Real R1R2_pow[5][5] = {6.76,5.57,8.94,4.88,7.62,5.57, 9.36,5.52,7.51,12.88,8.94,5.52,8.07,7.13,8.76,4.88,7.50,7.13,8.64,8.35,7.61,12.89,8.76,8.35,12.25+_interface_extra_punish/2};
                    Real clash=0;
                    for( Size iatom=0; iatom<5; ++iatom ) {
                        // Gly has no CB
                        if(sequence.at(ires-1) == 'G' && iatom == 4) continue;
                        
                        for( Size jatom=0; jatom<5; ++jatom ) {
                            // Gly has no CB
                            if(target_sequence.at(jres-1) == 'G' && jatom == 4) continue;

                            Vec  iatom_xyz = pose.xyz(ires, atom_types[iatom]);
                            Vec  jatom_xyz = target_op->xyz(jres, atom_types[jatom]);
                            Real dist = (iatom_xyz - jatom_xyz).squaredNorm();
                            Real s = (R1R2_pow[iatom][jatom]+_interface_extra_punish/2) - dist;
                            if(s>0.0) clash += s;
                        }
                    }
                    clash_total+=clash*_clash_multiplier;
                }
            }
        }
    }
    return clash_total;
}

} // namespace scoring


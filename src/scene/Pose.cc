#include "scene/Residue.hh"
#include "scene/Conformation.hh"
#include "scene/Pose.hh"
#include "basic/macros.hh"
#include "utils/random_util.hh"
#include "utils/utils.hh"
#include "utils/string_util.hh"
#include "utils/math_util.hh"
#include "basic/global.hh"

#include <fstream>
//#include <unordered_map>
#include <map>
#include <algorithm>
#include <iostream>

namespace scene {

using namespace basic;

Pose::Pose(Size nres, Size num_repeats, std::string symmetry) : 
    _nres(nres),
    _num_repeats(num_repeats),
    _symmetry(symmetry),
    _conformation(nres,symmetry),
    _energy_manager(nres,symmetry,num_repeats),
    _target_poseOP(nullptr),
    _voxel_clash_free_residues(nres,false),
    _freeze_root(false)
{
    initialize();
    _conformation.set_sequence(_sequence );
}
Pose::Pose(const Pose& pose,bool copy_enery,bool copy_trajctory)
{
    
    _nres = pose._nres;
    _num_repeats = pose._num_repeats;
    _num_chains = pose._num_chains;
    _symmetry = pose._symmetry;
    _sequence = pose._sequence;
    _ss = pose._ss;
    _conformation = Conformation(pose._conformation);
    // fix this op copy here when need
    if(copy_enery)_energy_manager=scoring::EnergyManager(pose._energy_manager);
    _fixed_segments = pose._fixed_segments;
    _target_poseOP = pose._target_poseOP==nullptr?nullptr:std::make_shared<Pose>(Pose(*pose._target_poseOP,false,false));
    if(copy_trajctory)_trajactory=std::vector<scene::Pose>(pose._trajactory);
}
Pose::Pose(std::string pdb_name,char chain ,std::string symmetry,Size repeat_num):
    _num_repeats(repeat_num)
{
    std::map<std::string,ATOM_TYPE> aa1toatom= {{"N",ATOM_N},{"C",ATOM_C},{"O",ATOM_O},{"S",ATOM_S},{"H",ATOM_H}};
    Size segment_len(0);
    std::string segment_seq("");

    // load pdb file
    // only need N, CA, C
    // O, H, CB, use ideal length and angle
    
    std::vector<std::string> pdb_lines;
    std::vector<std::vector<AtomInfo>> all_atom_info;
    std::ifstream pdb_file(pdb_name);
    {
        Size cur_residue_number=0;
        std::string cur_residue_origin_number;
        for(std::string line;std::getline(pdb_file,line);)
        {   
            if( line.substr(0,4) == "ATOM" &&(chain==' ' || line.at(21)==chain) )
            {
                if(cur_residue_origin_number!=line.substr(22,4)){
                    cur_residue_origin_number=line.substr(22,4);
                    cur_residue_number++;
                    all_atom_info.push_back(std::vector<AtomInfo>());
                }
                if(line.substr(13,2)=="N " || line.substr(13,3)=="CA " || line.substr(13,2)=="C ")
                {
                    pdb_lines.push_back(line);
                    if( line.substr(13,3)=="CA " ){
                        ++segment_len;
                        segment_seq += BASIC_AA3to1[line.substr(17,3)];
                    }
                }
                else{
                    AtomInfo atom_info;
                    if(line.substr(77,1)=="H")continue;
                    atom_info.atom_name = utils::trim(line.substr(12,4));
                    atom_info.res_name = line.substr(17,3);
                    atom_info.local_postion = Vec(std::stof(line.substr(30,8)),std::stof(line.substr(38,8)),std::stof(line.substr(46,8)));
                    atom_info.element = utils::trim(line.substr(77,2));
                    atom_info.atom_type = aa1toatom[line.substr(77,1)];
                    all_atom_info.back().push_back(atom_info);

                }
            }
        }
    }
    pdb_file.close();

    _nres = segment_len;
    // num * [N, CA, C]
    assert(pdb_lines.size() == 3*segment_len);
    assert(segment_len==all_atom_info.size());
    // attention here!! conformation and energy manager should init before the pose intialize

    _conformation = Conformation(_nres,symmetry);
    _sequence = segment_seq;
    _symmetry=symmetry;
    _energy_manager  = scoring::EnergyManager(_nres,symmetry,1);
    _voxel_clash_free_residues = std::vector<bool>(_nres,false);
    initialize();

    _conformation.set_sequence(_sequence );

    _fixed_segments.clear();
    snapshot();
    Size start_res = 1;
    Size end_res   = segment_len;
    // TODO: add repeat protein support. 2021-10-15

    // extract coordinates
    std::vector<Vec> N, CA, C;

    //
    for(Size i=0;i<pdb_lines.size();i++){
        switch (i % 3) {
            case 0:
                N.push_back(Vec(std::stof(pdb_lines[i].substr(30,8)),std::stof(pdb_lines[i].substr(38,8)),std::stof(pdb_lines[i].substr(46,8))));
                break;
            case 1:
                CA.push_back(Vec(std::stof(pdb_lines[i].substr(30,8)),std::stof(pdb_lines[i].substr(38,8)),std::stof(pdb_lines[i].substr(46,8))));
                break;
            case 2:
                C.push_back(Vec(std::stof(pdb_lines[i].substr(30,8)),std::stof(pdb_lines[i].substr(38,8)),std::stof(pdb_lines[i].substr(46,8))));
                break;
        }
    }
    for(Size ii=1; ii<segment_len; ++ii) {
        set_omega( ii+1, utils::get_dihedral( CA[ii-1], C[ii-1], N[ii],  CA[ii] ) );
        set_phi(   ii+1, utils::get_dihedral( C[ii-1],  N[ii],   CA[ii], C[ii]) );
        set_psi(  ii, utils::get_dihedral( N[ii-1], CA[ii-1], C[ii-1], N[ii]));
    }


    // set the global xform of that residue
    for(Size ii=1; ii<=segment_len; ++ii){
        EigenXform ca_stub = utils::xform_from_3points(N[ii-1],CA[ii-1],C[ii-1]);
        _conformation.set_stub(ii, ATOM_CA, ca_stub );
        for(auto & atom_info : all_atom_info[ii-1]){
            atom_info.local_postion = ca_stub.inverse()*atom_info.local_postion;
            _conformation.push_atom_info(ii,atom_info);
        }
    }
    for(Size ii=2; ii<=segment_len; ++ii){
        _conformation.set_stub(ii, ATOM_N, utils::xform_from_3points(C[ii-2],N[ii-1],CA[ii-1]));
    }
    for(Size ii=1; ii<segment_len; ++ii) {
        _conformation.set_stub(ii, ATOM_C,  utils::xform_from_3points(CA[ii-1],C[ii-1],N[ii]) );
    }
    // this end c stub should be correct evem without the below code but a coordinate update but in fact no
    EigenXform end_stub = EigenXform(EigenXform::Identity());
    end_stub.translation() = C[segment_len-1];
    _conformation.set_stub(segment_len, ATOM_C,  end_stub );
    end_stub.translation() = N[0];
    _conformation.set_stub(1, ATOM_N,  end_stub );
    for(Size ii=2; ii<=segment_len; ++ii) {
        _conformation.set_stub(ii, ATOM_CX,  utils::xform_from_3points(CA[ii-2],C[ii-2],N[ii-1]) );
    }
    for(Size ii=1; ii<=segment_len; ++ii){
        bool update_N_CA = ii!=1?true:false;
        bool update_C = ii!=segment_len?true:false;
        _conformation.update_local_xform_from_global_xform(ii,update_N_CA,update_N_CA,update_C);
    }
    set_root_index(_conformation.center_residue());
    _freeze_root = true;
    // update_coordinates();
}

Pose::~Pose(){}

void Pose::initialize()
{
    //
    assert(_nres % _num_repeats == 0);

    // initialize the sequence and dssp
    if(_sequence.size()!=_nres)_sequence = std::string(_nres,'G');
    _ss = std::string(_nres,'L');


    update_coordinates();
    snapshot();
    //with nomal initialize, also could set root at middle without load pdb
    set_root_index(_nres/2,false);
    _fixed_segments.clear();

    // get the number of chains
    _num_chains = _conformation.num_chains();

    // initialize energy table here
    _energy_manager.add_energy_twobody_intra(CLASH,_num_repeats);
    _energy_manager.add_energy_twobody_intra(RPX,_num_repeats);
    for(Size i =2;i<=_num_chains;i++){
        Real weight = _conformation.chain_weights()[i-2];
        if(weight==0)continue;
        _energy_manager.add_energy_twobody_inter_symmetry(CLASH,1,i,weight);
        _energy_manager.add_energy_twobody_inter_symmetry(RPX,1,i,weight);
    }
}

void Pose::reset_coords()
{
    for( Size ires=1; ires<=_nres; ++ires){
        if(_ss.at(ires-1) == 'H') {
            set_omega(ires, 3.14);
            set_phi(ires, -0.99483);
            set_psi(ires, -0.8203);
        } else if (_ss.at(ires-1) == 'E') {
            set_omega(ires, 3.14);
            set_phi(ires, -2.36);
            set_psi(ires, 2.36);
        } else {
            set_omega(ires, 3.14);
            set_phi(ires, 3.14);
            set_psi(ires, 3.14);
        }
    }

    if(_num_chains>1 && (!_freeze_root) ) random_root();
    _conformation.update_coordinates();
    _energy_manager.marker_changed_segment_moved(1, _nres, _conformation.root_index());
}

void Pose::extend_pose(Size Nter_extension, Size Cter_extension)
{
    // by default, the extended regions are set as loop and GGG
    _nres += Nter_extension + Cter_extension;
    _ss = std::string(Nter_extension, 'L') + _ss + std::string(Cter_extension, 'L');
    _sequence = std::string(Nter_extension, 'G') + _sequence + std::string(Cter_extension, 'G');

    // fixed segments
    // TODO: for repeat proteins
    // but for repeat proteins, what's the logic for pose extension
    std::set<std::pair<Size, Size> > old_segments = _fixed_segments;
    _fixed_segments.clear();
    for(std::pair<Size, Size> p : old_segments) {
        _fixed_segments.insert(std::pair<Size, Size>(p.first+Nter_extension, p.second+Nter_extension));
    }

    _conformation.extend_conformation(Nter_extension, Cter_extension);


    // energy manager
    // TODO: other energy tables
    // OneBody
    _energy_manager = scoring::EnergyManager(_nres,_symmetry,_num_repeats);
    _energy_manager.add_energy_twobody_intra(CLASH,_num_repeats);
    _energy_manager.add_energy_twobody_intra(RPX,_num_repeats);
    for(Size i =2;i<=_num_chains;i++){
        Real weight = _conformation.chain_weights()[i-2];
        if(weight==0)continue;
        _energy_manager.add_energy_twobody_inter_symmetry(CLASH,1,i,weight);
        _energy_manager.add_energy_twobody_inter_symmetry(RPX,1,i,weight);
    }
    _energy_manager.marker_changed_segment_moved(1, _nres, _conformation.root_index());

    update_coordinates();
    snapshot();
}

// torsions: (omega, phi, psi)*3
void
Pose::insert_fragment(Size ires, std::vector<Real> const & torsions, bool update_coordinates)
{
    Size frag_size = torsions.size()/3;
    assert(0 == torsions.size()%3);
    Size len_per_repeat = _nres / _num_repeats;
    for( Size ii=0; ii<frag_size; ++ii){
        Size insert_site = (ires+ii-1)%len_per_repeat+1;
        for(Size jj=0; jj<_num_repeats; ++jj) {
            set_omega(insert_site + jj*len_per_repeat, torsions[ii*3+0]);
            set_phi(insert_site + jj*len_per_repeat, torsions[ii*3+1]);
            set_psi(insert_site + jj*len_per_repeat, torsions[ii*3+2]);
        }
    }

    
    if(update_coordinates) { _conformation.update_coordinates(); }

    // deal with the fixed segments to speed up the score calculation
    std::vector<bool> temp_moved(_nres, false);
    Size start_res((ires-1)%len_per_repeat+1), end_res((ires+frag_size-1-1)%len_per_repeat+1);
    if(start_res <= end_res) {
        for(Size ii=0; ii<_num_repeats; ++ii) {
            std::fill(temp_moved.begin()+start_res+ii*len_per_repeat-1, temp_moved.begin()+ii*len_per_repeat+end_res, true);
        }
    } else {
        for(Size ii=0; ii<_num_repeats; ++ii) {
            std::fill(temp_moved.begin()+1+ii*len_per_repeat-1, temp_moved.begin()+end_res+ii*len_per_repeat, true);
            std::fill(temp_moved.begin()+start_res+ii*len_per_repeat-1, temp_moved.begin()+len_per_repeat+ii*len_per_repeat, true);
        }
    }

    // revert the fixed segments to false
    for(std::pair<Size, Size> const & s : _fixed_segments) {
        //maybe a bug : the two side residue's torsion angle of fixed segment could change, means this two residue also changed
        // std::fill(temp_moved.begin()+s.first-1, temp_moved.begin()+s.second, false);
        if(s.second<=s.first+1)continue;
        std::fill(temp_moved.begin()+s.first-1+1, temp_moved.begin()+s.second-1, false);
    }

    // scan the array to find moved segments
    bool previous = false;
    start_res = 0;
    end_res = 0;
    for(Size idx=0; idx<temp_moved.size(); ++idx) {
        if( !previous && temp_moved[idx] ){
            start_res  = idx+1;
        }
        if( previous && !temp_moved[idx] ){
            end_res = idx;
            _energy_manager.marker_changed_segment_moved(start_res, end_res, _conformation.root_index());
            start_res = 0;
            end_res = 0;
        }
        previous = temp_moved[idx];
    }
    if(previous)
        _energy_manager.marker_changed_segment_moved(start_res, temp_moved.size(), _conformation.root_index());   
}

void 
Pose::freeze_segment(Size start_pos, Size end_pos)
{
    _fixed_segments.insert(std::pair<Size, Size>(start_pos, end_pos));

    freeze_torsions(start_pos,false,false,true);
    for(Size i_tmp=start_pos+1; i_tmp<=end_pos-1; ++i_tmp) {
        freeze_torsions(i_tmp, true, true, true);
    }
    freeze_torsions(end_pos,true,true,false);
}

void
Pose::perturb_root(Real ang_mag, Real trans_mag, Size mode)
{
    _conformation.perturb_root(ang_mag, trans_mag, mode);
    // only if there are multiple chains

    _energy_manager.marker_changed_segment_moved(1,_nres,_conformation.root_index(),ONE_BODY);
    if(_num_chains>1){
        _energy_manager.marker_changed_segment_moved(1,_nres,_conformation.root_index(),TWO_BODY_INTER_CHAIN);
    }
}

void
Pose::random_root(bool perturb_ori, bool perturb_x, bool perturb_y, bool perturb_z)
{
    _conformation.random_root( perturb_ori, perturb_x, perturb_y, perturb_z);

    _energy_manager.marker_changed_segment_moved(1,_nres,_conformation.root_index(),ONE_BODY);

    if(_symmetry.at(0) == 'C') {
        //
        if(_num_chains>1 && (perturb_ori || perturb_x || perturb_y)) {
            _energy_manager.marker_changed_segment_moved(1,_nres,_conformation.root_index(),TWO_BODY_INTER_CHAIN);
        }
    } else if(_symmetry.at(0) == 'D') {
        //
        if((!perturb_ori) && (!perturb_x) && (!perturb_y) && perturb_z) {
            _energy_manager.marker_changed_segment_moved_by_chain(1,_nres,_conformation.root_index(),_num_chains/2+1,_num_chains);
        } else if(perturb_ori || perturb_x || perturb_y) {
            _energy_manager.marker_changed_segment_moved(1,_nres,_conformation.root_index(),TWO_BODY_INTER_CHAIN);
        } else {
            // WHAT???
            exit(0);
        }
    } else if (_symmetry.at(0) == 'H') {
        _energy_manager.marker_changed_segment_moved(1,_nres,_conformation.root_index(),TWO_BODY_INTER_CHAIN);
    } else {
        // NOT IMPLEMENTED ERROR!
        exit(0);
    }
}

void
Pose::random_jump(Real ang_mean, Real trans_mean)
{
    _conformation.random_jump(ang_mean, trans_mean);
    // only if there are multiple chains
    // _energy_manager.marker_changed_segment_moved(1,_nres,_conformation.root_index(),ONE_BODY);
    if(_num_chains>1){
        _energy_manager.marker_changed_segment_moved(1,_nres,_conformation.root_index(),TWO_BODY_INTER_CHAIN);
    }
}

void
Pose::perturb_jump(Real ang_mag, Real trans_mag)
{
    _conformation.perturb_jump(ang_mag, trans_mag);

    // only if there are multiple chains
    // _energy_manager.marker_changed_segment_moved(1,_nres,_conformation.root_index(),ONE_BODY);
    if(_num_chains>1){
        _energy_manager.marker_changed_segment_moved(1,_nres,_conformation.root_index(),TWO_BODY_INTER_CHAIN);
    }
}

void
Pose::load_pdb(std::string const & fname, Size insert_pos, bool freeze_cart, bool freeze_torsions, bool update,bool voxel_clash_free) // [a,b]
{
    std::map<std::string,ATOM_TYPE> aa1toatom= {{"N",ATOM_N},{"C",ATOM_C},{"O",ATOM_O},{"S",ATOM_S},{"H",ATOM_H}};
    Size segment_len(0);
    std::string segment_seq("");


    std::vector<std::string> raw_pdb_lines;
    bool is_gz;
    if   ( fname.substr(fname.size()-4,4)==".pdb" )         is_gz = false;
    else if ( fname.substr(fname.size()-7,7)==".pdb.gz" )   is_gz = true;
    else                                                  { std::cout << "Unknown Format!" << std::endl; exit(0); }

    if(is_gz) {
        igzstream pdb_file(fname.c_str(), std::ios_base::in);
        if(pdb_file.fail()) {
            std::cout << "Failed to load pdb: " << fname << std::endl;
            exit(0);
        }
        for(std::string line;std::getline(pdb_file,line);) {
            raw_pdb_lines.push_back(line);
        }
        pdb_file.close();
    } else {
        std::ifstream pdb_file(fname, std::ios_base::in);
        if(pdb_file.fail()) {
            std::cout << "Failed to load pdb: " << fname << std::endl;
            exit(0);
        }
        for(std::string line;std::getline(pdb_file,line);) {
            raw_pdb_lines.push_back(line);
        }
        pdb_file.close();
    }

    // load pdb file
    // only need N, CA, C
    // O, H, CB, use ideal length and angle
    std::vector<std::string> pdb_lines;

    Size cur_residue_number=0;
    std::string cur_residue_origin_number;
    std::vector<std::vector<AtomInfo>> all_atom_info;
    for(std::string line : raw_pdb_lines )
    {   

        // ATOM line
        if( line.substr(0,4) == "ATOM"  )
        {
            if(cur_residue_origin_number!=line.substr(22,4)){
                cur_residue_origin_number=line.substr(22,4);
                cur_residue_number++;
                all_atom_info.push_back(std::vector<AtomInfo>());
            }
            if(line.substr(13,2)=="N " || line.substr(13,3)=="CA " || line.substr(13,2)=="C ")
            {
                pdb_lines.push_back(line);
                if( line.substr(13,3)=="CA " ){
                    ++segment_len;
                    segment_seq += BASIC_AA3to1[line.substr(17,3)];
                }
            }
            else{
                AtomInfo atom_info;
                if(line.substr(77,1)=="H"||line.substr(13,2)=="O ")continue;
                atom_info.atom_name = utils::trim(line.substr(12,4));
                atom_info.res_name = line.substr(17,3);
                atom_info.local_postion = Vec(std::stof(line.substr(30,8)),std::stof(line.substr(38,8)),std::stof(line.substr(46,8)));
                atom_info.element = utils::trim(line.substr(77,3));
                atom_info.atom_type = aa1toatom[line.substr(77,1)];
                all_atom_info.back().push_back(atom_info);

            }
        }

        // REMARK line
        // ATOM line
        if( line.substr(0,21) == "REMARK PDBinfo-LABEL:"  )
        {
            //
            std::cout << "The extra pdb info of the loaded pdb: " << line << std::endl; 
            std::vector<std::string> temp_split = utils::string_split(line, ' ');
            for(auto s : temp_split) {
                if(s.substr(0,12) == "LIGAND_XFORM") {
                    _pdb_info.push_back(s);
                }
            }
        }
    }

    if(segment_len==0)
    {
        std::cout<<"Error: no CA atom in the pdb file"<<std::endl;
        exit(1);
    }
    // num * [N, CA, C]
    assert(pdb_lines.size() == 3*segment_len);
    assert(insert_pos+segment_len-1 <= _nres);
    // consider the special case for repeat proteins
    // label the fixed segments
    Size len_per_repeat = _nres / _num_repeats;
    // make sure the length of the loaded segment is shorter than the length of the repeat
    assert(segment_len <= len_per_repeat);     

    // update some information in the pose
    // update the sequence using the segment seq
       
    Size start_res = (insert_pos-1)%len_per_repeat+1;
    Size end_res   = (insert_pos+segment_len-1-1)%len_per_repeat+1;
    if(start_res <= end_res) {
        for(Size ii=0; ii<_num_repeats; ++ii) {
            _fixed_segments.insert(std::pair<Size, Size>(start_res+ii*len_per_repeat,end_res+ii*len_per_repeat));
            set_sequence(segment_seq, start_res+ii*len_per_repeat);
            
            // TODO: also consider repeat protein
            // only works for repeat==1
            if(voxel_clash_free) {
                for(Size itmp=start_res; itmp<=end_res; ++itmp) {
                    set_voxel_clash_free(itmp+ii*len_per_repeat, true);    
                }
            }
        }

    } else {
        // this is stupid.
        for(Size ii=0; ii<_num_repeats; ++ii) {
            _fixed_segments.insert(std::pair<Size, Size>(1+ii*len_per_repeat,end_res+ii*len_per_repeat));
            _fixed_segments.insert(std::pair<Size, Size>(start_res+ii*len_per_repeat,len_per_repeat+ii*len_per_repeat));
            set_sequence(segment_seq.substr(segment_len-end_res, end_res), 1+ii*len_per_repeat);
            set_sequence(segment_seq.substr(0, segment_len-end_res), start_res+ii*len_per_repeat);
        }
    }

    // TODO: add repeat protein support. 2021-10-15

    // extract coordinates
    std::vector<Vec> N, CA, C;

    //
    for(Size i=0;i<pdb_lines.size();i++){
        switch (i % 3) {
            case 0:
                N.push_back(Vec(std::stof(pdb_lines[i].substr(30,8)),std::stof(pdb_lines[i].substr(38,8)),std::stof(pdb_lines[i].substr(46,8))));
                break;
            case 1:
                CA.push_back(Vec(std::stof(pdb_lines[i].substr(30,8)),std::stof(pdb_lines[i].substr(38,8)),std::stof(pdb_lines[i].substr(46,8))));
                break;
            case 2:
                C.push_back(Vec(std::stof(pdb_lines[i].substr(30,8)),std::stof(pdb_lines[i].substr(38,8)),std::stof(pdb_lines[i].substr(46,8))));
                break;
        }
    }

    // TODO: code is ugly

    // set the torsion angles for that residue
    // omega, phi, psi

    // wraped the old code with a for loop and the residue index is updated by
    // (XXX-1)%len_per_repeat+1 + irepeat*len_per_repeat
    // not compatable with repeat,but nerver mind
    for(Size ii=1; ii<=segment_len; ++ii){
        EigenXform ca_stub = utils::xform_from_3points(N[ii-1],CA[ii-1],C[ii-1]);
        for(auto & atom_info : all_atom_info[ii-1]){
            atom_info.local_postion = ca_stub.inverse()*atom_info.local_postion;
            _conformation.push_atom_info(ii+insert_pos-1,atom_info);
        }
    }
    for(Size irepeat=0; irepeat<_num_repeats; ++irepeat) {
        set_omega((insert_pos-1)%len_per_repeat+1 + irepeat*len_per_repeat, 3.14);
        set_phi(  (insert_pos-1)%len_per_repeat+1 + irepeat*len_per_repeat, 3.14);
        for(Size ii=1; ii<segment_len; ++ii) {
            set_omega( (insert_pos+ii-1)%len_per_repeat+1 + irepeat*len_per_repeat, utils::get_dihedral( CA[ii-1], C[ii-1], N[ii],  CA[ii] ) );
            set_phi(   (insert_pos+ii-1)%len_per_repeat+1 + irepeat*len_per_repeat, utils::get_dihedral( C[ii-1],  N[ii],   CA[ii], C[ii]) );
            set_psi(   (insert_pos+ii-1-1)%len_per_repeat+1 + irepeat*len_per_repeat, utils::get_dihedral( N[ii-1], CA[ii-1], C[ii-1], N[ii]));
        }
        set_psi((insert_pos+segment_len-1-1)%len_per_repeat+1 + irepeat*len_per_repeat, 3.14);
        // set the global xform of that residue
        for(Size ii=0; ii<segment_len; ++ii){
            _conformation.set_stub((insert_pos+ii-1)%len_per_repeat+1 + irepeat*len_per_repeat, ATOM_CA, utils::xform_from_3points(N[ii],CA[ii],C[ii]) );
        }
        for(Size ii=1; ii<segment_len; ++ii){
            _conformation.set_stub((insert_pos+ii-1)%len_per_repeat+1 + irepeat*len_per_repeat, ATOM_N, utils::xform_from_3points(C[ii-1],N[ii],CA[ii]));
        }
        for(Size ii=0; ii<segment_len-1; ++ii) {
            EigenXform stub_temp = utils::xform_from_3points(CA[ii],C[ii],N[ii+1]);
            _conformation.set_stub((insert_pos+ii-1)%len_per_repeat+1 + irepeat*len_per_repeat, ATOM_C, stub_temp );
            // change the stub of ATOM_CX of the next residue
            if((insert_pos+ii-1)%len_per_repeat+1 + irepeat*len_per_repeat + 1 <= _nres)
                _conformation.set_stub((insert_pos+ii-1)%len_per_repeat+1 + irepeat*len_per_repeat + 1, ATOM_CX, stub_temp);
        }

        // update the local xform and freeze the torsions
        if( freeze_torsions ) {
            if(segment_len > 1) {
            	_conformation.update_local_xform_from_global_xform((insert_pos-1)%len_per_repeat+1 + irepeat*len_per_repeat, false, false, true);
            	_conformation.freeze_torsions( (insert_pos-1)%len_per_repeat+1 + irepeat*len_per_repeat, false, false, true );
            }
            for(Size ii=1; ii<segment_len-1; ++ii) {
                _conformation.update_local_xform_from_global_xform((insert_pos+ii-1)%len_per_repeat+1 + irepeat*len_per_repeat, true, true, true);
                _conformation.freeze_torsions( (insert_pos+ii-1)%len_per_repeat+1 + irepeat*len_per_repeat, true, true, true );
            }
            if(segment_len > 1) {
            	_conformation.update_local_xform_from_global_xform((insert_pos+segment_len-1-1)%len_per_repeat+1 + irepeat*len_per_repeat, true, true, false);
            	_conformation.freeze_torsions( (insert_pos+segment_len-1-1)%len_per_repeat+1 + irepeat*len_per_repeat, true, true, false );
        	}
        }

    }

    if( !freeze_torsions ) {
        _fixed_segments.clear();
    }

    // root
    if( freeze_cart ) {

        _freeze_root = true;

        //with update of repeat residues, idx should be in the middle of all repeat 
        // Size idx = insert_pos + (segment_len/2) + (_num_repeats/2)*len_per_repeat;
        //but after a thought, it seems no much optimization to do this, so just leave the root at first repeat for a clear define 
        // but with another thought, this way below also can make the root to second repeat, so i just put it at middle wahtever
        // Size idx = insert_pos + (segment_len/2);
        Size idx = insert_pos + (segment_len/2) + (_num_repeats/2)*len_per_repeat;
        set_root_index(idx, false);
        // this is used to update the multiple chain coord
        // not a good design
        if( _num_chains > 1)
            set_root_index(idx, false); // ??? True or False ???
    }

    // update coordinates
    if(update){
        update_coordinates();
        for(Size ires=1; ires<=_nres; ++ires){
            _conformation.update_twin_residue(ires);
        }
    }
}

void
Pose::dump_pdb(std::string const & fname,bool dump_target,bool dump_gzip,bool dump_trajactory,bool append, std::string const & extra_info)
{
    ogzstream f_gz;
    std::ofstream f;
    if(dump_gzip){
        if(!dump_target&&!append){
            f_gz.open(fname.c_str());
        }else{
            // gzstream can not append, but why ? maybe not as simple as normal format, might be fixed but not now
            // no append nor read/write mode
        // if ((mode & std::ios::ate) || (mode & std::ios::app)
        //     || ((mode & std::ios::in) && (mode & std::ios::out)))
        // return (gzstreambuf*)0;
            f_gz.open(fname.c_str(),std::ios::out);
        }
    }else{
        if(!dump_target&&!append){
            f.open(fname.c_str());
        }else{
            f.open(fname.c_str(),std::ios::app);
        }
    }
    if(dump_trajactory){
        if(dump_gzip){
            f_gz<<dump_pdb(false,false,true).str();
            f_gz.close();
        }else{
            f<<dump_pdb(false,false,true).str();
            f.close();
        }
        return;
    }
    std::vector<ATOM_TYPE> atom_types = {ATOM_N, ATOM_CA, ATOM_C, ATOM_O, ATOM_CB, ATOM_H};
    // std::vector<std::string> atom_names = {"N", "CA", "C", "O", "CB", "H"};
    std::vector<std::string> atom_names = {"N", "CA", "C", "O", "CB"};
    // std::vector<std::string> atom_type_names = {"N", "C", "C", "O", "C", "H"};
    std::vector<std::string> atom_type_names = {"N", "C", "C", "O", "C"};
    
    //chain identifiers
    std::vector<char> chain_identifiers = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};

    char buf[128];
    // std::cout<<traj<<std::endl;
    // _conformation.update_coordinates();
    if(dump_target) chain_identifiers={{'B'}};
    Size counter = 0;
    for(Size ichain=1; ichain<=_num_chains; ++ichain) {
        for(Size ires=1; ires<=_nres; ++ires) {
            char aa = _sequence.at(ires-1);
            for(Size ii = 0; ii<atom_names.size(); ++ii) {
                bool real_info=false;
                for(auto atom_info:_conformation.atom_info(ires)){
                    if(atom_info.atom_name==atom_names[ii]){
                        real_info=true;
                        break;
                    }
                }
                if(real_info)continue;
                Vec pos = _conformation.xyz(ires, atom_types[ii], ichain);
                //gly no CB
                if ('G' == aa && ii == 4) continue;
                // atom number
                ++counter;
                snprintf(buf, 128, "%6s%5d %*s%*s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n",
                    "ATOM  ",
                    counter,
                    CALC_CENTER_POSITION_PREV(4, atom_names[ii].c_str()),
                    atom_names[ii].c_str(),
                    CALC_CENTER_POSITION_POST(4, atom_names[ii].c_str()),
                    "",
                    ' ',
                    BASIC_AA1to3[aa].c_str(),
                    chain_identifiers[ichain-1],
                    ires,
                    ' ',
                    pos(0),
                    pos(1),
                    pos(2),
                    1.0,
                    0.0,
                    atom_type_names[ii].c_str(),
                    " "
                );

                if(dump_gzip)f_gz<<buf;
                else f<<buf;

            }
            EigenXform const & ca_stub = _conformation.stub(ires,ichain,ATOM_CA);
            for(auto & atom_info:_conformation.atom_info(ires)){
                Vec pos = atom_info.local_postion;
                pos = ca_stub*pos;
                ++counter;
                snprintf(buf, 128, "%6s%5d %*s%*s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n",
                    "ATOM  ",
                    counter,
                    CALC_CENTER_POSITION_PREV(4, atom_info.atom_name.c_str()),
                    atom_info.atom_name.c_str(),
                    CALC_CENTER_POSITION_POST(4, atom_info.atom_name.c_str()),
                    "",
                    ' ',
                    BASIC_AA1to3[aa].c_str(),
                    chain_identifiers[ichain-1],
                    ires,
                    ' ',
                    pos(0),
                    pos(1),
                    pos(2),
                    1.0,
                    0.0,
                    atom_info.element.c_str(),
                    " "
                );


                if(dump_gzip)f_gz<<buf;
                else f<<buf;
            }
            
        }
    }
    if(_fixed_segments.size()!=0){
        for(auto segment:_fixed_segments){
            if(segment.second-segment.first+1==_nres) continue;
            for(Size i = segment.first;i<=segment.second;i++){
                if(dump_gzip) f_gz<<"REMARK PDBinfo-LABEL:  "<<i<<" MOTIF"<<std::endl;
                else f<<"REMARK PDBinfo-LABEL:  "<<i<<" MOTIF"<<std::endl;
            }
        }
    }
    if(!dump_gzip && extra_info != "") {
        f << extra_info << std::endl;
    }
    f.close();
    if(_target_poseOP!=nullptr && !dump_gzip)_target_poseOP->dump_pdb(fname,true,dump_gzip);
    if(_target_poseOP!=nullptr && dump_gzip)f_gz<<_target_poseOP->dump_pdb(true).str();
    if(dump_gzip && extra_info != "") {
        f_gz << extra_info << std::endl;
    }
    f_gz.close();
    //f << "REMARK total_score: " << _energy.total_score() / ( _nres * _num_chains ) << std::endl;
    //f << "REMARK sc_neighbors: " << sidechain_neighbors(*this) << std::endl;

}
std::stringstream
Pose::dump_pdb(bool dump_target,bool dump_gzip,bool dump_trajactory,bool append)
{
    std::stringstream context;
    std::vector<ATOM_TYPE> atom_types = {ATOM_N, ATOM_CA, ATOM_C, ATOM_O, ATOM_CB, ATOM_H};
    // std::vector<std::string> atom_names = {"N", "CA", "C", "O", "CB", "H"};
    std::vector<std::string> atom_names = {"N", "CA", "C", "O", "CB"};
    // std::vector<std::string> atom_type_names = {"N", "C", "C", "O", "C", "H"};
    std::vector<std::string> atom_type_names = {"N", "C", "C", "O", "C"};
    
    //chain identifiers
    std::vector<char> chain_identifiers = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};

    char buf[128];
    for(Size traj=0;traj<(dump_trajactory?_trajactory.size():1);traj++){
        // std::cout<<traj<<std::endl;
        if(dump_trajactory){

            snprintf(buf, 128, "%6s%6s\n","MODEL ",std::to_string(traj+1).c_str());
            context<<buf;
            context<<(_trajactory[traj].dump_pdb(false,dump_gzip,false,true)).str();
            snprintf(buf, 128, "%6s\n","ENDMDL");
            context<<buf;
            continue;
        }
        // _conformation.update_coordinates();
        if(dump_target) chain_identifiers={{'B'}};
        Size counter = 0;
        for(Size ichain=1; ichain<=_num_chains; ++ichain) {
            for(Size ires=1; ires<=_nres; ++ires) {
                char aa = _sequence.at(ires-1);
                for(Size ii = 0; ii<atom_names.size(); ++ii) {
                    bool real_info=false;
                    for(auto atom_info:_conformation.atom_info(ires)){
                        if(atom_info.atom_name==atom_names[ii]){
                            real_info=true;
                            break;
                        }
                    }
                    if(real_info)continue;
                    Vec pos = _conformation.xyz(ires, atom_types[ii], ichain);
                    //gly no CB
                    if ('G' == aa && ii == 4) continue;
                    // atom number
                    ++counter;
                    snprintf(buf, 128, "%6s%5d %*s%*s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n",
                        "ATOM  ",
                        counter,
                        CALC_CENTER_POSITION_PREV(4, atom_names[ii].c_str()),
                        atom_names[ii].c_str(),
                        CALC_CENTER_POSITION_POST(4, atom_names[ii].c_str()),
                        "",
                        ' ',
                        BASIC_AA1to3[aa].c_str(),
                        chain_identifiers[ichain-1],
                        ires,
                        ' ',
                        pos(0),
                        pos(1),
                        pos(2),
                        1.0,
                        0.0,
                        atom_type_names[ii].c_str(),
                        " "
                    );

                    context<<buf;

                }
                EigenXform  ca_stub = EigenXform(_conformation.stub(ires,ichain,ATOM_CA));
                for(auto & atom_info:_conformation.atom_info(ires)){
                    Vec pos = atom_info.local_postion;
                    pos = ca_stub*pos;
                    ++counter;
                    snprintf(buf, 128, "%6s%5d %*s%*s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n",
                        "ATOM  ",
                        counter,
                        CALC_CENTER_POSITION_PREV(4, atom_info.atom_name.c_str()),
                        atom_info.atom_name.c_str(),
                        CALC_CENTER_POSITION_POST(4, atom_info.atom_name.c_str()),
                        "",
                        ' ',
                        BASIC_AA1to3[aa].c_str(),
                        chain_identifiers[ichain-1],
                        ires,
                        ' ',
                        pos(0),
                        pos(1),
                        pos(2),
                        1.0,
                        0.0,
                        atom_info.element.c_str(),
                        " "
                    );


                    context<<buf;
                }
                
            }
        }
        if(_target_poseOP!=nullptr)context<<(_target_poseOP->dump_pdb(true,dump_gzip)).str();
    }
    return context;
    //f << "REMARK total_score: " << _energy.total_score() / ( _nres * _num_chains ) << std::endl;
    //f << "REMARK sc_neighbors: " << sidechain_neighbors(*this) << std::endl;

}

}

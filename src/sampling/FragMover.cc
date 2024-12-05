#include "scene/Pose.hh"
#include "basic/macros.hh"
#include "sampling/FragMover.hh"
#include "basic/assert.hh"
#include "utils/random_util.hh"

#include <ctime>
#include <random>

#include <fstream>
#include <map>

namespace sampling {

using basic::Real;
using basic::Size;

FragMover::FragMover(Size protein_len, Size frag_len) :
    _protein_len(protein_len),
    _frag_len(frag_len),
    _runif(0.0,1.0),
    _distribution_ires(1,protein_len-frag_len+1)
{
}

FragMover::~FragMover() {}

Size
FragMover::frag_length() const {
    return _frag_len;
}
void
FragMover::frag_length(Size l) {

    //std::cout << "Why do you want to dynamically change the length of a fragment?
    //Undefined behavior will happen." << std::endl;
    //exit(0);

    _frag_len = l;
    _distribution_ires = std::uniform_int_distribution<Size>(1,_protein_len-_frag_len+1);
}
Size
FragMover::protein_length() const {
    return _protein_len;
}
void FragMover::protein_length(Size l) {
    _protein_len = l;
    _distribution_ires = std::uniform_int_distribution<Size>(1,_protein_len-_frag_len+1);
}

void FragMover::bias_loop_sampling_by_dssp( std::string const dssp, Real loop_sampling_weight, Size include_adjacent_residues )
{

    ALWAYS_ASSERT(dssp.length() == _protein_len);

    std::vector<char> dsspSS(_protein_len);
    for( Size ii=0; ii<_protein_len; ++ii) dsspSS[ii] = dssp[ii];

    /*
    if ( dsspSS[0]=='L' && dsspSS[1]!='L' ) { //often dssp mis-assignes residues at the beginning and end of chain.
        dsspSS[0]=dsspSS[1];
    }
    if ( protein_length>2 ) {
        if ( dsspSS[_protein_len-1]=='L' && dsspSS[_protein_len-2] != 'L' ) {
            dsspSS[_protein_len-1]=dsspSS[_protein_len-2];
        }
    }
    */

    std::vector<Size> adjacent_pos;
    for(Size ii=0; ii<dsspSS.size(); ++ii) {
        if(dssp[ii]=='L') {
            Size start = ii-include_adjacent_residues;
            if(start<0) start=0;
            Size end = ii+include_adjacent_residues;
            if(end>=_protein_len) end=_protein_len-1;
            for(Size jj=start;jj<=end;++jj) {
                adjacent_pos.push_back(jj);
            }
        }
    }
    for(Size ii : adjacent_pos) dsspSS[ii] = 'L';

    // generate weights for each position
    std::vector<Real> sampling_weights(_protein_len);
    for(Size ii=0; ii<_protein_len; ++ii) {
        if(dsspSS[ii]=='L')
            sampling_weights[ii] = loop_sampling_weight;
        else if (dsspSS[ii]=='X')
            sampling_weights[ii] = 0.0;
        else
            sampling_weights[ii] = 1.0;
    }

    // call the function
    bias_sampling_by_weight(sampling_weights);
}

void FragMover::bias_sampling_by_weight(std::vector<Real> & weights)
{
    ALWAYS_ASSERT(weights.size() == _protein_len)

    Size vec_len = _protein_len - _frag_len + 1;

    std::vector<Real> sampling_weights_on_frags(vec_len);
    for( Size ii = 0; ii < vec_len; ++ii ) {
        sampling_weights_on_frags[ii] = 0.0;
        for( Size i_shift=0; i_shift<_frag_len; ++i_shift ) {
            sampling_weights_on_frags[ii] += weights[ii+i_shift];
        }
    }
    Real weights_sum = 0.0;
    for( Size ii=0; ii<vec_len; ++ii ) {
        weights_sum += sampling_weights_on_frags[ii];
    }
    for( Size ii=0; ii<vec_len; ++ii ) {
        sampling_weights_on_frags[ii] = sampling_weights_on_frags[ii] * vec_len / weights_sum;
    }
    // generate outputs
    _sampling_alias.clear();
    _sampling_alias.resize(vec_len);
    _sampling_probs.clear();
    _sampling_probs.resize(vec_len);


    Size cur_small_block;
    Size cur_large_block;
    Size num_small_block = 0;
    Size num_large_block = 0;
    std::vector<Real> small_block(vec_len);
    std::vector<Real> large_block(vec_len);
    for( Size ires=vec_len-1; ires>=0; --ires ) {
        if(sampling_weights_on_frags[ires]<1) {
            small_block[num_small_block++] = ires;
        } else {
            large_block[num_large_block++] = ires;
        }
    }
    while ( num_small_block && num_large_block ) {
        cur_small_block = small_block[--num_small_block];
        cur_large_block = large_block[--num_large_block];
        _sampling_probs[cur_small_block] = sampling_weights_on_frags[cur_small_block];
        _sampling_alias[cur_small_block] = cur_large_block;
        sampling_weights_on_frags[cur_large_block] = sampling_weights_on_frags[cur_large_block] + sampling_weights_on_frags[cur_small_block] - 1;
        if(sampling_weights_on_frags[cur_large_block]<1){
            small_block[num_small_block++] = cur_large_block;
        } else {
            large_block[num_large_block++] = cur_large_block;
        }
    }

    while (num_large_block)
        _sampling_probs[large_block[--num_large_block]] = 1;

    while (num_small_block)
        _sampling_probs[small_block[--num_small_block]] = 1;
}

Size FragMover::pick_position(bool bias) 
{
    if( bias ) {
        Real rand1 = _runif(utils::rng()); // wtf
        Size rand2 = _distribution_ires(utils::rng())-1; // 0-index vs 1-index
        return (rand1 < _sampling_probs[rand2] ? rand2 : _sampling_alias[rand2])+1;
    } else {
        return _distribution_ires(utils::rng());
    }
}

VanillaFragMover::VanillaFragMover(Size protein_len, Size frag_len, Size frag_num) : 
    FragMover(protein_len, frag_len),
    _frag_num(frag_num),
    _distribution_ifrag(1,frag_num)
{}

VanillaFragMover::~VanillaFragMover() {}

void VanillaFragMover::load_fragments(std::string const & fname)
{
    _frags.resize(_protein_len - _frag_len + 1);
    for(Size ii=0; ii<_frags.size(); ++ii)
    {
        _frags[ii].resize(_frag_num);
        for(Size jj=0; jj<_frag_num; ++jj) {
            _frags[ii][jj].resize(_frag_len*3);
        }
    }

    Size count(0);
    std::ifstream f(fname);
    Real omega, phi, psi;
    while(!f.eof()){
        f >> omega >> phi >> psi;

        ++count;
        Size x = (count-1) / (_frag_num*_frag_len);
        Size residue = (count-1) % (_frag_num*_frag_len);
        Size y = residue / _frag_len;
        Size z = residue % _frag_len;

        _frags[x][y][z*3+0] = omega;
        _frags[x][y][z*3+1] = phi;
        _frags[x][y][z*3+2] = psi;
    }

    assert(count == (_protein_len - _frag_len + 1 ) * _frag_num * _frag_len );
}

void VanillaFragMover::apply(scene::Pose & pose, Size ires)
{
    Size jfrag = _distribution_ifrag(utils::rng());

    pose.insert_fragment(ires, _frags[ires-1][jfrag-1]);

}

FragMapMover::FragMapMover(Size protein_len, Size frag_len) :
    FragMover(protein_len, frag_len)
{}
FragMapMover::~FragMapMover()
{}

void FragMapMover::apply(scene::Pose & pose, Size ires)
{
    // get dssp and hash
    uint64_t hash = hash_ss(pose.dssp().substr(ires-1, _frag_len));

    if(_frag_map[hash].size()==0){
        std::cerr<<"get a strange dssp frag which have no frag in library: "<<pose.dssp().substr(ires-1, _frag_len)<<std::endl;
        return;
    }
    
    Size ifrag = std::uniform_int_distribution<Size>(0,_frag_map[hash].size()-1)(utils::rng());
    pose.insert_fragment(ires, _frag_map[hash][ifrag]);
}

bool FragMapMover::apply_with_ifrag(scene::Pose & pose, Size ires, Size ifrag)
{
    uint64_t hash = hash_ss(pose.dssp().substr(ires-1, _frag_len));

    if(_frag_map[hash].size()==0){
        std::cerr<<"get a strange dssp frag which have no frag in library: "<<pose.dssp().substr(ires-1, _frag_len)<<std::endl;
        return false;
    }

    if( ifrag != -1 && ifrag > _frag_map[hash].size() ) {
        return false;
    }

    if(ifrag == -1) {
        ifrag = std::uniform_int_distribution<Size>(1,_frag_map[hash].size())(utils::rng());
    }

    pose.insert_fragment(ires, _frag_map[hash][ifrag-1]);

    return true;
}


uint64_t FragMapMover::hash_ss(std::string const & ss){

    uint64_t hash=0;
    for(Size i =0 ;i<ss.size();i++){
        hash *= 10;
        hash += ss.at(i)=='H'?0:(ss.at(i)=='L'?1:2);
    }
    return hash;

}

// std::string FragMapMover::reverse_ss_hash(Size hash,Size len){
//     std::string ss;
//     for(Size i=0,k=pow(10,len-1);i<len;i++,k/=10){
//         ss.push_back(hash/k ==0?'H':(hash/k==1?'L':'E'));
//         if(hash/k!=0)hash =hash%k;
//     }
//     return ss;
// }

void FragMapMover::load_fragments(std::string const & fname)
{

    // this is because of the hashing procedure
    // overflow
    // not too big
    // What't the maximum length of a fragment to be hashed as int??
    // it might be OK to be larger
    assert(_frag_len <= 7);
    
    std::ifstream f(fname, std::ios::binary);
    
    _frag_map.clear();
    while(!f.eof()){
        
        int32_t counts = 0;
        uint64_t hash = 0;
        f.read((char*)&hash,sizeof(int64_t)); // to be compatible with the python script
        f.read((char*)&counts,sizeof(int32_t)); // to be compatible with the python script
        // std::cout<<reverse_ss_hash(hash,_frag_len)+":"<<counts<<std::endl;
        if(counts==0)continue;
        //std::cout << "Hash: " << hash << " Counts: " << counts << std::endl;
        _frag_map[hash].resize(counts);
        for(int32_t i =0;i<counts;i++){
            _frag_map[hash][i].resize(_frag_len*3);
            for(Size j=0;j<3*_frag_len;j++){
                float torsion;
                f.read((char*)&torsion, sizeof(float));
                _frag_map[hash][i][j] = Real(torsion);
            }
        }
    }
    assert(_frag_map.size()>0);
}

}

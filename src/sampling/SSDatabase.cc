#include "scene/Pose.hh"
#include "basic/macros.hh"
#include "sampling/SSDatabase.hh"
#include "basic/assert.hh"
#include "utils/random_util.hh"

#include <ctime>
#include <random>

#include <fstream>
#include <map>

namespace sampling {

using basic::Real;
using basic::Size;

SSDatabase::SSDatabase(Size protein_len) :
    _protein_len(protein_len),
    _num(0)
{
}

SSDatabase::~SSDatabase() {}

void
SSDatabase::load_database(std::string const & fname) {

    Size count(0);
    std::ifstream f(fname);
    std::string prefix, dssp;
    Size pos;
    while(!f.eof()){
        f >> prefix >> dssp >> pos;

        if( dssp.length() != _protein_len ) {
            std::cout << "The length of dssp str doesn't not match the protein length!!!!" << std::endl;
            exit(0);
        }

        ++count;
        
        _dssp_prefix.push_back(prefix);
        _dssp.push_back(dssp);
        _motif_pos.push_back(pos);
    }
    _num = count;

    _distribution_iss = std::uniform_int_distribution<Size>(0,count-1);
}

void
SSDatabase::fetch_dssp_config(std::string & prefix, std::string & dssp, Size & motif_pos)
{
    Size which_iss = _distribution_iss(utils::rng());
    prefix = _dssp_prefix[which_iss];
    dssp = _dssp[which_iss];
    motif_pos = _motif_pos[which_iss];
}

}
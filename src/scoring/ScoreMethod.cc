#include "scoring/ScoreMethod.hh"

#include <cmath>
#include <string>
#include <fstream>
#include <Eigen/Geometry>
#include <parallel_hashmap/phmap_dump.h>
#include <map>
namespace scoring {

using namespace basic;
using namespace scene;

// basic score method

void BaseScoreMethod::set_designable_res(const std::string & designable_res_str) {

    std::vector< std::string > parts;
    Size i(0), j(0);
    while ( j != std::string::npos ) {
        j = designable_res_str.find( ',', i );
        std::string const part = designable_res_str.substr(i,j-i);
        parts.push_back( part );
        i = j+1;
    }

    _designable_res.clear();
    for(const std::string s : parts) {
        if( s == "") continue;
        _designable_res.insert(stoi(s));
    }
}

void BaseScoreMethod::set_designable_res(std::vector<Size> designable_res) {
    _designable_res.clear();
    for(Size i : designable_res ) {
        _designable_res.insert(i);
    }
}


} // namespace scoring


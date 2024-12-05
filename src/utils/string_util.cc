#include "basic/types.hh"
#include "utils/string_util.hh"

#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>

namespace utils {

using namespace basic;

std::vector<std::string> string_split(const std::string & in, char splitchar)
{
    std::vector< std::string > parts;
    Size i(0), j(0);
    while ( j != std::string::npos ) {
        j = in.find( splitchar, i );
        std::string const part = in.substr(i,j-i);
        parts.push_back( part );
        i = j+1;
    }
    return parts;
}


}

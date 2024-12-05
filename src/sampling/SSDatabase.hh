#ifndef INCLUDED_sampling_SSDatabase_hh
#define INCLUDED_sampling_SSDatabase_hh

#include "basic/types.hh"
#include "scene/Pose.hh"

#include <memory>
#include <random>

namespace sampling {

using basic::Real;
using basic::Size;

class SSDatabase
{
    public:

        SSDatabase(Size protein_len);
        ~SSDatabase();

        void load_database(std::string const & fname);

        void fetch_dssp_config(std::string & prefix, std::string & dssp, Size & moitf_pos);

    protected:
        Size _protein_len;
        Size _num;
        std::vector<std::string> _dssp_prefix;
        std::vector<std::string> _dssp;
        std::vector<Size> _motif_pos;

        std::uniform_int_distribution<Size> _distribution_iss;
};

}

#endif

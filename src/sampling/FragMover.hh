#ifndef INCLUDED_sampling_FragMover_hh
#define INCLUDED_sampling_FragMover_hh

#include "basic/types.hh"
#include "scene/Pose.hh"

#include <memory>
#include <random>

namespace sampling {

using basic::Real;
using basic::Size;

class FragMover
{
    public:

        FragMover(Size protein_len, Size frag_len=7);
        ~FragMover();

        virtual void load_fragments(std::string const & fname) = 0;

        Size frag_length() const;
        void frag_length(Size l);
        Size protein_length() const;
        void protein_length(Size l);

        void bias_loop_sampling_by_dssp( std::string const dssp, Real loop_sampling_weight, Size include_adjacent_residues=0 /*0 means not including the adjacent residues*/);
        void bias_sampling_by_weight(std::vector<Real> & weights);
        Size pick_position(bool bias);

        virtual void apply(scene::Pose & pose, Size ires) = 0;

    protected:
        Size _protein_len;
        Size _frag_len;
        std::uniform_real_distribution<Real> _runif;
	    std::vector<Size> _sampling_alias;
	    std::vector<Real> _sampling_probs;

        std::uniform_int_distribution<Size> _distribution_ires;
};

class VanillaFragMover : public FragMover
{
    public:

        VanillaFragMover(Size protein_len, Size frag_len, Size frag_num=200);
        ~VanillaFragMover();

        void load_fragments(std::string const & fname) override;

        void apply(scene::Pose & pose, Size ires);
    private:
        
        Size _frag_num;
        
        std::vector<std::vector<std::vector<Real> > > _frags;

        std::uniform_int_distribution<Size> _distribution_ifrag;
};

class FragMapMover : public FragMover
{
    public:
        FragMapMover(Size protein_len, Size frag_len=7);
        ~FragMapMover();

        void load_fragments(std::string const & fname) override;

        void apply(scene::Pose & pose, Size ires);

        bool apply_with_ifrag(scene::Pose & pose, Size ires, Size ifrag=-1);

    private:
        uint64_t hash_ss(std::string const & ss);
        //std::string reverse_ss_hash(Size hash,Size len);

    private:
        std::unordered_map<uint64_t,std::vector<std::vector<Real>>> _frag_map;
};

}

#endif

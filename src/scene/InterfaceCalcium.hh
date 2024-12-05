#ifndef INCLUDED_scene_InterfaceCalcium_hh
#define INCLUDED_scene_InterfaceCalcium_hh

#include "basic/types.hh"

#include "scene/InterfaceMetal.hh"

#include <set>
#include <vector>
#include <memory>

namespace scene {

using namespace basic;


class InterfaceCalcium : public InterfaceMetal
{
    public:
        InterfaceCalcium();
        ~InterfaceCalcium();

    private:
        void load_chain1_metal_configs(std::string json_fname) override;
        void load_chain2_metal_configs(std::string json_fname) override;
        void compute_pairwise_distance() override;
};

typedef std::shared_ptr<InterfaceCalcium> InterfaceCalciumOP;

}

#endif

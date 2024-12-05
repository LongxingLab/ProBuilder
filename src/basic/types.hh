#ifndef INCLUDED_basic_types_hh
#define INCLUDED_basic_types_hh

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <unordered_map>
#include <parallel_hashmap/phmap.h>

#include <xbin/xbin.hpp>

namespace basic {

typedef int Size;
typedef float Real;
typedef Eigen::Transform<Real, 3, Eigen::AffineCompact > EigenXform;
typedef Eigen::AngleAxis<Real> AngleAxis;
typedef Eigen::Matrix<Real,3,1> Vec;
typedef Eigen::Matrix<Real, 3, 3> Mat;



typedef xbin::XformHash_bt24_BCC6<EigenXform, uint64_t> wefoldHasher;


template <typename K, typename V>
using hash_table = phmap::flat_hash_map<K, V >;

template <typename K, typename V>
using hash_set_table = phmap::flat_hash_map<K, phmap::flat_hash_set<V> >;



enum ATOM_TYPE { ATOM_CX, ATOM_N, ATOM_CA, ATOM_C , ATOM_O, ATOM_CB, ATOM_H,ATOM_S,ATOM_P};

// enum ENERGY_NAME {};
// RIF,SIDECHAIN_TARGET,SIDECHAIN,CLASH_WITH_VOXEL,RPX_TARGET,RPX_INTRA,RPX_INTER
// enum ENERGY_TYPE {RESIDUE,INTRA_CHAIN_PAIR,INTER_CHAIN_PAIR};
enum ENERGY_TYPE {ENERGY_TYPE_UNDEFINED,WHOLE_BODY,ONE_BODY,TWO_BODY_INTRA_CHAIN,TWO_BODY_INTER_CHAIN};

// enum PAIR_TYPE {PAIR_UNDEFINED,INTRA_CHAIN_PAIR,INTER_CHAIN_PAIR};

enum SCORE_TYPE {SCORE_TYPE_UNDEFINED,
                 SIDECHAIN_TARGET,
                 SIDECHAIN,
                 RPX,RPX1SIDE,
                 PREBUILD_RPX1SIDE,
                 CLASH,
                 DISULFIDE,
                 PRIVILEGED_MOTIF,
                 PRIVILEGED_INTERFACE_MOTIF,
                 PRIVILEGED_PAIR,
                 TARGET_RIF,
                 TARGET_RIF2,
                 TARGET_CONTEXT_CLASH,
                 TARGET_CONTEXT_CLASH2, 
                 METAL_COORDINATION};

enum CONSTRAINT_TYPE {ATOMPAIR_CONSTRAINT};

enum SS_TYPE {HELIX,BETASTRAND,LOOP,UNDEFINED_SS};

typedef std::pair <Size,ATOM_TYPE> AtomID;

static std::unordered_map<std::string,ATOM_TYPE> atom_name_map{
    {"CA",ATOM_CA},{"N",ATOM_N},{"C",ATOM_C},{"CB",ATOM_CB},{"O",ATOM_O},{"H",ATOM_H}};

}

#endif


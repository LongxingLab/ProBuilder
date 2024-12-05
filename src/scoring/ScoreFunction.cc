#include "scoring/ScoreFunction.hh"
#include "utils/hash_util.hh"

#include <cmath>
#include <string>
#include <fstream>
#include <Eigen/Geometry>
#include <parallel_hashmap/phmap_dump.h>
#include <map>
namespace scoring {

using namespace basic;
using namespace scene;

void
ScoreFunction::regist_method(std::string method_name,ScoreMethodOP methodOP){
    if(_method_map.count(method_name)!=0){
        std::cout<<"override an exist score method named : "<< method_name << std::endl;
    }
    _method_map[method_name] = methodOP;
}

void
ScoreFunction::delete_method(std::string method_name){
    if(_method_map.count(method_name)==0){
        std::cout<<"delete an unexist score method named : "<< method_name << std::endl;
        std::exit(0);
    }
    _method_map.erase(method_name);
}

ScoreMethodOP
ScoreFunction::get_method(std::string method_name){
    if(_method_map.count(method_name)==0){
        std::cout<<"cant find score method :" <<method_name<<std::endl;
        exit(0);
    }
    return _method_map.at(method_name);
}

std::vector<ScoreMethodOP>
ScoreFunction::get_all_methods () const {
    std::vector<ScoreMethodOP> all_score_methods;
    for(auto pair: _method_map){
        all_score_methods.push_back(pair.second);
    }
    return all_score_methods;
}

Real
ScoreFunction::score(scene::Pose & pose){
    Real total_score=0;
    for(auto pair: _method_map){
        total_score += pair.second->weighted_score(pose);
    }
    return total_score;
}

Real
ScoreFunction::score(scene::Pose & pose,std::vector<std::string> exclude_methods){
    Real total_score=0;
    for(auto pair: _method_map){
        if(std::find(exclude_methods.begin(),exclude_methods.end(),pair.first)!=exclude_methods.end())continue;
        total_score += pair.second->weighted_score(pose);
    }
    return total_score;
}

} // namespace scoring


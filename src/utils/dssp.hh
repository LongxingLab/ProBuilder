#ifndef INCLUDED_utils_dssp_hh
#define INCLUDED_utils_dssp_hh

#include "apps/args.hh"
#include "basic/types.hh"
#include "utils/random_util.hh"
#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>

namespace utils {

using basic::Size;

class Dssp{
    public:
    Dssp(Options opts){
        std::cout<<"generate dssp "<<std::endl;
        if(opts.segment_num!=-1){
            segments.resize(1);
            for(Size i =0;i<opts.segment_num;i++){segments[0].push_back('H');segments[0].push_back('L');}
            segments[0].pop_back();
        }else{
            segments.resize(4);
            for(Size i =3;i<=6;i++){
                for(Size j =0;j<i;j++){
                    segments[i].push_back('H');segments[i].push_back('L');
                }
                segments[i].pop_back();
            }
        }
        for(auto segment:segments){
            f1("",segment,0,opts);
        }
        std::cout<<"dssp sample num: "<<all_dssp.size()<<std::endl;
        assert(all_dssp.size()>0);
    };
    void search_dssp(std::string segments,Options opts);
    void f1 (std::string ss_now,std::string & segments,Size idx,Options & opts);
    std::string pick_dssp(){return all_dssp[random_int(0,all_dssp.size()-1)];}
    std::vector<std::string> all_dssp;
    std::vector<std::string> segments;
};

void Dssp::f1 (std::string ss_now,std::string & segments,Size idx,Options  & opts){
        
        std::string ss_append;
        if((opts.max_len<ss_now.size()))return;
        if(idx==segments.size()){
            all_dssp.push_back(ss_now);
            std::cout<<ss_now<<std::endl;
            return;
        }
        char this_ss = segments.at(idx);
        // std::cout<<ss_now<<std::endl;
        // std::cout<<this_ss<<std::endl;
        if(this_ss=='L'){
            for(Size i =opts.min_loop_len;i++;i<=opts.max_loop_len){
                f1(ss_now+(std::string(i,'L')),segments,idx+1,opts);
                // std::cout<<ss_now.append(std::string(i,'L'))<<std::endl;
            }
        }
        if(this_ss=='H'){
            for(Size i =opts.min_helix_len;i++;i<=opts.max_helix_len){
                f1(ss_now+(std::string(i,'H')),segments,idx+1,opts);
            }
        }
        else{
            std::cout<<"unkonw segment"<<std::endl;
            exit(0);

        }
}

}

#endif

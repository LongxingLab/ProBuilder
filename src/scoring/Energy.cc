#include "scoring/Energy.hh"

#include <numeric>
#include <iostream>

namespace scoring {

using namespace basic;

void OneBodyEnergy::marker_changed_segment_moved(Size start_res, Size end_res, Size root_res){
    if(end_res<root_res)std::fill(_element_moved.begin(),_element_moved.begin()+end_res,true);
    else if(start_res>root_res)std::fill(_element_moved.begin()+start_res-1,_element_moved.end(),true);
    else std::fill(_element_moved.begin(),_element_moved.end(),true);
}

void IntraTwobodyEnergy::marker_changed_segment_moved(Size start_res, Size end_res, Size root_res) 
{

    Size repeat_len = _nres/_num_repeats;
    start_res = (start_res-1)%(_nres/_num_repeats)+1;
    end_res = end_res>repeat_len?repeat_len:end_res;//make start and end residue to first repeat
    root_res = (root_res-1)%(_nres/_num_repeats)+1;
    
    //first block cal
    for(Size ii=1; ii<start_res;++ii)
        for(Size jj=start_res; jj<=repeat_len;++jj)
            _element_moved[(ii-1)*_nres+jj-1] = true;
    for(Size ii=start_res; ii<=end_res;++ii)
        for(Size jj=ii+1;jj<=repeat_len;++jj)
            _element_moved[(ii-1)*_nres+jj-1] = true;
    //second block have some unchanged residue pair but third block and etc have no unchanged pair, so i calculate
    //whole block of 2,3,4.. for code to read easily
    for(Size ii=1;ii<=repeat_len;ii++){
        for(Size jj=repeat_len+1;jj<=_nres;jj++){
            _element_moved[(ii-1)*_nres+jj-1] = true;
        }
    }
}

Real IntraTwobodyEnergy::total_score(){
    _total_score = 0;
    if(_num_repeats==1)_total_score = std::accumulate(_element_score.begin(), _element_score.end(), 0.0);
    else{
        Size len_per_repeat = _nres/_num_repeats;
        for(Size ii = 0; ii<len_per_repeat;ii++){
            for(Size jj=0; jj<_num_repeats;jj++){
                _total_score += std::accumulate(_element_score.begin()+ii*_nres, _element_score.begin()+(ii+1)*_nres-jj*len_per_repeat, 0.0);
            }
        }
    }
    return _total_score;
}

void SymetryInterTwobodyEnergy::marker_changed_segment_moved(Size start_res, Size end_res, Size root_res){
    if( start_res == end_res && start_res == root_res ) {
        //
        std::fill(_element_moved.begin(), _element_moved.end(), true);
        // early stop
        return;
    } else {
        //

        if( start_res >= root_res ) {
            // not very efficient??
            for(Size ii=1; ii<=_nres; ++ii)
                for(Size jj=1; jj<=_nres; ++jj)
                    if(ii>=start_res || jj>=start_res)
                        _element_moved[(ii-1)*_nres+jj-1] = true;
        } 
        if( end_res <= root_res) {
            //
            for(Size ii=1; ii<=_nres;++ii)
                for(Size jj=1; jj<=_nres; ++jj)
                    if(ii<=end_res || jj<=start_res)
                        _element_moved[(ii-1)*_nres+jj-1] = true;
        } 
        if( start_res <= root_res && end_res >= root_res ){
            std::fill(_element_moved.begin(), _element_moved.end(), true);
        }
    }
}

}

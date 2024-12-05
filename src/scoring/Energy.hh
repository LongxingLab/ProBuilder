#ifndef INCLUDED_scoring_Energy_hh
#define INCLUDED_scoring_Energy_hh

#include "basic/types.hh"

#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <memory>
#include <string>

namespace scoring {

using namespace basic;

class BaseEnergy
{
    public:
        BaseEnergy(Size nres,Real weight=1):_nres(nres),_weight(weight){};

        // moved position mark function should rewrite based on diffrent energy type 
        virtual void marker_changed_segment_moved(Size start_res, Size end_res, Size root_res) = 0 ; // [a,b]

        // score funtion should rewrite for diffrent energy type
        //rewrite method for one body enery
        ENERGY_TYPE energy_type(){return _energy_type;}
        void set_energy_type(ENERGY_TYPE energy_type){_energy_type=energy_type;}
        SCORE_TYPE score_type(){return _score_type;}
        void set_score_type(SCORE_TYPE score_type){_score_type=score_type;}
        // PAIR_TYPE pair_type(){return _pair_type;}

        // rewrite method for two body energy
        void set_weight(Real weight){_weight=weight;}

        // auto rewrite for diffrent enerygy type
        virtual Real total_score(){return std::accumulate(_element_score.begin(),_element_score.end(),0.0);};
        Real weighted_score(){return _weight*total_score();};
        virtual void snapshot()=0;
        virtual void rollback()=0;

    protected:
        std::vector<bool> _element_moved;
        std::vector<bool> _element_moved_snapshot;
        std::vector<Real> _element_score;
        std::vector<Real> _element_score_snapshot;
        Real _total_score;
        Real _total_score_snapshot;
        Real _weight;
        Size _nres;

        SCORE_TYPE _score_type;
        ENERGY_TYPE _energy_type;
        // onebody energy type must be intra and actually pair_type belongs to twobody energy,but for conveniece 
        // to use diffrent intra and inter moved mark funtion
        // PAIR_TYPE _pair_type;
};

class OneBodyEnergy: public BaseEnergy{
    public:
    OneBodyEnergy(Size nres,SCORE_TYPE score_type,Size ichain=1,Real weight=1):BaseEnergy(nres,weight){
        // _energy_type = ONE_BODY;
        //_pair_type = PAIR_UNDEFINED;
        set_score_type(score_type);
        set_energy_type(ONE_BODY);
        _ichain = ichain;
        //_score_type = score_type;
        _element_moved.resize(_nres,true);
        _element_moved_snapshot.resize(_nres,true);
        _element_score.resize(_nres,0);
        _element_score_snapshot.resize(_nres,0);
    };
    bool is_changed_residue(Size ires){
        return _element_moved[ires-1];
    };
    void set_score_per_residue(Size ires,Real value){
        _element_score[ires-1] = value;
        _element_moved[ires-1] = false;
    };
    Real get_score_per_residue(Size ires){
        return _element_score[ires-1];
    };
    Size ichain(){return _ichain;};
    virtual void snapshot() override{
        std::copy(_element_moved.begin(), _element_moved.end(), _element_moved_snapshot.begin());
        std::copy(_element_score.begin(), _element_score.end(), _element_score_snapshot.begin());
        _total_score_snapshot = _total_score;

    };
    virtual void rollback() override{
        std::copy(_element_moved_snapshot.begin(), _element_moved_snapshot.end(), _element_moved.begin());
        std::copy(_element_score_snapshot.begin(), _element_score_snapshot.end(), _element_score.begin());
        _total_score = _total_score_snapshot;

    };
    virtual void marker_changed_segment_moved(Size start_res, Size end_res, Size root_res) override;
    protected:
    Size _ichain;

};

class TwoBodyEnergy: public BaseEnergy{

    public:

    TwoBodyEnergy(Size nres,SCORE_TYPE score_type,ENERGY_TYPE energy_type,Real weight=1):BaseEnergy(nres,weight){
        // _pair_type = pair_type;
        // _energy_type = TWO_BODY;
        set_score_type(score_type);
        set_energy_type(energy_type);
        _element_moved.resize(_nres*_nres,true);
        _element_score_snapshot.resize(_nres*_nres,0);
        _element_score.resize(_nres*_nres,0);
        _element_moved_snapshot.resize(_nres*_nres,true);
        _weight = weight;
    };
    
    bool is_changed_pair(Size ires, Size jres){
        return _element_moved[(ires-1)*_nres+jres-1];
    };
    void set_pair_score(Size ires, Size jres, Real value){
        _element_score[(ires-1)*_nres+jres-1] = value;
        _element_moved[(ires-1)*_nres+jres-1] = false;
    };
    Real get_pair_score(Size ires, Size jres){
        return _element_score[(ires-1)*_nres+jres-1];
    };
    Size ichain() const {return _ichain;};
    Size jchain() const {return _jchain;};
    virtual void snapshot() override{
        std::copy(_element_moved.begin(), _element_moved.end(), _element_moved_snapshot.begin());
        std::copy(_element_score.begin(), _element_score.end(), _element_score_snapshot.begin());
        _total_score_snapshot = _total_score;

    };
    virtual void rollback() override{
        std::copy(_element_moved_snapshot.begin(), _element_moved_snapshot.end(), _element_moved.begin());
        std::copy(_element_score_snapshot.begin(), _element_score_snapshot.end(), _element_score.begin());
        _total_score = _total_score_snapshot;

    };


    protected:

    Size _ichain,_jchain;

};

class InterTwobodyEnergy:public TwoBodyEnergy{
    public:

    InterTwobodyEnergy(Size nres,Size ichain,Size jchain,SCORE_TYPE score_type,Real weight=1):TwoBodyEnergy(nres,score_type,TWO_BODY_INTER_CHAIN,weight)
    {
        _ichain=ichain;
        _jchain=jchain;
    };
    

    protected:
    
};

class SymetryInterTwobodyEnergy:public InterTwobodyEnergy{
    public:
    using InterTwobodyEnergy::InterTwobodyEnergy;

    virtual void marker_changed_segment_moved(Size start_res, Size end_res, Size root_res) override;

    protected:
    
};

class IntraTwobodyEnergy:public TwoBodyEnergy{
    public:

    IntraTwobodyEnergy(Size nres,SCORE_TYPE score_type,Size num_repeats=1,Real weight=1):TwoBodyEnergy(nres,score_type,TWO_BODY_INTRA_CHAIN,weight),_num_repeats(num_repeats){
        _ichain=1;
        _jchain=1;
    };
    virtual void marker_changed_segment_moved(Size start_res, Size end_res, Size root_res) override;
    virtual Real total_score() override;

    private:
    Size _num_repeats;

};

typedef std::shared_ptr<BaseEnergy> EnergyOP;
typedef std::shared_ptr<TwoBodyEnergy> TwoBodyEnergyOP;
typedef std::shared_ptr<OneBodyEnergy> OneBodyEnergyOP;

class EnergyManager{
    
    public:
    EnergyManager(){};
    EnergyManager(Size nres, std::string symmetry, Size num_repeats):_nres(nres),_symmetry(symmetry),_num_repeats(num_repeats){};
    void add_energy_onebody(SCORE_TYPE score_type,Real weight=1.0){

        if(energy_map.count(score_type)==0){
            energy_map[score_type] = std::vector<EnergyOP>(1,EnergyOP(new OneBodyEnergy(_nres,score_type,weight)));
        }else{
            energy_map[score_type].push_back(EnergyOP(new OneBodyEnergy(_nres,score_type,weight)));
        }

    }
    void add_energy_twobody_intra(SCORE_TYPE score_type,Size num_repeats =1,Real weight=1.0){

        if(energy_map.count(score_type)==0){
            energy_map[score_type] = std::vector<EnergyOP>(1,EnergyOP(new IntraTwobodyEnergy(_nres,score_type,weight)));

        }else{
            energy_map[score_type].push_back(EnergyOP(new IntraTwobodyEnergy(_nres,score_type,num_repeats,weight)));
        }

    }
    void add_energy_twobody_inter_symmetry(SCORE_TYPE score_type,Size ichain,Size jchain,Real weight = 1.0){

        if(energy_map.count(score_type)==0){
            energy_map[score_type] = std::vector<EnergyOP>(1,EnergyOP(new SymetryInterTwobodyEnergy(_nres,ichain,jchain,score_type,weight)));

        }else{
            energy_map[score_type].push_back(EnergyOP(new SymetryInterTwobodyEnergy(_nres,ichain,jchain,score_type,weight)));
        }
    }

    void marker_changed_segment_moved(Size start_res, Size end_res, Size root_res,ENERGY_TYPE energy_type = ENERGY_TYPE_UNDEFINED){
        for(auto pair :energy_map){
            for (auto energyOP : pair.second){
                // when mode is intra actually inter also need update, so only check if is inter
                if(energy_type!=ENERGY_TYPE_UNDEFINED && energy_type!=energyOP->energy_type()) continue;
                energyOP->marker_changed_segment_moved(start_res,end_res,root_res);
            }
        }
    }

    void marker_changed_segment_moved_by_chain(Size start_res, Size end_res, Size root_res,Size start_chain,Size end_chain){
        for(auto pair :energy_map){
            for (auto energyOP : pair.second){
                if(energyOP->energy_type()!=TWO_BODY_INTER_CHAIN)continue;
                TwoBodyEnergyOP twobody_energyOP = std::static_pointer_cast<TwoBodyEnergy>(energyOP);
                if(twobody_energyOP == nullptr){
                    std::cout<<"energy op from base to twobody fail!! check energy initialize is right";
                    exit(0);
                }
                if(twobody_energyOP->jchain()<start_chain || twobody_energyOP->jchain()>end_chain) continue;
                twobody_energyOP->marker_changed_segment_moved(start_res,end_res,root_res);
            }
        }
    }
    void snapshot(ENERGY_TYPE energy_type=ENERGY_TYPE_UNDEFINED){
        for(auto pair :energy_map){
            for (auto energyOP : pair.second){
                if(energy_type!=ENERGY_TYPE_UNDEFINED && energy_type!=energyOP->energy_type()) continue;
                energyOP->snapshot();
            }
        }
    }
    void rollback(ENERGY_TYPE energy_type=ENERGY_TYPE_UNDEFINED){
        for(auto pair :energy_map){
            for (auto energyOP : pair.second){
                if(energy_type!=ENERGY_TYPE_UNDEFINED && energy_type!=energyOP->energy_type()) continue;
                energyOP->rollback();
            }
        }
    }
    std::vector<EnergyOP> get_energy_tables(SCORE_TYPE score_type){return energy_map[score_type];}
    private:
    Size _nres;
    std::string _symmetry;
    Size _num_repeats;
    std::unordered_map <SCORE_TYPE,std::vector<EnergyOP>> energy_map;
};


}

#endif

//
//  local_structure_index.cpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 7/19/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#include "local_structure_index.hpp"
#define MAXNEIGHBORS 25

//use threads

namespace trajectoryAnalysis {
    LocalStructureIndex::LocalStructureIndex(Trajectory& traj, double rcutoff):OrderParameter(traj)
    {
        _rcutoff = rcutoff;
        _calculated = false;
        _condition.first=_rcutoff*_rcutoff;
        _condition.second=std::numeric_limits<int>::max();
        _lsi.resize(_snap->_center_of_mass_list.size());
        _LSIHist.setToNormalize(true);
    }
    
    LocalStructureIndex::~LocalStructureIndex(){
        
    }
    
    void LocalStructureIndex::setBinSize(double binsize){
        _LSIHist.setBinSize(binsize);
        _LSIHist.setToNormalize(true);
    }
    
    void LocalStructureIndex::print(const char* filename){
        _LSIHist.insert(0.,0.);
        if (!_calculated) compute();
        std::ofstream file(filename);
        file << _LSIHist;
        std::cout << "LSI successfully computed!\n";
        
    }
    
    void LocalStructureIndex::compute(){
        trajectory_t* traj = &(_trajectory->_trajectory);
        
        for (unsigned int f=0; f<traj->size(); f++) {
            _snap = &(*traj)[f];
            _computeNearestNeighbors(); //get information about nearest neighbors
            coord_list_t* com = &(_snap->_center_of_mass_list);
            for (unsigned int i=0; i<com->size(); i++) {
                _computeLSI(i);
            }
            _updateLSIframe();
        }
        _calculated = true;
    }
    
    void LocalStructureIndex::_computeLSI(unsigned int i){
        //double lsi=0;
        if (_neighbor_count[i] == 0) return;
        unsigned long m = std::count_if(_nearest_neighbors[i].begin(), _nearest_neighbors[i].begin()+_neighbor_count[i],std::bind2nd(std::less_equal<double_unsigned_pair_t>(),_condition));
        _lsi[i].resize(m);
        for (unsigned int j=0; j < m; j++) {
            _lsi[i][j] = sqrt(_nearest_neighbors[i][j+1].first) - sqrt(_nearest_neighbors[i][j].first);
            //can be made faster by computing sqrt outside
        }
        
    }
    
    void LocalStructureIndex::_updateLSIframe(){
        for (unsigned int j=0; j<_lsi.size(); j++) {
            double avg = 0;
            
            double ni = (double)_lsi[j].size();
            if (ni==0) continue;
            
            for (auto& i : _lsi[j]) avg += i;
            avg /= ni; _I = 0.;
            for (auto& i : _lsi[j]) _I += (i-avg)*(i-avg);
            
            _I /= ni;
            _LSIHist.insert(_I);
        }
    }
}
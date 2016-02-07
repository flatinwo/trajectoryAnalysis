//
//  order_parameter.cpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 10/13/15.
//  Copyright © 2015 Folarin Latinwo. All rights reserved.
//

#include "order_parameter.hpp"
#include "io.hpp"
#include "spatial.hpp"
#include <cassert>
#include <algorithm>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#include <boost/bind.hpp>
#include <boost/ref.hpp>

//try order-N algorithm in Allen & Tildesey or Frenkel & Smit book for autocorrelation
//it may be worth it to data in consistent vector format, that is data[0] represents a given set
//data[1] represents a given set, and so on....
//an alternative would be to load using multidimensional library in boost
//use R as library to CPP package

using namespace boost::accumulators;

namespace trajectoryAnalysis {
    
#define MAX_NUMBER_OF_NEIGHBORS 36
    
    OrderParameter::OrderParameter(const char* filename){
        load(filename, _data);
        std::cout << "Size of data is\t" << _data.size() << "\t" << _data[0].size() << std::endl;
        _restructureData();
    }
    
    OrderParameter::OrderParameter(Trajectory& traj):_trajectory(&traj),_mode(GLOBAL){
        _rcutoff = 0.74;
        _useMaxNumberOfNeighbors = false;
        _max_number_of_neighbors = 0;
        _snap = &_trajectory->_trajectory[0];
        _nmolecules = (unsigned int)_snap->_center_of_mass_list.size();
        _nearest_neighbors.resize(_snap->_center_of_mass_list.size(),
                                  double_unsigned_pair1d_t (MAX_NUMBER_OF_NEIGHBORS, std::pair<double, unsigned int>(0.,0)));
        _number_of_neighbors.resize(_snap->_center_of_mass_list.size(),0.);
    }
    
    OrderParameter::~OrderParameter(){
    }
    
    
#pragma mark SETS
    void OrderParameter::setRcutOff(double rcutoff){
        assert(rcutoff > 0.);
        _rcutoff = rcutoff;
    }
    
    void OrderParameter::setMaxNumberOfNearestNeighbors(unsigned int n_nghbrs){
        assert(n_nghbrs > 0);
        _max_number_of_neighbors = n_nghbrs;
        _useMaxNumberOfNeighbors = true;
        
    }
    
#pragma mark COMPUTES
    
    void OrderParameter::_computeNearestNeighbors(){
        coord_list_t* com = &(_snap->_center_of_mass_list);
        double rcutsqd = _rcutoff*_rcutoff;
        
        for (unsigned int i=0; i<com->size(); i++) {
            unsigned int k=0;
            for (unsigned int j=0; j<com->size(); j++) {
                if (i==j) continue;
                double rsq = distancesq((*com)[i], (*com)[j], _snap->box);
                if (rsq < rcutsqd) {
                    _nearest_neighbors[i][k].first = rsq;
                    _nearest_neighbors[i][k].second = j;
                    k++;
                    assert(k < MAX_NUMBER_OF_NEIGHBORS);
                }
            }
            std::sort(_nearest_neighbors[i].begin(), _nearest_neighbors[i].begin()+k); //implement sort up to
            assert(k >= _max_number_of_neighbors);
            assert(k > 0);  //Number of nearest numbers too many
        }
    }
    
    void OrderParameter::_refreshNeighbors(){
        //number of nearest neighbors
        for (auto& n :_number_of_neighbors) n = 0;
    }
    
    void OrderParameter::compute(){
        int virtual_function_overriden=0;
        assert(virtual_function_overriden);
    }
    
    
    
    
    void OrderParameter::_restructureData(){
        coord_list_t new_data;
        
        for (unsigned int i=0; i<_data[0].size(); i++) {
            coord_t x;
            for (unsigned int j=0; j<_data.size(); j++) {
                x.push_back(_data[j][i]);
            }
            new_data.push_back(x);
        }
        
        _data = new_data;
        std::cout << _data.size() << "\t" << _data[0].size() << "\n";
    }
}

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

//use try catch block for cases where the number of neighbors is large

using namespace boost::accumulators;

namespace trajectoryAnalysis {
    
#define MAX_NUMBER_OF_NEIGHBORS 100
    
    /*OrderParameter::OrderParameter(const char* filename){
        //load(filename, _data);
        std::cout << "I am not going to do anything\n";
        //std::cout << "Size of data is\t" << _data.size() << "\t" << _data[0].size() << std::endl;
        //_restructureData();
    }*/
    
    OrderParameter::OrderParameter(Trajectory& traj):_trajectory(&traj),_mode(GLOBAL){
        _initialize();
    }
    
    OrderParameter::~OrderParameter(){
    }
    
    
#pragma mark SETS
    void OrderParameter::setRcutOff(double rcutoff){
        assert(rcutoff > 0.);
        double minperiod =
        *(std::min_element(_snap->box.box_period.begin(), _snap->box.box_period.end()));
        if (rcutoff > 0.5*minperiod){
            std::cerr << "Your choice for rcutoff is larger than half of minimum period!\n";
            std::cout << "I will fix this by choosing maximum allowable cutoff\n";
            _rcutoff = 0.5*minperiod;
        }
        else _rcutoff = rcutoff;
    }
    
    void OrderParameter::setMaxNumberOfNearestNeighbors(unsigned int n_nghbrs){
        assert(n_nghbrs > 0);
        _max_number_of_neighbors = n_nghbrs;
        _useMaxNumberOfNeighbors = true;
        
    }
    
    void OrderParameter::setRmin(double min){
        assert(min > 0);
        assert(min < _rcutoff);
        _rminsq = min*min;
    }
    
#pragma mark COMPUTES
    
    //calculate neighbor distribution
    void OrderParameter::_computeNearestNeighbors(){
        coord_list_t* com = &(_snap->_center_of_mass_list);
        double rcutsqd = _rcutoff*_rcutoff;
        
        for (unsigned int i=0; i<com->size(); i++) {
            unsigned int k=0;
            for (unsigned int j=0; j<com->size(); j++) {
                if (i==j) continue;
                double rsq = distancesq((*com)[i], (*com)[j], _snap->box);
                if (rsq < rcutsqd && rsq > _rminsq) {
                    _nearest_neighbors[i][k].first = rsq;
                    _nearest_neighbors[i][k].second = j;
                    k++;
                    assert(k < MAX_NUMBER_OF_NEIGHBORS);
                }
            }
            _neighbor_count[i]=k;
            if (_neighborHist!=nullptr) _neighborHist->insert(k);
            std::sort(_nearest_neighbors[i].begin(), _nearest_neighbors[i].begin()+k); //implement sort up to
            assert(k >= _max_number_of_neighbors);
        }
    }
    
    void OrderParameter::printNeighborDistribution(const char* filename){
        _neighborHist = new stats_utils::HistogramDynamic<unsigned int>(1,true);
        for (auto& i : _trajectory->_trajectory){
            _snap = &i;
            _computeNearestNeighbors();
        }
        std::cout << "Neighbor distribution for r between " << sqrt(_rminsq)
        << " and " << _rcutoff << std::endl;
        std::ofstream file(filename);
        file << *_neighborHist;
        file.close();
        delete _neighborHist;
    }
    
    void OrderParameter::_refreshNeighbors(){
        //number of nearest neighbors
        for (auto& n :_number_of_neighbors) n = 0;
        _neighbor_count = _number_of_neighbors;
    }
    
    /*void OrderParameter::_refresh(unsigned int i){
        assert(false);
    }*/
    
    void OrderParameter::_refresh(){
        assert(false);
    }
    
    void OrderParameter::compute(){
        int virtual_function_overriden=0;
        assert(virtual_function_overriden);
    }
    
    void OrderParameter::_initialize(){
        _rcutoff = 0.74;
        _useMaxNumberOfNeighbors = false;
        _max_number_of_neighbors = 0;
        _snap = &_trajectory->_trajectory[0];
        _nmolecules = (unsigned int)_snap->_center_of_mass_list.size();
        _nearest_neighbors.resize(_snap->_center_of_mass_list.size(),
                                  double_unsigned_pair1d_t (MAX_NUMBER_OF_NEIGHBORS, std::pair<double, unsigned int>(0.,0)));
        _number_of_neighbors.resize(_snap->_center_of_mass_list.size(),0.);
        _neighbor_count = _number_of_neighbors;
        _localflag = false;
        _rminsq = 0;
        _neighborHist = nullptr;
    }
    
    void OrderParameter::setCalcType(Calc_t mode){
        _mode = mode;
        if (_mode==LOCAL || _mode==GLOBALANDLOCAL) _localflag = true;
        _update();
    }
    
    void OrderParameter::_update(){
        _refresh();
        assert(false);
    }
    
    /*
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
    */
}

//
//  tetrahedral_order_parameter.cpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 2/2/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#include "tetrahedral_order_parameter.hpp"
#include "spatial.hpp"
#include <sstream>
#include <iostream>
#include <iomanip>


//consider migrating to Boost Histograms

namespace trajectoryAnalysis {
    
#define THREEOVEREIGHT 0.3750000000
#define ONEOVERTHREE   0.3333333333

#define HISTMAX 1.0
#define HISTMIN -3.0
    
    TetrahedralOrderParameter::TetrahedralOrderParameter(Trajectory& traj):OrderParameter(traj){
        _max_number_of_neighbors = 4;
        _useMaxNumberOfNeighbors = true;
        _requireBinQvalues = true;
        if (_requireBinQvalues){
            _QHist.setBinSize(0.01);
            _QHist.setToNormalize(true);
            
            #ifndef USEHISTMAX
            _QHist.insert(HISTMAX,0);
            _QHist.insert(HISTMIN,0);
            #endif
        }
    }
    
    TetrahedralOrderParameter::~TetrahedralOrderParameter(){
        
    }
    
    
#pragma mark SETS
    void TetrahedralOrderParameter::setMaxNumberOfNearestNeighbors(unsigned int n_nghbrs){
        assert(n_nghbrs == 4);
        _max_number_of_neighbors = 4;
        _useMaxNumberOfNeighbors = true;
    }
    
#pragma mark COMPUTES
    void TetrahedralOrderParameter::compute(){
        _computeWithMaxNeighbors();
    }
    
    
    void TetrahedralOrderParameter::_computeWithMaxNeighbors(){
        trajectory_t* traj = &(_trajectory->_trajectory);
        
        for (unsigned int f=0; f<traj->size(); f++) {
            _snap = &(*traj)[f];
            _computeNearestNeighbors(); //get information about nearest neighbors
            coord_list_t* com = &(_snap->_center_of_mass_list);
            _Qs.resize(com->size());
            for (unsigned int i=0; i<com->size(); i++) _computeQ(i);
            _updateQframe();
        }
    }
    
    
    void TetrahedralOrderParameter::_computeQ(unsigned int i){
        double sum=0.;
        coord_list_t* com = &(_snap->_center_of_mass_list);
        for (unsigned int j=0; j < _max_number_of_neighbors - 1; j++){
            for (unsigned int k=j+1; k < _max_number_of_neighbors; k++) {
                unsigned int indj = _nearest_neighbors[i][j].second;
                unsigned int indk = _nearest_neighbors[i][k].second;
                double x = cosine_angle((*com)[i], (*com)[indj], (*com)[indk], _snap->box);
                sum += (x+ONEOVERTHREE)*(x+ONEOVERTHREE);
            }
        }
        _Qs[i] = (1. - THREEOVEREIGHT*sum);
        if (_requireBinQvalues) _QHist.insert(_Qs[i]);
    }
    
    void TetrahedralOrderParameter::_updateQframe(){
        _Q = 0.;
        for (double& n: _Qs) _Q += n;
        _Q /= (double) _snap->_center_of_mass_list.size();
        _Qframe.push_back(_Q);
    }
    
    void TetrahedralOrderParameter::_refresh(){
        _refreshNeighbors();
    }
    
    
#pragma mark IOs
    void TetrahedralOrderParameter::_openFiles(){
        if (_requireBinQvalues) _ofiles.resize(2);
        else _ofiles.resize(1);
        
        std::ostringstream os;
        os << std::setprecision(4);
        os << "q_tetrahedral.txt";
        _ofiles[0].reset(new std::ofstream(os.str().c_str()));
        
        if (_requireBinQvalues) {
            os.str("");
            os.clear();
            os << "q_Histogram.txt";
            _ofiles[1].reset(new std::ofstream(os.str().c_str()));
        }
        
    }
    
    void TetrahedralOrderParameter::_closeFiles(){
        for (auto& i : _ofiles) i->close();
    }
    
    void TetrahedralOrderParameter::print(){
        assert(_Qframe.size()>0);
        double timestep = _trajectory->_time_step;
        
        _openFiles();
        for (unsigned int i=0; i<_Qframe.size(); i++) {
            double time = ((double) i)*timestep;
            *_ofiles[0] << time << "\t" << _Qframe[i] << std::endl;
        }
        
        if (_requireBinQvalues) *_ofiles[1] << _QHist << std::endl;
        _closeFiles();
    }
}

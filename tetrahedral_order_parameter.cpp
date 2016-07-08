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
#include <functional>
#include <algorithm>


//consider migrating to Boost Histograms

namespace trajectoryAnalysis {
    
#define THREEOVEREIGHT 0.3750000000
#define ONEOVERTHREE   0.3333333333
#define NINEOVERFOUR   2.2500000000

#define HISTMAX 1.0
#define HISTMIN -3.0
    
    TetrahedralOrderParameter::TetrahedralOrderParameter(Trajectory& traj):OrderParameter(traj){
        _max_number_of_neighbors = 4;
        _useMaxNumberOfNeighbors = true;
        _requireBinQvalues = true;
        _requirePositionValues = false;
        _tHQs = nullptr;
        _condition.first=0.;
        _condition.second=std::numeric_limits<int>::max();
        
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
        if (_tHQs!= nullptr) delete _tHQs;
    }
    
    
#pragma mark SETS
    void TetrahedralOrderParameter::setMaxNumberOfNearestNeighbors(unsigned int n_nghbrs){
        assert(n_nghbrs == 4);
        _max_number_of_neighbors = 4;
        _useMaxNumberOfNeighbors = true;
    }
    
    void TetrahedralOrderParameter::setPositionTetrahedrality(bool flag){
        _requirePositionValues = flag;
        if (_requirePositionValues) _tHQs = new std::vector<tHQs>();
        else{
            if (_tHQs != nullptr) delete _tHQs;
        }
    }
    
    void TetrahedralOrderParameter::addRmax(double rmax){
        assert(_tHQs != nullptr);
        assert(rmax > 0.);
        
        (*_tHQs).push_back(tHQs(rmax*rmax));
        _counts.resize(_tHQs->size());
        
        if (_requireBinQvalues) {
            tHQs* HQs = &((*_tHQs)[_tHQs->size()-1]);            
#ifndef USEHISTMAX
            HQs->_QHist.insert(HISTMAX,0);
            HQs->_QHist.insert(HISTMIN,0);
#endif
            
            
            
        }
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
            if (_tHQs==nullptr) _Qs.resize(com->size());
            else _resize();
            for (unsigned int i=0; i<com->size(); i++) {
                if (_tHQs==nullptr)_computeQ(i);
                else _computeQR(i);
            }
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
        if (_tHQs == nullptr){
            _Q = 0.;
            for (double& n: _Qs) _Q += n;
            _Q /= (double) _snap->_center_of_mass_list.size();
            _Qframe.push_back(_Q);
        }
        else{
            for (unsigned int i=0; i<_tHQs->size(); i++) {
                tHQs* hq = &(*_tHQs)[i];
                _Q=0;
                for (double& n: hq->_Qs) _Q+=n;
                _Q /= (double) _snap->_center_of_mass_list.size();
                hq->_Qframe.push_back(_Q);
            }
        }
    }
    
    void TetrahedralOrderParameter::_computeQR(unsigned int i){
        double sum=0.;
        double neighbors;
        coord_list_t* com = &(_snap->_center_of_mass_list);
        
        for (unsigned int j=0; j<_tHQs->size(); j++) {
            _condition.first = (*_tHQs)[j]._rmaxsqd ;
            _counts[j] = std::count_if(_nearest_neighbors[i].begin(), _nearest_neighbors[i].begin()+4,
                                      std::bind2nd(std::less_equal<double_unsigned_pair_t>(),_condition));
            assert(_counts[j]<=4);
            //assert(_counts[j]>=2);
            
            //update stats for each _tHQ
            (*_tHQs)[j]._numberstats.first++;
            (*_tHQs)[j]._numberstats.second+=_counts[j];
            
            if (_counts[j]<2) continue;
            
            neighbors=sum=0.;
            for (unsigned int jj=0; jj < _counts[j] - 1; jj++){
                for (unsigned int k=jj+1; k < _counts[j]; k++) {
                    neighbors++;
                    unsigned int indj = _nearest_neighbors[i][jj].second;
                    unsigned int indk = _nearest_neighbors[i][k].second;
                    double x = cosine_angle((*com)[i], (*com)[indj], (*com)[indk], _snap->box);
                    sum += (x+ONEOVERTHREE)*(x+ONEOVERTHREE);
                }
            }
             double q = (1. - NINEOVERFOUR*sum/neighbors);
            (*_tHQs)[j]._Qs[i] = q;
            if (_requireBinQvalues) (*_tHQs)[j]._QHist.insert(q);
            
        }
    }
    
    void TetrahedralOrderParameter::_refresh(){
        _refreshNeighbors();
    }
    
    void TetrahedralOrderParameter::_resize(){
        for (auto& i : *_tHQs) i._Qs.resize(_snap->_center_of_mass_list.size());
    }
    
#pragma mark IOs
    void TetrahedralOrderParameter::_openFiles(){
        if (_tHQs==nullptr) {
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
        else{
            
            if (_requireBinQvalues) _ofiles.resize(2*_tHQs->size());
            else _ofiles.resize(_tHQs->size());
            
            std::ostringstream os;
            os << std::setprecision(4);
            int fc=0;
            unsigned long fci = _tHQs->size();
            for (unsigned int j=0; j<_tHQs->size(); j++) {
                os.str("");
                os.clear();
                tHQs* hq = &(*_tHQs)[j];
                double rcut = sqrt(hq->_rmaxsqd);
                os << "q_tetrahedral_"<< std::to_string(rcut)<<".txt";
                _ofiles[fc++].reset(new std::ofstream(os.str().c_str()));
                
                if (_requireBinQvalues) {
                    os.str("");
                    os.clear();
                    os << "q_Histogram_"<<std::to_string(rcut)<<".txt";
                    _ofiles[fci++].reset(new std::ofstream(os.str().c_str()));
                }
            }
        }
    }
    
    void TetrahedralOrderParameter::_closeFiles(){
        for (auto& i : _ofiles) i->close();
    }
    
    void TetrahedralOrderParameter::print(){
        double timestep = _trajectory->_time_step;
        
        _openFiles();
        if (_tHQs==nullptr){
            assert(_Qframe.size()>0);
            for (unsigned int i=0; i<_Qframe.size(); i++) {
                double time = ((double) i)*timestep;
                *_ofiles[0] << time << "\t" << _Qframe[i] << std::endl;
            }
            
            if (_requireBinQvalues) *_ofiles[1] << _QHist << std::endl;
        }
        else{
            for (unsigned int i=0; i<(*_tHQs)[0]._Qframe.size();i++) {
                double time = ((double) i)*timestep;
                unsigned int fc=0;
                for (unsigned int j=0; j<_tHQs->size(); j++) {
                    tHQs* hq = &(*_tHQs)[j];
                    assert(hq->_Qframe.size()==(*_tHQs)[0]._Qframe.size());
                    *_ofiles[fc++] << time << "\t" << hq->_Qframe[i] << std::endl;
                }
            }
            
            if (_requireBinQvalues) {
                unsigned long fc = _tHQs->size();
                for (unsigned int j=0; j<_tHQs->size(); j++) {
                    tHQs* hq = &(*_tHQs)[j];
                    *_ofiles[fc++] << hq->_QHist << std::endl;
                }
            }
            
            std::ofstream os("number_info.txt");
            for (unsigned int j=0; j<_tHQs->size();j++) {
                tHQs& i = (*_tHQs)[j];
                os << sqrt(i._rmaxsqd) << "\t" << i._numberstats.second/((double)i._numberstats.first) << std::endl;
            }
            os.close();
        }
        
        _closeFiles();
    }
}

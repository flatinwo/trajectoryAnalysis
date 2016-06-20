//
//  averaged_bond_order_parameter.cpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 6/13/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#include "averaged_bond_order_parameter.hpp"
#include "trajectory.hpp"
#include "struct_def.h"

//try to thread this

namespace trajectoryAnalysis {
    AveragedBondOrderParameter::AveragedBondOrderParameter(Trajectory& traj, int l):
    BondOrderParameter(traj,l){
        
    }
    
    AveragedBondOrderParameter::~AveragedBondOrderParameter(){
        
    }
    
    
    void AveragedBondOrderParameter::_computeWithRcutOff(){
        trajectory_t* traj = &(_trajectory->_trajectory);
        
        for (unsigned int f=0; f<traj->size(); f++) {
            _snap = &(*traj)[f];
            _max_number_of_neighbors=1;
            _computeNearestNeighbors(); //get information about nearest neighbors
            coord_list_t* com = &(_snap->_center_of_mass_list);
            for (unsigned int i=0; i<com->size(); i++) {
                for (unsigned int j=0; j<com->size(); j++) {
                    if (i==j) continue;
                    _computeHarmonics(i, j);
                    
                }
                for (unsigned int m=0; m <_qlm_i[i].size(); m++) {
                    assert(_number_of_neighbors[i] > 0);
                    _qlm_i[i][m] /= (double) _number_of_neighbors[i];
                    _Qlm[m] += _qlm_i[i][m];
                }
            }
            for (unsigned int i=0; i<com->size(); i++) {
                _max_number_of_neighbors = _number_of_neighbors[i];
                _computeql_i(i);
                _computeWl_i(i);
            }

            _computeQl();
            if (_requireThirdOrderInvaraints) _computeWl();
        }
        if (_localflag) _tallyLocals();
    }

    
    void AveragedBondOrderParameter::_computeWithMaxNeighbors(){
        trajectory_t* traj = &(_trajectory->_trajectory);
        
        for (unsigned int f=0; f<traj->size(); f++) {
            _snap = &(*traj)[f];
            _computeNearestNeighbors(); //get information about nearest neighbors
            coord_list_t* com = &(_snap->_center_of_mass_list);
            
            
            for (unsigned int i=0; i<com->size(); i++) {
                for (unsigned int j=0; j < _max_number_of_neighbors; j++) {
                    unsigned int k = _nearest_neighbors[i][j].second;
                    _computeHarmonics(i, k);
                }
                for (unsigned int m=0; m <_qlm_i[i].size(); m++) {
                    assert(_number_of_neighbors[i] > 0);
                    assert(_number_of_neighbors[i] == _max_number_of_neighbors);
                    _qlm_i[i][m] /= (double) _number_of_neighbors[i];
                    _Qlm[m] += _qlm_i[i][m];
                }
            }
            
            for (unsigned int i=0; i<com->size(); i++) {
                _computeql_i(i);
                _computeWl_i(i);
            }
            _computeQl();
           
            if (_requireThirdOrderInvaraints) _computeWl();
        }
        if (_localflag) _tallyLocals();
    }
    
    void AveragedBondOrderParameter::_computeql_i(unsigned int i){
        _average(i);
    }
    
    void AveragedBondOrderParameter::_average(unsigned int i){
        //component_list_t temp(_qlm_i[i].size(),component_t(0,0));//=  _qlm_i[i];
        component_list_t temp =  _qlm_i[i];
        double count=0;
        
        for (unsigned int m=0; m<_max_number_of_neighbors;m++){
            count++;
            unsigned int k = _nearest_neighbors[i][m].second; //neighbor index
            for (unsigned int j=0; j<_qlm_i[i].size();j++) temp[j] += _qlm_i[k][j];
        }
        
        double qli = 0;
        for (unsigned int j=0; j<temp.size(); j++) {
            temp[j] = temp[j]/count;
            qli += std::norm(temp[j]);
        }
        qli *= (4.*M_PI/(double) (2*_l + 1));
        qli=sqrt(qli);
        _qli[i] += qli;
        if (_requireBinValues) _qHist.insert(qli);
    }
}
//
//  averaged_tetrahedral_order_parameter.cpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 6/14/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#include "averaged_tetrahedral_order_parameter.hpp"
#include "spatial.hpp"

namespace trajectoryAnalysis {
    
#define THREEOVEREIGHT 0.3750000000
#define ONEOVERTHREE   0.3333333333
#define COUNT 1
    
    AveragedTetrahedralOrderParameter::AveragedTetrahedralOrderParameter(Trajectory& traj):
    TetrahedralOrderParameter(traj){
         _localqs.resize(_snap->_center_of_mass_list.size(),coord_t(6,0));
    }
    
    void AveragedTetrahedralOrderParameter::_computeQ(unsigned int i){
        coord_list_t* com = &(_snap->_center_of_mass_list);
        for (unsigned int j=0; j < _max_number_of_neighbors - 1; j++){
            for (unsigned int k=j+1; k < _max_number_of_neighbors; k++) {
                unsigned int indj = _nearest_neighbors[i][j].second;
                unsigned int indk = _nearest_neighbors[i][k].second;
                double x = cosine_angle((*com)[i], (*com)[indj], (*com)[indk], _snap->box);
                _localqs[i][k] = x + ONEOVERTHREE;
            }
        }
    }
    
    void AveragedTetrahedralOrderParameter::_average(unsigned int i){
        _sum = 0;
        _dummyclt = _localqs[i];
        for (unsigned int j=0; j < _max_number_of_neighbors; j++){
            unsigned int indj = _nearest_neighbors[i][j].second;
            for (unsigned int k=0; k<_dummyclt.size(); k++)_dummyclt[k] += _localqs[indj][k];

        }
        double maxnumsq = (double) (_max_number_of_neighbors+COUNT);
        maxnumsq *= maxnumsq;
        for (auto& j : _localqs[i] ) _sum += (j*j)/maxnumsq;
    }
    

    void AveragedTetrahedralOrderParameter::_updateQframe(){
        coord_list_t* com = &(_snap->_center_of_mass_list);
        for (unsigned int i=0; i<com->size(); i++) {
            _average(i);
            _Qs[i] = (1. - THREEOVEREIGHT*_sum);
            if (_requireBinQvalues) _QHist.insert(_Qs[i]);
        }
        _Q = 0.;
        for (double& n: _Qs) _Q += n;
        _Q /= (double) _snap->_center_of_mass_list.size();
        _Qframe.push_back(_Q);
    }
}
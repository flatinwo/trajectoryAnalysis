//
//  averaged_bond_order_parameter.cpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 6/13/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#include "averaged_bond_order_parameter.hpp"
#include "struct_def.h"

namespace trajectoryAnalysis {
    AveragedBondOrderParameter::AveragedBondOrderParameter(Trajectory& traj, int l):
    BondOrderParameter(traj,l){
        
    }
    
    AveragedBondOrderParameter::~AveragedBondOrderParameter(){
        
    }
    
    void AveragedBondOrderParameter::_computeql_i(unsigned int i){
        _average(i);
        BondOrderParameter::_computeql_i(i);
    }
    
    void AveragedBondOrderParameter::_average(unsigned int i){
        component_list_t temp = _qlm_i[i];
        double count=0;
        for (unsigned int m=0; m<_number_of_neighbors[i];m++){
            count++;
            unsigned int k = _nearest_neighbors[i][m].second; //neighbor index
            for (unsigned int j=0; j<_qlm_i[i].size();j++) temp[j] += _qlm_i[k][j];
        }
        
        for (unsigned int j=0; j<temp.size(); j++) {
            _qlm_i[i][j] = temp[j]/count;
        }
    }
}
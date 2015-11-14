//
//  order_parameter.cpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 10/13/15.
//  Copyright © 2015 Folarin Latinwo. All rights reserved.
//

#include "order_parameter.hpp"
#include "io.hpp"
#include <cassert>

//try order-N algorithm in Allen & Tildesey or Frenkel & Smit book for autocorrelation

namespace trajectoryAnalysis {
    
    OrderParameter::OrderParameter(const char* filename){
        load(filename, _data);
    }
    
    OrderParameter::OrderParameter(Trajectory& traj):_trajectory(&traj){
        
    }
    
    OrderParameter::~OrderParameter(){
    }
    
#pragma mark COMPUTES
    
    //brute force auto-correlation function
    void OrderParameter::computeAutoCorrelation(){
        assert(_data.size()>100);
        assert(_data[0].size()>2);
        
        int dataSize = (int) _data.size();
        
        corr_point_t zeros;
        zeros.first = zeros.second = 0;
        corr_point_list_t correlation((int) floor((double) dataSize/3), zeros);
        
        for (unsigned int i=0; i<_data.size(); i++) {
            for (unsigned int j=0; j<correlation.size(); j++) {
                if (i+j >= _data.size())
                    continue;
                correlation[j].second += _data[i][2]*_data[i+j][2];
                correlation[j].first++;
            }
        }
        
        _correlation.clear();
        double normalization = (double) correlation[0].first / correlation[0].second;
        
        for (auto it=correlation.begin(); it != correlation.end(); ++it) {
            _correlation.push_back((it->second / (double) it->first)*normalization);
        }
        
        
        assert(_correlation.size() == correlation.size());
    }
}

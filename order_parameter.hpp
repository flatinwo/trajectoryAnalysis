//
//  order_parameter.hpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 10/13/15.
//  Copyright © 2015 Folarin Latinwo. All rights reserved.
//

#ifndef order_parameter_hpp
#define order_parameter_hpp

#include <stdio.h>
#include "struct_def.h"
#include "trajectory.hpp"

//set all calculation flags then output all results
//have default flags

namespace trajectoryAnalysis {
    
    /**
     \brief A class to perform specific calculations on/for order parameters
     */
    
    class OrderParameter{
    public:
        OrderParameter(Trajectory&);
        OrderParameter(const char*);
        ~OrderParameter();
        
        void computeAutoCorrelation(int index=2);
        void computeCrossCorrelation();
        void computeAverageAndVariance(int j=1);
        
        //virtual void compute();
    
        void printCorrelation();
        
    protected:
        Trajectory* _trajectory;
        coord_list_t _data;
        coord_t _correlation;
        double _average;
        double _variance;
        
        void _restructureData();
        
    };
}

#endif /* order_parameter_hpp */

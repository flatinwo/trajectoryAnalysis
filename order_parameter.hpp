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

namespace trajectoryAnalysis {
    
    class OrderParameter{
    public:
        OrderParameter(Trajectory&);
        OrderParameter(const char*);
        ~OrderParameter();
        
        void computeAutoCorrelation();
        void computeCrossCorrelation();
        
        void computeBondOrderParameter();
        
    protected:
        Trajectory* _trajectory;
        coord_list_t _data;
        coord_t _correlation;
        
    };
}

#endif /* order_parameter_hpp */

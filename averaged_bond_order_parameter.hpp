//
//  averaged_bond_order_parameter.hpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 6/13/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#ifndef averaged_bond_order_parameter_hpp
#define averaged_bond_order_parameter_hpp

#include <stdio.h>
#include "bond_order_parameter.hpp"

namespace trajectoryAnalysis {
    class AveragedBondOrderParameter : public BondOrderParameter{
    public:
        AveragedBondOrderParameter(Trajectory&, int l=6);
        ~AveragedBondOrderParameter();
    protected:
        void _computeql_i(unsigned int);
        void _average(unsigned int);

    };
}

#endif /* averaged_bond_order_parameter_hpp */

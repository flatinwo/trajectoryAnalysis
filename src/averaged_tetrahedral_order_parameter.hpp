//
//  averaged_tetrahedral_order_parameter.hpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 6/14/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#ifndef averaged_tetrahedral_order_parameter_hpp
#define averaged_tetrahedral_order_parameter_hpp

#include <stdio.h>
#include "tetrahedral_order_parameter.hpp"

namespace trajectoryAnalysis{
    /**
     \brief A derived class to perform specific calculations for Averaged Tetrahedral order parameters
     Value should range between -3 and +1.
     */
    class AveragedTetrahedralOrderParameter : public TetrahedralOrderParameter{
    public:
        AveragedTetrahedralOrderParameter(Trajectory&);
       // ~AveragedTetrahedralOrderParameter();
        
    private:
        coord_list_t _localqs;
        coord_t _dummyclt;
        
    protected:
        double _sum;
        void _average(unsigned int);
        void _computeQ(unsigned int);
        void _updateQframe();
    };
}

#endif /* averaged_tetrahedral_order_parameter_hpp */

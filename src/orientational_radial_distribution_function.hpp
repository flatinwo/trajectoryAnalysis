//
//  orientational_radial_distribution_function.hpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 6/21/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#ifndef orientational_radial_distribution_function_hpp
#define orientational_radial_distribution_function_hpp

#include <stdio.h>
#include "radial_distribution_function.hpp"

namespace trajectoryAnalysis {
    
    typedef RadialDistributionFunction RDFt;
    
    class OrientationalRDF : public RDFt{
    public:
        OrientationalRDF(xyztrajectory_t&, Box&,double binsize=0.01);
        
    protected:
        
        
    };
}

#endif /* orientational_radial_distribution_function_hpp */

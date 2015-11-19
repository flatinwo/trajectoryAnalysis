//
//  visualizer.hpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 11/19/15.
//  Copyright © 2015 Folarin Latinwo. All rights reserved.
//

#ifndef visualizer_hpp
#define visualizer_hpp

#include <stdio.h>
#include "io.hpp"

namespace trajectoryAnalysis {
    class Visualizer{
    public:
        Visualizer();
        ~Visualizer();
        
        virtual void visualize(trajectory_t&);
        virtual void visualize(Snap&);
        
    protected:
        
        
    };
}
#endif /* visualizer_hpp */

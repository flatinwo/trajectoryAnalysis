//
//  visualizer_xyz.hpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 11/19/15.
//  Copyright © 2015 Folarin Latinwo. All rights reserved.
//

#ifndef visualizer_xyz_hpp
#define visualizer_xyz_hpp

#include <stdio.h>
#include "visualizer.hpp"
#include <fstream>

namespace trajectoryAnalysis {
    class VisualizerXYZ : public Visualizer{
    public:
        VisualizerXYZ(const char*);
        ~VisualizerXYZ();
        
        void setTypeMax(typelog_t& typemax);
        void visualize(xyztrajectory_t& traj);
        
    protected:
        std::ofstream _os;
        typelog_t _typemax;
        
    };
}

#endif /* visualizer_xyz_hpp */

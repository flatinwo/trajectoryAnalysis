//
//  trajectory.hpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 10/12/15.
//  Copyright © 2015 Folarin Latinwo. All rights reserved.
//

#ifndef trajectory_hpp
#define trajectory_hpp

#include <stdio.h>
#include "struct_def.h"

//to do write static_cast, dynamic_cast or something similar for trajectoryAnalysis::xyztrajectory_t to trajectoryAnalysis::trajectory_t

namespace trajectoryAnalysis {
    
    class Trajectory{
        
        friend class OrderParameter;
        friend class BondOrderParameter;
        friend class TetrahedralOrderParameter;
        
    public:
        Trajectory();
        Trajectory(const char* filename, bool=false, unsigned int i=1, unsigned int=5);
        Trajectory(const char* filename, FILETYPE=GRO, unsigned int i=1, unsigned int=4);
        ~Trajectory();
        
        unsigned int maxCorrelationLength;
        
        int getNumberOfSnaps();
        double getTimeStep();
        trajectory_t& getTrajectory();
        void unfold();
        
        void computeMeanSquaredDisplacement();
        void computeFskt();
        void computeGofR();
        
        friend std::ostream& operator << (std::ostream&, const Trajectory&);
        
    protected:
        trajectory_t _trajectory;
        double _time_step;
        bool _unfolded;
        
        void computeTimeStep();
        void computeMaxCorrelationLength();
        coord_t _correlation;
    };
}

#endif /* trajectory_hpp */

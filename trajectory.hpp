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
#include "histogram_dynamic.h"

//to do write static_cast, dynamic_cast or something similar for trajectoryAnalysis::xyztrajectory_t to trajectoryAnalysis::trajectory_t

namespace trajectoryAnalysis {
    
    typedef std::vector< stats_utils::HistogramDynamic < double > > Hists1Dt;
    typedef std::vector<Hists1Dt> Hists2Dt;
    
    class Trajectory{
        
        friend class OrderParameter;
        friend class BondOrderParameter;
        friend class TetrahedralOrderParameter;
        friend class AveragedBondOrderParameter;
        
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
        
        void setComputeFskt(bool=true,double=4.0);
        void setComputeFrequency(double);
        void setVanHoveBinSize(double=0.01);
        void setUseVanHove(bool);
        
        void computeMeanSquaredDisplacement();
        void computeGofR();
        
        friend std::ostream& operator << (std::ostream&, const Trajectory&);
        
    protected:
        trajectory_t _trajectory;
        double _time_step,_k;
        bool _useVanHove;
        bool _unfolded;
        bool _computeFskt;
        
        
        void computeTimeStep();
        void computeMaxCorrelationLength();
        void _printFskt();
        
        coord_t _correlation;
        
        coord_list_t* _Fskts;
        coord_t* _ks;
        
        Hists2Dt* _vanHovefxn;      //self-part of Van Hove function
    };
}

#endif /* trajectory_hpp */

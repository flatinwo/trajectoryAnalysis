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
#include <utility>

//to do write static_cast, dynamic_cast or something similar for trajectoryAnalysis::xyztrajectory_t to trajectoryAnalysis::trajectory_t

namespace trajectoryAnalysis {
    
    typedef stats_utils::HistogramDynamic<double> Hist1Dt;
    typedef std::vector< Hist1Dt > Hists1Dt;
    typedef std::vector<Hists1Dt> Hists2Dt;
    typedef std::pair<bool,unsigned int> Skipt;
    
    class Trajectory{
        
        friend class OrderParameter;
        friend class BondOrderParameter;
        friend class TetrahedralOrderParameter;
        friend class AveragedBondOrderParameter;
        friend class LocalStructureIndex;
        
    public:
        Trajectory();
        Trajectory(const char* filename, bool=false, unsigned int i=1, unsigned int=5);
        Trajectory(const char* filename, FILETYPE=GRO, unsigned int i=1, unsigned int=4, unsigned int=50000);
        Trajectory(trajectory_t&);
        ~Trajectory();
        
        friend Trajectory operator+(const Trajectory& rhs1, const Trajectory& rhs2);
        
        
        int getNumberOfSnaps();
        int getNumberOfFrames();
        double getTimeStep();
        trajectory_t& getTrajectory();
        void unfold();
        
        void setTimeStep(double);
        void setComputeFskt(bool=true,double=4.0);
        void setComputeFrequency(double);
        void setVanHoveBinSize(double=0.01);
        void setUseVanHove(bool);
        void setMaxNumberOfFrames(int n);
        void setUnFolded(bool);
        void setSkipInfo(Skipt = std::make_pair(true,100),short=2);
        void setMaxCorrelationLength(unsigned int mcl);
        
        void computeMeanSquaredDisplacement();
        function1d_t computeGofR(double binsize=0.01);
        Hist1Dt computeNeighborDistribution(short ind, double binsize=0.05);
        
        void printFskt(double=1.);
        
        friend std::ostream& operator << (std::ostream&, const Trajectory&);
        
    protected:
        trajectory_t _trajectory;
        double _time_step,_k;
        unsigned int _maxframes;
        bool _useVanHove;
        bool _unfolded;
        bool _computeFskt;
        short _skipfactor;
        unsigned int maxCorrelationLength;
        
        
        void computeTimeStep();
        void computeMaxCorrelationLength();
        
        coord_t _correlation;
        
        coord_list_t* _Fskts;
        coord_t* _ks;
        
        FILETYPE _type;
        
        Hists2Dt* _vanHovefxn;      //self-part of Van Hove function
        Skipt* _skipinfo;
    };
}

#endif /* trajectory_hpp */

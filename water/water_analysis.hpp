//
//  water_analysis.hpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 5/17/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#ifndef water_analysis_hpp
#define water_analysis_hpp

#include "struct_def.h"
#include "trajectory.hpp"
#include <stdio.h>
#include <memory>

namespace  trajectoryAnalysis {
#define COS30 0.8660254
#define COS50 0.6427876
#define MOLSIZE 3
#define BREAKINGEVENTS 40000
    
    class WaterAnalysis{
        
        typedef std::vector<unsigned_list_t> unsigned_list2d_t;
        
        struct HydrogenBondBreakInfo{
            std::vector<unsigned int> frames;
            std::vector<std::pair<short,short>> oxygen_ids;   //first is initial O, second is final O
            std::vector < std::vector< std::pair<unsigned long,double> > > trajA; //typedef this
            std::vector < std::vector< std::pair<unsigned long,double> > > trajB;
        };
    public:
        WaterAnalysis(Trajectory&);
        
        void setRCutOff(double);
        void setAngleCutOff(double);
        void setSkipFrame(int);
        void setTimeStep(double);
        
        int  getNumberOfBondBreakingEvents();
        
        void printBonded();
        void printNumberHydrogenBonds();
        void printROO();
        void print();
        
        
        void compute();
        
    protected:
        std::vector<unsigned_list2d_t> _bondedList;
        std::vector<HydrogenBondBreakInfo> _bondBreakInfos;
        
        //file operations
        std::vector< std::unique_ptr<std::ofstream> > _ofiles;
        void _openFiles();
        void _closeFiles();
        
        int _skip_frame;
        double _time_step;
        bool _requireThetaValues;
        
        void _computeHydrogenBonds(int i=0);
        void _computeAllHydrogenBonds();
        void _analyzeHydrogenBonds();
        void _computeROO();
        void _normalize();
        
        Trajectory* _traj;
        trajectory_t* _trajectory;
        
        corr_point_list_t trajAavg;
        corr_point_list_t trajBavg;
        
        double _rcutoff;
        double _anglecutoff;
        
    };
}
#endif /* water_analysis_hpp */

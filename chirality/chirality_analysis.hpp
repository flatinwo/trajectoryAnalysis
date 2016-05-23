//
//  chirality_analysis.hpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 5/23/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#ifndef chirality_analysis_hpp
#define chirality_analysis_hpp

#include <stdio.h>
#include "trajectory.hpp"
#include "io.hpp"
#include "visualizer_xyz.hpp"

namespace trajectoryAnalysis {
    class ChiralityAnalysis{
    public:
        ChiralityAnalysis(const char* filename, Box& box);
        ChiralityAnalysis(int argc, const char* argv[]);
        
        void analyze();
        void visualize();
        
    protected:
        double _averageC;
        std::string _filein;
        std::string _fileout;
        std::vector<unsigned int> _typecount, _typecountmax;
        std::vector<std::string> _typematch;
        
        xyztrajectory_t _traj;
        typelog_t   _logtypes;
        Box         _box;
        
        void _initialize();
        void _analyzeChiralityXYZ(xyzfile& snap);
        std::string _chiralityunwrap(coord_list_t& x, double& zetad);
        
    };
}

#endif /* chirality_analysis_hpp */

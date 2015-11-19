//
//  visualizer_xyz.cpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 11/19/15.
//  Copyright © 2015 Folarin Latinwo. All rights reserved.
//

#include "visualizer_xyz.hpp"
#include <algorithm>
#include <cassert>

namespace trajectoryAnalysis {
    VisualizerXYZ::VisualizerXYZ(const char* filename):_os(filename){
        
    }
    
    VisualizerXYZ::~VisualizerXYZ(){
        
    }
    
    void VisualizerXYZ::setTypeMax(typelog_t& logtype){
        _typemax = logtype;
    }
    
    
    void VisualizerXYZ::visualize(xyztrajectory_t& xyztraj){
        assert(_typemax.size() > 0);
        unsigned int n = 0; //total number of types
        
        for (typelog_t::iterator it=_typemax.begin(); it!=_typemax.end(); ++it)
            n += it->second;
    
        std::cout << "Total number of types is\t" << n << "\n";
        
        
        
        for (unsigned int i=0; i<xyztraj.size(); i++) {
            xyz_info info;
            info.type = xyztraj[i].type;
            _os << n << "\n\n";
            
            for (typelog_t::iterator it=_typemax.begin(); it!=_typemax.end(); ++it) {
                for (unsigned int j=0; j< it->second; j++) {
                    if (info.type[j] == it->first) {
                        _os << info.type[j] << "\t" << xyztraj[i].x[j] << "\n";
                    }
                    else
                        _os << info.type[j] << "\t0\t0\t-1\n";
                }
            }

        }
        
        _os.close();
        std::cout << "Total Number of Frames written is:\t" << xyztraj.size() << std::endl;
        
        
    }
}


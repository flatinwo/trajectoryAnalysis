//
//  quick_test_water.cpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 5/27/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include <cassert>

#include "trajectoryAnalysis.h"

using namespace trajectoryAnalysis;

int qtest_water(int argc, const char * argv[]) {
    // insert code here...
    //std::cout << "Hello, World!\n";
    
    //const char* filename = "/Users/Folarin/Documents/vmd_views/water/patchy_colloids/test_snaps/dump22f_2b.xyz";
    const char* filename = "/Users/Folarin/Documents/dump22f_2b.xyz";
    
    Trajectory traj(filename,true,1,1);
    traj.unfold();
    WaterAnalysis analyze(traj);
    analyze.compute();
    //analyze.printNumberHydrogenBonds();
    //analyze.printBonded();
    //std::cout << analyze.getNumberOfBondBreakingEvents() << std::endl;
    analyze.printROO();
    
    
    return 0;
}
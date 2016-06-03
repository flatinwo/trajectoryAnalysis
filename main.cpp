//
//  main.cpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 10/12/15.
//  Copyright © 2015 Folarin Latinwo. All rights reserved.
//

#include <iostream>
#include <cassert>

#include "trajectoryAnalysis.h"

using namespace trajectoryAnalysis;

int main(int argc, const char * argv[]) {
    // insert code here...
    //std::cout << "Hello, World!\n";

    //const char* filename = "/Users/Folarin/Documents/vmd_views/water/patchy_colloids/test_snaps/dump22f_2b.xyz";
    const char* filename = "/Users/Folarin/Documents/Tests/testme5.xyz";
    Box mybox;
    mybox.box_lo = coord_t(3,0.);
    mybox.box_hi[0] = 16.176239;
    mybox.box_hi[1] = 16.362369;
    mybox.box_hi[2] = 16.401357;
    
    StructureChiralityAnalysis mychiral(filename,mybox);
    mychiral.computeRDF();
    
    
    return 0;
}

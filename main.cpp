
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

    //const char* filename = "/Users/Folarin/Documents/Tests/gro_files/test_ql/em_conf_240.gro";
    //const char* filename = "/Users/Folarin/Documents/Tests/LJ/dump-npt.gro";
    
    const char* filename = "/Users/Folarin/Documents/Tests/xyz_files/testFskt.xyz";
    
    Trajectory traj(filename,true,1,5);
    traj.setComputeFskt(true,4.0);
    traj.computeMeanSquaredDisplacement();
    traj.printFskt();
    std::ofstream ofile("MSD.dat");
    ofile << traj;
    ofile.close();
    
    /*BondOrderParameter bop(traj,6);
    bop.setCalcType(OrderParameter::Calc_t::LOCAL);
    bop.setRcutOff(5.0);
    bop.setMaxNumberOfNearestNeighbors(4);
    bop.compute();
    bop.print();
    
    AveragedTetrahedralOrderParameter analyzeq(traj);
    analyzeq.setRcutOff(5.0);
    analyzeq.setRmin(2.00);
    analyzeq.compute();
    analyzeq.print();*/
    
    /*Trajectory Ensemble(filename,FILETYPE::GRO,1,4);
    AveragedBondOrderParameter BOP(Ensemble,6);
    //BondOrderParameter BOP(Ensemble,6);

    BOP.setCalcType(OrderParameter::Calc_t::LOCAL);
    //BOP.setRcutOff(1.40);
    BOP.setRcutOff(0.45);
    BOP.setMaxNumberOfNearestNeighbors(4);
    BOP.compute();
    BOP.print();*/
    
    std::cout << "Hello world\n";
    
    
    /*const char* filename = "/Users/Folarin/Documents/Tests/testme5.xyz";
    
    
    
    Box mybox;
    mybox.box_lo = coord_t(3,0.);
    mybox.box_hi[0] = 16.176239;
    mybox.box_hi[1] = 16.362369;
    mybox.box_hi[2] = 16.401357;
    
    StructureChiralityAnalysis mychiral(filename,mybox);
    mychiral.computeRDF();*/
    
    
    return 0;
}

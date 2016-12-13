
//
//  main.cpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 10/12/15.
//  Copyright © 2015 Folarin Latinwo. All rights reserved.
//

#include <iostream>
#include <cassert>
#include <time.h>

#include "trajectoryAnalysis.h"

using namespace trajectoryAnalysis;

void printTimeTaken(clock_t& tStart){
    printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
}

int main(int argc, const char * argv[]) {
    // insert code here...
    //std::cout << "Hello, World!\n";

    //const char* filename = "/Users/Folarin/Documents/Tests/gro_files/test_ql/em_conf_240.gro";
    //const char* filename = "/Users/Folarin/Documents/Tests/LJ/dump-npt.gro";
   // const char* filename = "/Users/Folarin/Documents/Tests/xyz_files/em_conf.gro";
    
    //const char* filename = "/Users/Folarin/Documents/Tests/xyz_files/filenames";
    
    const char* filename = "/Users/Folarin/Documents/Tests/xdr_files/otrj.xtc";
    
    //trajectory_t system;
    //xdr_info info;
    
    //loadxtc(filename, system, info);
    
    
    //test(filename);
    clock_t tStart = clock();
    std::unique_ptr<Trajectory> traj;
    
    //traj.reset(new Trajectory(getCombinedTrajectory(filename)));
    
    traj.reset(new Trajectory(filename,FILETYPE::XTC,1,1,1000));
    
    traj->setComputeFskt(true,2.*M_PI/traj->getTrajectory()[0].box.box_period[0]);
    traj->computeMeanSquaredDisplacement();
    traj->setTimeStep(0.002); //ps
    traj->printFskt();
    
    printTimeTaken(tStart);
    
   /* function1d_t gr = traj->computeGofR();//need to write copy
    std::cout << gr;*/
    
    /*LSI lsi(*traj,2.3);
    lsi.setRcutOff(3.5);
    lsi.setBinSize(0.002);
    lsi.print();*/
    
    
    /*traj.setComputeFskt(true,4.0);
    traj.computeMeanSquaredDisplacement();
    traj.printFskt();
    std::ofstream ofile("MSD.dat");
    ofile << traj;
    ofile.close();*/
    
    /*BondOrderParameter bop(traj,6);
    bop.setCalcType(OrderParameter::Calc_t::LOCAL);
    bop.setRcutOff(5.0);
    bop.setMaxNumberOfNearestNeighbors(4);
    bop.compute();
    bop.print();*/
    
    
    /*OrderParameter op(*traj);
    op.setRcutOff(1.251);
    //op.setRmin(0.99);
    op.printNeighborDistribution();*/
    
    /*TetrahedralOrderParameter analyzeq(*traj);
    analyzeq.setRcutOff(5.0);
    //analyzeq.setRmin(2.00);
    
    analyzeq.setPositionTetrahedrality(true);
    //analyzeq.addRmax(1.25);
    //analyzeq.addRmax(1.50);
    //analyzeq.addRmax(1.75);
    analyzeq.addRmax(2.00);
    //analyzeq.addRmax(2.25);
    
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

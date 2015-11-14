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


//idea... load trajectory
//set flag of outputs and then out data
//out directory provides all info about the trajectory, gofr, MSD, FskT

using namespace trajectoryAnalysis;

std::string chirality(coord_list_t x, const Box& box){
    
    Box _box = box;
    assert(x.size() == 4);
    double_coord_t dr1 = distancesqandvec(x[1], x[0], _box);
    double_coord_t dr2 = distancesqandvec(x[2], x[1], _box);
    double_coord_t dr3 = distancesqandvec(x[3], x[2], _box);
    
    //compute zeta
    double c1 = (dr2.second[1]*dr3.second[2] - dr2.second[2]*dr3.second[1]);
    double c2 = (dr2.second[0]*dr3.second[2] - dr2.second[2]*dr3.second[0]);
    double c3 = (dr2.second[0]*dr3.second[1] - dr2.second[1]*dr3.second[0]);
    
    double zeta = dr1.second[0]*c1 - dr1.second[1]*c2 + dr1.second[2]*c3;
    zeta /= (sqrt(dr1.first*dr2.first*dr3.first));
    
    //std::cout << zeta <<"\n";
    
    if (zeta > 0.33) {
        return "1";
    }
    else if (zeta < -0.33){
        return "2";
    }
    else
        return "3";

}

void analyzeChiralityXYZ(xyzfile& snap, const Box& box){
    
    int natoms = snap.n;
    assert(natoms%4 == 0);
    int nmolecules = natoms/4;
    unsigned int k ,l;
    k=l= 0;
    
    for (int i=0; i<nmolecules; i++) {
        coord_list_t x;
        for (unsigned int j=0; j<4; j++) {
            x.push_back(snap.x[k++]);
        }
        std::string zeta = chirality(x,box);
        for (unsigned int j=0; j<4; j++) {
            snap.type[l++] = zeta;
        }
    }
    
}


int main(int argc, const char * argv[]) {
    // insert code here...
    std::cout << "Hello, World!\n";
    
    xyzfile data;
    Box box;
    box.box_hi[0] = 30.6910;
    box.box_hi[1] = 30.6910;
    box.box_hi[2] = 30.6910;
    box.updatePeriod();

    
    xyztrajectory_t traj;
    loadxyz("/Users/Folarin/Documents/vmd_views/chirality/dumpfinal3_4.xyz", traj);
    savexyz("/Users/Folarin/Documents/vmd_views/chirality/dumpfinal3_4_done.xyz", traj);
    
    
    loadxyz("//Users/Folarin/Documents/vmd_views/chirality/dumpfinal4.xyz", data);
    analyzeChiralityXYZ(data, box);
    savexyz("//Users/Folarin/Documents/vmd_views/chirality/dumpfinal4process.xyz", data);
    
    
    
    /*
    trajectoryAnalysis::trajectory_t A;
    trajectoryAnalysis::loadxyz("/Users/Folarin/Documents/vmd_views/water/patchy_colloids/test_snaps/confined/NVT_ensemble_T=0.09_density=0.1542.xyz", A);
     */
    
    //trajectoryAnalysis::Trajectory A("/Users/Folarin/Documents/vmd_views/water/patchy_colloids/test_snaps/confined/NVT_ensemble_T=0.09_density=0.1542.xyz");
    
    //trajectoryAnalysis::OrderParameter B("/Users/Folarin/Library/Developer/Xcode/DerivedData/Build/Products/Debug/log_energy_NVT_T=0.07099_density0.07984.txt");
    
    //B.computeAutoCorrelation();
    
    /*trajectoryAnalysis::Trajectory trajA("/Users/Folarin/Documents/vmd_views/water/patchy_colloids/test_snaps/confined/NVT_ensemble_T=0.09_density=0.1542.xyz");
    trajA.computeMeanSquaredDisplacement();
    std::cout << trajA << "\n";*/
    
    
    
    
    return 0;
}

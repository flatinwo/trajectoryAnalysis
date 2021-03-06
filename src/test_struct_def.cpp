//
//  test_struct_def.cpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 10/12/15.
//  Copyright © 2015 Folarin Latinwo. All rights reserved.
//

#include "struct_def.h"
#include "trajectoryAnalysis.h"
#include <iostream>

using namespace trajectoryAnalysis;


std::string chirality(coord_list_t& x, const Box& box){
    
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

std::string chiralityunwrap(coord_list_t& x, const Box& box){
    
    Box _box = box;
    assert(x.size() == 4);
    double_coord_t dr1 = distancesqandvec(x[1], x[0], _box);
    x[1] = x[0] + dr1.second;
    
    
    double_coord_t dr2 = distancesqandvec(x[2], x[1], _box);
    x[2] = x[1] + dr2.second;
    
    
    
    double_coord_t dr3 = distancesqandvec(x[3], x[2], _box);
    x[3] = x[2] + dr3.second;
    
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
    unsigned int k,l,m;
    k=l=m=0;
    
    for (int i=0; i<nmolecules; i++) {
        coord_list_t x;
        for (unsigned int j=0; j<4; j++) {
            x.push_back(snap.x[k++]);
        }
        std::string zeta = chiralityunwrap(x,box);
        //this is very very inefficent
        for (unsigned int j=0; j<4; j++) {
            snap.type[l++] = zeta;
            snap.x[m++] = x[j];
        }
    }
    
}


int main(int argc, const char * argv[]) {
    // insert code here...
    std::cout << "Hello, World!\n";
    
    /*const char* filename1 = "/Users/Folarin/Library/Developer/Xcode/DerivedData/Build/Products/Debug/logme.txt";
     coord_list_t _data;
     load(filename1, _data);
     
     std::vector<double> v;
     for (unsigned int i=0; i<200; i++) {
     for (unsigned int j=0; j < _data[i].size(); j++)
     if (j==1) v.push_back(_data[i][j]);
     }
     
     stats_utils::Correlations<double> corr(v);
     corr.computeAuto();
     corr.print();*/
    
    //const char* filename = "/Users/Folarin/Documents/vmd_views/water/initial_configt.xyz";
    const char* filename = "/Users/Folarin/Documents/vmd_views/water/patchy_colloids/test_snaps/dumpTraj2.xyz";
    
    Trajectory traj(filename,true,1,3);
    
    BondOrderParameter bop(traj,6);
    bop.setRcutOff(5.75);
    bop.setMaxNumberOfNearestNeighbors(4);
    bop.setThirdOrderInvariants(true);
    bop.compute();
    std::cout << bop.getQl() << "\t" << bop.getWl() << std::endl;
    bop.print();
    
    
    TetrahedralOrderParameter analyzeq(traj);
    analyzeq.setRcutOff(5.75);
    analyzeq.compute();
    analyzeq.print();
    
    //test boost
    /*accumulator_set<double, stats<tag::mean, tag::moment<2> > > acc;
     
     //push in some data
     acc(1.2);
     acc(2.3);
     acc(3.4);
     acc(4.5);
     
     //Display the results
     std::cout << "Mean:\t" << mean(acc) << std::endl;
     std::cout << "Moment:\t" << boost::accumulators::moment<2> (acc) << std::endl;
     
     
     //Now test average and variance
     OrderParameter o1("/Users/Folarin/Library/Developer/Xcode/DerivedData/Build/Products/Debug/logme.txt");
     o1.computeAverageAndVariance(1);
     
     
     std::cout <<"\n\nCorrelations\n\n";
     o1.computeAutoCorrelation(1);
     o1.printCorrelation();
     */
    
    /*
     xyzfile data;
     Box box;
     box.box_hi[0] = atof(argv[3]);
     box.box_hi[1] = atof(argv[4]);
     box.box_hi[2] = atof(argv[5]);
     box.updatePeriod();
     
     
     const char* filename0 = argv[1];
     const char* filename1 = argv[2];
     
     
     xyztrajectory_t traj;
     loadxyz(filename0, traj);
     
     std::vector<unsigned int> typecount(3,0), typecountmax(3,0);
     std::vector<std::string> typematch(3,"1");
     typematch[1] = "2"; typematch[2] = "3";
     
     for (unsigned int i=0; i<traj.size(); i++){
     analyzeChiralityXYZ(traj[i],box);
     for (unsigned int j=0; j<3; j++){
     typecount[j] = (int) std::count(traj[i].type.begin(),traj[i].type.end(),typematch[j]);
     typecountmax[j] = std::max(typecount[j],typecountmax[j]);
     }
     }
     typelog_t logtypes;
     for (unsigned int i=0; i<typematch.size(); i++) logtypes[typematch[i]] = typecountmax[i];
     VisualizerXYZ chiral(filename1);
     chiral.setTypeMax(logtypes);
     chiral.visualize(traj);
     
     
     
     //    savexyz(filename1, traj);
     std::cout << typecountmax[0] << "\t" << typecountmax[1] << "\t" << typecountmax[2] << "\n";
     
     */   
    
    
    
    return 0;
}

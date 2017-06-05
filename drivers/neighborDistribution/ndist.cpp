//
//  ndist.cpp
//  trajectoryAnalysis
//
//  Created by DrFlo on 6/5/17.
//  Copyright © 2017 Folarin Latinwo. All rights reserved.
//

//
//  main.cpp
//  neighborDistribution
//
//  Created by DrFlo on 6/5/17.
//  Copyright © 2017 Folarin Latinwo. All rights reserved.
//

#include <iostream>
#include <fstream>
#include "trajectoryAnalysis.h"


using namespace trajectoryAnalysis;

//compute 4,6,8,10 at l.

void log_call(int argc, const char* argv[]){
    // write log of command that calls functino
    std::ofstream runlog("run.log",std::ofstream::out);
    for (unsigned int i=0;i<argc;i++) runlog << argv[i] << " ";
    runlog << "\n";
    runlog.close();
}

int main(int argc, const char * argv[]) {
    
    if (!( argc == 3 | argc == 4)){
        std::cout << "Usage:\n" << argv[0] << "\t[string] filename [short] ref. index [opt., double=0.05] binsize\n";
        exit(-1);
    }
    
    log_call(argc, argv);
    
    std::string str(argv[1]);
    int m = atoi(argv[2]);
    double binsize=0.05;
    if (argc == 4) binsize=atof(argv[3]);
    
    Trajectory traj(str.c_str(),true,1,1);
    Hist1Dt ht = traj.computeNeighborDistribution(m,binsize);
    std::ofstream output("neighbor_distribution.txt");
    output << ht;
    
    output.close();
    
    return 0;
}



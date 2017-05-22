//
//  main.cpp
//  chiralityBOP
//
//  Created by Folarin Latinwo on 4/7/17.
//  Copyright © 2017 Folarin Latinwo. All rights reserved.
//

#include <iostream>
#include <fstream>
#include "trajectoryAnalysis.h"


using namespace trajectoryAnalysis;

//compute 4,6,8,10 at l.

int main(int argc, const char * argv[]) {

    if (argc != 3){
        std::cout << "Usage:\n" << argv[0] << "\tfilename [char]\tl-value [int]\n";
        exit(-1);
    }
    std::string str(argv[1]);
    int l = atoi(argv[2]);
    
    
    Box tbox;
    tbox.box_hi = coord_t(3,100.);
    tbox.updatePeriod();
    
    std::ofstream output("chirality_bop.txt",std::ofstream::app | std::ofstream::out);
    
    BondOrderChiralityAnalysis boca(str.c_str(),tbox,l,true);
    boca.computeBOP(output);
    
    output.close();
    
    return 0;
}

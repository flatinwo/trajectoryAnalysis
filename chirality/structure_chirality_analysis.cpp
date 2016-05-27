//
//  structure_chirality_analysis.cpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 5/23/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#include "structure_chirality_analysis.hpp"
namespace trajectoryAnalysis {
    
#define myEPS 0.01
    bool testchiral(unsigned int i, unsigned int j){
        if (floor((i+ myEPS)/4) == floor((j+myEPS)/4)) return true;
        else return false;
    }

    StructureChiralityAnalysis::StructureChiralityAnalysis(int argc, const char * argv[]):
    ChiralityAnalysis(argc, argv){
    }
    
    StructureChiralityAnalysis::StructureChiralityAnalysis(const char* filename, Box& box):
    ChiralityAnalysis(filename,box){
    }
    
    StructureChiralityAnalysis::~StructureChiralityAnalysis(){
    }
    
    void StructureChiralityAnalysis::computeRDF(){
        analyze();
        _rdf = new RadialDistributionFunction(_traj,_box);
        _rdf->addSimpleExcludeRule(testchiral);
        
        
        _rdf->compute();
        _rdf->print();
        
        delete _rdf;
    }
    
}
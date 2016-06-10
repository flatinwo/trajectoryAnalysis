//
//  structure_chirality_analysis.hpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 5/23/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#ifndef structure_chirality_analysis_hpp
#define structure_chirality_analysis_hpp

#include <stdio.h>
#include "chirality_analysis.hpp"
#include "radial_distribution_function.hpp"

namespace trajectoryAnalysis {
    class StructureChiralityAnalysis : public ChiralityAnalysis{
    public:
        StructureChiralityAnalysis(int argc, const char* argv[]);
        StructureChiralityAnalysis(const char* filename, Box& box);
        ~StructureChiralityAnalysis();
        
        void refresh();
        void computeRDF();
        
        void setBinSize(double);
        
    protected:
        double _binsize;
        RadialDistributionFunction* _rdf;
        
    };
}


#endif /* structure_chirality_analysis_hpp */

//
//  bond_order_chirality_analysis.hpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 3/16/17.
//  Copyright © 2017 Folarin Latinwo. All rights reserved.
//

#ifndef bond_order_chirality_analysis_hpp
#define bond_order_chirality_analysis_hpp

#include <stdio.h>
#include "chirality_analysis.hpp"
#include "bond_order_parameter.hpp"
#include <array>
#include <fstream>

namespace trajectoryAnalysis{
    class BondOrderChiralityAnalysis : public ChiralityAnalysis{
    public:
        BondOrderChiralityAnalysis(int argc, const char* argv[], short, bool=true);
        BondOrderChiralityAnalysis(const char* filename, Box& box, short, bool=true);
        ~BondOrderChiralityAnalysis();
        
        void refresh();
        void computeBOP(std::ostream& = std::cout);
        
        void setMaxNumberOfNeighbors(int);
        void setRcutOff(double);
        
    protected:
        
        bool useCOM;
        short _idx;
        std::array<double, 3> ref_point={0.,0.,0};
        
        BondOrderParameter* _bop;
        Trajectory* _trajobj;
        
    };
}



#endif /* bond_order_chirality_analysis_hpp */

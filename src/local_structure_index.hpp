//
//  local_structure_index.hpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 7/19/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#ifndef local_structure_index_hpp
#define local_structure_index_hpp

#include <stdio.h>
#include "order_parameter.hpp"

namespace trajectoryAnalysis {
    class LocalStructureIndex : public OrderParameter{
    public:
        LocalStructureIndex(Trajectory&, double rcutoff=3.7);
        ~LocalStructureIndex();
        
        void setBinSize(double);
        
        void compute();
        void clear();
        void print(const char* filename="LSI.dat");
        
        double getLSI();
        
    protected:
        bool _calculated;
        coord_list_t _lsi;
        double _I;
        double_unsigned_pair_t _condition;
        
        void _computeLSI(unsigned int);
        void _updateLSIframe();
        
        
        stats_utils::HistogramDynamic<double> _LSIHist;
        
    };
    
    typedef LocalStructureIndex LSI;
}

#endif /* local_structure_index_hpp */

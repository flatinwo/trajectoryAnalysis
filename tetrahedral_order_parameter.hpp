//
//  tetrahedral_order_parameter.hpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 2/2/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#ifndef tetrahedral_order_parameter_hpp
#define tetrahedral_order_parameter_hpp

#include <stdio.h>
#include "order_parameter.hpp"
#include "histogram_dynamic.h"
#include <fstream>
#include <memory>

namespace trajectoryAnalysis {
    
    /**
     \brief A derived class to perform specific calculations on/for Tetrahedral order parameters
            Value should range between -3 and +1.
     */
    
    class TetrahedralOrderParameter : public OrderParameter{
    public:
        TetrahedralOrderParameter(Trajectory&);
        ~TetrahedralOrderParameter();
        
        void setMaxNumberOfNearestNeighbors(unsigned int);
        
        double getQ();
        
        void print();
        void compute();
        
    protected:
        double _Q;
        bool _requireBinQvalues;
        std::vector<double> _Qs;
        std::vector<double> _Qframe;
        
        void _computeWithMaxNeighbors();            //also considered as maximum number of bonds
        void _computeQ(unsigned int);
        
        void _updateQframe();
        
        void _resize();
        
        void _refresh();
        void _refresh(unsigned int);
        
        void _openFiles();
        void _closeFiles();
        
        stats_utils::HistogramDynamic<double> _QHist;
        
    };
}
#endif /* tetrahedral_order_parameter_hpp */

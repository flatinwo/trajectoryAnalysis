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
    
    struct tHQs {
        tHQs(double rmaxsqd):_rmaxsqd(rmaxsqd){
            _numberstats = std::make_pair(0.,0.);
            _QHist = stats_utils::HistogramDynamic<double>(0.01,true);
        }
        
        double _rmaxsqd;
        corr_point_t _numberstats;
        std::vector<double> _Qs;
        std::vector<double> _Qframe;
        stats_utils::HistogramDynamic<double> _QHist;
    };
    
    /**
     \brief A derived class to perform specific calculations on/for Tetrahedral order parameters
            Value should range between -3 and +1.
     */
    
    class TetrahedralOrderParameter : public OrderParameter{
    public:
        TetrahedralOrderParameter(Trajectory&);
        ~TetrahedralOrderParameter();
        
        void setMaxNumberOfNearestNeighbors(unsigned int);
        void setPositionTetrahedrality(bool);
        void setMinimumCount(unsigned int);
        void addRmax(double);
        
        double getQ();
        
        void print();
        void compute();
        
    protected:
        double _Q;
        double_unsigned_pair_t _condition;
        unsigned int _minimum_count;
        
        std::vector<tHQs>* _tHQs;
        std::vector<short> _counts; 
        
        bool _requireBinQvalues;
        bool _requirePositionValues;
        
        
        std::vector<double> _Qs;
        std::vector<double> _Qframe;
        
        void _computeWithMaxNeighbors();            //also considered as maximum number of bonds
        
        virtual void _computeQ(unsigned int);
        virtual void _updateQframe();
        
        void _computeQR(unsigned int);
        
        void _resize();
        
        void _refresh();
        void _refresh(unsigned int);
        
        void _openFiles();
        void _closeFiles();
        
        stats_utils::HistogramDynamic<double> _QHist;
        
    };
}
#endif /* tetrahedral_order_parameter_hpp */

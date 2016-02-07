//
//  order_parameter.hpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 10/13/15.
//  Copyright © 2015 Folarin Latinwo. All rights reserved.
//

#ifndef order_parameter_hpp
#define order_parameter_hpp

#include <stdio.h>
#include "struct_def.h"
#include "trajectory.hpp"
#include <fstream>
#include <memory>


//set all calculation flags then output all results
//have default flags

namespace trajectoryAnalysis {
    
    /**
     \brief A class to perform specific calculations for order parameters
     */
    
    class OrderParameter{
    public:
        OrderParameter(Trajectory&);
        OrderParameter(const char*);
        ~OrderParameter();
        
        enum Calc_t {GLOBAL, LOCAL, GLOBALANDLOCAL};
        
        void setRcutOff(double);
        virtual void setMaxNumberOfNearestNeighbors(unsigned int);
        
        virtual void compute();
    
        void printCorrelation();
        
    protected:
        double _rcutoff;                            //maximum distance for neighbors
        
        Trajectory* _trajectory;
        
        coord_list_t _data;
        coord_t _correlation;
        double _average;
        double _variance;
        
        unsigned_list_t _number_of_neighbors;       //count of number of neighbors
        unsigned int _max_number_of_neighbors;
        unsigned int _nmolecules;
        bool _useMaxNumberOfNeighbors;
        std::vector<double_unsigned_pair1d_t> _nearest_neighbors;
        
        Snap* _snap;
        Calc_t _mode;
        coord_t _coord;
        
        void _computeNearestNeighbors();
        void _refreshNeighbors();
        void _refresh(unsigned int);
        
        //file operations
        std::vector< std::unique_ptr<std::ofstream> > _ofiles;
        virtual void _openFiles(){};
        virtual void _closeFiles(){};
        
        //virtual void _computeWithRcutOff(){};
        //virtual void _computeWithMaxNeighbors(){};            //also considered as maximum number of bonds
        
        void _restructureData();


        
    };
}

#endif /* order_parameter_hpp */

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
//also considering adding an exclude rule for neighbors

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
        void setRmin(double);
        void setCalcType(Calc_t mode);
        virtual void setMaxNumberOfNearestNeighbors(unsigned int);
        
        virtual void compute();
        
        void addSimpleIncludeRule(bool (*funcp)(unsigned, unsigned));
        void addSimpleExcludeRule(bool (*funcp)(unsigned, unsigned));   //usually for intermolecule
    
        void printCorrelation();
        void printNeighborDistribution(const char* filename="NeighborDistribution.dat");
        
    protected:
        Trajectory* _trajectory;
        double _rcutoff;                            //maximum distance for neighbors
        double _rminsq;                             //set _rminsq
        coord_list_t _data;
        
        unsigned_list_t _number_of_neighbors;       //count of number of neighbors
        unsigned_list_t _neighbor_count;            //2nd counter of neighbors
        
        unsigned int _max_number_of_neighbors;
        unsigned int _nmolecules;
        bool _useMaxNumberOfNeighbors;
        std::vector<double_unsigned_pair1d_t> _nearest_neighbors;
        
        Snap* _snap;
        stats_utils::HistogramDynamic<unsigned int>* _neighborHist;
        Calc_t _mode;
        coord_t _coord;
        bool _localflag;
        
        
        std::vector<bool (*)(unsigned, unsigned)> _exclude_funcptrs;
        std::vector<bool (*)(unsigned, unsigned)> _include_funcptrs;
        
        void _initialize();
        void _computeNearestNeighbors();
        void _computeNeighborDistribution();
        
        void _refreshNeighbors();
        //virtual void _refresh(unsigned int);
        virtual void _refresh();
        virtual void _update();
        
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

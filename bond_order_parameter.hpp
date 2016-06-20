//
//  bond_order_parameter.hpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 1/2/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#ifndef bond_order_parameter_hpp
#define bond_order_parameter_hpp

#include <stdio.h>
#include "order_parameter.hpp"
#include "histogram_dynamic.h"
#include <fstream>
#include <memory>

namespace trajectoryAnalysis {
    
    /**
     \brief A derived class to perform specific calculations on/for Steinhardt bond order parameters
     */
    
    class BondOrderParameter : public OrderParameter{
    public:
        BondOrderParameter(Trajectory&, int l=6);
        ~BondOrderParameter();
        
        void setLvalue(int);
        void setThirdOrderInvariants(bool=true);
        
        double getQl();
        double getWl();
        void print();
        
        void addLvalue(int);
        
        void compute();
        
    protected:
        int _l;                                     //make this a vector also to accommodate several L values
        bool _requireThirdOrderInvaraints;          //should I compute third order invariants
        bool _requireBinValues;
        component_list_t _Qlm;                      //make this a vector to accommodate several L values
        component_list_t _Wl_i;                     //ditto
        
        std::vector<component_list_t> _qlm_i;       //make this a vector also to accommodate several LM values
        
        double _Ql,_Wl;
        std::vector<double> _Qls;                   //Qls..
        std::vector<double> _Wls;                   //Wls..
        std::vector<double> _qli;                   //local q'is
        std::vector<double> _wli;                   //local w's
        
        virtual void _computeWithRcutOff();
        virtual void _computeWithMaxNeighbors();            //also considered as maximum number of bonds
        
        void _computeHarmonics(unsigned int, unsigned int);
        double _ThreeJSymbol(int, int, int, int, int, int); //well-defined Wigner-3j symbol
        
        void _computeQl();
        virtual void _computeql_i(unsigned int);
        void _computeql();
        
        void _computeWl();
        void _computeWl_i(unsigned int);
        virtual void _computewl_i(unsigned int);
        
        void _tallyLocals();
        
        
        void _resize();
        void _refresh();
        void _refresh(unsigned int);
        void _update();
        
        void _openFiles();
        void _closeFiles();
        
        stats_utils::HistogramDynamic<double> _qHist; //can make this a vector of vectors
        stats_utils::HistogramDynamic<double> _wHist; //can make this a vector of vectors
        
    };
}

#endif /* bond_order_parameter_hpp */

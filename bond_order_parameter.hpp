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

namespace trajectoryAnalysis {
    
    /**
     \brief A derived class to perform specific calculations on/for bond order parameters
     */
    
    class BondOrderParameter : public OrderParameter{
    public:
        BondOrderParameter(Trajectory&, int l=6);
        ~BondOrderParameter();
        
        enum Calc_t {GLOBAL, LOCAL, GLOBALANDLOCAL};
        
        void setLvalue(int);
        
        double getQl();
        double getWl();
        
        void addLvalue(int);
        
        void compute();
        
    protected:
        int _l;                                     //make this a vector also to accommodate several L values
        double _rcutoff;                            //maximum distance for neighbors
        bool _requireThirdOrderInvaraints;          //should I compute third order invariants
        unsigned_list_t _number_of_neighbors;       //count of number of neighbors
        component_list_t _Qlm;                      //make this a vector to accommodate several L values
        component_list_t _Wl_i;                     //ditto
        
        std::vector<component_list_t> _qlm_i;       //make this a vector also to accommodate several LM values
        Snap* _snap;
        Calc_t _mode;
        coord_t _coord;
        
        double _Ql,_Wl;
        
        void _computeHarmonics(unsigned int, unsigned int);
        double _ThreeJSymbol(int, int, int, int, int, int); //well-defined Wigner-3j symbol
        
        void _computeQl();
        void _computeWl();
        void _computeWl_i(unsigned int);
        void _resize();
        
    };
}

#endif /* bond_order_parameter_hpp */

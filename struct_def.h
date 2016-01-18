//
//  struct_def.h
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 10/12/15.
//  Copyright © 2015 Folarin Latinwo. All rights reserved.
//

#ifndef struct_def_h
#define struct_def_h

#include <vector>
#include <string>
#include <complex>
#include <map>
#include <cassert>

namespace trajectoryAnalysis {
    
    typedef std::vector<double> coord_t;
    typedef std::vector<bool> bool_list_t;
    typedef std::vector<coord_t> coord_list_t;
    typedef std::vector<unsigned int>  unsigned_list_t;
    
    typedef std::map<double,double> function1d_t;
    
    typedef std::map<int,coord_list_t> vector1d_t;
    
    typedef std::pair<double,coord_t> double_coord_t;
    typedef std::pair<double, double>  double_pair_t;
    typedef std::pair<double, unsigned int> double_unsigned_pair_t;
    typedef std::vector<double_unsigned_pair_t> double_unsigned_pair1d_t;
    
    typedef std::complex<double> component_t;
    typedef std::vector<component_t> shpdesc_t;
    
    typedef std::pair<unsigned long long, double> corr_point_t;
    typedef std::vector< corr_point_t > corr_point_list_t;
    
    typedef std::vector< std::complex<double> > component_list_t;
    typedef void* arg_t;
    
    
    typedef std::map<std::string, unsigned int> typelog_t;
    
    
    //overloaded operators
    
    inline coord_t operator+(const coord_t& x, const coord_t& y){
        assert(x.size() == y.size());
        coord_t z(x.size());
        for (unsigned int i=0; i<x.size(); i++)
            z[i] = x[i] + y[i];
        return z;
    }
    
    
    struct Box{
        
        Box():
        box_lo(3,0.0),
        box_hi(3,10.0),
        box_period(3,10.0),
        periodic(3,true)
        {}
        
        coord_t box_lo;
        coord_t box_hi;
        coord_t box_period;
        
        bool_list_t periodic;
        
        void clear(){
            box_lo.clear();
            box_hi.clear();
            box_period.clear();
        }
        
        void updatePeriod(){
            for (unsigned int i=0; i< box_lo.size(); i++)
                box_period[i] = box_hi[i] - box_lo[i];
        }
        
        
    };
    
    struct Snap{
        Snap():
        timestep(0)
        {}
        
        int timestep;
        Box box;
        coord_list_t _center_of_mass_list;
    };
    
    typedef std::vector<Snap> trajectory_t;
}


#endif /* struct_def_h */

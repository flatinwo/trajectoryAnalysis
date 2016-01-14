//
//  spatial.hpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 10/12/15.
//  Copyright © 2015 Folarin Latinwo. All rights reserved.
//

#ifndef spatial_hpp
#define spatial_hpp

#include <stdio.h>
#include "struct_def.h"


namespace trajectoryAnalysis {
    
    struct orient_info{
        orient_info():
        cutoffsq(4.0),
        iswithincutoffsq(true){
        }
        
        double cutoffsq;
        bool iswithincutoffsq;
    };
    
    //the squared distance between two particles
    double distancesq(coord_t& x1, coord_t& x2);
    
    //the squared distance between two particles
    double distancesq(coord_t& x1, coord_t& x2, coord_t&, bool_list_t&);
    
    //the squared distance between two particles in a box
    double distancesq(coord_t& x1, coord_t& x2, Box&);
    
    //the squared distance and vector distance between two particles in a box
    double_coord_t distancesqandvec(coord_t& x1, coord_t& x2, Box&);
    
    //the r, angular theta and phi between two particles
    coord_t spherical_orientation(coord_t& x1, coord_t& x2, Box&);
    
    //angular theta and phi between two particles
    coord_t orientation(coord_t& x1, coord_t& x2, Box&);
    
    //the angular theta and phi between two particles within a given distance and determine if within distance
    coord_t orientation(coord_t& x1, coord_t& x2, Box&, orient_info&);
    
    
    //the distance between two particles
    double distance(coord_t& x1, coord_t&x2);
    
    //the distance between two particles
    double distance(coord_t& x1, coord_t& x2, Box&);
    
    //vector separation between two particles
    coord_t distancevec(coord_t& x1, coord_t& x2, Box&);
    
    //are particles in the box
    bool are_particles_in_box(coord_t& x, Box& box);
    
    //are particles in the box
    bool are_particles_in_box(coord_list_t& x, Box& box);
    
    void pbc(coord_t& x, const coord_t& period, const bool_list_t&);
    
    void pbc(double& x, double period, bool periodic=true);
    void pbcwithfloor(double& x, double period, bool periodic=true);
}


#endif /* spatial_hpp */

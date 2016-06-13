//
//  io.hpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 10/12/15.
//  Copyright © 2015 Folarin Latinwo. All rights reserved.
//

#ifndef io_hpp
#define io_hpp

#include <stdio.h>
#include "struct_def.h"

#include <iostream>
#include <string>


namespace trajectoryAnalysis{
    
    void load(const char*, coord_list_t&, arg_t=0x0);
    void loadxyz(const char* filename, coord_list_t& x, arg_t=0x0);
    
    struct xyz_info{
        xyz_info() : outstream(0x0),instream(0x0){
            xres = coord_t(3, 0.0);
        }
        
        Box box;
        std::vector<std::string> type;
        std::ostream* outstream;
        std::istream* instream;
        std::map<std::string, int> reservoir;
        coord_t xres;
    };
    
    void loadxyz(const char* filename, trajectory_t& system, arg_t arg=0x0);
    void loadxyzfancy(const char* filename, trajectory_t& system, unsigned int=1, unsigned int=1);
    void loadgrofancy(const char* filename, trajectory_t& system, unsigned int=1, unsigned int=1);
    
    struct xyzfile{
        unsigned int n;
        coord_list_t x;
        std::string commentstr;
        std::vector<std::string> type;
    };
    
    
    // list xyzfile_t
    typedef std::vector<xyzfile> xyztrajectory_t;
    
    void loadxyz(const char*, xyzfile&);
    void loadxyz(std::istream&, xyzfile&);
    void loadxyz(const char*, xyztrajectory_t&);
    
    void savexyz(const char* filename, coord_list_t& x, xyz_info&);
    void savevarxyz(const char* filename, coord_list_t& x, xyz_info&);
    void savexyz(const char* filename, xyzfile&);
    void savexyz(const char*, xyztrajectory_t&);
    
    std::ostream& operator << (std::ostream&, Box&);
    std::ostream& operator << (std::ostream&, function1d_t&);
    std::ostream& operator << (std::ostream&, coord_list_t&);
    std::ostream& operator << (std::ostream&, coord_t&);
    
    //std::ostream& operator << (std::ostream&, shpdesc_t&);
    
    std::istream& operator >> (std::istream&, coord_list_t&);
    std::istream& operator >> (std::istream&, coord_t&);
    //std::istream& operator >> (std::istream&, shpdesc_t&);
    
    void deltype(coord_list_t& x, std::vector<std::string>& types, std::string);
}


#endif /* io_hpp */

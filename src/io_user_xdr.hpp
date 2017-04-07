//
//  io_user_xdr.hpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 8/30/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//
//  Downloaded on 8/30/16 from http://www.gromacs.org/Developer_Zone/Programming_Guide/XTC_Library
//  References http://manual.gromacs.org/online/xtc.html


#ifndef io_user_xdr_hpp
#define io_user_xdr_hpp

#ifdef IO_XDR

#include <stdio.h>
#include "io.hpp"
#include "xdrfile.h"


namespace trajectoryAnalysis {
    struct xdr_info{
        short nskipframe;
        short every;
        short com_id;
        int max_frame;
        unsigned long nframes_used;
        unsigned long count;
        
        xdr_info():
        nskipframe(1),
        every(1),
        com_id(1),
        max_frame(50000),
        nframes_used(0),
        count(0){
            
        }
        
    };
    
    void loadxtc(const char* filename, trajectory_t& system, xdr_info& );
    XDRFILE* myread_xtc(const char* filename);
    void test(const char* filename);

}

#endif

#endif /* io_user_xdr_hpp */

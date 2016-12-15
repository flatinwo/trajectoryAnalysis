//
//  io_user_xdr.cpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 8/30/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#ifdef IO_XDR

#include "io_user_xdr.hpp"
#include "xdrfile_xtc.h"

namespace trajectoryAnalysis{
    
    void loadxtc(const char* filename, trajectory_t& system, xdr_info& settings){
        int natoms=0, step= 0;
        float prec, time;
        matrix box;
        
        XDRFILE* xd = myread_xtc(filename);
        
        //sets up filename to be used by read_xtc
        std::string str(filename);
        const int SDIM=200;
        char cfilename[SDIM]="";
        assert(str.size()<SDIM);
        for (unsigned int i=0;i<str.size();i++) cfilename[i] = str[i];
        
        //get number of atoms and assume number of atoms is fixed
        //implement error handling at some point
        read_xtc_natoms(cfilename, &natoms);
        rvec x[natoms];
        
        
        assert(natoms>0);
        
        settings.count=0;
        int nsnap=0;
        int index = settings.com_id;
        int every = settings.every;
        
        coord_t com_idx(DIM,0.);
        
        while (read_xtc(xd, natoms, &step, &time, box, x, &prec) == 0) {
  
            
            if ((settings.count%settings.nskipframe) == 0) {
                //settings.nframes_used++;
                system.push_back(Snap());
                system[nsnap].timestep = step;
                
                for (unsigned int k=0; k<DIM; k++)system[nsnap].box.box_hi[k]=10.*box[k][k]; // convert to angstroms
                
                system[nsnap].box.updatePeriod();
                
                
                for (unsigned int i=0; i<natoms; i++) {
                    if ((i%every) == index-1){
                        for (unsigned int k=0;k<DIM;k++) com_idx[k]= 10.*x[i][k]; //convert to angstroms
                        system[nsnap]._center_of_mass_list.push_back(com_idx);
                        
                    }
                }
                nsnap++;
                settings.nframes_used++;
            }
            settings.count++;
            
            if (settings.nframes_used > settings.max_frame){
                std::cerr << "Reached maximum number of frames allowable " << settings.nframes_used << "\n";
                break;
            }
        }
    }
   
    XDRFILE* myread_xtc(const char* filename){
        return xdrfile_open(filename, "r");
    }
    
    void test(const char* filename){
        return;
        
    }
    
    
}

#endif

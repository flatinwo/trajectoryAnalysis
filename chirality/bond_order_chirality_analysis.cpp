//
//  bond_order_chirality_analysis.cpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 3/16/17.
//  Copyright © 2017 Folarin Latinwo. All rights reserved.
//

#include "bond_order_chirality_analysis.hpp"


namespace trajectoryAnalysis{
    BondOrderChiralityAnalysis::BondOrderChiralityAnalysis(const char* filename, Box& box, short l,bool usecom):
    ChiralityAnalysis(filename,box),useCOM(usecom){
        
        trajectory_t mytraj(_traj.size(),Snap());
        coord_t com(3,0.);
        
        for (unsigned i=0; i<_traj.size();i++){
            if (useCOM){
                _idx = _computeMolecularCOM(_traj[i],&com);
                mytraj[i]._center_of_mass_list = _molecular_com;
            }
            else mytraj[i]._center_of_mass_list = _traj[i].x;
            mytraj[i].box = _box;
            mytraj[i]._type_list = _traj[i].type;
        }
        
        visualizeCOM();
        _trajobj = new Trajectory(mytraj);
        _bop = new BondOrderParameter(*_trajobj,l);
        _bop->setCalcType(OrderParameter::Calc_t::LOCAL);
        _bop->setMaxNumberOfNearestNeighbors(6);
        _bop->setRcutOff(4.5);

        
    }
    
    BondOrderChiralityAnalysis::~BondOrderChiralityAnalysis(){
        delete _bop;
        delete _trajobj;
    }
    
    void BondOrderChiralityAnalysis::computeBOP(){
        _bop->compute();
        _bop->print();
        
        std::cout << _idx << "\t"
                << _bop->getqli(_idx) <<"\t"
                << _bop->getwli(_idx) << "\n";
    }
    
    
}

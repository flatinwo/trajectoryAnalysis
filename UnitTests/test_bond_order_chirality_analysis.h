//
//  test_bond_order_chirality_analysis.h
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 3/16/17.
//  Copyright © 2017 Folarin Latinwo. All rights reserved.
//

#ifndef test_bond_order_chirality_analysis_h
#define test_bond_order_chirality_analysis_h

#include "bond_order_chirality_analysis.hpp"

using namespace trajectoryAnalysis;

struct TestBondOrderChiralityAnalysis : public testing::Test{
    std::unique_ptr<BondOrderChiralityAnalysis> boca;
    Box tbox;
    
    void compute(){
        boca->computeBOP();
        //return tBOP(bop.getqli(ref),bop.getwli(ref),i);
        
    }
};

/*TEST_F(TestBondOrderChiralityAnalysis,SetUp){
    const char* fn = "/Users/Folarin/Documents/cplusplustutor/ladapo_test/Latinwo_tests/trajectoryAnalysis/trajectoryAnalysis/UnitTests/data/racemic_lattice.xyz";
    
    tbox.box_hi = coord_t(3,100.);
    tbox.updatePeriod();
    
    boca.reset(new BondOrderChiralityAnalysis(fn,tbox,4,true));
    compute();
    
    boca.reset(new BondOrderChiralityAnalysis(fn,tbox,6,true));
    compute();
    
    
    boca.reset(new BondOrderChiralityAnalysis(fn,tbox,8,true));
    compute();
    
    boca.reset(new BondOrderChiralityAnalysis(fn,tbox,10,true));
    compute();
}*/

TEST_F(TestBondOrderChiralityAnalysis,SetUp){
    const char* fn = "/Users/Folarin/Documents/cplusplustutor/ladapo_test/Latinwo_tests/trajectoryAnalysis/trajectoryAnalysis/UnitTests/data/racemic_best22.xyz";
    
    tbox.box_hi = coord_t(3,100.);
    tbox.updatePeriod();
    
    boca.reset(new BondOrderChiralityAnalysis(fn,tbox,4,true));
    compute();
    
    boca.reset(new BondOrderChiralityAnalysis(fn,tbox,6,true));
    compute();
    
    
    boca.reset(new BondOrderChiralityAnalysis(fn,tbox,8,true));
    compute();
    
    boca.reset(new BondOrderChiralityAnalysis(fn,tbox,10,true));
    compute();
}

#endif /* test_bond_order_chirality_analysis_h */

//
//  test_bond_order_parameter.h
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 3/14/17.
//  Copyright © 2017 Folarin Latinwo. All rights reserved.
//

#ifndef test_bond_order_parameter_h
#define test_bond_order_parameter_h


#include "bond_order_parameter.hpp"
#include "gtest/gtest.h"
#include "trajectoryAnalysis.h"

using namespace trajectoryAnalysis;

struct tBOP{
    double Q;
    double W;
    short idx;
    
    tBOP(double q=0.,double w=0.,short ix=0):Q(q),W(w),idx(ix){
        
    }
};

struct TestBondOrderParameter : public testing::Test{
    std::unique_ptr<Trajectory> traj;
    
     tBOP compute(short i=6, short ref=0, short max_neighbor=12){
         BondOrderParameter bop(*traj,i);
         bop.setCalcType(OrderParameter::Calc_t::LOCAL);
         bop.setMaxNumberOfNearestNeighbors(max_neighbor);
         bop.setRcutOff(10.);
         bop.compute();
         return tBOP(bop.getqli(ref),bop.getwli(ref),i);
        
    }
};

TEST_F(TestBondOrderParameter,FCC){
    const char* fn = "/Users/Folarin/Documents/cplusplustutor/ladapo_test/Latinwo_tests/trajectoryAnalysis/trajectoryAnalysis/UnitTests/data/fcc_13.xyz";
    
    traj.reset(new Trajectory(fn,true,1,1));
    
    //testing 6th order bops
    tBOP bop = compute(6,0,12);
    EXPECT_NEAR(bop.Q,0.57452,0.0001);
    EXPECT_NEAR(bop.W,-0.013161,0.0001);
    
    //testing 4th order bops
    bop = compute(4,0,12);
    EXPECT_NEAR(bop.Q,0.19094,0.0001);
    EXPECT_NEAR(bop.W,-0.159317,0.0001);
}



TEST_F(TestBondOrderParameter,HCP){
    const char* fn = "/Users/Folarin/Documents/cplusplustutor/ladapo_test/Latinwo_tests/trajectoryAnalysis/trajectoryAnalysis/UnitTests/data/hcp_13.xyz";
    
    traj.reset(new Trajectory(fn,true,1,1));
    
    //testing 6th order bops
    tBOP bop = compute(6,0,12);
    EXPECT_NEAR(bop.Q,0.48476,0.0001);
    EXPECT_NEAR(bop.W,-0.012442,0.0001);
    
    //testing 4th order bops
    bop = compute(4,0,12);
    EXPECT_NEAR(bop.Q,0.09722,0.0001);
    EXPECT_NEAR(bop.W,0.134097,0.0001);
}

TEST_F(TestBondOrderParameter,ICO){
    const char* fn = "/Users/Folarin/Documents/cplusplustutor/ladapo_test/Latinwo_tests/trajectoryAnalysis/trajectoryAnalysis/UnitTests/data/ico_13.xyz";
    
    //note that the file above's reference particles is on row 13
    //and has index (or ref.) 12
    
    traj.reset(new Trajectory(fn,true,1,1));
    
    //testing 6th order bops
    tBOP bop = compute(6,12);
    EXPECT_NEAR(bop.Q,0.66332,0.0001);
    EXPECT_NEAR(bop.W,-0.169754,0.0001);
    
    //testing 4th order bops
    bop = compute(4,12);
    EXPECT_NEAR(bop.Q,0.,0.0001);
    //EXPECT_NEAR(bop.W,0.,0.0001); // it appears W4 is not zero
}

TEST_F(TestBondOrderParameter,SC){
    const char* fn = "/Users/Folarin/Documents/cplusplustutor/ladapo_test/Latinwo_tests/trajectoryAnalysis/trajectoryAnalysis/UnitTests/data/sc_7.xyz";
    
    //note that the file above's reference particles is on row 13
    //and has index (or ref.) 12
    
    traj.reset(new Trajectory(fn,true,1,1));
    
    //testing 6th order bops
    tBOP bop = compute(6,0,6);
    //EXPECT_NEAR(bop.Q,0.66332,0.0001);
    EXPECT_NEAR(bop.W,0.013161,0.0001);
    
    //testing 4th order bops
    bop = compute(4,0,6);
    //EXPECT_NEAR(bop.Q,0.,0.0001);
    EXPECT_NEAR(bop.W,0.159317,0.0001); // it appears W4 is not zero
    
}


/*TEST_F(TestBondOrderParameter,OTHER){
    const char* fn = "/Users/Folarin/Documents/cplusplustutor/ladapo_test/Latinwo_tests/trajectoryAnalysis/trajectoryAnalysis/UnitTests/data/tg_7.xyz";
    
    //note that the file above's reference particles is on row 13
    //and has index (or ref.) 12
    
    traj.reset(new Trajectory(fn,true,1,1));
    
    //testing 6th order bops
    tBOP bop = compute(6,0,6);
    
    
    //EXPECT_NEAR(bop.Q,0.66332,0.0001);
//EXPECT_NEAR(bop.W,0.013161,0.0001);
    
    //testing 4th order bops
    bop = compute(4,0,6);
    //EXPECT_NEAR(bop.Q,0.,0.0001);
    //EXPECT_NEAR(bop.W,0.159317,0.0001); // it appears W4 is not zero
    
}*/



#endif /* test_bond_order_parameter_h */

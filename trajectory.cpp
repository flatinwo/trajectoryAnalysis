//
//  trajectory.cpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 10/12/15.
//  Copyright © 2015 Folarin Latinwo. All rights reserved.
//

#include "trajectory.hpp"
#include "io.hpp"
#include "spatial.hpp"
#include <cassert>

namespace trajectoryAnalysis {
    
#define TRAJMIN 1
    
    //constructor 1
    Trajectory::Trajectory():maxCorrelationLength(0),_time_step(100), _unfolded(false){
    }
    
    //constructor 2
    Trajectory::Trajectory(const char* filename, bool fancy, unsigned int index, unsigned int every){

	_time_step = 100; _unfolded = false; maxCorrelationLength=0;

        //load trajectory
        if (fancy) loadxyzfancy(filename, _trajectory, index, every);
        else loadxyz(filename, _trajectory);
        
        //compute time step and maximum correlation length
        computeTimeStep();
        computeMaxCorrelationLength();

    }
    
    
    //destructor
    Trajectory::~Trajectory(){
        
    }
    
#pragma mark GETS
    
    //return number of snaps
    int Trajectory::getNumberOfSnaps(){
        return (int) _trajectory.size();
    }
    
    //return time step
    double Trajectory::getTimeStep(){
        return _time_step;
    }
    
    //return trajectory
    trajectory_t& Trajectory::getTrajectory(){
        return _trajectory;
    }
    
    
#pragma mark OTHERS
    // write operator for trajectory
    void Trajectory::unfold(){
        //make sure you have required number of particles
        assert(_trajectory.size()>1);
        
        //coord_list_t positions;
        coord_t dr(3,0.);
        
        for (unsigned int i=0; i <_trajectory.size()-1; i++) {
            assert(_trajectory[i]._center_of_mass_list.size() == _trajectory[i+1]._center_of_mass_list.size());
            for (unsigned int j=0; j< _trajectory[i]._center_of_mass_list.size(); j++) {
                for (unsigned int k=0; k<_trajectory[i]._center_of_mass_list[j].size(); k++) { 
                    dr[k] = _trajectory[i+1]._center_of_mass_list[j][k] - _trajectory[i]._center_of_mass_list[j][k];
                    pbc(dr[k], _trajectory[i].box.box_period[k]);
                    _trajectory[i+1]._center_of_mass_list[j][k] =  _trajectory[i]._center_of_mass_list[j][k] + dr[k];
                }
                
            }
            
        }
        _unfolded = true;
    }
    
    
    void Trajectory::computeMeanSquaredDisplacement(){
        
        if (_unfolded != true)
            unfold();
        
        assert(maxCorrelationLength > 10);
        
        corr_point_t zeros;
        zeros.first = zeros.second = 0;
        corr_point_list_t correlation(maxCorrelationLength, zeros);
        
        for (unsigned int i=0; i < _trajectory.size(); i++) {
            for (unsigned int j=0; j < maxCorrelationLength; j++) {
                if (i+j >= _trajectory.size())
                    continue;
                else{
                    coord_t dr(3,0.);
                    for (unsigned int l=0; l < _trajectory[i]._center_of_mass_list.size(); l++){
                        for (unsigned int k=0; k<_trajectory[i]._center_of_mass_list[l].size(); k++) {
                            dr[k] = _trajectory[i+j]._center_of_mass_list[l][k] - _trajectory[i]._center_of_mass_list[l][k]; //write operator for this
                            correlation[j].second += dr[k]*dr[k];
                        }
                    }
                    correlation[j].first++;
                }
            }
        }
        
        _correlation.clear();
        double normalization = 1./(double) _trajectory[0]._center_of_mass_list.size();
        
        for (auto it=correlation.begin(); it != correlation.end(); ++it) {
            _correlation.push_back((it->second / (double) it->first)*normalization);
        }
        
    }
    
    void Trajectory::computeTimeStep(){
        
        if (_trajectory.size() > TRAJMIN) {
            _time_step = (double) _trajectory[1].timestep - _trajectory[0].timestep;
            assert(_time_step > 0);
        }
        else
            _time_step = 0;
    }
    
    void Trajectory::computeMaxCorrelationLength(){
        
        if (_trajectory.size() > 10) {
            maxCorrelationLength = (unsigned int) (floor((double)_trajectory.size()/3));
        }
        else
            maxCorrelationLength = (unsigned int)_trajectory.size() - 1;
    }
    
    void Trajectory::computeGofR(){
        
    }
    
    
#pragma mark OPERATORS
    
    std::ostream& operator << (std::ostream& _os, const Trajectory& _trajectory){
        
        for (unsigned int i=0; i < _trajectory._correlation.size(); i++) {
            _os << _trajectory._time_step * i << "\t" << _trajectory._correlation[i] << std::endl;
        }
        
        return _os;
        
    }
    
    
    
}

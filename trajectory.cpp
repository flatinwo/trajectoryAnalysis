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
#include <algorithm>
#include <fstream>
#include <sstream>

namespace trajectoryAnalysis {
    
#define TRAJMIN 1
#define SMALL 0.0001
    
    //constructor 1
    Trajectory::Trajectory():maxCorrelationLength(0),_time_step(100), _unfolded(false){
        _Fskts = nullptr;_ks=nullptr;_vanHovefxn=nullptr;_useVanHove=false;
    }
    
    //constructor 2
    //overlap true
    Trajectory::Trajectory(const char* filename, bool fancy, unsigned int index, unsigned int every){

        _time_step = 100; _unfolded = false; maxCorrelationLength=0;_computeFskt=false;_k=0;
        
        _Fskts = nullptr;_ks=nullptr;_vanHovefxn=nullptr;_useVanHove=false;

        //load trajectory
        if (fancy) loadxyzfancy(filename, _trajectory, index, every);
        else loadxyz(filename, _trajectory);
        
        //compute time step and maximum correlation length
        computeTimeStep();
        computeMaxCorrelationLength();

    }
    
    Trajectory::Trajectory(const char* filename, FILETYPE type, unsigned int index, unsigned int every){
        _time_step = 100; _unfolded = false; maxCorrelationLength=0;_computeFskt=false;_k=0;
        
        _Fskts=nullptr;_ks=nullptr;_vanHovefxn=nullptr;_useVanHove=false;
        
        
        if (type==GRO) loadgrofancy(filename,_trajectory,index,every);
        else { std::cerr << "Unknown filetype\n"; exit(-1);}
        
        //compute time step and maximum correlation length
        computeTimeStep();
        computeMaxCorrelationLength();
    }
    
    
    //destructor
    Trajectory::~Trajectory(){
        if (_Fskts!=nullptr && _computeFskt) _printFskt();
        if (_Fskts!=nullptr) delete _Fskts;
        if (_ks !=nullptr) delete _ks;
        if (_vanHovefxn != nullptr) delete _vanHovefxn;
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
    
#pragma mark SETS
    void Trajectory::setUseVanHove(bool flag){
        _useVanHove = flag;
    }
    
    void Trajectory::setComputeFskt(bool shouldi, double k){
        _computeFskt = shouldi;
        
        if (_computeFskt) {
            assert(k>0);
            _k = k;
        }
        else{
            if (_Fskts!=nullptr){
                delete _Fskts;
                delete _ks;
            }
        }
        
        if (_Fskts==nullptr && _computeFskt){
            _Fskts = new coord_list_t();
            _ks = new coord_t();
            
            _ks->push_back(_k);
            _Fskts->push_back(coord_t(maxCorrelationLength,0.));
            
            if (_vanHovefxn==nullptr && _useVanHove) _vanHovefxn = new Hists2Dt(_trajectory[0]._center_of_mass_list.size(),
                                                                 Hists1Dt(maxCorrelationLength)); //this assumes the number of particles remains fixed
        }
        
    }
    
    void Trajectory::setVanHoveBinSize(double binsize){
        if (_vanHovefxn==nullptr) {
            std::cerr << "Cannot set Van Hove Function without setting compute Fskt\n";
            exit(-1);
        }
        else{
            assert(binsize>0);
            for (auto& i :*_vanHovefxn) {
                for (auto& it : i) {
                    it.setBinSize(binsize);
                    it.insert(0.,0.);
                }
            }
        }
    }
    
    void Trajectory::setComputeFrequency(double k){
        if (!_computeFskt || _Fskts == nullptr) {
            std::cerr << "Cannot set compute frequency...\n"
                      << "Need to setComputeFskt first\n";
            exit(-1);
        }
        else{
            assert(k>0);
            if (std::find_if(_ks->begin(), _ks->end(), [k](double b){return std::abs(k-b) < SMALL;})==_ks->end()) {
                std::cout << "k value, " << k <<  "already exists in computation\n";
            }
            else{
                _ks->push_back(k);
                _Fskts->push_back(coord_t(maxCorrelationLength,0.));
            }
            
        }
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


        coord_t dr(3,0.);
        double deltar=0.;
        for (unsigned int i=0; i < _trajectory.size(); i++) {
            for (unsigned int j=0; j < maxCorrelationLength; j++) {
                if (i+j >= _trajectory.size())
                    continue;
                else{
                    for (unsigned int l=0; l < _trajectory[i]._center_of_mass_list.size(); l++){
                        deltar=0.;
                        for (unsigned int k=0; k<_trajectory[i]._center_of_mass_list[l].size(); k++) {
                            dr[k] = _trajectory[i+j]._center_of_mass_list[l][k] - _trajectory[i]._center_of_mass_list[l][k]; //write operator for this
                            correlation[j].second += dr[k]*dr[k];
                            if (_computeFskt) deltar += dr[k]*dr[k];
                        }
                        if (_computeFskt) {
                            assert(_Fskts->size()==_ks->size());
                            for (unsigned int f=0; f<_Fskts->size();f++) {
                                double kdeltar = (*_ks)[f]*sqrt(deltar);
                                (*_Fskts)[f][j] += sin(kdeltar)/kdeltar;
                            }
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
        
        //normalize Fskt
        if (_computeFskt) {
            for (unsigned int f=0; f<_Fskts->size();f++) {
                for (unsigned int j=1; j<(*_Fskts)[f].size(); j++) {
                     (*_Fskts)[f][j] /= (normalization*correlation[j].first);
                }
            }
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
    
    void Trajectory::_printFskt(){
        
        for (unsigned int i=0; i<_Fskts->size(); i++) {
            std::string str="Fskt"+ std::to_string((*_ks)[i])+".dat";
            std::ofstream myfile(str.c_str());
            for (unsigned int j=1; j< (*_Fskts)[i].size(); j++) {
                myfile << _time_step*j << "\t" << (*_Fskts)[i][j]/((*_Fskts)[i][1]) << std::endl;
            }
            myfile.close();
        }
    }
    
    
    
}

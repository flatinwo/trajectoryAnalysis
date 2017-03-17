//
//  trajectory.cpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 10/12/15.
//  Copyright © 2015 Folarin Latinwo. All rights reserved.
//

#include "trajectory.hpp"
#include "radial_distribution_function.hpp"
#include "io.hpp"
#include "spatial.hpp"
#include <cassert>
#include <algorithm>
#include <fstream>
#include <sstream>
#include "io_user_xdr.hpp"

namespace trajectoryAnalysis {
    
#define TRAJMIN 1
#define SMALL 0.0000001
    
    //constructor 1
    Trajectory::Trajectory():maxCorrelationLength(1),_time_step(100), _unfolded(false){
        _Fskts = nullptr;_ks=nullptr;_vanHovefxn=nullptr;_useVanHove=false;_maxframes=50000;
        _skipinfo=_skipinfo=new Skipt(false,1);_skipfactor=2;
    }
    
    //constructor 2
    //overlap true
    Trajectory::Trajectory(const char* filename, bool fancy, unsigned int index, unsigned int every){

        _time_step = 100; _unfolded = false; maxCorrelationLength=0;_computeFskt=false;_k=0;
        _maxframes=50000;
        
        _Fskts = nullptr;_ks=nullptr;_vanHovefxn=nullptr;_useVanHove=false;_skipinfo=_skipinfo=new Skipt(false,1);
        _skipfactor=2;

        //load trajectory
        if (fancy) loadxyzfancy(filename, _trajectory, index, every);
        else loadxyz(filename, _trajectory);
        
        //compute time step and maximum correlation length
        computeTimeStep();
        computeMaxCorrelationLength();
        _type = XYZ;

    }
    
    Trajectory::Trajectory(trajectory_t& traj){
        _type = UNSPECIFIED;
        _time_step = 100; _unfolded = false; maxCorrelationLength=0;_computeFskt=false;_k=0;
        _maxframes = 50000;_skipfactor=2;
        
        
        _Fskts=nullptr;_ks=nullptr;_vanHovefxn=nullptr;_useVanHove=false;_skipinfo=new Skipt(false,1);
        //compute time step and maximum correlation length
        
        _trajectory = traj;
        
        
        computeTimeStep();
        computeMaxCorrelationLength();
        
        
    }
    
    Trajectory::Trajectory(const char* filename, FILETYPE type, unsigned int index, unsigned int every, unsigned int maxframes){
        _time_step = 100; _unfolded = false; maxCorrelationLength=0;_computeFskt=false;_k=0;
        _maxframes = maxframes;_skipfactor=2;
        
        _Fskts=nullptr;_ks=nullptr;_vanHovefxn=nullptr;_useVanHove=false;_skipinfo=new Skipt(false,1);
        
        
        if (type==GRO) {
            loadgrofancy(filename,_trajectory,index,every);
            _type = GRO;
        }
#ifdef IO_XDR
        else if (type==XTC){
            xdr_info settings;
            settings.com_id = index;
            settings.every = every;
            settings.max_frame = _maxframes;
            loadxtc(filename,_trajectory,settings);
            _type = XTC;
        }
#endif
        else { std::cerr << "Unknown filetype\n"; exit(-1);}
        
        //compute time step and maximum correlation length
        computeTimeStep();
        computeMaxCorrelationLength();
    }
    
    
    //destructor
    Trajectory::~Trajectory(){
        std::cout << "Number of frames processed is " << getNumberOfSnaps() << std::endl;
        std::cout << "Trajectory is destructing...\n";
        if (_Fskts!=nullptr) delete _Fskts;
        if (_ks !=nullptr) delete _ks;
        if (_vanHovefxn != nullptr) delete _vanHovefxn;
        if (_skipinfo != nullptr) delete _skipinfo;
    }
    
#pragma mark GETS
    
    //return number of snaps
    int Trajectory::getNumberOfSnaps(){
        return (int) _trajectory.size();
    }
    
    //return number of frames, same as number of snaps
    int Trajectory::getNumberOfFrames(){
        return getNumberOfSnaps();
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
    void Trajectory::setMaxNumberOfFrames(int n){
        assert(n < INT_MAX);
        if (_trajectory.size() < n) n = (int)_trajectory.size();
        _maxframes = n;
        _trajectory.resize(_maxframes);
        computeMaxCorrelationLength();
        
    }
    
    void Trajectory::setMaxCorrelationLength(unsigned int mcl){
        assert(mcl > 10);
        if (3*mcl < _trajectory.size()) maxCorrelationLength = mcl;
        else std::cerr << "Specified max correlation length is out of range.\n"
                        << "I will use the max possible value of " << maxCorrelationLength << " .\n"
                        << "To override choose length <  " << _trajectory.size()/3 << ".\n";
    }
    
    void Trajectory::setTimeStep(double dt){
        _time_step = dt;
    }
    
    
    void Trajectory::setUseVanHove(bool flag){
        _useVanHove = flag;
    }
    
    void Trajectory::setSkipInfo(Skipt _skip,short sf){
        if (_skipinfo == nullptr) _skipinfo = new Skipt();
        _skipinfo->first = _skip.first;
        _skipinfo->second = _skip.second;
        _skipfactor = sf;
        
        assert(_skipfactor >= 1);
        assert(_skipinfo->second >= 1);
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
    
    void Trajectory::setUnFolded(bool flag){
        _unfolded = flag;
    }
    
    //not yet done
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
        //efficient to skip a few frames...
        //implement at a latter time... log type spacing or add flag for this
        for (unsigned int i=0; i < _trajectory.size(); i++) {
            for (unsigned int j=0; j < maxCorrelationLength; j++) {
                if (i+j >= _trajectory.size()) continue;
                else if (_skipinfo->first && j > _skipinfo->second && j%_skipfactor != 0) continue; //maybe make skipinfo or struct in case i want a different skip_info
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
                                if (kdeltar < SMALL) (*_Fskts)[f][j] += 1.0;
                                else (*_Fskts)[f][j] += sin(kdeltar)/kdeltar;
                            }
                        }

                    }
                    correlation[j].first++;
                }
            }
        }
        
        _correlation.clear();
        double normalization = 1./(double) _trajectory[0]._center_of_mass_list.size();
        
        _correlation.push_back(0.); //origin
        for (auto it=correlation.begin()+1; it != correlation.end(); ++it) {
            _correlation.push_back((it->second / (double) it->first)*normalization);
        }
        
        //normalize Fskt
        if (_computeFskt) {
            for (unsigned int f=0; f<_Fskts->size();f++) {
                for (unsigned int j=0; j<(*_Fskts)[f].size(); j++) {
                     (*_Fskts)[f][j] *= (normalization/correlation[j].first);
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
    
    
    function1d_t Trajectory::computeGofR(double binsize){
        
        stats_utils::HistogramDynamic < double > _hist;
        _hist.setBinSize(binsize);
        _hist.insert(0,0);
        
        for (auto& f : _trajectory) {
            coord_list_t* com = &(f._center_of_mass_list);
            double rcutsqd = *(std::min_element(f.box.box_period.begin(), f.box.box_period.end()));
            rcutsqd *= (0.5-SMALL);
            rcutsqd *= rcutsqd;
            
            
            for (unsigned int i=0; i<com->size(); i++) {
                for (unsigned int j=i+1; j<com->size(); j++) {
                    double rsq = distancesq((*com)[i], (*com)[j], f.box);
                    if (rsq < rcutsqd ) _hist.insert(sqrt(rsq),2.0);
                }
            }

        }
        
        function1d_t _gofr, normalization, probability;
        double _vol = 1.;
        for (auto& i : _trajectory[0].box.box_period) _vol *= i;
        
        double npart = _trajectory[0]._center_of_mass_list.size();
        double density = npart/_vol;
        double delta = 0.5*_hist.getBinSize();
        
        for (unsigned int i=0; i<_hist.getNumberOfBins(); i++) {
            std::pair<double,double> bin = _hist.getBin(i);
            double r = bin.first;
            double count = bin.second;
            
            double vol_shell = 4.0/3.0*M_PI*(pow(r+delta,3) - pow(r-delta,3));
            double np_ideal_gas = density*vol_shell;
            
            if (normalization.find(r) == normalization.end()) {
                normalization[r] = 0.0;
                probability[r] = 0.0;
            }
            
            normalization[r] += np_ideal_gas*npart;
            probability[r] += count;
        }
        
        for (function1d_t::iterator i=normalization.begin(); i!=normalization.end(); ++i)
            _gofr[i->first] = probability[i->first]/(i->second*_trajectory.size()); //needs work
        
        return _gofr;
    }
    
   /* void Trajectory::computeGofR(){
        
    }*/
    
    
#pragma mark OPERATORS
    
    std::ostream& operator << (std::ostream& _os, const Trajectory& _trajectory){
        
        for (unsigned int i=0; i < _trajectory._correlation.size(); i++) {
            _os << _trajectory._time_step * i << "\t" << _trajectory._correlation[i] << std::endl;
        }
        
        return _os;
        
    }
    
    Trajectory operator+(const Trajectory& rhs1, const Trajectory& rhs2){
        Trajectory trajr = rhs1;
        
        if (rhs1._trajectory[0]._center_of_mass_list.size() != rhs2._trajectory[0]._center_of_mass_list.size()) {
            std::cerr << "Warning! Number of atoms changed on combining trajectories\n";
        }
        
        trajr._trajectory.insert(trajr._trajectory.end(), rhs2._trajectory.begin(), rhs2._trajectory.end());
        
        trajr.computeTimeStep();
        trajr.computeMaxCorrelationLength();
        return trajr;
    }
    
    void Trajectory::printFskt(double dstep){
        if (_Fskts==nullptr){
            std::cerr << "Fskt has not been computed\n";
            exit(-1);
        }
        
        if (dstep == 1.|| _type == XTC) dstep = _trajectory[1].timestep - _trajectory[0].timestep;
        if (dstep <= 0.) dstep=1.;
        
        for (unsigned int i=0; i<_Fskts->size(); i++) {
            std::string str="Fskt"+ std::to_string((*_ks)[i])+".dat";
            std::ofstream myfile(str.c_str());
            for (unsigned int j=0; j< (*_Fskts)[i].size() && j < maxCorrelationLength; j++) {
                if (_skipinfo->first && j > _skipinfo->second && j%_skipfactor != 0) continue;
                myfile << dstep*_time_step*(j) << "\t" << (*_Fskts)[i][j] << std::endl;
            }
            myfile.close();
        }
    }
    
    
    
}

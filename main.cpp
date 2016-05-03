//
//  main.cpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 10/12/15.
//  Copyright © 2015 Folarin Latinwo. All rights reserved.
//

#include <iostream>
#include <cassert>
#include <algorithm>
#include <numeric>

#include "trajectoryAnalysis.h"


#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>


//idea... load trajectory
//set flag of outputs and then out data
//out directory provides all info about the trajectory, gofr, MSD, FskT
//check with rakesh q6 of snap

using namespace trajectoryAnalysis;
using namespace boost::accumulators;

#define COS30 0.8660254
#define COS50 0.6427876
#define MOLSIZE 3
#define BREAKINGEVENTS 40000

class WaterAnalysis{
    typedef std::vector<unsigned_list_t> unsigned_list2d_t;
    
    struct HydrogenBondBreakInfo{
        std::vector<unsigned int> frames;
        std::vector<std::pair<short,short>> oxygen_ids;   //first is initial O, second is final O
        std::vector < std::vector< std::pair<unsigned long,double> > > trajA; //typedef this
        std::vector < std::vector< std::pair<unsigned long,double> > > trajB;
    };
public:
     WaterAnalysis(Trajectory&);
    
    void setRCutOff(double);
    void setAngleCutOff(double);
    void setSkipFrame(int);
    void setTimeStep(double);
    
    int  getNumberOfBondBreakingEvents();
    
    void printBonded();
    void printNumberHydrogenBonds();
    void printROO();
    
    
    void compute();
    
protected:
    std::vector<unsigned_list2d_t> _bondedList;
    std::vector<HydrogenBondBreakInfo> _bondBreakInfos;
    
    int _skip_frame;
    double _time_step;
    
    void _computeHydrogenBonds(int i=0);
    void _computeAllHydrogenBonds();
    void _analyzeHydrogenBonds();
    void _computeROO();
    void _normalize();

    Trajectory* _traj;
    trajectory_t* _trajectory;
    
    corr_point_list_t trajAavg;
    corr_point_list_t trajBavg;
    
    double _rcutoff;
    double _anglecutoff;
    
};


WaterAnalysis::WaterAnalysis(Trajectory& traj):_traj(&traj),_skip_frame(50),_time_step(20.),_rcutoff(3.5),_anglecutoff(COS30){
    _trajectory = &_traj->getTrajectory();
}

void WaterAnalysis::setRCutOff(double cutoff){
    assert(cutoff > 0);
    _rcutoff = cutoff;
}

void WaterAnalysis::setAngleCutOff(double cutoff){
    assert(cutoff > 0);
    _anglecutoff = cutoff;
}

void WaterAnalysis::setSkipFrame(int skip){
    assert(skip > 0);
    _skip_frame = skip;
}

void WaterAnalysis::setTimeStep(double step){
    _time_step = step;
}

int WaterAnalysis::getNumberOfBondBreakingEvents(){
    int sum = 0;
    for (auto& i : _bondBreakInfos)
        sum += i.frames.size();
    return sum;
}

void WaterAnalysis::_computeROO(){
    double thetaa,thetab = 0;
    std::pair<int,double> info;
    long total_count = 0;
    short ostar = -MOLSIZE;
    
    trajAavg.resize(2*_skip_frame,corr_point_t(0,0));
    trajBavg.resize(2*_skip_frame,corr_point_t(0,0));

    
    for (unsigned int i=0; i < _bondBreakInfos.size(); i++) {
        _bondBreakInfos[i].trajA.resize(_bondBreakInfos[i].frames.size());
        _bondBreakInfos[i].trajB.resize(_bondBreakInfos[i].frames.size());
        if ( i%2 == 0 ){
            ostar += MOLSIZE;
        }
        for (unsigned int j=0; j <_bondBreakInfos[i].frames.size(); j++) {
            int frame = _bondBreakInfos[i].frames[j];
            
            short oa = _bondBreakInfos[i].oxygen_ids[j].first;
            short ob = _bondBreakInfos[i].oxygen_ids[j].second;
  
            
            for (int k = -_skip_frame; k < (int)_skip_frame; k++) {
                int l = frame + k;
                if (l < 0 || l >= _trajectory->size()) continue;
                Snap* _snap = &(*_trajectory)[l];
                coord_list_t* com = &(_snap->_center_of_mass_list);
                
                double ra = distance((*com)[ostar],
                                     (*com)[oa], _snap->box);
                double rb = distance((*com)[ostar],
                                     (*com)[ob], _snap->box);
                //std::cout << k << "\t" << ra << "\t" << rb << std::endl;
                if (ra > 4.00 && rb > 4.00) continue;
                thetaa = cosine_angle((*com)[ostar], (*com)[oa], (*com)[ostar+(i%2)+1], _snap->box);
                thetab = cosine_angle((*com)[ostar], (*com)[ob], (*com)[ostar+(i%2)+1], _snap->box);
                if ((ra < 4.00 && thetaa > COS50 && k <= 0) || (rb < 4.00 && thetab > COS50 && k > 0)){
                    //info.first = k;
                    //info.second = ra;_bondBreakInfos[i].trajA[j].push_back(info);
                    //info.second = rb; _bondBreakInfos[i].trajB[j].push_back(info);
                    
                    assert(k+_skip_frame < trajAavg.size());
                    trajAavg[k+_skip_frame].first++;
                    trajAavg[k+_skip_frame].second += ra;
                    
                    trajBavg[k+_skip_frame].first++;
                    trajBavg[k+_skip_frame].second += rb;
                    
                    total_count++;
                }
                        
            }
            
            
        }
        if (total_count > BREAKINGEVENTS*2*_skip_frame){
            std::cout << "I am going to break this\n";
            break;
        }
        
    }
    
    std::cout << (double) total_count /(16000*2*_skip_frame) << std::endl;
    
}

void WaterAnalysis::_analyzeHydrogenBonds(){
    assert(_bondedList.size() > 0);
    _bondBreakInfos.clear();
    _bondBreakInfos.resize(_bondedList[0].size());
    
    std::pair<short, short> ids;
    
    for (unsigned int i=0; i < _bondBreakInfos.size() ; i++){
        for (unsigned int j=_skip_frame; j < _bondedList.size() - 1; j++) {
            if ( _bondedList[j][i].size() == 0 || _bondedList[j+1][i].size() == 0) continue;
            if (std::find(_bondedList[j+1][i].begin(), _bondedList[j+1][i].end(), _bondedList[j][i][0]) == _bondedList[j+1][i].end()){
                _bondBreakInfos[i].frames.push_back(j);
                ids.first = _bondedList[j][i][0];
                ids.second = _bondedList[j+1][i][0]; //just pick the first oxygen
                _bondBreakInfos[i].oxygen_ids.push_back(ids);
            }
        }
    }
    
}

void WaterAnalysis::_computeHydrogenBonds(int i){

    Snap* _snap = &(*_trajectory)[i];
    coord_list_t* com = &(_snap->_center_of_mass_list);
    double rcutsqd = _rcutoff*_rcutoff;
    
    unsigned nmolecules = (unsigned) (ceil((double) com->size()/MOLSIZE));
    
    unsigned_list2d_t hbond(2*nmolecules);

    
    for (unsigned int i=0; i<com->size(); i+= MOLSIZE) {
        unsigned twomolid = 2*((unsigned) (floor ((double) i / (double) MOLSIZE))) ;
        //std::cout << twomolid << "\t" << twomolid+1-1 << "\t" << twomolid+2-1 << "\n";
        for (unsigned int j=0; j<com->size(); j+= MOLSIZE) {
            if (i==j) continue;
            double rsq = distancesq((*com)[i], (*com)[j], _snap->box);
            if (rsq < rcutsqd) {
                for (unsigned k=1;k<MOLSIZE;k++){
                    double x = cosine_angle((*com)[i], (*com)[i+k], (*com)[j], _snap->box);
                    if (x > _anglecutoff) {
                        hbond[twomolid+k-1].push_back(j);
                    }
                }
            }
        }
        
    }
    
    _bondedList.push_back(hbond);
    
}

void WaterAnalysis::_computeAllHydrogenBonds(){
    _bondedList.clear();
    for (unsigned int i=0; i<_trajectory->size(); i++) {
        _computeHydrogenBonds(i);
    }
}

void WaterAnalysis::compute(){
    //_computeHydrogenBonds(0);
    _computeAllHydrogenBonds();
    _analyzeHydrogenBonds();
    _computeROO();
    _normalize();
}

void WaterAnalysis::_normalize(){
    for (unsigned int i=0; i < trajBavg.size(); i++) {
        if (trajBavg[i].first != 0) trajBavg[i].second /= (double) trajBavg[i].first;
        if (trajAavg[i].first != 0) trajAavg[i].second /= (double) trajAavg[i].first;
    }
}

void WaterAnalysis::printROO(){
    for (int i=0; i<trajBavg.size(); i++) {
        std::cout << _time_step*(i - _skip_frame) << "\t" << trajAavg[i].second << "\t" << trajBavg[i].second <<
        "\t" << trajBavg[i].first << std::endl;
    }
}


void WaterAnalysis::printBonded(){
    for (unsigned int i=0; i<_bondedList.size(); i++) {
        std::cout << i ;
        for (unsigned int j=0; j<_bondedList[i][4].size(); j++) {
            std::cout << "\t" << _bondedList[i][4][j];
        }
        std::cout << "\n";
    }
}

void WaterAnalysis::printNumberHydrogenBonds(){
    for (unsigned int i=0; i<_bondedList.size(); i++) {
        unsigned nbond = 0;
        for (auto& j : _bondedList[i]) nbond += j.size();
        std::cout << i  << "\t" << nbond << "\n";
    }
}

int main(int argc, const char * argv[]) {
    // insert code here...
    //std::cout << "Hello, World!\n";

    const char* filename = "/Users/Folarin/Documents/vmd_views/water/patchy_colloids/test_snaps/dump22f_2b.xyz";
    
    Trajectory traj(filename,true,1,1);
    WaterAnalysis analyze(traj);
    analyze.compute();
    //analyze.printNumberHydrogenBonds();
    //analyze.printBonded();
    //std::cout << analyze.getNumberOfBondBreakingEvents() << std::endl;
    analyze.printROO();
    
    
    return 0;
}

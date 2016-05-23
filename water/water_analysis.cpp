//
//  water_analysis.cpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 5/17/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#include "water_analysis.hpp"
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <iostream>
#include "spatial.hpp"
#include <algorithm>
#include <numeric>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <memory>

using namespace boost::accumulators;

namespace trajectoryAnalysis {
    WaterAnalysis::WaterAnalysis(Trajectory& traj):_traj(&traj),_skip_frame(50),_time_step(20.),_rcutoff(3.5),_anglecutoff(COS30){
        _trajectory = &_traj->getTrajectory();
        _requireThetaValues = false;
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
        _openFiles();
        for (int i=0; i<trajBavg.size(); i++) {
            *_ofiles[0] << _time_step*(i - _skip_frame) << "\t" << trajAavg[i].second << "\t" << trajBavg[i].second <<
            "\t" << trajBavg[i].first << std::endl;
        }
        _closeFiles();
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

#pragma mark IOs
    void WaterAnalysis::_openFiles(){
        if (_requireThetaValues) _ofiles.resize(2);
        else _ofiles.resize(1);
        
        std::ostringstream os;
        os << std::setprecision(4);
        os << "ROO.txt";
        _ofiles[0].reset(new std::ofstream(os.str().c_str()));
        
        if (_requireThetaValues) {
            os.str("");
            os.clear();
            os << "Theta.txt";
            _ofiles[1].reset(new std::ofstream(os.str().c_str()));
        }
        
    }
    
    void WaterAnalysis::_closeFiles(){
        for (auto& i : _ofiles) i->close();
    }
    
    void WaterAnalysis::print(){
        printROO();
       /* assert(_Qframe.size()>0);
        double timestep = _trajectory->_time_step;
        
        _openFiles();
        for (unsigned int i=0; i<_Qframe.size(); i++) {
            double time = ((double) i)*timestep;
            *_ofiles[0] << time << "\t" << _Qframe[i] << std::endl;
        }
        
        if (_requireBinQvalues) *_ofiles[1] << _QHist << std::endl;
        _closeFiles();*/
    }
}

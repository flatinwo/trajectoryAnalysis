//
//  bond_order_parameter.cpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 1/2/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#include "bond_order_parameter.hpp"
#include "spatial.hpp"
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <iostream>
#include <cassert>
#include <gsl/gsl_sf_coupling.h>
#include <algorithm>
#include <sstream>
#include <iostream>

//when porting, will need to utilize boost and gsl
//now to compute local order and also include neighbor lists
//load trajectory with standard txt for trajectory
//useful references (1) http://www.pas.rochester.edu/~wangyt/algorithms/bop/
//                  (2) Y. Wang, S. Teitel and C. Dellago, Journal of Chemical Physics 122, 214722 (2005)
//

using namespace boost::math;

namespace trajectoryAnalysis {
    
#define MAX_NUMBER_OF_NEIGHBORS 50
    
    BondOrderParameter::BondOrderParameter(Trajectory& traj, int l):
    OrderParameter(traj),
    _l(l){
        assert(_l>1);
        _requireThirdOrderInvaraints = true;
        _Wl = 0.;
        _Wl_i.resize(_snap->_center_of_mass_list.size(), std::complex<double>(0.,0.));
        _qlm_i.resize(_snap->_center_of_mass_list.size());
        for (auto& m : _qlm_i) m.resize(2*_l+1, std::complex<double>(0,0));
        _Qlm.resize(2*_l+1,std::complex<double>(0,0));
    }
    
    BondOrderParameter::~BondOrderParameter(){
        
    }
    
    void BondOrderParameter::setLvalue(int l){
        assert(l>1);
        _l = l;
    }
    
    
    void BondOrderParameter::setThirdOrderInvariants(bool flag){
        _requireThirdOrderInvaraints = flag;
    }
    
    double BondOrderParameter::getQl(){
        return _Ql;
    }
    
    double BondOrderParameter::getWl(){
        if (!_requireThirdOrderInvaraints) std::cerr << "ERROR: BondOrderParameter::getWl()\n" <<
                                                        "Cannot get Wl without requiring it\n";
        return _Wl;
            
    }

    
    void BondOrderParameter::compute(){
        if (!_useMaxNumberOfNeighbors) _computeWithRcutOff();
        else _computeWithMaxNeighbors();
    }
    
    void BondOrderParameter::_computeWithMaxNeighbors(){
        trajectory_t* traj = &(_trajectory->_trajectory);
        
        for (unsigned int f=0; f<traj->size(); f++) {
            _snap = &(*traj)[f];
            _computeNearestNeighbors(); //get information about nearest neighbors
            coord_list_t* com = &(_snap->_center_of_mass_list);
            
            
            for (unsigned int i=0; i<com->size(); i++) {
                for (unsigned int j=0; j < _max_number_of_neighbors; j++) {
                    unsigned int k = _nearest_neighbors[i][j].second;
                    _computeHarmonics(i, k);
                }
                for (unsigned int m=0; m <_qlm_i[i].size(); m++) {
                    assert(_number_of_neighbors[i] > 0);
                    assert(_number_of_neighbors[i] == _max_number_of_neighbors);
                    _qlm_i[i][m] /= (double) _number_of_neighbors[i];
                    _Qlm[m] += _qlm_i[i][m];
                }
                if (_requireThirdOrderInvaraints) _computeWl_i(i);
            }
            _computeQl();
            if (_requireThirdOrderInvaraints) _computeWl();
        }
    }
    
    void BondOrderParameter::_computeWithRcutOff(){
        trajectory_t* traj = &(_trajectory->_trajectory);
        
        for (unsigned int f=0; f<traj->size(); f++) {
            _snap = &(*traj)[f];
            coord_list_t* com = &(_snap->_center_of_mass_list);
            for (unsigned int i=0; i<com->size(); i++) {
                for (unsigned int j=0; j<com->size(); j++) {
                    if (i==j) continue;
                    _computeHarmonics(i, j);
                    
                }
                for (unsigned int m=0; m <_qlm_i[i].size(); m++) {
                    assert(_number_of_neighbors[i] > 0);
                    _qlm_i[i][m] /= (double) _number_of_neighbors[i];
                    _Qlm[m] += _qlm_i[i][m];
                }
                if (_requireThirdOrderInvaraints) _computeWl_i(i);
            }
            _computeQl();
            if (_requireThirdOrderInvaraints) _computeWl();
        }
    }
    
    //compute Harmonics for system
    void BondOrderParameter::_computeHarmonics(unsigned int i, unsigned int j){
        
        _coord = spherical_orientation(_snap->_center_of_mass_list[i],
                                       _snap->_center_of_mass_list[j],
                                       _snap->box);
        if (_coord[0] > _rcutoff) return;
        else{
            double theta = _coord[1];
            double phi = _coord[2];
            assert(theta <= M_PI);
            assert(phi < 2*M_PI);
            for (int m = -_l; m < _l+1 ; m++) {
                _qlm_i[i][m + _l] += spherical_harmonic(_l, m, theta, phi);
            }
            _number_of_neighbors[i]++;
        }
    }
    
    double BondOrderParameter::_ThreeJSymbol(int ja, int jb, int jc, int ma, int mb, int mc){
        return gsl_sf_coupling_3j(2*ja, 2*jb, 2*jc, 2*ma, 2*mb, 2*mc);
    }
    
    //compute Ql of system
    void BondOrderParameter::_computeQl(){
        double sum = 0;
        for (unsigned int i=0; i<_Qlm.size();i++) sum += std::norm(_Qlm[i]);
        sum *= (4.*M_PI/(double) (2*_l + 1));
        _Ql = sqrt(sum)/(double) _snap->_center_of_mass_list.size();
        _Qls.push_back(_Ql);
        
        if (!_requireThirdOrderInvaraints) _refresh();
    }
    
    //compute Wl of system
    void BondOrderParameter::_computeWl_i(unsigned int i){
        
        for (int m1=-_l; m1 < _l+1; m1++) {
            for (int m2=-_l; m2 < _l+1; m2++) {
                for (int m3=-_l; m3<_l+1; m3++) {
                    if (m1 + m2 + m3 != 0) continue;
                    double w3j = _ThreeJSymbol(_l, _l, _l, m1, m2, m3);
                    _Wl_i[i] += w3j*(_qlm_i[i][m1+_l]*_qlm_i[i][m2+_l]*_qlm_i[i][m3+_l]);
                }
            }
        }
        double normalization = 0.;
        for (auto& m : _qlm_i[i]) normalization += std::norm(m);
        _Wl_i[i] /= std::pow(normalization, 1.5);
        
    }
    
    void BondOrderParameter::_computeWl(){
        _Wl = 0.;
        for (auto& i : _Wl_i) _Wl += i.real();
        _Wl /= (double) _Wl_i.size();
        _Wls.push_back(_Wl);
        _refresh();
    }
    
    void BondOrderParameter::print(){
        assert(_Qls.size()>0);
        if (_requireThirdOrderInvaraints){
            assert(_Wls.size() == _Qls.size());
        }
        double timestep = _trajectory->_time_step;
        
        _openFiles();
        for (unsigned int i=0; i<_Qls.size(); i++) {
            double time = ((double) i)*timestep;
            *_ofiles[0] << time << "\t" << _Qls[i] << std::endl;
            if (_requireThirdOrderInvaraints) *_ofiles[1] << time << "\t" << _Wls[i] << std::endl;
        }
        _closeFiles();
    }
    
    void BondOrderParameter::_refresh(){
        //refresh _qlmi
        for (auto& m : _qlm_i)
            for (auto& l : m) l = std::complex<double>(0.,0.);
        
        //refresh _Qlm
        for (auto& i: _Qlm) i = std::complex<double>(0.,0.);
        
        _refreshNeighbors();
        
        //refresh Wl_i
        if (_requireThirdOrderInvaraints)
            for (auto& i : _Wl_i) i = std::complex<double>(0.,0.);
        
    }
    
    void BondOrderParameter::_openFiles(){
        if (_requireThirdOrderInvaraints) _ofiles.resize(2);
        else _ofiles.resize(1);
        
        std::ostringstream os;
        os << std::setprecision(4);
        os << "Q_" << _l << ".txt";
        _ofiles[0].reset(new std::ofstream(os.str().c_str()));
        
        if (_requireThirdOrderInvaraints) {
            os.str("");
            os.clear();
            os << "W_" << _l << ".txt";
            _ofiles[1].reset(new std::ofstream(os.str().c_str()));
        }

    }
    
    void BondOrderParameter::_closeFiles(){
        for (auto& i : _ofiles) i->close();
    }
}

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

//when porting, will need to utilize boost and gsl
//now to compute local order and also include neighbor lists
//load trajectory with standard txt for trajectory
//useful references (1) http://www.pas.rochester.edu/~wangyt/algorithms/bop/
//                  (2) Y. Wang, S. Teitel and C. Dellago, Journal of Chemical Physics 122, 214722 (2005)
//

using namespace boost::math;

namespace trajectoryAnalysis {
    
#define MAX_NUMBER_OF_NEIGHBORS 120
    
    BondOrderParameter::BondOrderParameter(Trajectory& traj, int l):
    OrderParameter(traj),
    _l(l),
    _mode(GLOBAL){
        assert(_l>1);
        _rcutoff = 0.74;
        _requireThirdOrderInvaraints = true;
        _useMaxNumberOfNeighbors = false;
        _max_number_of_neighbors = 0;
        _snap = &_trajectory->_trajectory[0];
        _nearest_neighbors.resize(_snap->_center_of_mass_list.size(),
                                  double_unsigned_pair1d_t (MAX_NUMBER_OF_NEIGHBORS, std::pair<double, unsigned int>(0.,0)));
        _number_of_neighbors.resize(_snap->_center_of_mass_list.size(),0.);
        
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
    
    void BondOrderParameter::setRcutOff(double rcutoff){
        assert(rcutoff > 0.);
        _rcutoff = rcutoff;
    }
    
    void BondOrderParameter::setMaxNumberOfNearestNeighbors(unsigned int n_nghbrs){
        assert(n_nghbrs > 0);
        _max_number_of_neighbors = n_nghbrs;
        _useMaxNumberOfNeighbors = true;
        
    }
    
    double BondOrderParameter::getQl(){
        return _Ql;
    }
    
    double BondOrderParameter::getWl(){
        if (!_requireThirdOrderInvaraints) std::cerr << "ERROR: BondOrderParameter::getWl()\n" <<
                                                        "Cannot Wl without requiring it\n";
        return _Wl;
            
    }
    
    void BondOrderParameter::_computeNearestNeighbors(){
        coord_list_t* com = &(_snap->_center_of_mass_list);
        double rcutsqd = _rcutoff*_rcutoff;
        
        for (unsigned int i=0; i<com->size(); i++) {
            unsigned int k=0;
            for (unsigned int j=0; j<com->size(); j++) {
                if (i==j) continue;
                double rsq = distancesq((*com)[i], (*com)[j], _snap->box);
                if (rsq < rcutsqd) {
                    _nearest_neighbors[i][k].first = rsq;
                    _nearest_neighbors[i][k].second = j;
                    k++;
                    assert(k < MAX_NUMBER_OF_NEIGHBORS);
                }
            }
            std::sort(_nearest_neighbors[i].begin(), _nearest_neighbors[i].begin()+k); //implement sort up to
            assert(_nearest_neighbors[i].size() >= _max_number_of_neighbors);
        }
    }
    
    void BondOrderParameter::compute(){
        if (!_useMaxNumberOfNeighbors) _computeWithRcutOff();
        else _computeWithMaxNeighbors();
    }
    
    void BondOrderParameter::_computeWithMaxNeighbors(){
        _computeNearestNeighbors(); //get information about nearest neighbors
        coord_list_t* com = &(_snap->_center_of_mass_list);
        
        
        for (unsigned int i=0; i<com->size(); i++) {
            for (unsigned int j=0; j < _max_number_of_neighbors; j++) {
                unsigned int k = _nearest_neighbors[i][j].second;
                _computeHarmonics(i, k);
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
        
        for (auto& i : _number_of_neighbors)
            std::cout << i << std::endl;
    }
    
    void BondOrderParameter::_computeWithRcutOff(){
        coord_list_t* com = &(_snap->_center_of_mass_list);
        for (unsigned int i=0; i<com->size(); i++) {
            for (unsigned int j=0; j<com->size(); j++) {
                if (i==j) continue;
                _computeHarmonics(i, j);
                
            }
            //for (auto& z : _qlm_i[i]) z /= (double)_number_of_neighbors[i];  //for latter write operators
            //_Qlm += _qlm_i[i];
            for (unsigned int m=0; m <_qlm_i[i].size(); m++) {
                assert(_number_of_neighbors[i] > 0);
                _qlm_i[i][m] /= (double) _number_of_neighbors[i];
                _Qlm[m] += _qlm_i[i][m];
            }
            if (_requireThirdOrderInvaraints) _computeWl_i(i);
        }
        _computeQl();
        if (_requireThirdOrderInvaraints) _computeWl();
        
        for (auto& i : _number_of_neighbors)
            std::cout << i << std::endl;
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
        for (auto& i : _Wl_i) _Wl += i.real();
        _Wl /= (double) _Wl_i.size();
    }
}
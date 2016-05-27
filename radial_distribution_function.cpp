//
//  radial_distribution_function.cpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 5/25/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#include "radial_distribution_function.hpp"
#include <boost/math/special_functions/factorials.hpp>
#include "spatial.hpp"
#include "io.hpp"
#include <set>


//write a type converter from xyztrajectory_t and trajectory_t
//set only allows distinct containers
//maybe typedef boost

namespace trajectoryAnalysis {
    RadialDistributionFunction::RadialDistributionFunction(xyztrajectory_t& xyztraj, Box& box, double binsize):_trajx(&xyztraj),_box(&box),_bin_size(binsize){
        _initialize();
    }
    
    /* determine the number of unique types in snapshot */
    void RadialDistributionFunction::_analyzeFrames(){
        _distinct_container.clear();
        for (auto curr_it=_trajx->begin(), end = _trajx->end(); curr_it != end; ++curr_it) {
            for (auto curr_t = curr_it->type.begin(), end_t=curr_it->type.end(); curr_t != end_t; ++curr_t) {
                _distinct_container.insert(*(curr_t));
            }
        }
        _nunique_types = (int)_distinct_container.size();
    }
    
    void RadialDistributionFunction::addSimpleExcludeRule(bool (*funcp)(unsigned, unsigned)){
        _exclude_funcptrs.push_back(funcp);
    }
    
    void RadialDistributionFunction::_updateHist(int l, double r){
        _Hists[l].insert(r,2.0);
        _hist_infos[l]._ntypei++;
        _hist_infos[l]._ntypej++;
    }
    
    //might need orientational dependence here or do some preprocessing
    //this is dirty, can reduce by a factor of two by being smart about counting
    void RadialDistributionFunction::compute(){
        _number_of_frames = 0;
        for (unsigned int i=0; i<_trajx->size(); i++) {//skip_frames;
            xyzfile* xyzi = &(*_trajx)[i];
            coord_list_t* x = &(xyzi->x);
            _number_of_frames++;
            unsigned int natoms = (unsigned int) x->size();
            
            for (unsigned int j=0; j< natoms; j++) {
                for (unsigned int k=j+1; k<natoms; k++) {
                    bool testex = false;
                    //skip all coordinates that satisfy exclude rule
                    for (auto& ex : _exclude_funcptrs) {
                        if (ex(j,k)) testex = true;
                        
                        
                        if (testex) break;
                    }
                    if (testex) continue;
                    
                    double rsq = distancesq((*x)[j],(*x)[k],*_box);
                    if (rsq < _rmax) {
                        double r = sqrt(rsq);
                        for (unsigned int l=0; l<_hist_infos.size(); l++) {
                            if (xyzi->type[j].compare(xyzi->type[k])==0) { //two equal strings
                                if (xyzi->type[j].compare(_hist_infos[l]._typei) == 0) {
                                    _updateHist(l, r);
                                    break;
                                }
                            }
                            else{
                                if (xyzi->type[j]==_hist_infos[l]._typei
                                    && xyzi->type[k] ==_hist_infos[l]._typej) {
                                    _updateHist(l, r);
                                    break;
                                }
                                else if (xyzi->type[j]==_hist_infos[l]._typej
                                         && xyzi->type[k] ==_hist_infos[l]._typei){
                                    _updateHist(l, r);
                                    break;
                                }
                                
                            }
                        }
                    }
                }
            }
        }
    }
    
    void RadialDistributionFunction::_normalize(){
        assert(_Hists.size() == _hist_infos.size());
        _vol = _box->getVolume();
        for (unsigned int i=0; i<_Hists.size(); i++) _GofRs[i] = _normalize(_Hists[i], _hist_infos[i]);
    }
    
    function1d_t RadialDistributionFunction::_normalize(const stats_utils::HistogramDynamic<double>& _hist, _normalize_info& info){
        function1d_t _gofr, normalization, probability;
        stats_utils::HistogramDynamic<double> temp(_hist);
        
        double density = info._ntypei/_vol; //this density times number of frames
        double delta = 0.5*temp.getBinSize();
        
        for (unsigned int i=0; i<temp.getNumberOfBins(); i++) {
            std::pair<double,double> bin = temp.getBin(i);
            double r = bin.first;
            double count = bin.second;
            
            double vol_shell = 4.0/3.0*M_PI*(pow(r+delta,3) - pow(r-delta,3));
            double np_ideal_gas = density*vol_shell;
            
            if (normalization.find(r) == normalization.end()) {
                normalization[r] = 0.0;
                probability[r] = 0.0;
            }
            
            normalization[r] += np_ideal_gas*info._ntypej;
            probability[r] += count;
        }
        
        for (function1d_t::iterator i=normalization.begin(); i!=normalization.end(); ++i)
            _gofr[i->first] = probability[i->first]/(4*i->second/(info._ntypei)); //needs work
        return _gofr;
        
    }
    
    void RadialDistributionFunction::update(){
        compute();
    }
    
    void RadialDistributionFunction::print(){
        _normalize();
        std::cout << "After " << _number_of_frames << " frames the rdfs are in gr files\n";
        _openFiles();
        for (unsigned int i=0; i< _ofiles.size(); i++) *(_ofiles[i]) << _GofRs[i];
        _closeFiles();
    }
    
    void RadialDistributionFunction::_openFiles(){
        _ofiles.resize(_hist_infos.size());
        for (unsigned int i=0; i<_ofiles.size(); i++) {
            std::string str = "gr"+_hist_infos[i]._typei+_hist_infos[i]._typej+".txt";
            _ofiles[i].reset(new std::ofstream(str.c_str()));
        }
    }
    
    void RadialDistributionFunction::_closeFiles(){
        for (auto& i :_ofiles) i->close();
    }
    
    void RadialDistributionFunction::_initialize(){
        _analyzeFrames();
        assert(_nunique_types > 0);
        assert(_distinct_container.size()>0);
        unsigned int size=0;
        //number of pairs
        if (_nunique_types == 1) size = 1;
        else
            size = _nunique_types + (unsigned int)
            boost::math::factorial<double>(_nunique_types)/(2*(boost::math::factorial<double>(_nunique_types-2))); //because we have pair correlation function
        
        
        _hist_infos.resize(size);
        
        unsigned int k=0; //quick counter
        for (auto it=_distinct_container.begin(),end=_distinct_container.end(); it != end; ++it) {
            _hist_infos[k]._typei = _hist_infos[k]._typej = *it;
            _hist_infos[k]._same_type = true;
            _hist_infos[k]._ntypei = _hist_infos[k]._ntypej = 0.;
            k++;
        }
        
        for (unsigned int i=0; i<_nunique_types-1; i++) {
            for (unsigned int j=i+1; j<_nunique_types; j++) {
                _hist_infos[k]._typei = _hist_infos[i]._typei;
                _hist_infos[k]._typej = _hist_infos[j]._typej;
                _hist_infos[k]._same_type = false;
                _hist_infos[k]._ntypei = _hist_infos[k]._ntypej = 0.;
                k++;
            }
        }
        assert(k==_hist_infos.size());
        
        //set up dimension to be square of (half of box period in smallest dimension)
        _rmax = *(std::min_element(_box->box_period.begin(), _box->box_period.end()));
        _rmax *= 0.5;
        _rmax *= _rmax;
        
        //set-up histogram
        _Hists.resize(_hist_infos.size());
        _GofRs.resize(_hist_infos.size());
        for (auto& i : _Hists){
            i.setBinSize(_bin_size);
            i.insert(0,0);
        }

    }
    
}
//
//  correlations.cpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 2/4/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#include <stdio.h>
#include "correlations.h"

namespace stats_utils{
    
    template <class T>
    void Correlations<T>::computeAuto(){
        if (!_avgandvarcomputed) _computeAvgandVar();
        for (unsigned long i=0; i < _data_sets.size(); i++){
            corr_point_t zeros; zeros.first = zeros.second = 0.;
            corr_point_list_t correlation((int) floor(((double) _data_sets[i]->size())/3.),zeros);
            std::transform(_data_sets[i]->begin(),_data_sets[i]->end(), _data_sets[i]->begin(),
                           std::bind2nd(std::minus<double>(),_average[i]));
            
            for (unsigned int ik=0; ik < _data_sets[i]->size(); ik++) {
                for (unsigned int j=ik; j< _data_sets[i]->size(); j++) {
                    if (j-ik >= correlation.size())
                        continue;
                    correlation[j-ik].second += (*_data_sets[i])[ik]*(*_data_sets[i])[j];
                    correlation[j-ik].first++;
                }
            }
            _autocorrelations[i].clear();
            for (auto it=correlation.begin(); it != correlation.end(); ++it) {
                double corr = (it->second/(double) it->first)/_variance[i];
                _autocorrelations[i].push_back(corr);
            }
        }
        
    }
    
    template <class T>
    void Correlations<T>::_computeAvgandVar(){
        for (unsigned long i=0; i<_data_sets.size(); i++) {
            accumulator_set<T, stats< tag::mean, tag::variance > > acc; //typedef this;
            std::for_each(_data_sets[i]->begin(), _data_sets[i]->end(), boost::bind<void>(boost::ref(acc), _1));
            _average.push_back(mean(acc));
            _variance.push_back(variance(acc));
        }
        _avgandvarcomputed = true;
    }
    
    template class Correlations<int>;
    template class Correlations<double>;
    
}
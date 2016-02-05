//
//  correlations.h
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 2/3/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#ifndef correlations_h
#define correlations_h

#include <vector>
#include <cassert>
#include <map>
#include <algorithm>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#include <boost/bind.hpp>
#include <boost/ref.hpp>

using namespace boost::accumulators;

namespace stats {
    typedef std::pair<unsigned long long, double> corr_point_t;
    typedef std::vector< corr_point_t > corr_point_list_t;
    
    template <class T>
    class Correlations{
    public:
        Correlations(std::vector<T>& data){
            _data_sets.push_back(&data);
        }
        
        void addData(std::vector<T>& data){
            _data_sets.push_back(&data);
        }
        
        void computeAuto();
        
        void computeCross();
        
        void print();
        
    protected:
        std::vector< std::vector <T>* > _data_sets;
        std::vector<double> _average, _variance;
        
        std::vector< std::vector < double > > _autocorrelations;
        std::vector< std::vector < double > > _crosscorrelations;
        
        void _computeAvgandVar(){
            for (unsigned long i=0; i<_data_sets.size(); i++) {
                accumulator_set<T, stats< tag::mean, tag::variance > > acc; //typedef this;
                std::for_each((*_data_sets)[i].begin(), (*_data_sets)[i].end(), boost::bind<void>(boost::ref(acc), _1));
                _average.push_back(mean(acc));
                _variance.push_back(variance(acc));
            }
        }
        
        
    };
}


#endif /* correlations_h */

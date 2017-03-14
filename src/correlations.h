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
#include <iostream>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#include <boost/bind.hpp>
#include <boost/ref.hpp>

using namespace boost::accumulators;

namespace stats_utils {
    typedef std::pair<unsigned long long, double> corr_point_t;
    typedef std::vector< corr_point_t > corr_point_list_t;
    
    template <class T>
    class Correlations{
    public:
        Correlations(std::vector<T>& data):_avgandvarcomputed(false){
            _data_sets.push_back(&data);
            _autocorrelations.resize(_data_sets.size());
            _crosscorrelations.resize(_data_sets.size());
        }
        
        void addData(std::vector<T>& data){
            _data_sets.push_back(&data);
        }
        
        void computeAuto();
        void computeAuto(unsigned int i);
        void computeCross(unsigned int i=0, unsigned int j=1);
        
        void print(){
            unsigned int k=0;
            for (auto& i: _autocorrelations)
                for (auto& j : i)
                    std::cout << k++ << "\t" << j << std::endl;
        }
        
    protected:
        bool _avgandvarcomputed;
        
        std::vector< std::vector <T>* > _data_sets;
        std::vector<double> _average, _variance;
        
        std::vector< std::vector < double > > _autocorrelations;
        std::vector< std::vector < double > > _crosscorrelations;
        void _computeAvgandVar();
        
        
    };
}


#endif /* correlations_h */

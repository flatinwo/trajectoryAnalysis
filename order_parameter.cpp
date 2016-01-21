//
//  order_parameter.cpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 10/13/15.
//  Copyright © 2015 Folarin Latinwo. All rights reserved.
//

#include "order_parameter.hpp"
#include "io.hpp"
#include <cassert>
#include <algorithm>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#include <boost/bind.hpp>
#include <boost/ref.hpp>

//try order-N algorithm in Allen & Tildesey or Frenkel & Smit book for autocorrelation
//it may be worth it to data in consistent vector format, that is data[0] represents a given set
//data[1] represents a given set, and so on....
//an alternative would be to load using multidimensional library in boost
//use R as library to CPP package

using namespace boost::accumulators;

namespace trajectoryAnalysis {
    
    OrderParameter::OrderParameter(const char* filename){
        load(filename, _data);
        std::cout << "Size of data is\t" << _data.size() << "\t" << _data[0].size() << std::endl;
        _restructureData();
    }
    
    OrderParameter::OrderParameter(Trajectory& traj):_trajectory(&traj){
        
    }
    
    OrderParameter::~OrderParameter(){
    }
    
#pragma mark COMPUTES
    
    //brute force auto-correlation function
    void OrderParameter::computeAutoCorrelation(int index){
        assert(_data.size()>index);
        assert(_data[index].size()>100);
        
        int dataSize = (int) _data[index].size();
        
        corr_point_t zeros;
        zeros.first = zeros.second = 0.;
        corr_point_list_t correlation((int) floor((double) dataSize/3), zeros);
        
        //delete average from all data
        std::transform(_data[index].begin(), _data[index].end(), _data[index].begin(), std::bind2nd(std::minus<double>(), _average));

        for (unsigned int i=0; i<_data[index].size(); i++) {
            for (unsigned int j=i; j<_data[index].size(); j++) {
                if (j-i >= correlation.size())
                    continue;
                correlation[j-i].second += _data[index][i]*_data[index][j];
                correlation[j-i].first++;
            }
        }
        
        
        _correlation.clear();        
        
        for (auto it=correlation.begin(); it != correlation.end(); ++it) {
            double corr = (it->second/(double) it->first)/_variance;
            _correlation.push_back(corr);
        }
        
        
        assert(_correlation.size() == correlation.size());
    }
    
    void OrderParameter::computeAverageAndVariance(int j){
        accumulator_set<double, stats< tag::mean, tag::variance > > acc; //typede this
        assert(j<_data.size());
        
        
        std::for_each(_data[j].begin(), _data[j].end(), boost::bind<void>(boost::ref(acc), _1)); //not copying by value
        
        _average = mean(acc);
        _variance = boost::accumulators::variance(acc);
        
        std::cout << "Average is\t" << _average << std::endl;
        std::cout << "Variance is\t" << _variance << std::endl;

    }
    
    
    void OrderParameter::printCorrelation(){
        for (unsigned int i=0; i<_correlation.size(); i++)
            std::cout << i << "\t" << _correlation[i] << std::endl;
    }
    
    
    
    void OrderParameter::_restructureData(){
        coord_list_t new_data;
        
        for (unsigned int i=0; i<_data[0].size(); i++) {
            coord_t x;
            for (unsigned int j=0; j<_data.size(); j++) {
                x.push_back(_data[j][i]);
            }
            new_data.push_back(x);
        }
        
        _data = new_data;
        std::cout << _data.size() << "\t" << _data[0].size() << "\n";
    }
}

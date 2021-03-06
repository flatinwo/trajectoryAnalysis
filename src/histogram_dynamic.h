//
//  histogram_dynamic.h
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 2/2/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#ifndef histogram_dynamic_h
#define histogram_dynamic_h

#include <deque>
#include <iostream>

#include <cstdlib>

namespace stats_utils {
    template <class T>
    class HistogramDynamic{
    public:
        HistogramDynamic():_inc(0.01),_is_firsttime(true){
            
        }
        
        HistogramDynamic(T inc, bool to_normalize=false): _inc(inc),_is_firsttime(true),
        _to_normalize(to_normalize){
            
        }
        
        //copy constructor
        HistogramDynamic(const HistogramDynamic& rhs){
            _weight_of_bin = rhs._weight_of_bin;
            _min = rhs._min;
            _max = rhs._max;
            _inc = rhs._inc;
            _is_firsttime = rhs._is_firsttime;
            _to_normalize = rhs._to_normalize;
        }
        
        //get size of each bin
        T getBinSize() const{
            return _inc;
        }
        
        void combineWith(const HistogramDynamic& rhs){
            if (_inc != rhs._inc) {
                std::cerr << "Error:HistogramDynamic:\nAttempt to combine histograms with different bin size\n";
                exit(-1);
            }
            
            HistogramDynamic temp(rhs);
            if (_max > temp._max) temp.insert(_max - _inc/2.0,0);
            if (temp._max > _max) insert(temp._max - temp._inc/2.0,0);
            if (_min < temp._min) temp.insert(_min + _inc/2.0, 0);
            if (temp._min < _min) insert(temp._min + temp._inc/2.0,0);
            
            for (unsigned int i=0; i<temp._weight_of_bin.size(); i++)
                _weight_of_bin[i] += temp._weight_of_bin[i];
        }
        
        
        //set the bin size
        void setBinSize(T size){
            _inc = size;
            clear();
        }
        
        void setToNormalize(bool to_normalize){
            _to_normalize = to_normalize;
        }
        
        //indexing operator
        double & operator [] (int i){
            return _weight_of_bin[i];
        }
        
        //insert a data point into histogram
        void insert(T value){
            insert(value,1.0);
        }
        
        //insert data point into histogram with a given weight
        void insert(T value, double weight){
            if (_is_firsttime) {
                createFirstBin(value,weight);
                return;
            }
            
            int bin = computeBinContaining(value);
            if (bin == TOO_BIG) expandRight(value,weight);
            else if (bin == TOO_SMALL) expandLeft(value,weight);
            else _weight_of_bin[bin] += weight;
            
        }
        
        //clears the content of the histogram
        void clear(){
            _weight_of_bin.clear();
            _is_firsttime = true;
        }
        
        //clears the content of the histogram
        void reset(){
            clear();
        }
        
        //normalizes so that probabilities are normalized to one
        void normalize(){
            double sum = 0.;
            for (auto& i: _weight_of_bin) sum += i;
            for (auto& i :_weight_of_bin) i /= (sum * (double) _inc);
        }
        
        //prints histogram to standard out
        void print(){
            print(std::cout);
        }
        
        //prints histogram to stream
        void print(std::ostream &os){
            if (_to_normalize) normalize();
            for (unsigned int i=0; i<_weight_of_bin.size(); i++) {
                os << (T)(_min + ((_inc*i) + (_inc*(i+1)))*0.5) << "\t" << _weight_of_bin[i] << std::endl;
            }
        }
        
        unsigned int getNumberOfBins() const{
            return int(_weight_of_bin.size());
        }
        
        std::pair<T,T> getBin(unsigned int i){
            if (i >=_weight_of_bin.size()) {
                std::cerr << "Error:HistogramDynamic:\nAttempt to access out of bounds index\n";
                exit(-1);
            }
            std::pair<T,T> x;
            x.first = (T) (_min + ((_inc*i) + (_inc*(i+1)))*0.5);
            x.second = _weight_of_bin[i];
            return x;
        }
        
        std::pair<T,T> max(){
            std::pair<T,T> max;
            max.second = 1e-32;
            
            for (unsigned int i=0; i<_weight_of_bin.size(); i++) {
                if (_weight_of_bin[i] > max.second) {
                    max.first = (T) (_min + ((_inc*i) + (_inc*(i+1)))*0.5);
                    max.second = _weight_of_bin[i];
                }
            }
            return max;
        }

        int BinContaining(T value){
		return computeBinContaining(value);
        }
 

        
    protected:
        //expands the histogram in the positive x-direction
        void expandRight(T value, double weight){
            for (T i=_max; i<value; i+=_inc) {
                _max = i + _inc;
                _weight_of_bin.push_back(0);
            }
            _weight_of_bin[_weight_of_bin.size()-1] += weight;
        }
        
        //expands the histogram in the negative y-direction
        void expandLeft(T value, double weight){
            for (T i = _min; i > value; i-= _inc) {
                _min = i - _inc;
                _weight_of_bin.push_front(0);
            }
            _weight_of_bin[0] += weight;
        }
        
        //compute bin number containing a given data point
        int computeBinContaining(T value){
            if (value < _min) return TOO_SMALL;
            else if (value >= _max) return TOO_BIG;
            else return (int) ( ((double) value - _min)/(_max - _min)*_weight_of_bin.size() );
        }
        
        //creates the initial bin in the histogram
        void createFirstBin(T value, double weight){
            _weight_of_bin.push_back(weight);
            _min = value - _inc * 0.5;
            _max = value + _inc * 0.5;
            _is_firsttime = false;
        }
        
        
        std::deque<double> _weight_of_bin;
        
        T _min;
        T _max;
        T _inc;
        
        enum {TOO_SMALL = -2, TOO_BIG = -1};
        bool _is_firsttime;
        bool _to_normalize=false;
    };
    
    template <class TYPE>
    std::ostream& operator << (std::ostream& os, HistogramDynamic<TYPE>& h){
        h.print(os);
        return os;
    }
}


#endif /* histogram_dynamic_h */

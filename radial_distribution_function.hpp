//
//  radial_distribution_function.hpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 5/25/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#ifndef radial_distribution_function_hpp
#define radial_distribution_function_hpp

#include <stdio.h>
#include "struct_def.h"
#include "trajectory.hpp"
#include "histogram_dynamic.h"
#include "io.hpp"
#include <fstream>
#include <memory>
#include <set>

namespace trajectoryAnalysis{
    class RadialDistributionFunction{
    public:
        RadialDistributionFunction(trajectory_t&, double binsize=0.05);
        RadialDistributionFunction(xyztrajectory_t&, Box&,double binsize=0.05);
        
        void compute();
        void print();
        void update();
        
        void addSimpleExcludeRule(bool (*funcp)(unsigned, unsigned));   //usually for intermolecule
        void addSimpleIncludeRule(bool (*funcp)(unsigned, unsigned));
        
    protected:
        trajectory_t* _trajc;
        xyztrajectory_t* _trajx;
        Box* _box;
        
        double _bin_size;
        unsigned long _number_of_frames;
        double _rmax;
        double _vol;
        int _nunique_types;
        
        struct _normalize_info{
            _normalize_info():
            _same_type(false),
            _ntypei(0.),
            _ntypej(0.){
            }
            bool _same_type;
            double _ntypei;
            double _ntypej;
            std::string _typei;
            std::string _typej;
        };
        
        std::vector<function1d_t> _GofRs;
        std::vector< stats_utils::HistogramDynamic < double > > _Hists;
        std::vector< _normalize_info > _hist_infos;
        
        std::vector<bool (*)(unsigned, unsigned)> _exclude_funcptrs;
        std::vector<bool (*)(unsigned, unsigned)> _include_funcptrs;
        std::vector< std::shared_ptr<std::ofstream> > _ofiles;
        std::set<std::string> _distinct_container;
        
        void _openFiles();
        void _closeFiles();
        void _updateHist(int,double);
        
        void _initialize();
        void _normalize();
        
        void _analyzeFrames();
        
        
        function1d_t _normalize(const stats_utils::HistogramDynamic<double>&, _normalize_info&);
        
        
    };
}


#endif /* radial_distribution_function_hpp */

//
//  chirality_analysis.cpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 5/23/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#include "chirality_analysis.hpp"
#include "spatial.hpp"
#include <algorithm>
#include <limits>

namespace trajectoryAnalysis {
    ChiralityAnalysis::ChiralityAnalysis(int argc, const char* argv[]){
        
        _filein = argv[1];
        _fileout = argv[2];
        
        _box.box_hi[0] = atof(argv[3]);
        _box.box_hi[1] = atof(argv[4]);
        _box.box_hi[2] = atof(argv[5]);
        _box.updatePeriod();
        
        _initialize();
        
    }
    
    ChiralityAnalysis::ChiralityAnalysis(const char* filename, Box& box){
        _filein = filename;
        _fileout = "chiralized" + std::string(filename);
        _box = box;
        _box.updatePeriod();
        _initialize();
        
    }
    
    void ChiralityAnalysis::_initialize(){
        _averageC = 0.;
        loadxyz(_filein.c_str(), _traj);
        
        //type info
        _typecount.resize(3,0);
        _typecountmax = _typecount;
        
        _typematch.resize(3,"1");
        _typematch[1] = "2";
        _typematch[2] = "3";
    }
    
    void ChiralityAnalysis::analyze(){
        std::ofstream fileEE("postAnalysisEE.out");
        for (unsigned int i=0; i<_traj.size(); i++){
            _analyzeChiralityXYZ(_traj[i]);
            fileEE << i << "\t" << _averageC << std::endl;
            
            //allows for visualization with vmd by accounting for changing types on move)
            for (unsigned int j=0; j<3; j++){
                _typecount[j] = (int) std::count(_traj[i].type.begin(),_traj[i].type.end(),_typematch[j]);
                _typecountmax[j] = std::max(_typecount[j],_typecountmax[j]);
            }
        }
        for (unsigned int i=0; i<_typematch.size(); i++) _logtypes[_typematch[i]] = _typecountmax[i];
        fileEE.close();
    }
    
    void ChiralityAnalysis::visualize(){
        VisualizerXYZ chiral(_fileout.c_str());
        chiral.setTypeMax(_logtypes);
        chiral.visualize(_traj);
    }
    
    
    int ChiralityAnalysis::_computeMolecularCOM(xyzfile& snap, coord_t* nearest_x){
        int natoms = snap.n;
        assert(natoms%4 == 0);
        int nmolecules = natoms/4;
        unsigned int k,l,m;
        double xl,yl,zl;
        k=l=m=0;
        
        int nidx=-1;
        double min = std::numeric_limits<double>::max();
        
        coord_t x(3,0.);
        _molecular_com.resize(nmolecules,coord_t(3,0.));
        
        
        for (int i=0; i<nmolecules; i++){
            for (unsigned int j=0;j<4;j++){
                xl = snap.x[k][0]; yl = snap.x[k][1]; zl = snap.x[k][2];
                pbc(xl,_box.box_period[0]);
                pbc(yl,_box.box_period[1]);
                pbc(zl,_box.box_period[2]);
                
                x[0] += xl;
                x[1] += yl;
                x[2] += zl;
                k++;
            }
            for (unsigned int l=0;l<3;l++) {
                _molecular_com[i][l] = x[l]/4;
            }
            
            if (nearest_x != nullptr){
                double rsq = distancesq(_molecular_com[i], *nearest_x);
                if ( rsq < min) {
                    nidx = i;
                    min =rsq;
                }
            }
            
            x[0]=x[1]=x[2]=0;
            
        }
        assert(k = natoms+1);
        return nidx;
    }
    
    void ChiralityAnalysis::_analyzeChiralityXYZ(xyzfile& snap){
        int natoms = snap.n;
        double zeta_d;
        assert(natoms%4 == 0);
        int nmolecules = natoms/4;
        unsigned int k,l,m;
        k=l=m=0;
        _averageC = 0.;
        
        for (int i=0; i<nmolecules; i++) {
            coord_list_t x;
            for (unsigned int j=0; j<4; j++) {
                x.push_back(snap.x[k++]);
            }
            std::string zeta = _chiralityunwrap(x,zeta_d);
            _averageC += zeta_d;
            for (unsigned int j=0; j<4; j++) {
                /*snap.type[l++] = zeta;*/
                if (j==0 || j==3) snap.type[l++]=zeta+"1";
                else snap.type[l++]=zeta+"2";
                snap.x[m++] = x[j];
            }
        }
        _averageC /= (double) nmolecules;
    }
    
    std::string ChiralityAnalysis::_chiralityunwrap(coord_list_t& x, double& zetad){
        double_coord_t dr1 = distancesqandvec(x[1], x[0], _box);
        x[1] = x[0] + dr1.second;
        double_coord_t dr2 = distancesqandvec(x[2], x[1], _box);
        x[2] = x[1] + dr2.second;
        double_coord_t dr3 = distancesqandvec(x[3], x[2], _box);
        x[3] = x[2] + dr3.second;
        
        //compute zeta
        double c1 = (dr2.second[1]*dr3.second[2] - dr2.second[2]*dr3.second[1]);
        double c2 = (dr2.second[0]*dr3.second[2] - dr2.second[2]*dr3.second[0]);
        double c3 = (dr2.second[0]*dr3.second[1] - dr2.second[1]*dr3.second[0]);
        
        double zeta = dr1.second[0]*c1 - dr1.second[1]*c2 + dr1.second[2]*c3;
        zeta /= (sqrt(dr1.first*dr2.first*dr3.first));
        zetad = zeta;
        
        //std::cout << zeta <<"\n";
        
        if (zeta > 0.33) {
            return "1";
        }
        else if (zeta < -0.33){
            return "2";
        }
        else
            return "3";
    }
    
    
}


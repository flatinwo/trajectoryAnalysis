//
//  io.cpp
//  trajectoryAnalysis
//
//  Created by Folarin Latinwo on 10/12/15.
//  Copyright © 2015 Folarin Latinwo. All rights reserved.
//

#include "io.hpp"
#include "spatial.hpp"
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <stdexcept>
#include <cassert>
#include <map>
#include <assert.h>
#include <stack>

//implement strides for calculations

namespace trajectoryAnalysis {
    
    /**
     \brief Combines all trajectories in a filename
     */
    Trajectory getCombinedTrajectory(const char* filename){
        std::vector<std::string> filepaths;
        std::ifstream datafile(filename);
        
        if (datafile.is_open()) {
            std::cerr << filename << " is now open...\n";
            std::string str;
            while (datafile >> str) filepaths.push_back(str);
        }
        else{
            std::cerr << "Unable to open to file " << filename << std::endl;
            std::exit(-1);
        }
        
        assert(filepaths.size()>0);
        Trajectory traj(filepaths[0].c_str(),true,1,5);
        for (unsigned int i=1; i<filepaths.size(); i++) {
            Trajectory temp(filepaths[i].c_str(),true,1,5);
            traj = traj + temp;
        }
        
        return traj;
    }

    
    //Brief: Load raw coordinates
    void load(const char* filename, coord_list_t& x, arg_t arg){
        
        std::ifstream file(filename);
        if (file.fail()) {
            std::cerr << "Error: load: can't openfile "
            << filename << ".\n";
            exit(1);
        }
        
        std::string str;
        while (std::getline(file,str)){
            std::istringstream iss(str.c_str());
            double temp;
            coord_t xi;
            while (iss >> temp) {
                xi.push_back(temp);
            }
            if (xi.size() != 0) {
                x.push_back(xi);
            }
        }
    }
    
    
    // delete type from coordinate list
    void deltype(coord_list_t& x, std::vector<std::string>& types, std::string t){
        assert(x.size() == types.size());
        
        coord_list_t xtmp = x;
        std::vector<std::string> typetmp = types;
        
        types.clear();
        x.clear();
        
        for (unsigned int i=0; i<xtmp.size(); i++) {
            if (typetmp[i] != t) {
                x.push_back(xtmp[i]);
                types.push_back(typetmp[i]);
            }
        }
    }
    
    
    //converts file to ifstream and loads!
    
    void loadxyz(const char* filename, xyzfile& xyz){
        std::ifstream file(filename);
        assert(file.good());
        loadxyz(file,xyz);
        file.close();
    }
    
    void loadxyz(const char* filename, xyztrajectory_t& xyztraj){
        std::ifstream file(filename);
        assert(file.good());
        
        //clear vector
        xyztraj.clear();
        
        //load up xyztraj
        while (!file.eof()) {
            xyzfile xyz;
            loadxyz(file, xyz);
            xyztraj.push_back(xyz);
        }
        
        //typically, we have a dangling xyzfile
        xyztraj.pop_back();
        
        file.close();
        std::cout << "Total Number of Frames read is:\t" << xyztraj.size() << std::endl;
        
    }
    
    
    void loadxyz(std::istream& is, xyzfile& xyz){
        int N=0;
        is >> N;
        xyz.n = N;
        std::getline(is,xyz.commentstr); // will get rest of first line
        std::getline(is,xyz.commentstr);
        
        xyz.x.resize(N,std::vector<double>(3,0.0));
        xyz.x.resize(N);
        xyz.type.resize(N);
        
        for (int i=0; i<N; i++) {
            double x,y,z;
            std::string typestr;
            is >> typestr >> x >> y >> z;
            xyz.type[i] = typestr;
            xyz.x[i][0] = x;
            xyz.x[i][1] = y;
            xyz.x[i][2] = z;
        }
        
    }
    
    
    void loadxyz(const char* filename, coord_list_t& x, arg_t arg){
        xyz_info* info = (xyz_info*) arg;
        xyz_info temp_info;
        
        if (arg == 0x0){
            info = &temp_info;
        }
        
        std::istream* is;
        std::ifstream ifs;
        
        if (info->instream != 0x0){
            is = info->instream;
        }
        else{
            ifs.open(filename);
            is = &ifs;
        }
        
        assert(is->good());
        
        int n = 0;
        std::string str;
        std::getline(*is, str);
        std::istringstream iss(str.c_str());
        iss >> n;
        std::getline(*is,str);
        x.clear();
        info->type.clear();
        
        for (int i=0; i<n; i++) {
            std::string typei;
            double xi, yi, zi;
            *is >> typei >> xi >> yi >> zi;
            info->type.push_back(typei);
            coord_t xyz(3);
            xyz[0] = xi; xyz[1] = yi; xyz[2] = zi;
            x.push_back(xyz);
        }
        std::getline(*is,str);
        
        if (info->instream == 0x0) {
            ifs.close();
        }
        
    }
    
    void loadxyzfancy(const char* filename, trajectory_t& system, unsigned int index, unsigned int every){
        xyz_info* info = new xyz_info();
        
        std::istream* is;
        std::ifstream ifs;
        
        if (info->instream != 0x0){
            is = info->instream;
        }
        else{
            ifs.open(filename);
            is = &ifs;
        }
        
        assert(is->good());
        assert(index > 0);
        assert(every > 0);
        
        int nsnap = 0;
        std::string str, word;
        
        //set exceptions
        //is->exceptions(std::ifstream::failbit|std::ifstream::badbit);
        is->exceptions(std::ifstream::badbit);
        
        
        try{
            while (std::getline(*is, str)) {
                
                system.push_back(Snap());
                
                //box reference
                Box& box = system[nsnap].box;
                
                //get total number of particles
                int n = 0;
                std::istringstream iss(str.c_str());
                iss >> n;
                
                //std::cout << str << "\t" << n << std::endl;
                assert(n>0);
                
                //now attempt to get box size
                //a quick long winded approach and assumes that box_lo = zeroes
                std::getline(*is,str);
                std::stack<std::string> text;
                std::istringstream iss1(str);
                while (iss1.good()) {
                    iss1 >> word;
                    text.push(word);
                }
                
                //it is iterator
                for (auto it=box.box_hi.rbegin(); it != box.box_hi.rend(); ++it) {
                    *it = std::stod(text.top());
                    text.pop();
                }
                //update to compute period
                box.updatePeriod();
                text.pop(); text.pop();
                system[nsnap].timestep = std::stoi(text.top());
                
                coord_list_t x;
                info->type.clear();
                
                for (int i=0; i<n; i++) {
                    std::string typei;
                    double xi, yi, zi;
                    std::getline(*is,str);
                    if ((i%every) == index-1 ){
                        std::istringstream is1(str);
                        is1 >> typei >> xi >> yi >> zi;
                        info->type.push_back(typei);
                        coord_t xyz(3), rxyz(3);
                        xyz[0] = xi; xyz[1] = yi; xyz[2] = zi;
                        x.push_back(xyz);
                        system[nsnap]._center_of_mass_list.push_back(x[0]);
                        system[nsnap]._type_list.push_back(typei);
                        x.clear();
                    }
                    if (str.empty()){
                        system.pop_back();
                        delete info;
                        return;
                    }
                }
                //std::getline(*is,str);
                nsnap++;
            }
            
            if (info->instream == 0x0) {
                ifs.close();
            }
            
            delete info;

        }
        catch(std::ifstream::failure e){
            std::cerr << "Exception happened: " << e.what() << "\n"
            << "Error bits are: "
            << "\nfailbit: " << is->fail()
            << "\neofbit: " << is->eof()
            << "\nbadbit: " << is->bad()
            << "\nlikely provided incomplete file" << std::endl;
            system.pop_back();
        }
    }
    
    
    void loadgrofancy(const char* filename, trajectory_t& system, unsigned int index, unsigned int every){
        xyz_info* info = new xyz_info();
        
        std::istream* is;
        std::ifstream ifs;
        
        if (info->instream != 0x0){
            is = info->instream;
        }
        else{
            ifs.open(filename);
            is = &ifs;
        }
        
        assert(is->good());
        assert(index > 0);
        assert(every > 0);
        
        int nsnap = 0;
        std::string str, word;
        
        while (std::getline(*is, str)) {
            std::getline(*is,str); //get first line
            system.push_back(Snap());
            
            //box reference
            Box& box = system[nsnap].box;
            
            //get total number of particles
            int n = 0;
            std::istringstream iss(str.c_str());
            iss >> n;
            
            //std::cout << str << "\t" << n << std::endl;
            assert(n>0);
            
            box.updatePeriod();
            
            system[nsnap].timestep = nsnap;
            
            coord_list_t x;
            info->type.clear();
            
            double xi, yi, zi;
            for (int i=0; i<n; i++) {
                std::string typei,junk;
                std::getline(*is,str);
                if ((i%every) == index-1 ){
                    std::istringstream is1(str);
                    is1 >> junk >>  typei >> junk >> xi >> yi >> zi;
                    info->type.push_back(typei);
                    coord_t xyz(3), rxyz(3);
                    xyz[0] = xi*10; xyz[1] = yi*10; xyz[2] = zi*10;//assumes output in nm
                    x.push_back(xyz);
                    system[nsnap]._center_of_mass_list.push_back(x[0]);
                    system[nsnap]._type_list.push_back(typei);
                    x.clear();
                }
                
                if (str.empty()){
                    system.pop_back();
                    delete info;
                    return;
                }

            }
            std::getline(*is,str);
            if (str.empty()){
                system.pop_back();
                delete info;
                return;
            }

            std::istringstream is1(str);
            is1 >> box.box_hi[0] >> box.box_hi[1] >> box.box_hi[2];
            for (auto &i : box.box_hi) i*=10.;//assumes output is in nm
            box.updatePeriod();
            nsnap++;
        }
        
        if (info->instream == 0x0) {
            ifs.close();
        }
        
        delete info;
        
    }


    //to be deprecated, same as loadxyzfancy with specific parameters
    void loadxyz(const char* filename, trajectory_t& system, arg_t arg){
        xyz_info* info = (xyz_info*) arg;
        xyz_info temp_info;
        
        if (arg == 0x0){
            info = &temp_info;
        }
        
        std::istream* is;
        std::ifstream ifs;
        
        if (info->instream != 0x0){
            is = info->instream;
        }
        else{
            ifs.open(filename);
            is = &ifs;
        }
        
        assert(is->good());
        
        int nsnap = 0;
        std::string str, word;
        while (std::getline(*is, str)) {
            
            system.push_back(Snap());
            
            //box reference
            Box& box = system[nsnap].box;
            
            //get total number of particles
            int n = 0;
            std::istringstream iss(str.c_str());
            iss >> n;
            
            //std::cout << str << "\t" << n << std::endl;
            assert(n>0);
            
            //now attempt to get box size
            //a quick long winded approach and assumes that box_lo = zeroes
            std::getline(*is,str);
            std::stack<std::string> text;
            std::istringstream iss1(str);
            while (iss1.good()) {
                iss1 >> word;
                text.push(word);
            }
            
            //it is iterator
            for (auto it=box.box_hi.rbegin(); it != box.box_hi.rend(); ++it) {
                *it = std::stod(text.top());
                text.pop();
            }
            //update to compute period
            box.updatePeriod();
            
            
            text.pop(); text.pop();
            
            system[nsnap].timestep = std::stoi(text.top());
            
            
            coord_list_t x,r;
            info->type.clear();
            
            for (int i=0; i<n; i++) {
                std::string typei;
                double xi, yi, zi, rxi, ryi, rzi, wxi = 0.0;
                
                if (i%5 != 0) *is >> typei >> xi >> yi >> zi >> rxi >> ryi >> rzi;
                else *is >> typei >> xi >> yi >> zi >> wxi >> rxi >> ryi >> rzi;
                
                info->type.push_back(typei);
                coord_t xyz(3), rxyz(3);
                xyz[0] = xi; xyz[1] = yi; xyz[2] = zi;
                
                if (i%5 != 0){
                    rxyz[0] = rxi; rxyz[1] = ryi; rxyz[2] = rzi;
                }
                else{
                    rxyz.resize(4);
                    rxyz[0] = wxi; rxyz[1] = rxi; rxyz[2] = ryi; rxyz[3] = rzi;
                }
                
                x.push_back(xyz);
                
                r.push_back(rxyz);
                if (x.size() == 5) {
                    system[nsnap]._center_of_mass_list.push_back(x[0]);
                    x.clear();
                    r.clear();
                }
            }
            std::getline(*is,str);
            nsnap++;
        }
        
        if (info->instream == 0x0) {
            ifs.close();
        }
        
    }
    
    
    //saves data with variate types
    
    void savevarxyz(const char* filename, coord_list_t& x, xyz_info& info){
        assert(info.type.size()>0);
        assert(info.reservoir.size()>0);
        std::ofstream xyzfile(filename);
        
        std::ostream *os;
        
        if (info.outstream == NULL) {
            assert(xyzfile.good());
            os = &xyzfile;
        }
        else{
            assert(info.outstream->good());
            os = info.outstream;
        }
        
        int n=0;
        typedef std::map<std::string, int> reservoir_t;
        reservoir_t::iterator r;
        
        for (r = info.reservoir.begin(); r!=info.reservoir.end(); ++r){
            n += r->second;
        }
        
        *os << n << "\n\n";
        
        for (r=info.reservoir.begin(); r!=info.reservoir.end();++r){
            std::string type = r->first;
            int npad = r->second;
            for (unsigned int i=0; i<x.size(); i++) {
                if (info.type[i] == type) {
                    *os << type << "\t" << x[i][0] << "\t" << x[i][1] << "\t" << x[i][2] << "\n";
                    npad--;
                }
            }
            for (int i=0; i<npad; i++) {
                *os << type << "\t" << info.xres[0] <<"\t"<< info.xres[1] <<"\t"<< info.xres[2] <<"\n";
            }
        }
        xyzfile.close();
    }
    
    
    void savexyz(const char* filename, xyzfile& data){
        xyz_info info;
        info.type = data.type;
        savexyz(filename, data.x, info);
    }
    
    void savexyz(const char* filename, xyztrajectory_t& xyztraj){
        assert(xyztraj.size() > 0);
        
        std::ofstream xyzfile(filename);
        std::ostream *os;
        
        for (unsigned int i=0; i<xyztraj.size(); i++) {
            unsigned long n = xyztraj[i].x.size();
            xyz_info info;
            info.type = xyztraj[i].type;
            
            if (info.outstream == NULL){
                assert(xyzfile.good());
                os = &xyzfile;
            }
            else{
                assert(info.outstream->good());
                os = info.outstream;
            }
            *os << n << "\n\n";
            if (info.type.size() > 0) {
                for (unsigned int j=0; j<n; j++){
                    *os << info.type[j] << "\t" << xyztraj[i].x[j] << "\n";
                }
            }
            else{
                for (unsigned int j=0; j<n; j++) *os << "H\t" << xyztraj[i].x[j] << "\n";
                
            }
        }
        
        xyzfile.close();
        std::cout << "Total Number of Frames written is:\t" << xyztraj.size() << std::endl;
        
    }

    
    
    //have a different version with xyz type
    void savexyz(const char* filename, coord_list_t& x, xyz_info& info){
        unsigned long n = x.size();
        
        std::ofstream xyzfile(filename);
        
        std::ostream *os;
        
        if (info.outstream == NULL){
            assert(xyzfile.good());
            os = &xyzfile;
        }
        else{
            assert(info.outstream->good());
            os = info.outstream;
        }
        
        *os << n << "\n\n";
        if (info.type.size() > 0) {
            for (unsigned int i=0; i<n; i++){
                *os << info.type[i] << "\t" << x[i][0] << "\t" << x[i][1] << "\t" << x[i][2] << "\n";
            }
        }
        else{
            for (unsigned int i=0; i<n; i++){
                *os << "H\t" << x[i][0] << "\t" << x[i][1] << "\t" << x[i][2] << "\n";
            }
        }
        xyzfile.close();
    }
    
    std::ostream& operator << (std::ostream& os, Box& box){
        os << box.box_lo << "\n" << box.box_hi << "\n" << box.box_period << "\n";
        return os;
    }
    
    std::ostream& operator << (std::ostream& os, coord_list_t& x){
        for (unsigned int i=0; i<x.size(); i++){
            os << x[i] << "\n";
        }
        return os;
    }
    
    std::ostream& operator << (std::ostream& os, coord_t& x){
        unsigned long sm1 = x.size()-1;
        for (unsigned long k=0; k<x.size();k++){
            os << x[k];
            if (k!=sm1) {
                os << "\t";
            }
        }
        return os;
    }
    
    std::ostream& operator << (std::ostream& os, shpdesc_t& sd){
        unsigned long sm1 = sd.size() - 1;
        for (unsigned long k=0; k<sd.size(); k++) {
            os << sd[k];
            if (k!=sm1) {
                os << "\t";
            }
        }
        return os;
    }
    
    std::ostream& operator << (std::ostream& os, function1d_t& f){
        for (function1d_t::iterator i=f.begin(); i!=f.end(); ++i){
            os << i->first << "\t" << i->second << "\n";
        }
        return os;
    }
    
    std::istream& operator >> (std::istream& is, coord_list_t& x){
        return is;
    }
    
    std::istream& operator >> (std::istream& is, coord_t& x){
        return is;
    }
    
    std::istream& operator >> (std::istream& is, shpdesc_t& sd){
        return is;
    }
}



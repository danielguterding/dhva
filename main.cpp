/*
* Copyright (c) 2013, Daniel Guterding <guterding@itp.uni-frankfurt.de>
*
* This file is part of dhva.
*
* dhva is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* dhva is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with dhva. If not, see <http://www.gnu.org/licenses/>.
*/

//main.cpp
#include <iostream>
#include <string>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

#include "files.hpp"
#include "settings.hpp"
#include "ruc.hpp"
#include "sc.hpp"
#include "orbit.hpp"
#include "eval.hpp"
using namespace std;

int main(int argc, char* argv[]){
  
  GlobalSettings settings;
  boost::filesystem::path filepath;
  
  if(argc != 12){
    cout << "Using precompiled settings." << endl;
    filepath = "sphere.bxsf";
    settings.inputinev = 0;
    settings.nksc = 60;
    settings.nsc = 4.0;
    settings.phi = 90.0/180*M_PI;
    settings.theta = 180.0/180*M_PI;
    settings.maxkdiff = 0.15;
    settings.maxfreqdiff = 0.10;
    settings.minimumfreq = 50;
    settings.ip = 0;
    settings.go = 0;
  }
  else {
    cout << "Using command line settings." << endl;
    filepath = argv[1];
    settings.inputinev = atoi(argv[2]);
    settings.nksc = atoi(argv[3]);
    settings.nsc = atoi(argv[4]);
    settings.phi = atof(argv[5])/180*M_PI;
    settings.theta = atof(argv[6])/180*M_PI;
    settings.maxkdiff = atof(argv[7]);
    settings.maxfreqdiff = atof(argv[8]);
    settings.minimumfreq = atof(argv[9]);
    settings.ip = atoi(argv[10]);
    settings.go = atoi(argv[11]);
  }
  
  string datadirstr = "data/";
  boost::filesystem::path datadir(datadirstr);
  mkdir(datadir);
  
  if(!boost::filesystem::exists(filepath)){
    printf("Error. Input file does not exist.\n");
  }
  else{
    cout << "Started reading input file." << endl;
    bxsf file(filepath, settings.inputinev);
    cout << "Finished reading input file." << endl;
    
    cout << "Started reconstruction of reciprocal unit cell." << endl;
    ReciprocalUnitCell ruc(file.get_nkpoints(), file.get_h(), file.get_energies());
    cout << "Finished reconstruction of reciprocal unit cell." << endl;
    
    cout << "Started populating super cell." << endl;
    SuperCell sc(settings, ruc);
    cout << "Finished populating super cell." << endl;
    
    cout << "Started orbit detection." << endl;
    OrbitFinder orbit(settings, sc);
    cout << "Finished orbit detection." << endl;
    
    cout << "Started evaluating orbits." << endl;
    OrbitEvaluator eval(orbit.get_orbits_pointer(), sc.get_sc_length(), settings.nsc);
    cout << "Finished evaluating orbits." << endl;
    
    cout << "Started matching fermi surface sheets." << endl;
    SheetMatcher match(eval.get_evaluated_orbits());
    cout << "Finished matching fermi surface sheets." << endl;
    
    cout << "Started singling out extremal frequencies." << endl;
    FrequencyCalculator freqcalc(settings, match.get_sheets(), ruc.get_h());
    cout << "Finished singling out extremal frequencies." << endl;
    
    cout << "Starting to write output file." << endl;
    string filenamestr = boost::lexical_cast<string>(filepath.filename());
    filenamestr.erase(0, 1);
    filenamestr.erase(filenamestr.size()-1);
    boost::filesystem::path outfilepath = datadirstr + boost::lexical_cast<string>(
					      boost::format("%s.%i_%i_%3.1f_%3.1f_%1.3f_%1.3f_%i_%i.out") 
					      % filenamestr % settings.nksc % settings.nsc 
					      % atof(argv[3]) % atof(argv[4]) % settings.maxkdiff 
					      % settings.maxfreqdiff % settings.minimumfreq % settings.ip);
    write_output(settings, outfilepath, freqcalc.get_properties());
    cout << "Finished writing output file." << endl;
    
    if(settings.go == 1){
      cout << "Started writing graphical output." << endl;
      
      boost::multi_array<fptype,3> energies;
      energies.resize(boost::extents[settings.nksc][settings.nksc][settings.nksc]);
      energies = sc.get_energies();
    
      for(int k=1;k<settings.nksc-1;k++){
        boost::filesystem::path outfilepath3(boost::lexical_cast<string>(boost::format("data/graphical%03i.txt") % k));
        boost::filesystem::ofstream outfilehandle3(outfilepath3);
   
        for(int i=settings.nksc-1;i>-1;i--){
          string line = boost::lexical_cast<string>(boost::format("%3i ") % i);
          for(int j=0;j<settings.nksc;j++){
	     if(energies[i][j][k] > 0){
	       line += "1";
	     }
	     else{
	       line += "0";
	     }
          }
          outfilehandle3 << line << endl;
        }
        outfilehandle3.close();
      }
      cout << "Finished writing graphical output." << endl;
    }
    cout << "Program finished." << endl;
  }
  
  return 0;
}

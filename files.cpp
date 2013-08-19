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

//files.cpp
#include "files.hpp"

// Functions for bxsf class

bxsf::bxsf(boost::filesystem::path path){
  
  h.resize(boost::extents[3][3]); 
  filepath = path;
  read();
}

void bxsf::read(){
  
  int BANDNUMBER_PASSED = 0;
  int END_BANDGRID_3D_PASSED = 0;
  
  boost::filesystem::ifstream filehandle(filepath);
  string line;
  while(getline(filehandle, line)){
    if(line.find("Fermi Energy:") != string::npos){
      vector<string> linestr;
      line = trim_all(line); //erase trailing and leading spaces, reduce intermediate spaces to one space
      boost::split(linestr, line, boost::is_any_of("\t "));
      fermi = boost::lexical_cast<fptype>(linestr[2]);
    }
    
    if(line.find("BANDGRID_3D_BANDS") != string::npos){
      vector<string> linestr;
      getline(filehandle, line); //skip useless line
      getline(filehandle, line); //this line contains nkpoints
      line = trim_all(line);
      boost::split(linestr, line, boost::is_any_of("\t "));
      for(int i=0;i<3;i++){
	nkpoints[i] = boost::lexical_cast<int>(linestr[i]);
      }
      getline(filehandle, line); //grid origin line
      for(int j=0;j<3;j++){
        getline(filehandle, line); //lines of the reciprocal unit cell matrix
        line = trim_all(line);
        boost::split(linestr, line, boost::is_any_of("\t "));
        for(int i=0;i<3;i++){
	  h[j][i] = 2*M_PI*INVBOHR2INVANGSTROM*boost::lexical_cast<fptype>(linestr[i]); //reciprocal lattice vectors are given in bohr^-1 and without 2PI prefactor
        }
      }
    }
    
    if(line.find("BAND:") != string::npos){
      vector<string> linestr;
      line = trim_all(line);
      boost::split(linestr, line, boost::is_any_of("\t "));
      bandnumber = boost::lexical_cast<int>(linestr[1]);
      BANDNUMBER_PASSED = 1;
    }
    
    if(line.find("END_BANDGRID_3D") != string::npos){
      END_BANDGRID_3D_PASSED = 1;
    }
    
    if((1 == BANDNUMBER_PASSED) && (0 == END_BANDGRID_3D_PASSED) && (line.find("BAND:") == string::npos)){
      vector<string> linestr;
      line = trim_all(line);
      boost::split(linestr, line, boost::is_any_of("\t "));
      for(int i=0;i<int(linestr.size());i++){
	energies_list.push_back(boost::lexical_cast<fptype>(linestr[i]) - fermi); //fermi energy is substracted
      }
    }
  }
  filehandle.close();
  
  energies.resize(boost::extents[nkpoints[0]][nkpoints[1]][nkpoints[2]]); //holds energies on the k-point grid
  
  if((nkpoints[0] * nkpoints[1] * nkpoints[2]) == energies_list.size()){
    for(int i=0;i<nkpoints[0];i++){
      for(int j=0;j<nkpoints[1];j++){
        for(int k=0;k<nkpoints[2];k++){
	  energies[i][j][k] = RYDBERG2EV*energies_list[i*nkpoints[1]*nkpoints[2] + j*nkpoints[2] + k];
        }
      }
    }
  }
  else{
    printf("Error: Number of energy values does not match gridsize.\n");
  }
}

boost::array<int, 3> bxsf::get_nkpoints(){
  
  return nkpoints;
}

boost::multi_array<fptype, 2> bxsf::get_h(){
  
  return h;
}

int bxsf::get_bandnumber(){
  
  return bandnumber;
}

vector<fptype> bxsf::get_energies_list(){
  
  return energies_list;
}

boost::multi_array<fptype, 3> bxsf::get_energies(){
  
  return energies;
}

// Main Functions

string trim_all(const std::string &str){  //with a more recent version of boost boost::trim_all() can be used instead of this function

  return boost::algorithm::find_format_all_copy(
    boost::trim_copy(str),
    boost::algorithm::token_finder (boost::is_space(),boost::algorithm::token_compress_on),
    boost::algorithm::const_formatter(" "));
}

void mkdir(boost::filesystem::path dir){
  
  if(!(boost::filesystem::exists(dir))){
    boost::filesystem::create_directory(dir);
  }
}

void write_output(GlobalSettings settings, boost::filesystem::path outfilepath, vector<AveragedOrbit> ao){
  
  boost::filesystem::ofstream outfilehandle(outfilepath);
  int naverages = ao.size();
  
  outfilehandle << boost::lexical_cast<string>(boost::format("# nksc : %i") % settings.nksc) << endl;
  outfilehandle << boost::lexical_cast<string>(boost::format("# phi : %f") % (settings.phi*180.0/M_PI)) << endl;
  outfilehandle << boost::lexical_cast<string>(boost::format("# theta : %f") % (settings.theta*180.0/M_PI)) << endl;
  outfilehandle << boost::lexical_cast<string>(boost::format("# maxkdiff : %f") % settings.maxkdiff) << endl;
  outfilehandle << boost::lexical_cast<string>(boost::format("# maxfreqdiff : %f") % settings.maxfreqdiff) << endl;
  outfilehandle << boost::lexical_cast<string>(boost::format("# minimumfreq : %f") % settings.minimumfreq) << endl;
  outfilehandle << boost::lexical_cast<string>(boost::format("# interpolation algo : %i") % settings.ip) << endl << endl;
  
  outfilehandle << "#Freq [T], SdevFreq [T], M [m_e], SdevM [m_e], X [0...1], SdevX, Y, SdevY, Z, SdevZ, Number of copies" << endl;
  for(int i=0;i<naverages;i++){
    outfilehandle << boost::lexical_cast<string>(boost::format("%5.1f %5.2f %f %f %f %f %f %f %f %f %i") 
                     % ao[i].f % ao[i].fsdev % ao[i].m % ao[i].msdev % ao[i].x % ao[i].xsdev % ao[i].y 
                     % ao[i].ysdev % ao[i].z % ao[i].zsdev % ao[i].n) << endl;
  }
  outfilehandle.close();
}
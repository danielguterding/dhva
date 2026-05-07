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
#include <stdexcept>

// Functions for bxsf class

BXSFKVectorUnit parse_kvector_unit(string unitstr){

  unitstr = boost::to_lower_copy(unitstr);
  if((unitstr == "legacy") || (unitstr == "bohr-no-2pi") || (unitstr == "bohr_no_2pi") || (unitstr == "crystal-bohr") || (unitstr == "crystal_bohr")){
    return BXSF_KVECTOR_LEGACY;
  }
  if((unitstr == "angstrom") || (unitstr == "angstrom^-1") || (unitstr == "a^-1") || (unitstr == "invang") || (unitstr == "vaspkit")){
    return BXSF_KVECTOR_ANGSTROM;
  }
  if((unitstr == "bohr") || (unitstr == "bohr^-1") || (unitstr == "invbohr")){
    return BXSF_KVECTOR_BOHR;
  }
  if((unitstr == "angstrom-no-2pi") || (unitstr == "angstrom_no_2pi") || (unitstr == "crystal-angstrom") || (unitstr == "crystal_angstrom")){
    return BXSF_KVECTOR_ANGSTROM_NO_2PI;
  }

  throw runtime_error(string("Unknown reciprocal-vector unit '") + unitstr + "'. Use legacy, angstrom, bohr, or angstrom-no-2pi.");
}

string kvector_unit_name(BXSFKVectorUnit unit){

  if(unit == BXSF_KVECTOR_LEGACY){
    return "legacy";
  }
  if(unit == BXSF_KVECTOR_ANGSTROM){
    return "angstrom";
  }
  if(unit == BXSF_KVECTOR_BOHR){
    return "bohr";
  }
  if(unit == BXSF_KVECTOR_ANGSTROM_NO_2PI){
    return "angstrom-no-2pi";
  }

  return "unknown";
}

bxsf::bxsf(boost::filesystem::path path, int inputinev_in, int requested_bandnumber_in, BXSFKVectorUnit kvector_unit_in){
  
  inputinev = inputinev_in;
  requested_bandnumber = requested_bandnumber_in;
  kvector_unit = kvector_unit_in;
  bandnumber = -1;
  fermi = 0;
  for(int i=0;i<3;i++){
    nkpoints[i] = 0;
  }
  h.resize(boost::extents[3][3]); 
  filepath = path;
  read();
}

fptype bxsf::kvector_unit_conversion_factor(){

  if(kvector_unit == BXSF_KVECTOR_LEGACY){
    return 2*M_PI*INVBOHR2INVANGSTROM;
  }
  if(kvector_unit == BXSF_KVECTOR_ANGSTROM){
    return 1.0;
  }
  if(kvector_unit == BXSF_KVECTOR_BOHR){
    return INVBOHR2INVANGSTROM;
  }
  if(kvector_unit == BXSF_KVECTOR_ANGSTROM_NO_2PI){
    return 2*M_PI;
  }

  throw runtime_error("Unknown reciprocal-vector unit.");
}

void bxsf::read(){
  
  bool grid_found = false;
  bool fermi_found = false;
  bool in_band = false;
  int current_bandnumber = -1;
  vector<fptype> current_band;
  vector<int> bandnumbers;
  vector<vector<fptype> > bands;

  boost::filesystem::ifstream filehandle(filepath);
  if(!filehandle){
    throw runtime_error("Could not open input file.");
  }

  string line;
  while(getline(filehandle, line)){
    if(line.find("Fermi Energy:") != string::npos){
      vector<string> linestr;
      line = trim_all(line); //erase trailing and leading spaces, reduce intermediate spaces to one space
      boost::split(linestr, line, boost::is_any_of("\t "));
      fermi = boost::lexical_cast<fptype>(linestr[2]);
      fermi_found = true;
    }
    
    if((line.find("BANDGRID_3D_BANDS") != string::npos) || (line.find("BEGIN_BANDGRID_3D") != string::npos)){
      vector<string> linestr;
      grid_found = true;
      getline(filehandle, line); //number of bands
      getline(filehandle, line); //this line contains nkpoints
      line = trim_all(line);
      boost::split(linestr, line, boost::is_any_of("\t "));
      for(int i=0;i<3;i++){
	nkpoints[i] = boost::lexical_cast<int>(linestr[i]);
      }
      getline(filehandle, line); //grid origin line
      fptype kvector_factor = kvector_unit_conversion_factor();
      for(int j=0;j<3;j++){
        getline(filehandle, line); //lines of the reciprocal unit cell matrix
        line = trim_all(line);
        boost::split(linestr, line, boost::is_any_of("\t "));
        for(int i=0;i<3;i++){
          h[j][i] = kvector_factor*boost::lexical_cast<fptype>(linestr[i]);
        }
      }
      continue;
    }
    
    if(line.find("BAND:") != string::npos){
      if(in_band){
        bandnumbers.push_back(current_bandnumber);
        bands.push_back(current_band);
        current_band.clear();
      }
      vector<string> linestr;
      line = trim_all(line);
      boost::split(linestr, line, boost::is_any_of("\t "));
      current_bandnumber = boost::lexical_cast<int>(linestr[1]);
      in_band = true;
      continue;
    }
    
    if(line.find("END_BANDGRID_3D") != string::npos){
      if(in_band){
        bandnumbers.push_back(current_bandnumber);
        bands.push_back(current_band);
        current_band.clear();
        in_band = false;
      }
      break;
    }
    
    if(in_band){
      vector<string> linestr;
      line = trim_all(line);
      if(line.empty()){
        continue;
      }
      boost::split(linestr, line, boost::is_any_of("\t "));
      for(int i=0;i<int(linestr.size());i++){
	current_band.push_back(boost::lexical_cast<fptype>(linestr[i]) - fermi); //fermi energy is substracted
      }
    }
  }
  filehandle.close();
  
  if(!fermi_found){
    throw runtime_error("Fermi energy not found in input file.");
  }
  if(!grid_found){
    throw runtime_error("BXSF band grid not found. Expected BANDGRID_3D_BANDS or BEGIN_BANDGRID_3D.");
  }
  if(bands.empty()){
    throw runtime_error("No BAND blocks found in input file.");
  }

  size_t expected_size = size_t(nkpoints[0]) * size_t(nkpoints[1]) * size_t(nkpoints[2]);
  for(int i=0;i<int(bands.size());i++){
    if(bands[i].size() != expected_size){
      throw runtime_error(string("Number of energy values in band ") + boost::lexical_cast<string>(bandnumbers[i]) + " does not match gridsize.");
    }
  }

  int selected_band = -1;
  if(requested_bandnumber > 0){
    for(int i=0;i<int(bandnumbers.size());i++){
      if(bandnumbers[i] == requested_bandnumber){
        selected_band = i;
      }
    }
    if(selected_band < 0){
      throw runtime_error(string("Requested band ") + boost::lexical_cast<string>(requested_bandnumber) + " not found in input file.");
    }
  }
  else if(bands.size() == 1){
    selected_band = 0;
  }
  else{
    vector<int> crossing_bands;
    for(int i=0;i<int(bands.size());i++){
      fptype min_energy = bands[i][0];
      fptype max_energy = bands[i][0];
      for(int j=1;j<int(bands[i].size());j++){
        if(bands[i][j] < min_energy){
          min_energy = bands[i][j];
        }
        if(bands[i][j] > max_energy){
          max_energy = bands[i][j];
        }
      }
      if((min_energy <= 0) && (max_energy >= 0)){
        crossing_bands.push_back(i);
      }
    }

    if(crossing_bands.size() == 1){
      selected_band = crossing_bands[0];
    }
    else if(crossing_bands.empty()){
      throw runtime_error("No band crosses the Fermi energy. Supply an explicit band number as the final command line argument if this is intended.");
    }
    else{
      string msg = "Multiple bands cross the Fermi energy:";
      for(int i=0;i<int(crossing_bands.size());i++){
        msg += " ";
        msg += boost::lexical_cast<string>(bandnumbers[crossing_bands[i]]);
      }
      msg += ". Supply one of these band numbers as the final command line argument.";
      throw runtime_error(msg);
    }
  }

  bandnumber = bandnumbers[selected_band];
  energies_list = bands[selected_band];
  energies.resize(boost::extents[nkpoints[0]][nkpoints[1]][nkpoints[2]]); //holds energies on the k-point grid
  for(int i=0;i<nkpoints[0];i++){
    for(int j=0;j<nkpoints[1];j++){
      for(int k=0;k<nkpoints[2];k++){
	energies[i][j][k] = energies_list[i*nkpoints[1]*nkpoints[2] + j*nkpoints[2] + k];
	if(!inputinev){
	  energies[i][j][k] *= RYDBERG2EV;
	}
      }
    }
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

void write_output(GlobalSettings settings, boost::filesystem::path outfilepath, vector<AveragedOrbit> ao, int bandnumber, BXSFKVectorUnit kvector_unit){
  
  boost::filesystem::ofstream outfilehandle(outfilepath);
  int naverages = ao.size();
  
  outfilehandle << boost::lexical_cast<string>(boost::format("# nksc : %i") % settings.nksc) << endl;
  if(bandnumber > 0){
    outfilehandle << boost::lexical_cast<string>(boost::format("# band : %i") % bandnumber) << endl;
  }
  outfilehandle << boost::lexical_cast<string>(boost::format("# kvector unit : %s") % kvector_unit_name(kvector_unit)) << endl;
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

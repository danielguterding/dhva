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

//orbit.cpp
#include "orbit.hpp"

OrbitFinder::OrbitFinder(GlobalSettings& settings, SuperCell& sc){
  
  nksc = settings.nksc;
  energies.resize(boost::extents[nksc][nksc][nksc]);
  kpoints.resize(boost::extents[nksc][nksc][nksc][3]);
  unchecked.resize(boost::extents[nksc][nksc][nksc]);
  orbitcont.set_slicecount(nksc);
  
  energies = *sc.get_energies_pointer();
  kpoints = *sc.get_kpoints_pointer();
  
  for(int i=0;i<nksc;i++){
    for(int j=0;j<nksc;j++){
      for(int k=0;k<nksc;k++){
	unchecked[i][j][k] = true;
      }
    }
  }
  
  start();
  orbitcont.delete_empty_and_open_orbits();
}

void OrbitFinder::start(){
  
  for(int k=1;k<nksc;k++){
    int i = 1, j = 1;
    while(i < nksc - 1){ //do if we are not finished with this slice
      while(j < nksc - 1){ //do if we are not at the end of a row
	if(unchecked[i][j][k]){
	  unchecked[i][j][k] = false;
	  if(energies[i][j][k] <= 0){
	    stepper(i, j, k);
	  }
	  else{
	    j++; //step right
	  }
	}
	else{
	  j++; //step right
	}
      } //end inner while
      i++; //go to start of next row
      j=1;
    } //end outer while
  } //end for
}

void OrbitFinder::stepper(const int i_in, const int j_in, const int k_in){
  
  //cout << "Stepper started." << endl;
  i = i_in;
  j = j_in;
  k = k_in;
  orbitcont.new_orbit(k);
  i_bef = i;
  j_bef = j-1;
  glance_east(); //first glance east
  stepper_done = false;
  for(int m=0;m<6;m++){
    ip[m] = -m-1;
    jp[m] = -m-1;
  }
  
  while(!stepper_done){
    set_checked();
    //cout << "I am at " << i << " " << j << " " << k << endl;
    if(glanced_outside_fs()){
      record_fs();
      if(orbit_closed()){
	stepper_done = true;
      }
      else{
	glance1();
	//cout << "glance1" << endl;
      }
    }
    else{
      step_to_glanced_point();
      if(point_on_sc_border() || circle_detected()){
	stepper_done = true;
      }
      else{
	glance2();
	//cout << "glance2" << endl;
      }  
    }
  }
}

void OrbitFinder::glance1(){
  
  if(!point_on_sc_border()){
    dir = new_glance_direction_for_glance1();
  
    switch(dir){
      case 0: glance_north(); break;
      case 1: glance_east(); break;
      case 2: glance_south(); break;
      case 3: glance_west(); break;
    } 
  }
  else{
    stepper_done = true;
  }
}

void OrbitFinder::glance2(){

  dir = new_glance_direction_for_glance2();
  switch(dir){
    case 0: glance_north(); break;
    case 1: glance_east(); break;
    case 2: glance_south(); break;
    case 3: glance_west(); break;
  } 
}

inline int OrbitFinder::new_glance_direction_for_glance1(){
  
  int dir = -1; //0==north, 1==east, 2==south, 3 == west
  if(north_of(i, j, ig, jg)){
    dir = 1;
  }
  else if(east_of(i, j, ig, jg)){
    dir = 2;
  }
  else if(south_of(i, j, ig, jg)){
    dir = 3;
  }
  else if(west_of(i, j, ig, jg)){
    dir = 0;
  }
  else{
    cout << "Error." << endl;
  }
  
  return dir;
}

inline int OrbitFinder::new_glance_direction_for_glance2(){
  
  int dir = -1; //0==north, 1==east, 2==south, 3 == west
  if(north_of(i, j, i_bef, j_bef)){
    dir = 1;
  }
  else if(east_of(i, j, i_bef, j_bef)){
    dir = 2;
  }
  else if(south_of(i, j, i_bef, j_bef)){
    dir = 3;
  }
  else if(west_of(i, j, i_bef, j_bef)){
    dir = 0;
  }
  else{
    cout << "Error." << endl;
  }
  
  return dir;
}

inline bool OrbitFinder::north_of(int i1, int j1, int i2, int j2){
  
  return ((i1 == (i2 - 1)) && (j1 == j2));
}

inline bool OrbitFinder::east_of(int i1, int j1, int i2, int j2){
  
  return ((i1 == i2) && (j1 == (j2 - 1)));
}

inline bool OrbitFinder::south_of(int i1, int j1, int i2, int j2){
  
  return ((i1 == (i2 + 1)) && (j1 == j2));
}

inline bool OrbitFinder::west_of(int i1, int j1, int i2, int j2){
  
  return ((i1 == i2) && (j1 == (j2 + 1)));
}

void OrbitFinder::record_fs(){
  
  OrbitPoint p;
  fptype x = kpoints[i][j][k][0]; //these are actual sc k-space coordinates
  fptype y = kpoints[i][j][k][1];
  fptype x_bef = kpoints[i_bef][j_bef][k][0];
  fptype y_bef = kpoints[i_bef][j_bef][k][1];
  fptype xg = kpoints[ig][jg][k][0];
  fptype yg = kpoints[ig][jg][k][1];
  fptype E = energies[i][j][k];
  fptype Eg = energies[ig][jg][k];
  fptype E_bef = energies[i_bef][j_bef][k];
  
  p.i = i;
  p.j = j;
  p.ig = ig;
  p.jg = jg;
  p.x = x - E/(Eg - E)*(xg - x);
  p.y = y - E/(Eg - E)*(yg - y);
  p.dEparallel = (Eg - E)/((xg-x) + (yg-y));
  p.dEperpendicular = 0;
  p.dir = dir;
  
  if(!orbitcont.last_orbit_empty(k)){
    OrbitPoint plast = orbitcont.get_last_orbitpoint(k);
    plast.dEperpendicular = (E - E_bef)/((x-x_bef) + (y-y_bef));
    orbitcont.update_last_orbitpoint(k, plast);
  }
  //do not exchange updating plast and adding p
  orbitcont.add_orbitpoint(k, p);
}

inline void OrbitFinder::glance_north(){
  
  ig = i + 1;
  jg = j;
}

inline void OrbitFinder::glance_east(){
  
  ig = i;
  jg = j + 1;
}

inline void OrbitFinder::glance_south(){
  
  ig = i - 1;
  jg = j;
}

inline void OrbitFinder::glance_west(){
  
  ig = i;
  jg = j - 1;
}

void OrbitFinder::set_checked(){
  
  unchecked[ig][jg][k] = false;
}

bool OrbitFinder::orbit_closed(){
  
  return orbitcont.simple_orbit_closed(k);
}

bool OrbitFinder::glanced_outside_fs(){
 
  return (energies[ig][jg][k] > 0);
}

void OrbitFinder::step_to_glanced_point(){
  
  update_stephistory();
  
  i_bef = i;
  j_bef = j;
  i = ig;
  j = jg;
}

bool OrbitFinder::point_on_sc_border(){
  
  bool isonborder = false;
  if((i>=(nksc - 1)) || (j >= (nksc - 1)) || (k >= (nksc - 1)) || (i<=0) || (j<=0) || (k<=0)){
    isonborder = true;
  }
  return isonborder;
}

bool OrbitFinder::circle_detected(){
  
  return (((ip[0] == ip[4]) && (jp[0] == jp[4])) && ((ip[1] == ip[5]) && (jp[1] == jp[5])));
}

void OrbitFinder::update_stephistory(){
  
  for(int i=1;i<6;i++){
    ip[i-1] = ip[i];
    jp[i-1] = jp[i];
  }
  ip[5] = ig;
  jp[5] = jg;
}

OrbitContainer OrbitFinder::get_orbits(){
  
  return orbitcont;
}

OrbitContainer* OrbitFinder::get_orbits_pointer(){
 
  return &orbitcont;
}

OrbitContainer::OrbitContainer(){
  
}

void OrbitContainer::set_slicecount(int n){
  
  nslices = n;
  orbitdata.resize(nslices);
}

void OrbitContainer::new_orbit(int sliceindex){
  
  orbitdata[sliceindex].resize(orbitdata[sliceindex].size()+1);
}

void OrbitContainer::add_orbitpoint(int sliceindex, OrbitPoint p){
  
  orbitdata[sliceindex][orbitdata[sliceindex].size()-1].push_back(p);
  //cout << "Adding " << p.x << " " << p.y << "." << endl;
}

OrbitPoint OrbitContainer::get_first_orbitpoint(int sliceindex){
  
  int orbitnumber = orbitdata[sliceindex].size() - 1;
  return orbitdata[sliceindex][orbitnumber][0];
}

OrbitPoint OrbitContainer::get_last_orbitpoint(int sliceindex){
  
  int orbitnumber = orbitdata[sliceindex].size()-1;
  int pointnumber = orbitdata[sliceindex][orbitnumber].size()-1;
  return orbitdata[sliceindex][orbitnumber][pointnumber];
}

bool OrbitContainer::orbit_closed(int sliceindex, int orbitindex){
  
  int npoints = orbitdata[sliceindex][orbitindex].size();
  bool orbit_empty = (npoints < 3); //we need at least three points to calculate an area
  //cout << npoints << " " << orbit_empty << endl;
  bool orbit_open = false;
  if(!orbit_empty){
    orbit_open = (orbitdata[sliceindex][orbitindex][0].x != orbitdata[sliceindex][orbitindex][npoints-1].x) || 
                 (orbitdata[sliceindex][orbitindex][0].y != orbitdata[sliceindex][orbitindex][npoints-1].y);
  }
  bool closed = (!orbit_empty) && (!orbit_open);
  return closed;
}

bool OrbitContainer::simple_orbit_closed(int sliceindex){
  
  int orbitnumber = orbitdata[sliceindex].size()-1;
  int pointnumber = orbitdata[sliceindex][orbitnumber].size()-1;
  bool closed = false;
  if((get_first_orbitpoint(sliceindex).x == get_last_orbitpoint(sliceindex).x) &&
    (get_first_orbitpoint(sliceindex).y == get_last_orbitpoint(sliceindex).y) && (pointnumber > 0)){ //only orbits containing more than one point can be closed
    closed = true;
  }
  return closed;
}

void OrbitContainer::delete_empty_and_open_orbits(){
  
  vector< vector< vector< OrbitPoint > > > goodorb; //use a new vector because erasing from a vector is extremely slow
  
  for(int i=0;i<nslices;i++){
    vector<vector<OrbitPoint> > slice;
    for(uint j=0;uint(j<orbitdata[i].size());j++){
      bool orbit_ok = orbit_closed(i, j);
      if(orbit_ok){
	slice.push_back(orbitdata[i][j]);
      }
    }
    goodorb.push_back(slice);
  }
  orbitdata = goodorb;
}

void OrbitContainer::print_orbits(){
  
  for(uint i=0;i<orbitdata.size();i++){
    for(uint j=0;j<orbitdata[i].size();j++){
      if(orbitdata[i][j].size() > 0){
	for(uint k=0;k<orbitdata[i][j].size();k++){
	  cout << orbitdata[i][j][k].x << " " << orbitdata[i][j][k].y << endl;
	}
      }
    }
  }
}

void OrbitContainer::write_orbits(boost::filesystem::path filepath){
  
  boost::filesystem::ofstream filehandle(filepath);
  for(uint i=0;i<orbitdata.size();i++){ //loop over slices
    for(uint j=0;j<orbitdata[i].size();j++){ //loop over orbits
      if(orbitdata[i][j].size() > 0){
	for(uint k=0;k<orbitdata[i][j].size();k++){ //loop over points in orbit
	  filehandle << boost::lexical_cast<string>(boost::format("% .8f % .8f %4i %4i %4i %4i") 
	                                            % orbitdata[i][j][k].x % orbitdata[i][j][k].y % 
	                                            orbitdata[i][j][k].i % orbitdata[i][j][k].j %
	                                            orbitdata[i][j][k].ig % orbitdata[i][j][k].jg) << endl;
	}
	filehandle << endl;
      }
    }
    filehandle << endl;
  }
  filehandle.close();
}

int OrbitContainer::get_slicecount(){
  
  return orbitdata.size();
}

int OrbitContainer::get_orbitcount(int sliceindex){
  
  return orbitdata[sliceindex].size();
}

vector<OrbitPoint> OrbitContainer::get_orbit(int sliceindex, int orbitindex){
  
  return orbitdata[sliceindex][orbitindex];
}

void OrbitContainer::update_last_orbitpoint(int sliceindex, OrbitPoint p){
  
  int orbitnumber = orbitdata[sliceindex].size()-1;
  orbitdata[sliceindex][orbitnumber].pop_back();
  orbitdata[sliceindex][orbitnumber].push_back(p);
}

bool OrbitContainer::last_orbit_empty(int sliceindex){
  
  int orbitnumber = orbitdata[sliceindex].size()-1;
  int pcount = orbitdata[sliceindex][orbitnumber].size();
  
  return (pcount < 1);
}
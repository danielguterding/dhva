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

//eval.cpp
#include "eval.hpp"

OrbitEvaluator::OrbitEvaluator(OrbitContainer* orbits, fptype sc_length_in, fptype nsc_in){
  
  sc_length = sc_length_in;
  nsc = nsc_in;
  con = *orbits;
  nslices = con.get_slicecount();
  for(int i=0;i<nslices;i++){
    norbits.push_back(con.get_orbitcount(i));
  }
  
  calc_all();
}

void OrbitEvaluator::calc_all(){
  
  for(int i=0;i<nslices;i++){
    vector<EvaluatedOrbit> orbitsperslice;
    for(int j=0;j<norbits[i];j++){
      EvaluatedOrbit orb;
      orb.f = calc_frequency(i, j);
      orb.m = calc_mass(i, j);
      orb.z = calc_z(i);
      orb.cx = calc_center_x(i, j);
      orb.cy = calc_center_y(i, j);
      orb.sdevx = calc_standarddev_x(i, j, orb.cx);
      orb.sdevy = calc_standarddev_y(i, j, orb.cy);
      orb.minx = calc_min_x(i, j);
      orb.maxx = calc_max_x(i, j);
      orb.miny = calc_min_y(i, j);
      orb.maxy = calc_max_y(i, j);
      
      //cout << orb.f << endl;

      orbitsperslice.push_back(orb);
    }
    res.push_back(orbitsperslice);
  }
}

fptype OrbitEvaluator::calc_frequency(int sliceindex, int orbitindex){
  
  vector<OrbitPoint> orbit = con.get_orbit(sliceindex, orbitindex);
  int npoints = orbit.size();
  
  fptype area = 0;
  for(int i=0;i<(npoints-1);i++){
    area += orbit[i].x * orbit[i+1].y - orbit[i+1].x * orbit[i].y;
  }
  area*=0.5;
  area = fabs(area);
  
  fptype f = HBAR/(2*M_PI*ELCHARGE)*1e20 * area;//*1e20 is necessary because reciprocal lattice vectors are given in 1/Angstrom
  return f;
}

fptype OrbitEvaluator::calc_mass(int sliceindex, int orbitindex){
  
  vector<OrbitPoint> orbit = con.get_orbit(sliceindex, orbitindex);
  int npoints = orbit.size();
  
  fptype effmass = 0; //effective mass in units of the free electron mass
  for(int i=0;i<(npoints-1);i++){
    fptype den = get_energy_derivative(orbit[i], orbit[i+1]);
    fptype enu = sqrt(pow(orbit[i+1].x - orbit[i].x,2) + pow(orbit[i+1].y - orbit[i].y,2));
    //cout << enu/den << endl;
    effmass += enu/den;
  }
  //cout << effmass << endl;
  effmass *= HBAR/ELMASS/ELCHARGE*HBAR/2.0/M_PI*1e20; //1e20 compensates for Angstrom^-2
  
  return effmass;
}

fptype OrbitEvaluator::calc_z(int sliceindex){
  
  fptype z = sc_length*float(sliceindex)/nslices - sc_length/nsc; 
  return z;
}

fptype OrbitEvaluator::get_energy_derivative(OrbitPoint p1, OrbitPoint p2){ //points i and i+1
  
  fptype derivative;
  
  if(glance_direction_parallel(p1.dir, p2.dir)){
    derivative = sqrt(pow(p1.dEparallel,2) + pow(p1.dEperpendicular,2));
  }
  else{
    derivative = sqrt(pow(p1.dEparallel,2) + pow(p2.dEparallel,2));
  }
  
  return derivative;
}

bool OrbitEvaluator::glance_direction_parallel(int d1, int d2){
  
  bool parallel = false;
  
  if((d1 == 1) && (d2 == 1)){
    parallel = true;
  }
  else if((d1 == 1) && (d2 == 3)){
    parallel = true;
  }
  else if((d1 == 3) && (d2 == 1)){
    parallel = true;
  }
  else if((d1 == 3) && (d2 == 3)){
    parallel = true;
  }
  else if((d1 == 0) && (d2 == 0)){
    parallel = true;
  }
  else if((d1 == 0) && (d2 == 2)){
    parallel = true;
  }
  else if((d1 == 2) && (d2 == 0)){
    parallel = true;
  }
  else if((d1 == 2) && (d2 == 2)){
    parallel = true;
  }
  
  return parallel;
}

fptype OrbitEvaluator::calc_center_x(int sliceindex, int orbitindex){
  
  vector<OrbitPoint> orbit = con.get_orbit(sliceindex, orbitindex);
  int npoints = orbit.size();

  fptype centerx=0;
  for(int i=0;i<npoints;i++){
    centerx+=orbit[i].x;
  }
  centerx/=float(npoints);
  
  return centerx;
}

fptype OrbitEvaluator::calc_center_y(int sliceindex, int orbitindex){
  
  vector<OrbitPoint> orbit = con.get_orbit(sliceindex, orbitindex);
  int npoints = orbit.size();

  fptype centery=0;
  for(int i=0;i<npoints;i++){
    centery+=orbit[i].y;
  }
  centery/=float(npoints);
  
  return centery;
}

fptype OrbitEvaluator::calc_standarddev_x(int sliceindex, int orbitindex, fptype cx){
  
  vector<OrbitPoint> orbit = con.get_orbit(sliceindex, orbitindex);
  int npoints = orbit.size();
  
  vector<fptype> xvals;
  for(int i=0;i<npoints;i++){
    xvals.push_back(orbit[i].x);
  }
  
  return standarddev(xvals, cx);
}

fptype OrbitEvaluator::calc_standarddev_y(int sliceindex, int orbitindex, fptype cy){
  
  vector<OrbitPoint> orbit = con.get_orbit(sliceindex, orbitindex);
  int npoints = orbit.size();
  
  vector<fptype> yvals;
  for(int i=0;i<npoints;i++){
    yvals.push_back(orbit[i].y);
  }
  
  return standarddev(yvals, cy);
}

fptype OrbitEvaluator::calc_min_x(int sliceindex, int orbitindex){
  
  vector<OrbitPoint> orbit = con.get_orbit(sliceindex, orbitindex);
  sort(orbit.begin(), orbit.end(), xcomp);
  
  return orbit[0].x;
}

fptype OrbitEvaluator::calc_max_x(int sliceindex, int orbitindex){
  
  vector<OrbitPoint> orbit = con.get_orbit(sliceindex, orbitindex);
  sort(orbit.begin(), orbit.end(), xcomp);
  int nvals = orbit.size();
  
  return orbit[nvals-1].x;
}

fptype OrbitEvaluator::calc_min_y(int sliceindex, int orbitindex){
  
  vector<OrbitPoint> orbit = con.get_orbit(sliceindex, orbitindex);
  sort(orbit.begin(), orbit.end(), ycomp);
  
  return orbit[0].y;
}

fptype OrbitEvaluator::calc_max_y(int sliceindex, int orbitindex){
  
  vector<OrbitPoint> orbit = con.get_orbit(sliceindex, orbitindex);
  sort(orbit.begin(), orbit.end(), ycomp);
  int nvals = orbit.size();
  
  return orbit[nvals-1].y;
}

fptype standarddev(vector<fptype> values, fptype average){
  
  int nvalues = values.size();
  fptype variance=0;
  
  for(int i=0;i<nvalues;i++){
    variance += pow(average - values[i],2);
  }
  
  return sqrt(variance);
}

bool xcomp(OrbitPoint p1, OrbitPoint p2){ 
  
  return (p1.x<p2.x); 
}

bool ycomp(OrbitPoint p1, OrbitPoint p2){ 
  
  return (p1.y<p2.y); 
}

bool fcompaveragedorbit(AveragedOrbit o1, AveragedOrbit o2){
  
  return (o1.f<o2.f);
}

vector<vector<EvaluatedOrbit> > OrbitEvaluator::get_evaluated_orbits(){
  
  return res;
}

SheetMatcher::SheetMatcher(vector<vector<EvaluatedOrbit> > orbits_in){
  
  eorbits = orbits_in;
  nslices = orbits_in.size();
  
  for(int i=0;i<nslices;i++){
    int no = orbits_in[i].size();
    norbits.push_back(no);
    vector<bool> temp;
    for(int j=0;j<no;j++){
      temp.push_back(false);
    }
    matched.push_back(temp);
  }

  for(int i=1;i<nslices-1;++i)
  {
    for(int j=0;j<norbits[i];++j)
    {
      if(!matched[i][j])
      {
        find_sheet(i,j);
      }
    }
  }
}

void SheetMatcher::find_sheet(int sliceindex, int orbitindex){
  
  vector<EvaluatedOrbit> sh;
  sh.push_back(eorbits[sliceindex][orbitindex]);
  
  for(int i=sliceindex+1;i<nslices;i++){
    vector<PossibleMatch> pm;
    for(int j=0;j<norbits[i];j++){
      if(simple_matching_condition_fulfilled(sliceindex, orbitindex, i, j) && (!(matched[i][j]))){
	PossibleMatch m;
	m.o = eorbits[i][j];
	m.i = i;
	m.j = j;
	pm.push_back(m);
      }
    }
    int nfoundorbits = pm.size();
    if(nfoundorbits == 1){
      sh.push_back(pm[0].o);
      matched[pm[0].i][pm[0].j] = true;
    }
    else if(nfoundorbits > 1){
      PossibleMatch bestmatch = get_best_match(eorbits[sliceindex][orbitindex], pm);
      sh.push_back(bestmatch.o);
      matched[bestmatch.i][bestmatch.j] = true;
    }
  }
  
  sheets.push_back(sh);
}

bool SheetMatcher::simple_matching_condition_fulfilled(int i1, int j1, int i2, int j2){ //orbit1 is already matched, orbit2 is candidate
  
  EvaluatedOrbit orb1 = eorbits[i1][j1];
  EvaluatedOrbit orb2 = eorbits[i2][j2];
  
  bool sdevcx = (((orb1.cx - orb1.sdevx) < orb2.cx) && ((orb1.cx + orb1.sdevx) > orb2.cx)); //check standard deviations of centers
  bool sdevcy = (((orb1.cy - orb1.sdevy) < orb2.cy) && ((orb1.cy + orb1.sdevy) > orb2.cy));
  bool sdevc = sdevcx && sdevcy;
  
  bool sdevmaxx = (((orb1.maxx - orb1.sdevx) < orb2.maxx) && ((orb1.maxx + orb1.sdevx) > orb2.maxx)); //check standard deviations of maximum coordinates
  bool sdevmaxy = (((orb1.maxy - orb1.sdevy) < orb2.maxy) && ((orb1.maxy + orb1.sdevy) > orb2.maxy));
  bool sdevmax = sdevmaxx && sdevmaxy;
  
  bool sdevminx = (((orb1.minx - orb1.sdevx) < orb2.minx) && ((orb1.minx + orb1.sdevx) > orb2.minx)); //check standard deviations of minimum coordinates
  bool sdevminy = (((orb1.miny - orb1.sdevy) < orb2.miny) && ((orb1.miny + orb1.sdevy) > orb2.miny));
  bool sdevmin = sdevminx && sdevminy;
  
  return (sdevc && sdevmax && sdevmin);
}

PossibleMatch SheetMatcher::get_best_match(EvaluatedOrbit orbit1, vector<PossibleMatch> pm){
  
  vector<PossibleMatch> npm;
  npm = calc_matching_parameter(orbit1, pm);
  
  sort(npm.begin(), npm.end(), Bcomp);
  
  return npm[0]; //return value with lowest B-value
}

vector<PossibleMatch> SheetMatcher::calc_matching_parameter(EvaluatedOrbit orbit1, vector<PossibleMatch> pm){
  
  vector<PossibleMatch> npm;
  
  int nmatches = pm.size();
  for(int i=0;i<nmatches;i++){
    fptype B = pow(orbit1.cx - pm[i].o.cx,2) + pow(orbit1.cy - pm[i].o.cy,2) 
             + pow(orbit1.maxx - pm[i].o.maxx,2) + pow(orbit1.maxy - pm[i].o.maxy,2)
	     + pow(orbit1.minx - pm[i].o.minx,2) + pow(orbit1.miny - pm[i].o.miny,2);
    PossibleMatch np = pm[i];
    np.B = B;
    npm.push_back(np);
  }
  
  return npm;
}

vector<vector<EvaluatedOrbit> > SheetMatcher::get_sheets(){
  
  return sheets;
}

bool Bcomp(PossibleMatch o1, PossibleMatch o2){ 
  
  return (o1.B<o2.B); 
}

bool fcomp(ExtremalOrbitInRUC o1, ExtremalOrbitInRUC o2){ 
  
  return (o1.f<o2.f); 
}

FrequencyCalculator::FrequencyCalculator(GlobalSettings& settings, vector<vector<EvaluatedOrbit> > sheets_in, boost::multi_array<fptype, 2> hruc_in){
  
  sheets = sheets_in;
  nsheets = sheets.size();
  minimumfreq = settings.minimumfreq;
  
  Eigen::Matrix<fptype,3,3> h;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      h(i,j) = hruc_in[i][j];
    }
  }
  
  hinv = h.inverse();
  
  maxkdiff = settings.maxkdiff;
  maxfreqdiff = settings.maxfreqdiff;
  
  for(int i=0;i<nsheets;i++){
    int n = sheets[i].size();
    norbits.push_back(n);
  }
  
  //single out extremal orbits
  for(int i=0;i<nsheets;i++){
    for(int j=0;j<norbits[i];j++){
      if(orbit_extremal(i, j)){
	extremalorbits.push_back(sheets[i][j]);
      }
    }
  }
  
  //cout << extremalorbits.size() << endl;
  
  //transform coordinates back to ruc scaled x,y,z=0...1
  transform_to_ruc();
  
  //now that we have transformed back to scaled ruc coordinates, we can easily match copies of the same orbit
  //first group the orbits by maximum k-difference defined in settings.maxkdiff
  group_orbits();
  
  //now that we have groups of extremal orbits that are sorted by position and frequency, we can average their properties within the group
  average_properties();
}

bool FrequencyCalculator::orbit_extremal(int sheetindex, int orbitindex){
  
  EvaluatedOrbit o, ol, ou; //orbit to check and two neighbouring orbits on the sheet
  o = sheets[sheetindex][orbitindex];
  
  if(orbitindex == 0){
    ol = sheets[sheetindex][norbits[sheetindex]-1];
    ou = sheets[sheetindex][orbitindex+1];
  }
  else if(orbitindex == (norbits[sheetindex]-1)){
    ol = sheets[sheetindex][orbitindex-1];
    ou = sheets[sheetindex][0];
  }
  else{
    ol = sheets[sheetindex][orbitindex-1];
    ou = sheets[sheetindex][orbitindex+1];
  }
  
  bool minimum = ((o.f <= ou.f) && (o.f <= ol.f));
  bool maximum = ((o.f >= ou.f) && (o.f >= ol.f));
  bool extremal = (minimum || maximum);
  
  return extremal;
}

void FrequencyCalculator::transform_to_ruc(){
  
  Eigen::Matrix<fptype,3,1> vec;
  int nextremalorbits = extremalorbits.size();
  for(int i=0;i<nextremalorbits;i++){
    vec(0,0) = extremalorbits[i].cx;
    vec(1,0) = extremalorbits[i].cy;
    vec(2,0) = extremalorbits[i].z;
    
    vec = hinv*vec; //transform k-space coordinates to x,y,z = 0...1 ind ruc, fmod is for shift
    for(int j=0;j<3;j++){
      vec(j,0) = fmod(vec(j,0),1);
      if((vec(j,0)) < 0){
        vec(j,0) += 1.0;
      }
    }
    
    ExtremalOrbitInRUC orb;
    orb.f = extremalorbits[i].f;
    orb.m = extremalorbits[i].m;
    orb.x = vec(0,0);
    orb.y = vec(1,0);
    orb.z = vec(2,0);
    rucorbits.push_back(orb);
  }
}

void FrequencyCalculator::group_orbits(){
  
  int norbits = rucorbits.size();
  vector<vector<ExtremalOrbitInRUC> > sorted_by_k;
  vector<vector<ExtremalOrbitInRUC> > sorted_by_f;
  
  
  //first group by k-distance criterion
  vector<bool> grouped_by_k;
  for(int i=0;i<norbits;i++){
    grouped_by_k.push_back(false);
  }
  
  for(int i=0;i<norbits-1;i++){
    vector<ExtremalOrbitInRUC> group;
    if(!grouped_by_k[i]){
      group.push_back(rucorbits[i]);
      grouped_by_k[i] = true;
      for(int j=i+1;j<norbits;j++){
        if((!grouped_by_k[j]) && within_kdistance(i, j)){
	  group.push_back(rucorbits[j]);
	  grouped_by_k[j] = true;
        }
      }
      sorted_by_k.push_back(group);
    } 
  }
  
  //now build subgroups from frequency criterion
  int ngroups = sorted_by_k.size();
  for(int i=0;i<ngroups;i++){
    //sort this already k-grouped orbits by frequency
    vector<ExtremalOrbitInRUC> thisgroup = sorted_by_k[i];
    sort(thisgroup.begin(), thisgroup.end(), fcomp);
    
    int norbits = thisgroup.size();
    
    vector<bool> grouped_by_f;
    for(int j=0;j<norbits;j++){
      grouped_by_f.push_back(false);
    }
    
    for(int j=0;j<norbits-1;j++){
      vector<ExtremalOrbitInRUC> newgroup;
      if(!grouped_by_f[j]){
	newgroup.push_back(thisgroup[j]);
	grouped_by_f[j] = true;
	for(int k=j+1;k<norbits;k++){
	  if((!grouped_by_f[k]) && within_fdistance(thisgroup, j, k)){
	    newgroup.push_back(thisgroup[k]);
	    grouped_by_f[k] = true;
	  }
        }
        sorted_by_f.push_back(newgroup);
      }
    }
  }
  
  grouped_orbits = sorted_by_f;
}

bool FrequencyCalculator::within_kdistance(int i1, int i2){
  
  bool xok = (((rucorbits[i1].x - maxkdiff) < rucorbits[i2].x) && ((rucorbits[i1].x + maxkdiff) > rucorbits[i2].x));
  bool yok = (((rucorbits[i1].y - maxkdiff) < rucorbits[i2].y) && ((rucorbits[i1].y + maxkdiff) > rucorbits[i2].y));
  bool zok = (((rucorbits[i1].z - maxkdiff) < rucorbits[i2].z) && ((rucorbits[i1].z + maxkdiff) > rucorbits[i2].z));
  
  return (xok && yok && zok);
}

bool FrequencyCalculator::within_fdistance(vector<ExtremalOrbitInRUC>& orbits, int i, int j){
  
  
  bool fok = ((orbits[i].f*(1-maxfreqdiff) < orbits[j].f) && (orbits[i].f*(1+maxfreqdiff) > orbits[j].f));
  return fok;
}

void FrequencyCalculator::average_properties(){
  
  int ngroups = grouped_orbits.size();
  
  for(int i=0;i<ngroups;i++){
    vector<fptype> fvec, mvec, xvec, yvec, zvec;
    int norb=grouped_orbits[i].size();
    
    for(int j=0;j<norb;j++){
      fvec.push_back(grouped_orbits[i][j].f);
      mvec.push_back(grouped_orbits[i][j].m);
      xvec.push_back(grouped_orbits[i][j].x);
      yvec.push_back(grouped_orbits[i][j].y);
      zvec.push_back(grouped_orbits[i][j].z);
    }
    
    AveragedOrbit ao;
    ao.f = average(fvec);
    ao.m = average(mvec);
    ao.x = average(xvec);
    ao.y = average(yvec);
    ao.z = average(zvec);
    
    ao.fsdev = standarddev(fvec, ao.f);
    ao.msdev = standarddev(mvec, ao.m);
    ao.xsdev = standarddev(xvec, ao.x);
    ao.ysdev = standarddev(yvec, ao.y);
    ao.zsdev = standarddev(zvec, ao.z);
    
    ao.n = norb;
    
    if(ao.f > minimumfreq){
      averagevec.push_back(ao);
    }
  }
  
  //sort averaged groups by frequency
  sort(averagevec.begin(), averagevec.end(), fcompaveragedorbit);
}

vector<AveragedOrbit> FrequencyCalculator::get_properties(){
  
  return averagevec;
}

fptype average(vector<fptype> values){
  
  int nentries = values.size();
  fptype val = 0;
  for(int i=0;i<nentries;i++){
    val += values[i];
  }
  val /= nentries;
  return val;
}
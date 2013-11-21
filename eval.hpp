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

//eval.hpp
#include <iostream>
#include <vector>

#include "orbit.hpp"
#include "typedefs.hpp"

#ifndef EVAL_H
#define EVAL_H

struct EvaluatedOrbit{
  
  fptype f; //dhva frequency in T
  fptype m; //effective mass in m_e
  fptype z; //z coordinate
  fptype cx; //center of the x coordinates of this orbit
  fptype cy; //center of the y coordinates of this orbit
  fptype sdevx; //standard deviation of the x coordinates of this orbit
  fptype sdevy; //standard deviation of the y coordinates of this orbit
  fptype minx; //minimum x value of this orbit
  fptype maxx;
  fptype miny;
  fptype maxy;
};

class OrbitEvaluator{
  
  public:
    OrbitEvaluator(OrbitContainer* orbits, fptype sc_length_in, fptype nsc_in);
    vector<vector<EvaluatedOrbit> > get_evaluated_orbits();
  private:
    void calc_all();
    fptype sc_length, nsc;
    fptype calc_frequency(int sliceindex, int orbitindex);
    fptype calc_mass(int sliceindex, int orbitindex);
    fptype calc_z(int sliceindex);
    fptype get_energy_derivative(OrbitPoint p1, OrbitPoint p2);
    bool glance_direction_parallel(int d1, int d2);
    fptype calc_center_x(int sliceindex, int orbitindex);
    fptype calc_center_y(int sliceindex, int orbitindex);
    fptype calc_standarddev_x(int sliceindex, int orbitindex, fptype cx);
    fptype calc_standarddev_y(int sliceindex, int orbitindex, fptype cy);
    fptype calc_min_x(int sliceindex, int orbitindex);
    fptype calc_max_x(int sliceindex, int orbitindex);
    fptype calc_min_y(int sliceindex, int orbitindex);
    fptype calc_max_y(int sliceindex, int orbitindex);
    OrbitContainer con;
    int nslices;
    vector<int> norbits;
    vector<vector<EvaluatedOrbit> > res;
};

struct PossibleMatch{
  
  EvaluatedOrbit o;
  int i; //sliceindex
  int j; //orbitindex
  fptype B; //matching parameter
};

class SheetMatcher{
  
  public:
    SheetMatcher(vector<vector<EvaluatedOrbit> > orbits_in);
    vector<vector<EvaluatedOrbit> > get_sheets();
  private:
    void find_sheet(int sliceindex, int orbitindex);
    bool simple_matching_condition_fulfilled(int i1, int j1, int i2, int j2);
    PossibleMatch get_best_match(EvaluatedOrbit orbit1, vector<PossibleMatch> pm);
    vector<PossibleMatch> calc_matching_parameter(EvaluatedOrbit orbit1, vector<PossibleMatch> pm);
    vector<vector<EvaluatedOrbit> > eorbits;
    int nslices;
    vector<int> norbits;
    vector<vector<bool> > matched;
    vector<vector< EvaluatedOrbit> > sheets; 
};

struct ExtremalOrbitInRUC{
  
  fptype f;
  fptype m;
  fptype x;
  fptype y;
  fptype z;
};

struct AveragedOrbit{
  
  fptype f; //properties and standard deviations
  fptype fsdev;
  fptype m;
  fptype msdev;
  fptype x;
  fptype xsdev;
  fptype y;
  fptype ysdev;
  fptype z;
  fptype zsdev;
  int n;
};

class FrequencyCalculator{
  
  public:
    FrequencyCalculator(GlobalSettings& settings, vector<vector<EvaluatedOrbit> > sheets_in, boost::multi_array<fptype, 2> hruc_in);
    vector<AveragedOrbit> get_properties();
  private:
    bool orbit_extremal(int sheetindex, int orbitindex);
    bool within_kdistance(int i1, int i2);
    bool within_fdistance(vector<ExtremalOrbitInRUC>& orbits, int i, int j);
    void transform_to_ruc();
    void group_orbits();
    void average_properties();
    int nsheets;
    fptype minimumfreq;
    fptype maxkdiff;
    fptype maxfreqdiff;
    vector<int> norbits;
    vector<vector<EvaluatedOrbit> > sheets;
    vector<EvaluatedOrbit> extremalorbits;
    vector<ExtremalOrbitInRUC> rucorbits;
    vector<vector<ExtremalOrbitInRUC> > grouped_orbits;
    vector<AveragedOrbit> averagevec;
    Eigen::Matrix<fptype,3,3> hinv;
};

#endif

fptype average(vector<fptype> values);
fptype standarddev(vector<fptype> values, fptype average);
bool xcomp(OrbitPoint p1, OrbitPoint p2);
bool ycomp(OrbitPoint p1, OrbitPoint p2);
bool Bcomp(PossibleMatch o1, PossibleMatch o2);
bool fcomp(ExtremalOrbitInRUC o1, ExtremalOrbitInRUC o2);

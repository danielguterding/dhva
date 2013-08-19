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

//orbit.hpp
#include <cstdio>
#include <vector>
#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

#include "typedefs.hpp"
#include "sc.hpp"

using namespace std;

#ifndef ORBIT_H
#define ORBIT_H

struct OrbitPoint{
  int i;
  int j;
  int ig;
  int jg;
  fptype x;
  fptype y;
  fptype dEparallel;
  fptype dEperpendicular;
  int dir; //glancing direction
};

class OrbitContainer{
  public:
    OrbitContainer();
    void set_slicecount(int n);
    void new_orbit(int sliceindex);
    void add_orbitpoint(int sliceindex, OrbitPoint p);
    bool orbit_closed(int sliceindex, int orbitindex);
    bool simple_orbit_closed(int sliceindex);
    void delete_empty_and_open_orbits();
    void print_orbits();
    void write_orbits(boost::filesystem::path filepath);
    int get_slicecount();
    int get_orbitcount(int sliceindex);
    vector<OrbitPoint> get_orbit(int sliceindex, int orbitindex);
    OrbitPoint get_first_orbitpoint(int sliceindex);
    OrbitPoint get_last_orbitpoint(int sliceindex);
    void update_last_orbitpoint(int sliceindex, OrbitPoint p);
    bool last_orbit_empty(int sliceindex);
  private:
    int nslices;
    vector< vector< vector< OrbitPoint > > > orbitdata;
};

class OrbitFinder{
  public:
    OrbitFinder(GlobalSettings& settings, SuperCell& sc);
    OrbitContainer get_orbits();
    OrbitContainer* get_orbits_pointer();
  private:
    int nksc;
    boost::multi_array<fptype,3> energies;
    boost::multi_array<fptype,4> kpoints;
    boost::multi_array<bool,3> unchecked;
    OrbitContainer orbitcont;
    void start();
    void stepper(const int i_in, const int j_in, const int k_in);
    void record_fs();
    inline void glance_north();
    inline void glance_east();
    inline void glance_south();
    inline void glance_west();
    void set_checked();
    bool orbit_closed();
    bool glanced_outside_fs();
    void step_to_glanced_point();
    void glance1();
    void glance2();
    inline int new_glance_direction_for_glance1();
    inline int new_glance_direction_for_glance2();
    inline bool north_of(int i1, int j1, int i2, int j2); //is p2 north of p1 ?
    inline bool east_of(int i1, int j1, int i2, int j2);
    inline bool south_of(int i1, int j1, int i2, int j2);
    inline bool west_of(int i1, int j1, int i2, int j2);
    bool point_on_sc_border();
    bool circle_detected();
    void update_stephistory();
    //variables for stepper algorithm
    int i_bef, j_bef, i, j, k, ig, jg;
    int dir;
    bool stepper_done;
    //check last five points, if 1==5 and 2==6 we are in a loop
    int ip[6];
    int jp[6];  
};

#endif

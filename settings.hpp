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

//settings.hpp
#ifndef __SETTINGS_H_INCLUDED__ 
#define __SETTINGS_H_INCLUDED__ 

#include "typedefs.hpp"

struct GlobalSettings{
  int inputinev; //specifies the units in which energies are given in input files, 0==Rydberg, 1==eV
  int nksc; //number of k-points along one edge of the super cell
  fptype nsc; //number of reciprocal unit cells along on edge of the super cell
  fptype phi; //angle between of the magnetic field measured down from z_RUC toward xy_RUC
  fptype theta; //angle of the magnetic field measured from x_RUC toward y_RUC
  fptype maxkdiff; //maximum fraction of the reciprocal lattice vectors lengths extremal orbits can differ from each other to be taken as copies of each other
  fptype maxfreqdiff; //maximum fraction of the dhva frequency ...
  fptype minimumfreq; //minimum frequency, all smaller frequencies are neglected
  int ip; //interpolator type, 0=linear, 1=cubic
  int go; //graphical out put switch, 0=no, 1=yes
};

#endif

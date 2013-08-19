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

//ruc.cpp
#include "ruc.hpp"

ReciprocalUnitCell::ReciprocalUnitCell(const boost::array<int, 3>& nkpoints, const boost::multi_array<fptype, 2>& h_arr, const boost::multi_array<fptype, 3>& e_arr){
  
  nk = nkpoints;
  h.resize(boost::extents[3][3]);
  h = h_arr;
  energies.resize(boost::extents[nk[0]][nk[1]][nk[2]]);
  energies = e_arr;
}

boost::array<int,3> ReciprocalUnitCell::get_nk(){
  
  return nk;
}

boost::multi_array<fptype, 2> ReciprocalUnitCell::get_h(){
  
  return h;
}

boost::multi_array<fptype,3> ReciprocalUnitCell::get_energies(){
  
  return energies;
}
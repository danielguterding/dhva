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

//trilinear.cpp
#include "trilinear.hpp"

TriLinearInterpolator::TriLinearInterpolator(const boost::multi_array<fptype,3>& data_in, const boost::array<int,3>& nkpoints){
  
  n1 = nkpoints[0];
  n2 = nkpoints[1];
  n3 = nkpoints[2];
  data.resize(boost::extents[n1*n2*n3]);
  for(int k=0;k<n3;k++){
    for(int j=0;j<n2;j++){
      for(int i=0;i<n1;i++)
	data[index(i,j,k)] = data_in[i][j][k];
    }
  }
}

fptype TriLinearInterpolator::operator()(fptype x, fptype y, fptype z){
  
  fptype dx = fmod(x, 1), dy = fmod(y, 1), dz = fmod(z, 1); //determine the relative position in the box enclosed by nearest data points
  
  int xi = (int)floor(x); //calculate lower-bound grid indices
  int yi = (int)floor(y);
  int zi = (int)floor(z);
  
  fptype v000 = data[index(xi, yi, zi)];
  fptype v100 = data[index(xi+1, yi, zi)];
  fptype v010 = data[index(xi, yi+1, zi)];
  fptype v001 = data[index(xi, yi, zi+1)];
  fptype v101 = data[index(xi+1, yi, zi+1)];
  fptype v011 = data[index(xi, yi+1, zi+1)];
  fptype v110 = data[index(xi+1, yi+1, zi)];
  fptype v111 = data[index(xi+1, yi+1, zi+1)];
 
  fptype result = v000*(1-dx)*(1-dy)*(1-dz) + v100*dx*(1-dy)*(1-dz) + v010*(1-dx)*dy*(1-dz) 
                + v001*(1-dx)*(1-dy)*dz + v101*dx*(1-dy)*dz + v011*(1-dx)*dy*dz 
                + v110*dx*dy*(1-dz) + v111*dx*dy*dz;
  return result;
}

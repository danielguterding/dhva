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

//sc.cpp 
#include "sc.hpp"

SuperCell::SuperCell(GlobalSettings& settings, ReciprocalUnitCell& ruc){
  
  nksc = settings.nksc;
  phi = settings.phi;
  theta = settings.theta;
  nsc = settings.nsc;
  
  kpoints.resize(boost::extents[nksc][nksc][nksc][3]);
  kpoints_rucframe_reduced.resize(boost::extents[nksc][nksc][nksc][3]);
  kpoints_ip_indices.resize(boost::extents[nksc][nksc][nksc][3]);
  energies.resize(boost::extents[nksc][nksc][nksc]);
  
  calc_length_longest_ruc_vector(ruc.get_h());
  calc_sc_kgrid(nksc);
  calc_sc_kgrid_rucframe_reduced(ruc);
  calc_sc_kgrid_ip_indices(ruc);
  
  if(settings.ip == 0){
    TriLinearInterpolator ip(ruc.get_energies(), ruc.get_nk());
    calc_sc_energies_linear(ip);
  }
  else if(settings.ip == 1){
    fptype spacing = 1.0;
    TriCubicInterpolator ip(ruc.get_energies(), spacing, ruc.get_nk());
    calc_sc_energies_cubic(ip);
  }
  else{
    cout << "Error. Interpolation Method not present." << endl;
  }
}

void SuperCell::calc_anglematrix(){
  
  Eigen::Matrix<fptype,3,3> m_mat;
  
  fptype s = sin(phi), t = cos(phi), u = 1 - cos(phi), v = sin(theta), w = cos(theta);
  
  m_mat(0,0) = pow(v,2)*u + t;
  m_mat(0,1) = -v*w*u;
  m_mat(0,2) = -w*s;
  m_mat(1,0) = m_mat(0,1);
  m_mat(1,1) = pow(w,2)*u + t;
  m_mat(1,2) = -v*s;
  m_mat(2,0) = -m_mat(0,2);
  m_mat(2,1) = -m_mat(1,2);
  m_mat(2,2) = t;
  
  anglematrix = m_mat.inverse();
}

void SuperCell::calc_transformmatrix(const boost::multi_array<fptype, 2>& h){
  
  Eigen::Matrix<fptype,3,3> h_mat;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      h_mat(i,j) = h[i][j];
    }
  }
  
  transformmatrix = h_mat.inverse();
}

void SuperCell::calc_length_longest_ruc_vector(const boost::multi_array<fptype, 2>& h){
  
  boost::array<fptype,3> lengths;
  for(int i=0;i<3;i++){
    lengths[i] = 0;
    for(int j=0;j<3;j++){
      lengths[i] += pow(h[j][i],2);
    }
    lengths[i] = sqrt(lengths[i]);
  }
  
  fptype result;
  if((lengths[0] > lengths[1]) and (lengths[0] > lengths[2])){
    result = lengths[0];
  }
  else if((lengths[1] > lengths[0]) and (lengths[1] > lengths[2])){
    result = lengths[1];
  }
  else{
    result = lengths[2];
  }
  longest_rucvec_length = result;
}

void SuperCell::calc_sc_kgrid(int nksc){
  
  boost::multi_array<fptype,1> kvals;
  kvals.resize(boost::extents[nksc]);
  for(int i=0;i<nksc;i++){
    kvals[i] = float(i)/(nksc-1) * longest_rucvec_length *nsc - longest_rucvec_length; //these are super cell k-space coordinates
  }
  
  for(int i=0;i<nksc;i++){
    for(int j=0;j<nksc;j++){
      for(int k=0;k<nksc;k++){
	kpoints[i][j][k][0] = kvals[i];
	kpoints[i][j][k][1] = kvals[j];
	kpoints[i][j][k][2] = kvals[k];
      }
    }
  }
}

void SuperCell::calc_sc_kgrid_rucframe_reduced(ReciprocalUnitCell& ruc){
  
  Eigen::Matrix<fptype,3,1> vec;
  
  calc_anglematrix();
  calc_transformmatrix(ruc.get_h());
  
  for(int i=0;i<nksc;i++){
    for(int j=0;j<nksc;j++){
      for(int k=0;k<nksc;k++){
	for(int l=0;l<3;l++){
	  vec(l,0) = kpoints[i][j][k][l];
	}
	vec = anglematrix * vec;
	vec = shift_to_ruc(vec);
	for(int l=0;l<3;l++){
	  kpoints_rucframe_reduced[i][j][k][l] = vec(l,0); //these are reduced coordinates in ruc
	}
      }
    }
  }
}

void SuperCell::calc_sc_kgrid_ip_indices(ReciprocalUnitCell& ruc){
  
  boost::array<int,3> nk;
  nk = ruc.get_nk();
  
  for(int i=0;i<nksc;i++){
    for(int j=0;j<nksc;j++){
      for(int k=0;k<nksc;k++){
	for(int l=0;l<3;l++){
	  kpoints_ip_indices[i][j][k][l] = kpoints_rucframe_reduced[i][j][k][l] * (nk[l]-1);
	}
      }
    }
  }
}


void SuperCell::calc_sc_energies_linear(TriLinearInterpolator& ip){
  
  fptype x,y,z;
  
  for(int i=0;i<nksc;i++){
    for(int j=0;j<nksc;j++){
      for(int k=0;k<nksc;k++){
	x = kpoints_ip_indices[i][j][k][0];
	y = kpoints_ip_indices[i][j][k][1];
	z = kpoints_ip_indices[i][j][k][2];
	energies[i][j][k] = ip(x,y,z);
      }
    }
  }
}
  
void SuperCell::calc_sc_energies_cubic(TriCubicInterpolator& ip){
  
  fptype x,y,z;
  
  for(int i=0;i<nksc;i++){
    for(int j=0;j<nksc;j++){
      for(int k=0;k<nksc;k++){
	x = kpoints_ip_indices[i][j][k][0];
	y = kpoints_ip_indices[i][j][k][1];
	z = kpoints_ip_indices[i][j][k][2];
	energies[i][j][k] = ip(x,y,z);
      }
    }
  }
}

Eigen::Matrix<fptype,3,1> SuperCell::shift_to_ruc(Eigen::Matrix<fptype,3,1> vec){
  
  vec = transformmatrix * vec;
  for(int i=0;i<3;i++){
    vec(i,0) = fmod(vec(i,0),1);
    if(vec(i,0) < 0){
      vec(i,0) += 1.0;
    }
  } 
  return vec;
}

boost::multi_array<fptype,3> SuperCell::get_energies(){
  
  return energies;
}

boost::multi_array<fptype,3> * SuperCell::get_energies_pointer(){
  
  return &energies;
}

boost::multi_array<fptype,4> SuperCell::get_kpoints(){
  
  return kpoints;
}

boost::multi_array<fptype,4> * SuperCell::get_kpoints_pointer(){
  
  return &kpoints;
}

fptype SuperCell::get_sc_length(){
  
  return nsc*longest_rucvec_length;
}
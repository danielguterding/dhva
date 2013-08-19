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

//sc.hpp
#include <iostream>
#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <Eigen/Dense>

#include "typedefs.hpp"
#include "settings.hpp"
#include "ruc.hpp"
#include "tricubic.hpp"
#include "trilinear.hpp"

#ifndef SUPER_CELL_H
#define SUPER_CELL_H

using namespace std;

class SuperCell{
  public:
    SuperCell(GlobalSettings& settings, ReciprocalUnitCell& ruc);
    boost::multi_array<fptype,3> get_energies();
    boost::multi_array<fptype,3> * get_energies_pointer();
    boost::multi_array<fptype,4> get_kpoints();
    boost::multi_array<fptype,4> * get_kpoints_pointer();
    fptype get_sc_length();
  private:
    int nksc;
    float nsc;
    fptype phi, theta;
    fptype longest_rucvec_length;
    boost::multi_array<fptype,4> kpoints;
    boost::multi_array<fptype,4> kpoints_rucframe_reduced;
    boost::multi_array<fptype,4> kpoints_ip_indices;
    boost::multi_array<fptype,3> energies;
    Eigen::Matrix<fptype,3,3> anglematrix; //T^-1
    Eigen::Matrix<fptype,3,3> transformmatrix;  //M^-1
    void calc_anglematrix();
    void calc_transformmatrix(const boost::multi_array<fptype,2>& h);
    void calc_length_longest_ruc_vector(const boost::multi_array<fptype,2>& h);
    void calc_sc_kgrid(const int nksc);
    void calc_sc_kgrid_rucframe_reduced(ReciprocalUnitCell& ruc);
    void calc_sc_kgrid_ip_indices(ReciprocalUnitCell& ruc);
    void calc_sc_energies_linear(TriLinearInterpolator& ip);
    void calc_sc_energies_cubic(TriCubicInterpolator& ip);
    Eigen::Matrix<fptype,3,1> shift_to_ruc(Eigen::Matrix<fptype,3,1> vec);
};

#endif
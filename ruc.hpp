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

//ruc.hpp
#include <boost/array.hpp>
#include <boost/multi_array.hpp>

#include "typedefs.hpp"

#ifndef RECIPROCAL_UNIT_CELL_H
#define RECIPROCAL_UNIT_CELL_H

class ReciprocalUnitCell{
  public:
    ReciprocalUnitCell(const boost::array<int, 3>& nkpoints, const boost::multi_array<fptype, 2>& h_arr, const boost::multi_array<fptype, 3>& e_arr);
    boost::array<int,3> get_nk();
    boost::multi_array<fptype, 2> get_h();
    boost::multi_array<fptype,3> get_energies();
  private:
    boost::array<int, 3> nk; //number is number of entries
    boost::multi_array<fptype, 2> h; //number is number of dimensions, number of elements must be set in constructor
    boost::multi_array<fptype, 3> energies;
};

#endif

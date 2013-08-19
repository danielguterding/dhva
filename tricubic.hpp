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

//tricubic.hpp
#include <Eigen/Dense>
#include <boost/array.hpp>
#include <boost/multi_array.hpp>

#include "typedefs.hpp"

#ifndef TRI_CUBIC_INTERPOLATOR_H
#define TRI_CUBIC_INTERPOLATOR_H

//This code is adapted from https://github.com/deepzot/likely
class TriCubicInterpolator{
  // Performs tri-cubic interpolation within a 3D periodic grid.
  // Based on http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.89.7835
  public:
    TriCubicInterpolator(const boost::multi_array<fptype,3>& data, const fptype& spacing, const boost::array<int,3>& nkpoints);
    fptype operator()(fptype x, fptype y, fptype z);
  private:
    boost::multi_array<fptype,1> _data;
    fptype _spacing;
    int _n1, _n2, _n3;
    int _i1, _i2, _i3;
    bool _initialized;
    Eigen::Matrix<fptype,64,1> _coefs;
    Eigen::Matrix<fptype,64,64> _C;
    inline int _index(int i1, int i2, int i3) const {
        if((i1 %= _n1) < 0) i1 += _n1;
        if((i2 %= _n2) < 0) i2 += _n2;
        if((i3 %= _n3) < 0) i3 += _n3;
        return i1 + _n1*(i2 + _n2*i3);
	}
};

#endif
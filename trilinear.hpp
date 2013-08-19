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

//trilinear.hpp 
#include <iostream>
#include <boost/array.hpp>
#include <boost/multi_array.hpp>

#include "typedefs.hpp"

using namespace std;

#ifndef TRI_LINEAR_INTERPOLATOR_H
#define TRI_LINEAR_INTERPOLATOR_H

class TriLinearInterpolator{
  public:
    TriLinearInterpolator(const boost::multi_array<fptype,3>& data, const boost::array<int,3>& nkpoints);
    fptype operator()(fptype x, fptype y, fptype z);
  private:
    boost::multi_array<fptype,1> data;
    int n1, n2, n3;
    inline int index(int i1, int i2, int i3) const {
        if((i1 %= n1) < 0) i1 += n1;
        if((i2 %= n2) < 0) i2 += n2;
        if((i3 %= n3) < 0) i3 += n3;
        return i1 + n1*(i2 + n2*i3);
	}
};
#endif
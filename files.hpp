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

//files.hpp
#include <vector>
#include <string>
#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "eval.hpp"
#include "typedefs.hpp"
using namespace std;

#ifndef __BXSF_H_INCLUDED__
#define __BXSF_H_INCLUDED__

class bxsf{
  public:
    bxsf(boost::filesystem::path path);
    boost::array<int, 3> get_nkpoints();
    boost::multi_array<fptype, 2> get_h();
    int get_bandnumber();
    vector<fptype> get_energies_list();
    boost::multi_array<fptype, 3> get_energies();
  private:
    boost::filesystem::path filepath;
    boost::filesystem::ifstream filehandle;
    fptype fermi;
    boost::array<int, 3> nkpoints; //number is number of entries
    boost::multi_array<fptype, 2> h; //number is number of dimensions, number of elements must be set in constructor
    int bandnumber;
    vector<fptype> energies_list;
    boost::multi_array<fptype, 3> energies;
    void read();
};

#endif

string trim_all(const std::string &str);
void mkdir(boost::filesystem::path dir);
void write_output(GlobalSettings settings, boost::filesystem::path outfilepath, vector<AveragedOrbit> ao);
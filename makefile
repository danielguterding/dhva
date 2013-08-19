#
# Copyright (c) 2013, Daniel Guterding <guterding@itp.uni-frankfurt.de>
#
# This file is part of dhva.
#
# dhva is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# dhva is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with dhva. If not, see <http://www.gnu.org/licenses/>.
#

CXX      = g++
CXXFLAGS = -Wall -march=native -O3 -I/home/user/local/include/eigen3
LDFLAGS  = -lm -lboost_system -lboost_filesystem

OBJECTS = main.o files.o tricubic.o trilinear.o ruc.o sc.o orbit.o eval.o
DEFINES =

dhva : $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(DEFINES) $(OBJECTS) $(LDFLAGS) -o dhva

main.o : main.cpp files.hpp settings.hpp ruc.hpp sc.hpp orbit.hpp eval.hpp typedefs.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c main.cpp -o main.o

files.o : files.cpp files.hpp typedefs.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c files.cpp -o files.o
	
tricubic.o : tricubic.cpp tricubic.hpp typedefs.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c tricubic.cpp -o tricubic.o
	
trilinear.o : trilinear.cpp trilinear.hpp typedefs.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c trilinear.cpp -o trilinear.o
	
ruc.o : ruc.cpp ruc.hpp typedefs.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c ruc.cpp -o ruc.o

sc.o : sc.cpp sc.hpp typedefs.hpp settings.hpp ruc.hpp tricubic.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c sc.cpp -o sc.o
	
orbit.o : orbit.cpp orbit.hpp typedefs.hpp settings.hpp sc.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c orbit.cpp -o orbit.o
	
eval.o : eval.cpp eval.hpp typedefs.hpp settings.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c eval.cpp -o eval.o
	
clean:
	rm dhva $(OBJECTS)
#	rm -R data

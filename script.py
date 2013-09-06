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

#script for scanning a range of angles with dhva
import subprocess

def main():
  
  angles = range(-90,91,5) #sets the range of angles to scan
  jobfilename = "job.sh"
  
  filename = "data/input/example.bxsf"
  nksc = 400 #number of k-points along one side of the super cell
  nsc = 4 #number of reciprocal unit cells along one side of the super cell
  phi = 90 
  theta = 90
  maxkdiff = 1.0 #maximum k-space difference for center coordinates of orbits in one sheet 
  maxfdiff = 0.01 #maximum frequency difference among neighbouring orbits in one sheet
  minimumfreq = 50 #minimum frequency in tesla
  ip = 1 #interpolation method, 0==linear, 1==cubic
  go = 0 #graphical output, 0==no, 1==yes
  
  for an in angles:
    phi = an #we scan a range of phi angles with fixed theta
    command = './dhva %s %i %f %f %f %f %f %f %i %i' % (filename, nksc, nsc, phi, theta, maxkdiff, maxfdiff, minimumfreq, ip, go) 
    p = subprocess.Popen(command.split()) 
    p.wait()
    
  return 0
  
main()

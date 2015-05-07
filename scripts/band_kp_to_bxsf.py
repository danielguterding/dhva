#
# Copyright (c) 2015, Daniel Guterding <guterding@itp.uni-frankfurt.de>
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

import os
import sys
import numpy as np
from math import sqrt, cos, pi

def get_direct_and_reciprocal_cell_from_FPLO_in(infilename):
  angstroemtobohr = 1.88972612457
  infilehandle = open(infilename, 'r')
  lines = infilehandle.readlines()
  infilehandle.close()
  
  angstroem = False
  for l in lines:
    idx = l.find('real lattice_constants[3]={')
    if(idx>=0):
      a,b,c = l[l.find('{')+1:l.find('}')].strip().split(',')
    
    idx = l.find('real axis_angles[3]={')
    if(idx>=0):
      alpha, beta, gamma = l[l.find('{')+1:l.find('}')].strip().split(',')
      
    idx = l.find('Angstroem')
    if(idx>=0):
      angstroem = True
      
  a, b, c = float(a), float(b), float(c)
  alpha, beta, gamma = float(alpha)*pi/180.0, float(beta)*pi/180.0, float(gamma)*pi/180.0
      
  if(angstroem):
    a *= angstroemtobohr
    b *= angstroemtobohr
    c *= angstroemtobohr
      
  #construct lattice vectors 
  b3 = cos(alpha) * b;
  b2 = sqrt(b**2 - b3**2)
  a3 = cos(beta)*a
  a2 = (a*b*cos(gamma) - a3*b3) / b2
  a1 = sqrt(a**2 - a2**2 - a3**2)

  direct_cell = np.array([[a1, 0, 0], [a2, b2, 0], [a3, b3, c]])
  
  volume = direct_cell[:,0].dot(np.cross(direct_cell[:,1], direct_cell[:,2]))
  reciprocal_cell = np.zeros((3,3))
  reciprocal_cell[:,0]= 2*pi/volume * np.cross(direct_cell[:,1], direct_cell[:,2])
  reciprocal_cell[:,1] = 2*pi/volume * np.cross(direct_cell[:,2], direct_cell[:,0])
  reciprocal_cell[:,2] = 2*pi/volume * np.cross(direct_cell[:,0], direct_cell[:,1])
  
  return direct_cell, reciprocal_cell

def compose_energystring(band):
  #compose a neatly formatted string where all the energies are put in in groups of four
  out = ''
  for i in range(0, len(band),4):
    li = i
    ui = min([i+4, len(band)])
    vals = [e for e in band[li:ui]]
    temp = '      '
    for v in vals:
      temp += '% 1.6e ' % v
    temp += '\n'
    out += temp
    
  return out

def main():
  if(7 == len(sys.argv)):
    infilename_dotin = sys.argv[1]
    infilename = sys.argv[2]
    outfilename = sys.argv[3]
    nx = int(sys.argv[4])
    ny = int(sys.argv[5])
    nz = int(sys.argv[6])
    
    direct_cell, reciprocal_cell = get_direct_and_reciprocal_cell_from_FPLO_in(infilename_dotin)
    ax, bx, cx = reciprocal_cell[0,:]/2.0/pi
    ay, by, cy = reciprocal_cell[1,:]/2.0/pi
    az, bz, cz = reciprocal_cell[2,:]/2.0/pi
    
    infilehandle = open(infilename, 'r')
    lines = infilehandle.readlines()
    infilehandle.close()
  
    splitline = lines[0].strip().split()
    nbands = int(splitline[-1]) - int(splitline[-2])+1
    nkp = int(splitline[3])
    if(not(nx*ny*nz == nkp)):
      print 'Number of k-points in supplied input file is not equal to number of k-points calculated from input paramters.'
      print 'Aborted.'
      return 1
  
    lines = [l for l in lines if not('#' == l[0])]
    energies = np.zeros((nbands, nkp), dtype=float)
    for i,l in enumerate(lines):
      for j,e in enumerate(l.strip().split()[1:]):
        energies[j,i] = float(e)
        
    energies = [band for band in energies if ((np.min(band) <=0) and (np.max(band) >= 0))]
        
    for i,band in enumerate(energies):
      outfilehandle = open(outfilename + '.%03i' % (i+1), 'w')
      startstring = '''BEGIN_INFO\n  Fermi Energy: 0.00000\nEND_INFO\n'''
      outfilehandle.write(startstring)
    
      infostring = '''BEGIN_BLOCK_BANDGRID_3D\n  band_energies\n  BANDGRID_3D_BANDS\n     1\n     %i %i %i\n     0.0 0.0 0.0\n     % f % f % f \n     % f % f % f \n     % f % f % f\n  BAND:  1\n''' % (nx, ny, nz, ax, bx, cx, ay, by, cy, az, bz, cz)
      outfilehandle.write(infostring)
    
      energystring = compose_energystring(band)
      outfilehandle.write(energystring)
    
      endstring = '''  END_BANDGRID_3D\nEND_BLOCK_BANDGRID_3D'''
      outfilehandle.write(endstring)
    
      outfilehandle.close()
  else:
    print 'Wrong number of input arguments.'
    print 'Usage:   python band_kp_to_bxsf.py [string dotinfilename] [string infilename] [string outfilename] [int nx] [int ny] [int nz]'
    print 'Example: python band_kp_to_bxsf.py =.in +band_kp outfile.bxsf 30 30 10'
    
main()

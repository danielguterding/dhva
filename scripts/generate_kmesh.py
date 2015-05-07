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

def main():
  if(5 == len(sys.argv)):
    print 'Generating k-point file.'
    infilename_dotin = sys.argv[1]
    nx = int(sys.argv[2]) #number of points along x-axis
    ny = int(sys.argv[3]) #number of points along y-axis
    nz = int(sys.argv[4]) #number of points along z-axis
    
    direct_cell, reciprocal_cell = get_direct_and_reciprocal_cell_from_FPLO_in(infilename_dotin)
    
    #for now these variables are set manually here
    onlypartiallyoccupiedbands = 't'
    loweroffset = 0
    upperoffset = 0
    
    #construct cell scaling matrix
    cm = np.linalg.norm(direct_cell[:, 0])*reciprocal_cell/2.0/pi
    
    #determine number of kpoints in the plane
    nk = nx*ny*nz
  
    outfilename = '=.kp'
    outfilehandle = open(outfilename, 'w')
    outfilehandle.write('%i %s %i %i\n' % (nk, onlypartiallyoccupiedbands, loweroffset, upperoffset))
    for x in np.linspace(-0.5, 0.5, num=nx, endpoint=False):
      for y in np.linspace(-0.5, 0.5, num=ny, endpoint=False):
        for z in np.linspace(-0.5, 0.5, num=nz, endpoint=False):
          v = np.array([x, y, z])
          v = cm.dot(v)
          outfilehandle.write('% 1.14f % 1.14f % 1.14f\n' % (v[0], v[1], v[2]))
    outfilehandle.close()
    print 'k-point file successfully written to =.kp.\nProgram finished.'
  else:
    print 'Wrong number of input arguments.'
    print 'Usage:   python generate_kmesh.py [string dotinfile] [int nx] [int ny] [int nz]'
    print 'Example: python generate_kmesh.py =.in 30 30 10'
    
main()

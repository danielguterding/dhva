import os
import numpy as np
import sys
from mayavi import mlab

def main():
  if(2 <= len(sys.argv)):
    filenames = sys.argv[1:]
    mlab.figure(bgcolor=(1,1,1), fgcolor=(0,0,0))
    for filename in filenames:
      filehandle = open(str(filename), 'r')
      lines = filehandle.readlines()
      filehandle.close()
  
      fermiindex = 0
      infoindex = 0
      bandstart = 0
      bandend = 0
      for i in range(len(lines)):
        if(lines[i].find('Fermi Energy:') != -1):
          fermiindex = i
        if(lines[i].find('BANDGRID_3D_BANDS') != -1):
          infoindex = i+2
        if(lines[i].find('BAND:') != -1):
          bandstart = i+1
        if(lines[i].find('END_BANDGRID_3D') != -1):
          bandend = i
      
      fermi = float(lines[fermiindex].split()[2])
    
      energies = []
      for i in range(bandstart, bandend):
        line = lines[i].split()
        for val in line:
          energies.append(float(val))

      nkpoints = np.array(lines[infoindex].split(), dtype='int')
      energies = np.array(energies, dtype='float')
      energies = np.reshape(energies, (nkpoints[0], nkpoints[1], nkpoints[2]))
  
      a = np.array(lines[infoindex+2].split(), dtype='float')
      b = np.array(lines[infoindex+3].split(), dtype='float')
      c = np.array(lines[infoindex+4].split(), dtype='float')
    
      #be careful, this only works for orthorhombic volumes
      x,y,z = np.mgrid[0:a[0]:nkpoints[0]*1j,0:b[1]:nkpoints[1]*1j,0:c[2]:nkpoints[2]*1j]

      src = mlab.pipeline.scalar_field(x,y,z,energies)
      mlab.pipeline.iso_surface(src, contours=[fermi], color=(1,0,0))
    
    mlab.axes(ranges=[0,a[0],0,b[1],0,c[2]], x_axis_visibility=False, y_axis_visibility=False, z_axis_visibility=False, xlabel='k_x', ylabel='k_y', zlabel='k_z')
  
  
    mlab.outline() #draws box
    mlab.show() 
  else:
    print 'Wrong number of input arguments.'
    print 'Usage: python fs3d.py input.bxsf'
main()

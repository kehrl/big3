#Since matlab can't import shapefiles, this function converts all the shapefiles of ice front 
#position and rifts in "/Data/Shape_files/Ice_fronts (Rifts)/Helheim" to text files that can be read by Matlab.
#The code is called in matlab to streamline the process.
#
#LMK, UW, 4/15/2014

import shapefile
import math
import os
import sys
import glob
import numpy
sys.path.append(os.path.join(os.getenv("HOME"),"Dropbox/Code/Modules"))
import elmer_mesh as mesh

#Convert ice front files
D=os.path.join(os.getenv("HOME"),"Data/Shape_files/Ice_fronts/Helheim")
os.chdir(D)
files=glob.glob('*.shp')

for i in range(0,len(files)):
  fname = files[i]
  sf = shapefile.Reader(fname)
  shapes = sf.shapes()
  fid = open(fname[0:8]+".dat","w")
  for j in range(0,len(shapes[0].points)):
    x=shapes[0].points[j][0]
    y=shapes[0].points[j][1]
    fid.write('{} {} \n'.format(x,y))
  fid.close()

#Convert rift files
D=os.path.join(os.getenv("HOME"),"Data/Shape_files/Rifts/Helheim")
os.chdir(D)
files=glob.glob('*.shp')

for i in range(0,len(files)):
  fname = files[i]
  sf = shapefile.Reader(fname)
  shapes = sf.shapes()
  fid = open(fname[0:8]+".dat","w")
  for j in range(0,len(shapes[-1].points)):
    x=shapes[-1].points[j][0]
    y=shapes[-1].points[j][1]
    fid.write('{} {} \n'.format(x,y))
  fid.close()



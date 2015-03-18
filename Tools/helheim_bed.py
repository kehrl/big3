# This script loads different profiles for the Helheim bed.

# LMK, UW, 10/06/2014

import os
import sys
import numpy as np
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Helheim/Observations/Elevation"))
import coords, helheim_elevation
import scipy.interpolate
import geotiff

def cresis(year):

  if year == '2001':
    file=os.path.join(os.getenv("HOME"),"Data/Bed/Cresis/helheim_flightline_05212001_good_wgs84.dat")
    data=np.loadtxt(file,skiprows=1)
    y=data[3180:3297,0]
    x=data[3180:3297,1]
    H=data[3180:3297,3]
    x2,y2 = coords.convert(x,y,4326,3413)
    
    surf=helheim_elevation.atm('2001')
    zs = scipy.interpolate.griddata(surf['20010521'][:,0:2],surf['20010521'][:,2],np.column_stack([x2,y2]),method='nearest')
    bed = np.column_stack([x2,y2,zs-H,zs])
  
  else:
    file=os.path.join(os.getenv("HOME"),"Data/Bed/Cresis/Helheim_2006_2013_Composite/flightlines/Helheim_2006_2013_Composite_Flightlines.txt")
    fid = open(file,"r")
    x=[]
    y=[]
    zs=[]
    zb=[]
    years=[]
    lines = fid.readlines()
    for line in lines[1:-1]:
      fields = line.split()
      years.append(float(fields[8][0:4]))
      x.append(float(fields[5]))
      y.append(float(fields[4]))
      zb.append(float(fields[1]))
      zs.append(float(fields[0]))
    years=np.array(years)
    x=np.array(x)
    y=np.array(y)
    zb=np.array(zb)
    zs=np.array(zs)
    
    x,y = coords.convert(x,y,4326,3413)
    
    if year == '2008':
      ind=range(32654,33126)
    elif year == '2009a':
      ind=range(46758,46950)
    elif year == '2009b':
      ind=range(53233,53391)
    elif year == '2009c':
      ind=range(56248,56448)
    elif year == '2012':
      ind=range(63289,63409)  
    elif year == '2013':
      ind=range(68366,68966)
    else:
      print "Unrecognized CreSIS profile"
    
    bed = np.column_stack([x[ind],y[ind],zb[ind]]) 

  return bed

def morlighem(xpts,ypts):
  file = os.path.join(os.getenv("HOME"),"Data/Bed/Morlighem_2014/morlighem_helheim_bed.tif")
  [x,y,z]=geotiff.read(file)
  
  f = scipy.interpolate.RectBivariateSpline(y,x,z)
  bed = f.ev(ypts,xpts)
  bed[bed<-2000]='NaN'
      
  return bed
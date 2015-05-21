# This script loads different profiles for the Helheim bed.

# LMK, UW, 10/06/2014

import os
import sys
import numpy as np
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Helheim/Tools"))
import coords, helheim_elevation
import scipy.interpolate
import geotiff

def cresis(year,verticaldatum):

  if year == '2001':
    file=os.path.join(os.getenv("HOME"),"Data/Bed/Cresis/helheim_flightline_05212001_good_wgs84.dat")
    data=np.loadtxt(file,skiprows=1)
    ind=range(3180,3297)
    y=data[:,0]
    x=data[:,1]
    H=data[:,3]
    x2,y2 = coords.convert(x,y,4326,3413)
    
    surf = helheim_elevation.atm('2001','ellipsoid')
    zs = scipy.interpolate.griddata(surf['20010521'][:,0:2],surf['20010521'][:,2],np.column_stack([x2,y2]),method='nearest')
    zb = zs-H
  
  else:
    print "Using data set Helheim_2008_2012_Composite"
    # New data set
    file=os.path.join(os.getenv("HOME"),"Data/Bed/Cresis/Helheim_2006_2014_Composite/flightlines/Helheim_2006_2014_Composite_Flightlines.txt")
    # Old data set
    #file=os.path.join(os.getenv("HOME"),"Data/Bed/Cresis/Helheim_2008_2012_Composite/flightlines/Helheim_2008_2012_Composite_Flightlines.txt")
    fid = open(file,"r")
    x=[]
    y=[]
    zs=[]
    zb=[]
    years=[]
    lines = fid.readlines()
    for line in lines[1:-1]:
      fields = line.split()
      try:
        years.append(float(fields[8][0:4]))
      except:
        years.append(0)
      x.append(float(fields[5]))
      y.append(float(fields[4]))
      zb.append(float(fields[1]))
      zs.append(float(fields[0]))
    years=np.array(years)
    x=np.array(x)
    y=np.array(y)
    zb=np.array(zb)
    zs=np.array(zs)
    
    x2,y2 = coords.convert(x,y,4326,3413)
    
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
    elif year == 'all':
      ind=range(0,len(x))
    else:
      print "Unrecognized CreSIS profile"
  
  # Select what reference we want for the elevation  
  if verticaldatum == "geoid":
    geoidheight = coords.geoidheight(x2,y2)
    zb = zb-geoidheight
  elif verticaldatum == "ellipsoid":
    zb = zb
  else:
    print "Unknown vertical datum, defaulting to ellipsoid height"
    
  return np.column_stack([x2[ind],y2[ind],zb[ind]])

def morlighem_pts(xpts,ypts,verticaldatum):
  
  # Load bed DEM
  file = os.path.join(os.getenv("HOME"),"Data/Bed/Morlighem_2014/MCdataset-2015-04-27.tif")
  [x,y,z]=geotiff.read(file)
  
  f = scipy.interpolate.RegularGridInterpolator((y,x),z,method="linear")
  bed = f(np.column_stack([ypts,xpts]))
  bed[bed<-2000]='NaN'
  
  # Morlighem bed DEM is given as elevation above mean sea level (at geoid). So we need
  # to correct only if we want the ellipsoid height.
  if verticaldatum == "ellipsoid":
    geoidheight = coords.geoidheight(xpts,ypts)
    bed = bed+geoidheight
  elif verticaldatum == "geoid":
    bed = bed
  else:
    print "Unknown vertical datum, defaulting to geoid height"
      
  return bed
  
def morlighem_grid(xmin,xmax,ymin,ymax,verticaldatum):

  # Load Bed DEM
  file = os.path.join(os.getenv("HOME"),"Data/Bed/Morlighem_2014/MCdataset-2015-04-27.tif")
  [xb,yb,zb]=geotiff.read(file,xmin,xmax,ymin,ymax)
  
  # Load Geoid
  file = os.path.join(os.getenv("HOME"),"Data/Bed/Morlighem_2014/geoid.tif")
  [xg,yg,zg]=geotiff.read(file,xmin,xmax,ymin,ymax)
  
  # Morlighem bed DEM is given as elevation above mean sea level (at geoid). So we need
  # to correct only if we want the ellipsoid height.
  if verticaldatum == "ellipsoid":
    bed = zb+zg
  elif verticaldatum == "geoid":
    bed = zb
  else:
    print "Unknown vertical datum, defaulting to geoid height"
  
  return xb,yb,bed
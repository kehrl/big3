# This script loads different profiles for the Helheim bed.

# LMK, UW, 10/06/2014

import os
import sys
import numpy as np
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
sys.path.append(os.path.join(os.getenv("HOME"),"Code/BigThreeGlaciers/Tools"))
import coords, elevation
import scipy.interpolate
import geotiff

def cresis(year,glacier,verticaldatum):

  if (glacier == 'Helheim') and (year == '2001'):

    file = os.path.join(os.getenv("HOME"),"Data/Bed/Cresis/Helheim/helheim_20010521_good_wgs84.csv")
    ind = range(3180,3297)
    
    data=np.loadtxt(file,skiprows=1,delimiter=',')
    y=data[:,0]
    x=data[:,1]
    H=data[:,3]
    x2,y2 = coords.convert(x,y,4326,3413)
    
    surf = elevation.atm('2001','ellipsoid')
    zs = scipy.interpolate.griddata(surf['20010521'][:,0:2],surf['20010521'][:,2],np.column_stack([x2,y2]),method='nearest')
    zb = zs-H
  
  else:
    if glacier == 'Helheim':
      print "Using data set Helheim_2008_2014_Composite"
      file=os.path.join(os.getenv("HOME"),"Data/Bed/Cresis/Helheim/Helheim_2006_2014_Composite/flightlines/Helheim_2006_2014_Composite_Flightlines.txt")
    elif glacier == 'Kanger':
      print "Using data set Kanger_2008_2014_Composite"
      file=os.path.join(os.getenv("HOME"),"Data/Bed/Cresis/Kanger/Kangerdlugssuaq_2006_2014_Composite/flightlines/Kangerdlugssuaq_2006_2014_Composite_Flightlines.txt")

    else: 
      print "unknown glacier"
    
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
    
    if glacier == 'Helheim':
      if year == '2006a':
        ind = range(22285,22413)
      elif year == '2006b':
        ind = range(26856,26969)
      elif year == '2008a':
        ind = range(29200,29335)
      elif year == '2008b':
        ind = range(30121,30253)
      elif year == '2009':
        ind = range(53500,53595)
      elif year == '2011':
        ind = range(61240,61265)
      elif year == '2012':
        ind = range(63900,64013)
      elif year == '2013':
        ind = range(68400,68542)
      elif year == '2014':
        ind = range(77100,77285)
      elif year == 'all':
        ind=range(0,len(x))
    elif glacier == 'Kanger':
      if year == '2008':
        ind = range(23600,23853)
      elif year == '2009a':
        ind = range(30800,31032)
      elif year == '2009b':
        ind = range(29800,29927)
      elif year == '2012':
        ind = np.arange(36130,36065,-1)
      elif year == '2013':
        ind = np.arange(39370,39210,-1)
      elif year == '2014':
        ind = range(40445,40517)
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
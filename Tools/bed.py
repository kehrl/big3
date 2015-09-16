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

def cresis(year,glacier,verticaldatum='geoid',cleanup='True'):

  if (year == '2001'):
    if glacier == 'Helheim':
      file = os.path.join(os.getenv("HOME"),"Data/Bed/Cresis/Helheim/helheim_20010521_good_wgs84.csv")
      ind = range(3180,3297)
    elif glacier == 'Kanger':
      file = os.path.join(os.getenv("HOME"),"Data/Bed/Cresis/Kanger/Data_20010520_01.csv")
      ind=range(4387,4475)
      
    data = np.loadtxt(file,skiprows=1,delimiter=',')
    y=data[:,0]
    x=data[:,1]
    H=data[:,3]
    x2,y2 = coords.convert(x,y,4326,3413)
    
    surf = elevation.atm('2001','ellipsoid')
    if glacier == 'Helheim':
      zs = scipy.interpolate.griddata(surf['20010521'][:,0:2],surf['20010521'][:,2],np.column_stack([x2,y2]),method='nearest')
      dates = np.ones(len(x2))*20010521
    elif glacier == 'Kanger':
      zs = scipy.interpolate.griddata(surf['20010520'][:,0:2],surf['20010520'][:,2],np.column_stack([x2,y2]),method='nearest')
      dates = np.ones(len(x2))*20010520
    
    zb = zs-H
  
  else:
    if glacier == 'Helheim':
      print "Using data set Helheim_2006_2014_Composite"
      file=os.path.join(os.getenv("HOME"),"Data/Bed/Cresis/Helheim/Helheim_2006_2014_Composite/flightlines/Helheim_2006_2014_Composite_Flightlines.txt")
    elif glacier == 'Kanger':
      print "Using data set Kanger_2006_2014_Composite"
      file=os.path.join(os.getenv("HOME"),"Data/Bed/Cresis/Kanger/Kangerdlugssuaq_2006_2014_Composite/flightlines/Kangerdlugssuaq_2006_2014_Composite_Flightlines.txt")

    else: 
      print "unknown glacier"
    
    fid = open(file,"r")
    x=[]
    y=[]
    zs=[]
    zb=[]
    years=[]
    type=[]
    dates = []
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
      type.append(fields[3])
      dates.append(float(fields[8]))
    years=np.array(years)
    x=np.array(x)
    y=np.array(y)
    zb=np.array(zb)
    zs=np.array(zs)
    dates=np.array(dates)
    
    x2,y2 = coords.convert(x,y,4326,3413)
    
    if glacier == 'Helheim':
      badind = np.r_[range(22406,22513),range(26969,27020),range(29300,29336),range(30253,30408),range(33224,33440),range(33552,33620),range(37452,37531),\
      		range(37640,37675),range(37819,37836),range(44127,44207),range(46942,47030),range(53595,53663),\
      		range(53713,53793),range(53974,53987),range(56646,56726),range(64006,64013),range(61237,61529),range(61745,62000),range(68541,68810),\
      		range(69202,69475),range(75645,75904),range(77285,77538),range(77728,77970)] 
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
        ind = []
        for i in range(0,len(type)):
          if ("ICESat" not in type[i]) and (i not in badind):
            ind.append(i)
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
      elif year == 'all':
        ind = []
        for i in range(0,len(type)):
          if 'ICESat' not in type[i]:
            ind.append(i)
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
    
  return np.column_stack([x2[ind],y2[ind],zb[ind],dates[ind]])

def cresis_grid(glacier,verticaldatum='geoid'):

  # Read the ASCII grid. Why they use ASCII grids, no one knows. Ugh.
  if glacier == 'Helheim':
    print "Using data set Helheim_2006_2014_Composite"
    file = os.path.join(os.getenv("HOME"),"Data/Bed/Cresis/Helheim/Helheim_2006_2014_Composite/grids/helheim_2006_2014_composite_bottom.txt")
  elif glacier == 'Kanger':
    print "Using data set Kanger_2006_2014_Composite"
    file = os.path.join(os.getenv("HOME"),"Data/Bed/Cresis/Kanger/Kangerdlugssuaq_2006_2014_Composite/grids/kangerdlugssuaq_2006_2014_composite_bottom.txt")

  # Get file info
  fid = open(file)
  ny = int(fid.readline().split()[-1])
  nx = int(fid.readline().split()[-1])
  xllcorner = float(fid.readline().split()[-1])
  yllcorner = float(fid.readline().split()[-1])
  cellsize = float(fid.readline().split()[-1])
  nodata = float(fid.readline().split()[-1])
  fid.close()
  del fid
  
  # Set up x,y
  x = np.linspace(xllcorner+cellsize/2,xllcorner+cellsize*(nx-1),nx)
  y = np.linspace(yllcorner+cellsize/2,yllcorner+cellsize*(ny-1),ny)
  
  # Get data
  grid = np.flipud(np.loadtxt(file,skiprows=6))
  grid[grid==nodata] = 'NaN'
  
  if verticaldatum == 'geoid':
    x_grid,y_grid = np.meshgrid(x,y)
    geoidheight = coords.geoidheight(x_grid.flatten(),y_grid.flatten())
    grid = grid - np.reshape(geoidheight,(ny,nx))
  elif verticaldatum == 'ellipsoid':
    grid = grid
  else:
    print "Unknown datum, defaulting to ellipsoid"

  return x,y,grid

def cresis_grid_pts(xpts,ypts,glacier,verticaldatum='geoid'):

  # Get Cresis grid
  x,y,z = cresis_grid(glacier,verticaldatum)

  # Create interpolation function
  f = scipy.interpolate.RegularGridInterpolator((y,x),z,method='linear',bounds_error=False)
  
  # Get points
  bed = f(np.column_stack([ypts,xpts]))

  return bed

def morlighem_pts(xpts,ypts,glacier,verticaldatum='geoid'):
  
  # Load bed DEM
  file = os.path.join(os.getenv("HOME"),"Data/Bed/Morlighem_2014/MCdataset-2015-04-27.tif")
  [x,y,z]=geotiff.read(file)
  z[z==-9999] = 'NaN'
  
  f = scipy.interpolate.RegularGridInterpolator((y,x),z,method='linear',bounds_error=False)
  bed = f(np.column_stack([ypts,xpts]))
  
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
  zb[zb==-9999] = 'NaN'
  
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
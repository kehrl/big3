# This module manipulate ice front positions:

# distance_along_flowline: find terminus/rift position along flowline
# load_all: load all shapefiles into a big data array

import os
import sys
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
import shapefile
import fracyear
from shapely.geometry import LineString
import numpy as np

def distance_along_flowline(x,y,dists,glacier,type='icefront'):

  if type is 'icefront':
    DIRI=os.path.join(os.getenv("HOME"),"Data/ShapeFiles/IceFronts/"+glacier+"/")
  elif type is 'rift':
    DIRI=os.path.join(os.getenv("HOME"),"Data/ShapeFiles/Rifts/"+glacier+"/")

  files = os.listdir(DIRI)

  terminus_val = []
  terminus_time = []
  if glacier == "Helheim":
    lineflow = LineString(np.row_stack([np.column_stack([x,y]),[315715,-2577820]]))
  else:
    lineflow = LineString(np.column_stack([x,y]))
    
  n = 0
  for file in files:
    if file.endswith('.shp') and (not "moon" in file):
      sf = shapefile.Reader(DIRI+file)
      shapes = sf.shapes()
      termpts = np.array(shapes[0].points[:])
      lineterm = LineString(termpts)

      # Find intersection
      intersect = (lineflow.intersection(lineterm))
      if not(intersect.is_empty):
        flowind = (abs(x-intersect.x)).argmin()
        if x[flowind] > intersect.x: #make sure that we have the smaller value
          flowind=flowind-1;
    
        # Terminus position
        terminus_val.append( dists[flowind]+((intersect.x-x[flowind])**2+(intersect.y-y[flowind])**2)**(0.5) )
    
        # Time of that terminus position
        if ("TSX" in file) or ("moon" in file):
          terminus_time.append(fracyear.doy_to_fracyear(float(file[0:4]),float(file[5:8])))
        elif ("ASTER" in file) or ("Landsat" in file):
          terminus_time.append(fracyear.date_to_fracyear(float(file[0:4]),float(file[5:7]),float(file[8:10])))

  terminus_time = np.array(terminus_time)
  terminus_val = np.array(terminus_val)
  
  sortind=np.argsort(terminus_time,0)
  terminus_time = terminus_time[sortind]
  terminus_val = terminus_val[sortind]
  
  return terminus_val, terminus_time

def position(x,y,dists,glacier,year,month,day):

  # Get all terminus positions
  terminus_val,terminus_time = distance_along_flowline(x,y,dists,glacier,type='icefront')

  # Get fractional year
  time = fracyear.date_to_fracyear(float(year),float(month),float(day))
  
  # Interpolate terminus position for the time
  terminus = np.interp(time,terminus_time,terminus_val)

  return terminus

def load_all(time1,time2,type,glacier):
  
  if type is 'icefront':
    DIRI=os.path.join(os.getenv("HOME"),"Data/ShapeFiles/IceFronts/"+glacier+"/")
  elif type is 'rift':
    DIRI=os.path.join(os.getenv("HOME"),"Data/ShapeFiles/Rifts/"+glacier+"/")

  files = os.listdir(DIRI)

  shapefiles = []
  termt = []
  for file in files:
    if file.endswith('.shp') and (not "moon" in file):
      # Time of that terminus position
      if ("TSX" in file) or ("moon" in file):
        time = fracyear.doy_to_fracyear(float(file[0:4]),float(file[5:8]))
      elif ("ASTER" in file) or ("Landsat" in file):
        time = fracyear.date_to_fracyear(float(file[0:4]),float(file[5:7]),float(file[8:10]))
      if (time > time1) and (time < time2):
        termt.append(time)
        shapefiles.append(file)

  n = len(shapefiles)

  termx = np.zeros([300,n])
  termy = np.zeros([300,n])
  termx[:,:]= 'NaN'
  termy[:,:]= 'NaN'
  for i in range(0,n):
    file = shapefiles[i]
    # Load shapefile
    sf = shapefile.Reader(DIRI+file)
    shapes = sf.shapes()
    termpts = np.array(shapes[0].points[:])
    termx[0:len(termpts[:,0]),i] = termpts[:,0]
    termy[0:len(termpts[:,0]),i] = termpts[:,1]
      
  return termx, termy, termt
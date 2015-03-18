# This module manipulate ice front positions:

# distance_along_flowline: find terminus/rift position along flowline
# load_all: load all shapefiles into a big data array

import os
import shapefile
import jdcal
from shapely.geometry import LineString
import numpy as np

def distance_along_flowline(x,y,dists,type):

  if type is 'icefront':
    DIRI=os.path.join(os.getenv("HOME"),"Data/ShapeFiles/IceFronts/Helheim/")
  elif type is 'rift':
    DIRI=os.path.join(os.getenv("HOME"),"Data/ShapeFiles/Rifts/Helheim/")

  files = os.listdir(DIRI)

  terminus_val = []
  terminus_time = []

  lineflow = LineString(np.row_stack([np.column_stack([x,y]),[315715,-2577820]]))
  n = 0
  for file in files:
    if file.endswith('.shp') and (not "moon" in file):
      # Load shapefile
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
          year = float(file[0:4])
          day = float(file[5:8])
          day1 = jdcal.gcal2jd(year,1,1) # Get length of year
          day2 = jdcal.gcal2jd(year+1,1,1)
          terminus_time.append( year + day/(day2[1]+day2[0]-day1[0]-day1[1]))
        elif ("ASTER" in file) or ("Landsat" in file):
          year = float(file[0:4])
          day = jdcal.gcal2jd(year,float(file[5:7]),float(file[8:10]))
          day2 = jdcal.gcal2jd(year,12,31)
          day1 = jdcal.gcal2jd(year-1,12,31)
          doy = day[1]+day[0]-day1[1]-day1[0]
          terminus_time.append( year + doy/(day2[1]+day2[0]-day1[0]-day1[1]))

  terminus_time = np.array(terminus_time)
  terminus_val = np.array(terminus_val)
  
  sortind=np.argsort(terminus_time,0)
  terminus_time = terminus_time[sortind]
  terminus_val = terminus_val[sortind]
  
  return terminus_val, terminus_time

def load_all(time1,time2,type):
  
  if type is 'icefront':
    DIRI=os.path.join(os.getenv("HOME"),"Data/ShapeFiles/IceFronts/Helheim/")
  elif type is 'rift':
    DIRI=os.path.join(os.getenv("HOME"),"Data/ShapeFiles/Rifts/Helheim/")

  files = os.listdir(DIRI)

  shapefiles = []
  termt = []
  for file in files:
    if file.endswith('.shp') and (not "moon" in file):
      if ("TSX" in file) or ("moon" in file):
        year = float(file[0:4])
        day = float(file[5:8])
        day1 = jdcal.gcal2jd(year,1,1) # Get length of year
        day2 = jdcal.gcal2jd(year+1,1,1)
        time = ( year + day/(day2[1]+day2[0]-day1[0]-day1[1]))
      elif ("ASTER" in file) or ("Landsat" in file):
        year = float(file[0:4])
        day = jdcal.gcal2jd(year,float(file[5:7]),float(file[8:10]))
        day2 = jdcal.gcal2jd(year,12,31)
        day1 = jdcal.gcal2jd(year-1,12,31)
        doy = day[1]+day[0]-day1[1]-day1[0]
        time = ( year + doy/(day2[1]+day2[0]-day1[0]-day1[1]))
      if (time > time1) and (time < time2):
        termt.append(time)
        shapefiles.append(file)

  n = len(shapefiles)


  termx = np.zeros([150,n])
  termy = np.zeros([150,n])
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
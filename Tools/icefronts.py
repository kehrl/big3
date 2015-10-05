# This module manipulate ice front positions:

# distance_along_flowline: find terminus/rift position along flowline
# load_all: load all shapefiles into a big data array

import os
import sys
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"Util/Modules"))
import shapefile
import fracyear
from shapely.geometry import LineString
import numpy as np
from subprocess import call
import glob

def cleanup(glacier,type='icefront'):

  '''
  cleanup(glacier,type='icefront')
  
  QGIS stupidly doesn't cleanup the shapefiles after you edit/delete features,
  so we have to do it manually. Basically I save a temporary shapefile and then save it
  again with its original name. Saving the file automatically cleans the file. There is no 
  output.
  '''
  
  if type is 'icefront':
    DIRI=os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/IceFronts/"+glacier+"/")
  elif type is 'rift':
    DIRI=os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/Rifts/"+glacier+"/")

  files = os.listdir(DIRI)
  os.chdir(DIRI)
  for file in files:
    if file.endswith('shp'):
      # Create temporary shapefile that has been loaded by ogr, so now it's been cleaned up
      call(["ogr2ogr","temp.shp",file])
      shpfiles = glob.glob(file[0:-4]+"*")
      # Delete old shapefile
      for shpfile in shpfiles:
        os.remove(shpfile)
      tempfiles = glob.glob("temp*")
      # Save new shapefile as original filename
      call(["ogr2ogr",file,"temp.shp"])
      for tempfile in tempfiles:
        os.remove(tempfile)

  return "success"

def distance_along_flowline(x,y,dists,glacier,type='icefront'):

  '''
  terminus_val,terminus_time=distance_along_flowline(x,y,dists,glacier,type='icefront')
  
  Find terminus position (or rift position) along a flowline.
  
  Inputs:
  x,y,dists: x and y coordinates for flowline, and distance along the flowline
  glacier: glacier name
  type: icefront or rift positions
  
  Outputs:
  terminus_val: distance of terminus along flowline
  terminus_time: time for each terminus_val
  '''
  
  if type is 'icefront':
    DIRI=os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/IceFronts/"+glacier+"/")
  elif type is 'rift':
    DIRI=os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/Rifts/"+glacier+"/")

  files = os.listdir(DIRI)

  terminus_val = []
  terminus_time = []
  dates = []

  lineflow = LineString(np.column_stack([x,y]))
    
  n = 0

  for file in files:
    intersection = []
    if file.endswith('.shp') and (not "moon" in file):
      sf = shapefile.Reader(DIRI+file)
      shapes = sf.shapes()
      for shape in shapes:
        termpts = np.array(shape.points[:])
        # Only look for intersection if there are points in the shape
        if len(termpts) > 0:
          lineterm = LineString(termpts)
          # Find intersection
          intersect = (lineflow.intersection(lineterm))
          if not(intersect.is_empty):
            if len(intersection) > 0:
              print "More than one position along the flowline where the ice front intersects."
              print file
            else:
              intersection = np.array([intersect.x,intersect.y])

      if len(intersection) > 0:
        flowind = (abs(x-intersection[0])).argmin()
        if x[flowind] > intersection[0]: #make sure that we have the smaller value
          flowind=flowind-1;
    
        # Terminus position
        terminus_val.append( dists[flowind]+((intersection[0]-x[flowind])**2+(intersection[1]-y[flowind])**2)**(0.5) )
        dates.append(file[0:-4])
        
        # Time of that terminus position
        if ("TSX" in file) or ("moon" in file):
          terminus_time.append(fracyear.doy_to_fracyear(float(file[0:4]),float(file[5:8])))
        elif ("ASTER" in file) or ("Landsat" in file):
          terminus_time.append(fracyear.date_to_fracyear(float(file[0:4]),float(file[5:7]),float(file[8:10])))

  terminus_time = np.array(terminus_time)
  terminus_val = np.array(terminus_val)
  terminus_file = np.array(dates)
  
  sortind=np.argsort(terminus_time,0)
  terminus_time = terminus_time[sortind]
  terminus_val = terminus_val[sortind]
  terminus_file = terminus_file[sortind]
  
  return terminus_val, terminus_time

def position(x,y,dists,glacier,time):

  '''
  terminus = position(x,y,dists,glacier,time)
  
  Get the terminus position along the flowline for time "time."
  
  Inputs:
  x,y,dists: x,y coordinates and their distance along the flowline
  glacier: glacier name
  time: time when want the terminus position
  
  Outputs:
  terminus: distance of terminus position along flowline (using linear interpolation from
  		the nearby picked terminus positions)
  '''

  # Get all terminus positions
  terminus_val,terminus_time = distance_along_flowline(x,y,dists,glacier,type='icefront')
  
  if len(time) > 1:
    time = fracyear.date_to_fracyear(time[0],time[1],time[2])
  
  # Interpolate terminus position for the time
  terminus = np.interp(time,terminus_time,terminus_val)

  return terminus

def load_all(time1,time2,glacier,type='icefront'):
  
  '''
  termx,termy,termt = load_all(time1,time2,type,glacier)
  
  Load all terminus positions for the chosen glacier.
  
  Inputs:
  time1,time2: load terminus positions from time1 to time2
  glacier: glacier name
  type: icefront or rift
  
  Outputs:
  termx,termy: array of x,y coordinates of terminus positions for times in termt
  '''
  
  if type is 'icefront':
    DIRI=os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/IceFronts/"+glacier+"/")
  elif type is 'rift':
    DIRI=os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/Rifts/"+glacier+"/")

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
  numpoints = 0
  for i in range(0,n):
    file = shapefiles[i]
    # Load shapefile
    sf = shapefile.Reader(DIRI+file)
    shapes = sf.shapes()
    for shape in shapes:
      termpts = np.array(shape.points[:])
      if len(termpts[:,0]) > 0:
        termx[numpoints:numpoints+len(termpts[:,0]),i] = termpts[:,0]
        termy[numpoints:numpoints+len(termpts[:,0]),i] = termpts[:,1]
        numpoints = numpoints+len(termpts[:,0])
      
  return termx, termy, termt

def near_time(time,glacier):

  '''
  best_x,best_y,best_time=near_time(time,glacier)
  
  Inputs:
  time: time when we want terminus position
  glacier: glacier name
  
  Outputs:
  best_x,best_y: x,y coordinates for terminus position that is closest in time to "time"
  '''

  DIRI=os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/IceFronts/"+glacier+"/")

  files = os.listdir(DIRI)

  best_time = []
  min_diff = 1.0
  for file in files:
    if file.endswith('.shp') and (not "moon" in file):
      # Time of that terminus position
      if ("TSX" in file) or ("moon" in file):
        icetime = fracyear.doy_to_fracyear(float(file[0:4]),float(file[5:8]))
      elif ("ASTER" in file) or ("Landsat" in file):
        icetime = fracyear.date_to_fracyear(float(file[0:4]),float(file[5:7]),float(file[8:10]))
      if abs(icetime-time) < min_diff:
        best_x = np.zeros(0)
        best_y = np.zeros(0)
        min_diff = abs(icetime-time)
        sf = shapefile.Reader(DIRI+file)
        shapes = sf.shapes()
        for shape in shapes:
          try:
            termpts = np.array(shape.points[:])
            best_x = np.r_[best_x,termpts[:,0]]
            best_y = np.r_[best_y,termpts[:,1]]
          except:
            print ""
        best_time = icetime

  return best_x,best_y,best_time


def calving(glacier):
  
  '''
  behavior=calving(glacier)
  
  Load files that state the type of "calving" for each satellite image. 
  
  Inputs:
  glacier: glacier name
  
  Outputs:
  behavior: a column array of time, calvingstyle, satellite, file date
  '''
  
  DIR=os.path.join(os.getenv("DATA_HOME"),"CalvingStyle/"+glacier+"/")

  # Outputs
  time = []
  value = []
  type = []
  file = []

  fid = open(DIR+"Combined.dat")
  lines = fid.readlines()
  for line in lines:
    if not(line.startswith('#')):
      p = line.split()
      value.append(p[2])
      type.append(p[1])
      file.append(p[0])
      if p[1] == 'Landsat':
        time.append(fracyear.date_to_fracyear(float(p[0][0:4]),float(p[0][5:7]),float(p[0][8:10])))
      elif p[1] == 'TSX':
        time.append(fracyear.doy_to_fracyear(float(p[0][0:4]),float(p[0][5:8])))
      else:
        print "Not working, check ", p[0]
  
  # Make them arrays
  time = np.array(time)
  type = np.array(type)
  value = np.array(value)
  file = np.array(file)
 
  # Put it into a column array
  behavior = np.column_stack([time,value,type,file])

  return behavior


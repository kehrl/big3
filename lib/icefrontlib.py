# This module manipulate ice front positions:

# distance_along_flowline: find terminus/rift position along flowline
# load_all: load all shapefiles into a big data array

import os
import sys
import shapefile
import datelib
from shapely.geometry import LineString, Polygon
import numpy as np
from subprocess import call
import glob

def cleanup(glacier,type='Icefronts'):

  '''
  cleanup(glacier,type='icefront')
  
  QGIS stupidly doesn't (always) cleanup the shapefiles after you edit/delete features,
  so we have to do it manually. Basically I save a temporary shapefile and then save it
  again with its original name. Saving the file automatically cleans the file. There is no 
  output.
  '''
  
  if type is 'icefront' or type is 'IceFronts':
    DIRI=os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/IceFronts/"+glacier+"/")
  elif type is 'rift':
    DIRI=os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/Rifts/"+glacier+"/")
  else:
    DIRI=os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/"+type+"/"+glacier+"/")
    

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

def distance_along_flowline(x,y,dists,glacier,type='icefront',imagesource=False,
		time1=-np.inf, time2=np.inf):

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
  sensorname = []

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
        
        # Time of that terminus position
        if ("TSX" in file) or ("moon" in file):
          terminus_time.append(datelib.doy_to_fracyear(float(file[0:4]),float(file[5:8])))
          sensorname.append('TSX')
        elif ("ASTER" in file):
          terminus_time.append(datelib.date_to_fracyear(float(file[0:4]),float(file[5:7]),float(file[8:10])))
          sensorname.append('ASTER')
        elif ("Landsat7" in file):
          terminus_time.append(datelib.date_to_fracyear(float(file[0:4]),float(file[5:7]),float(file[8:10])))
          sensorname.append('Landsat7')
        elif ("Landsat" in file):
          terminus_time.append(datelib.date_to_fracyear(float(file[0:4]),float(file[5:7]),float(file[8:10])))
          sensorname.append('Landsat8')
        elif ("WV" in file):
          terminus_time.append(datelib.date_to_fracyear(float(file[0:4]),float(file[5:7]),float(file[8:10])))
          sensorname.append('WV')
        else:
          sys.exit("Don't know that date format for "+file)
  terminus_time = np.array(terminus_time)
  terminus_val = np.array(terminus_val)
  terminus_source = np.array(sensorname)
  
  # Need to double check that we imported the same number of times and 
  # terminus values. If we didn't something is horribly wrong.
  if len(terminus_time) == len(terminus_val):
    sortind=np.argsort(terminus_time,0)
    terminus_time = terminus_time[sortind]
    terminus_val = terminus_val[sortind]
    terminus_source = terminus_source[sortind]
  else:
    sys.exit("Length of terminus values and times are different. Something is very wrong")
  
  # Subset record to cover only a certain time period if time1,time2 are set.
  ind = np.where((terminus_time > time1) & (terminus_time < time2))[0]
  terminus_time = terminus_time[ind]
  terminus_val = terminus_val[ind]
  terminus_source = terminus_source[ind]
  
  # Sometimes we will want to know the image source for each digitized ice front position.
  if imagesource == False:
    return terminus_val, terminus_time
  else:
    return terminus_val, terminus_time, terminus_source

def box_method(glacier,imagesource=False,time1=-np.inf, time2=np.inf):
  
  DIRI=os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/IceFronts/"+glacier+"/")  
  files = os.listdir(DIRI)

  # Load shapefile for box
  sf_box = shapefile.Reader(os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/FluxGates/"+glacier+"/box_method.shp"))
  shapes_box = sf_box.shapes()
  
  # Lines that define back of box
  pts_north = np.array(shapes_box[0].points)
  line_north = LineString(pts_north)
  pts_south = np.array(shapes_box[1].points)
  line_south = LineString(pts_south)
  pts_west = np.array(shapes_box[2].points)

  # Box width
  width = np.sqrt((pts_west[0][0]-pts_west[1][0])**2+(pts_west[0][1]-pts_west[1][1])**2)

  # Set up variables for loading
  terminus_val = []
  terminus_time = []
  sensorname = []

  n = 0
  for file in files:
    intersect_north = []
    intersect_south = []
    if file.endswith('.shp') and (not "moon" in file) and (not "Landsat7" in file):
      sf = shapefile.Reader(DIRI+file)
      shapes = sf.shapes()
      if len(shapes) > 1:
        print "check ",file
      else:
        pts_terminus = np.array(shapes[0].points[:])
        if pts_terminus[0,1] > pts_terminus[-1,1]:
          pts_terminus = np.flipud(pts_terminus)
        # Only look for intersection if there are points in the shape
        if len(pts_terminus) > 0:
          line_terminus = LineString(pts_terminus)
          # Find intersection with sides of box
          intersect_north = (line_north.intersection(line_terminus))
          intersect_south = (line_south.intersection(line_terminus))
          if intersect_north.is_empty or intersect_south.is_empty:
            print "Ice front doesn't extend across entire domain ", file
          else:
            ind_in_box = np.where((pts_terminus[:,1] < intersect_north.y) & (pts_terminus[:,1] > intersect_south.y))[0]
            pts_box = np.row_stack([np.array(pts_west),np.r_[intersect_south.xy],pts_terminus[ind_in_box,:],np.r_[intersect_north.xy]])
            box = Polygon(pts_box)
            A = box.area
    
            # Terminus position
            terminus_val.append( A/width )
        
            # Time of that terminus position
            if ("TSX" in file) or ("moon" in file):
              terminus_time.append(datelib.doy_to_fracyear(float(file[0:4]),float(file[5:8])))
              sensorname.append('TSX')
            elif ("ASTER" in file):
              terminus_time.append(datelib.date_to_fracyear(float(file[0:4]),float(file[5:7]),float(file[8:10])))
              sensorname.append('ASTER')
            elif ("Landsat" in file):
              terminus_time.append(datelib.date_to_fracyear(float(file[0:4]),float(file[5:7]),float(file[8:10])))
              sensorname.append('Landsat8')
            elif ("WV" in file):
              terminus_time.append(datelib.date_to_fracyear(float(file[0:4]),float(file[5:7]),float(file[8:10])))
              sensorname.append('WV')
            else:
              sys.exit("Don't know that date format for "+file)
  
  terminus_time = np.array(terminus_time)
  terminus_val = np.array(terminus_val)
  terminus_source = np.array(sensorname)

  # Need to double check that we imported the same number of times and 
  # terminus values. If we didn't something is horribly wrong.
  if len(terminus_time) == len(terminus_val):
    sortind=np.argsort(terminus_time,0)
    terminus_time = terminus_time[sortind]
    terminus_val = terminus_val[sortind]
    terminus_source = terminus_source[sortind]
  else:
    sys.exit("Length of terminus values and times are different. Something is very wrong")

  # Subset record to cover only a certain time period if time1,time2 are set.
  ind = np.where((terminus_time > time1) & (terminus_time < time2))[0]
  terminus_time = terminus_time[ind]
  terminus_val = terminus_val[ind]
  terminus_source = terminus_source[ind]
  
  # Sometimes we will want to know the image source for each digitized ice front position.
  if imagesource == False:
    return terminus_val, terminus_time
  else:
    return terminus_val, terminus_time, terminus_source

def quantify_uncertainty(x,y,dists,glacier,time1,time2):
  
  # Get ice front positions
  terminus_val, terminus_time,terminus_source = distance_along_flowline(x,y,dists,
  		glacier,type='icefront',imagesource=True,time1=time1,time2=time2)
  
  # Figure out what times we have multiple ice front positions for comparison
  unique,unique_count = np.unique(terminus_time,return_counts=True)
  dup_times = unique[unique_count > 1]
  
  N = len(dup_times)
  delta = np.zeros(N)
  files = [] 
  for i in range(0,N):
    ind = []
    for j in range(0,len(terminus_time)):      
      if terminus_time[j] == dup_times[i]:
        ind.append(j)
    if len(ind) == 2:
      delta[i] = terminus_val[ind[0]] - terminus_val[ind[1]]
      files.append([terminus_source[ind[0]],terminus_source[ind[1]]])
    else:
      print "problem"  
  
  return delta  

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
    time = datelib.date_to_fracyear(time[0],time[1],time[2])
  
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
        time = datelib.doy_to_fracyear(float(file[0:4]),float(file[5:8]))
      elif ("ASTER" in file) or ("Landsat" in file) or ("WV" in file):
        time = datelib.date_to_fracyear(float(file[0:4]),float(file[5:7]),float(file[8:10]))
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
    numpoints=0
    file = shapefiles[i]
    # Load shapefile
    sf = shapefile.Reader(DIRI+file)
    shapes = sf.shapes()
    for shape in shapes:
      try:
        termpts = np.array(shape.points[:])
        if len(termpts[:,0]) > 0:
          termx[numpoints:numpoints+len(termpts[:,0]),i] = termpts[:,0]
          termy[numpoints:numpoints+len(termpts[:,0]),i] = termpts[:,1]
          numpoints = numpoints+len(termpts[:,0])
      except:
        pass
      
  return termx, termy, termt

def near_time(time,glacier,type='all'):

  '''
  best_x,best_y,best_time=near_time(time,glacier)
  
  Inputs:
  time: fractional year when we want terminus position
  glacier: glacier name
  type: data source for ice front position (TSX, WV, Landsat, or all)
  
  Outputs:
  best_x,best_y: x,y coordinates for terminus position that is closest in time to "time"
  '''

  # Ice front directory
  DIRI=os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/IceFronts/"+glacier+"/")
  files = os.listdir(DIRI)

  if type == 'all':
    type = ['TSX','Landsat','ASTER','WV']

  best_time = []
  min_diff = 1.0
  for file in files:
    if file.endswith('.shp') and (not "moon" in file):
      # Time of that terminus position
      if ("TSX" in file) and ("TSX" in type):
        icetime = datelib.doy_to_fracyear(float(file[0:4]),float(file[5:8]))
      elif ("ASTER" in file) and ("ASTER" in type):
        icetime = datelib.date_to_fracyear(float(file[0:4]),float(file[5:7]),float(file[8:10]))
      elif ("Landsat" in file) and ("Landsat" in type):
        icetime = datelib.date_to_fracyear(float(file[0:4]),float(file[5:7]),float(file[8:10]))
      elif ("WV" in file) and ("WV" in file):
        icetime = datelib.date_to_fracyear(float(file[0:4]),float(file[5:7]),float(file[8:10]))
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
            pass
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

  value = []
  type = []
  file = []
  time = []

  fid = open(DIR+"Combined.dat")
  lines = fid.readlines()
  for line in lines:
    if not(line.startswith('#')):
      p = line.split()
      value.append(p[2])
      type.append(p[1])
      file.append(p[0])
      if p[1] == 'Landsat' or p[1] == 'Worldview':
        time.append(datelib.date_to_fracyear(float(p[0][0:4]),float(p[0][5:7]),float(p[0][8:10])))
      elif p[1] == 'TSX':
        time.append(datelib.doy_to_fracyear(float(p[0][0:4]),float(p[0][5:8])))
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


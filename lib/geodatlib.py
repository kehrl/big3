# This module reads Ian's velocity data.
#
# Functions:
# readgeodat(filename) - reads the .geodat file, which stores the number of pixels
#         in the image, the pixel size and the location of the lower left corner
# readbinary(filename) - reads the binary velocity data given the information about the 
#         grid from "readgeodat"
# readvelocity(filename) - checks to see if there are geotiff velocity files in the 
#		  subdirectory; if not it loads the velocities from the binary files

import numpy as np
import struct
import math
import os
import jdcal
import datelib, geotifflib

def readgeodat(filename):
    geodatfile = open(filename,"r")
    xgeo = np.zeros((3,2))
    i = 0
    while True:
        line = geodatfile.readline().split()
        if len(line)==2:
            try:
                xgeo[i,0] = float(line[0])
                xgeo[i,1] = float(line[1])
                i = i+1
            except ValueError:
                i = i
        if len(line)==0:
            break
    geodatfile.close()
    xgeo[2,:] = xgeo[2,:]*1000.0
    
    return xgeo

def readtime(filename):
	metafile = open(filename+".meta","r")	
	lines = metafile.readlines()
	
	if lines[0][0] == 'N': # Using nominal date
	  time = float(lines[0][15:21])
	  interval = 'NaN'
	else: # Using Central Julian date
	  jdate = float(lines[0][36:47])	
	
	  # Get date
	  year,month,day,fracday=jdcal.jd2gcal(jdate,0)

	  time = datelib.date_to_fracyear(year,month,day+fracday)
	  
	  # Get time interval for velocity (time of second image - time of first image)
	  month1 = datelib.month(lines[1][32:35])
	  day1 = float(lines[1][36:38])
	  year1 = float(lines[1][39:43])
	  
	  month2 = datelib.month(lines[2][33:36])
	  day2 = float(lines[2][37:39])
	  year2 = float(lines[2][40:44])
	
	  interval = datelib.date_to_fracyear(year2,month2,day2)-datelib.date_to_fracyear(year1,month1,day1)
	
	return time,interval


def readbinary(filename,nodatavalue=float('nan')):
    
  # Get time
  if os.path.exists(filename+".meta"):
    time,interval = readtime(filename)
  else:
    time = 'NaN'
    interval = 'NaN'
    
  # Get the data about the grid from the .geodat file
  xgeo = readgeodat(filename+".vx.geodat")
  (nx, ny) = (int(xgeo[0,0]), int(xgeo[0,1]))
  (dx, dy) = (xgeo[1,0], xgeo[1,1])
  (xo, yo) = (xgeo[2,0], xgeo[2,1])
    
  x = np.linspace(xo,xo+(nx-1)*dx,nx)
  y = np.linspace(yo,yo+(ny-1)*dy,ny)

  # Open the x-velocity file and read in all the binary data
  vxfile = open(filename+".vx","rb")
  vx_raw_data = vxfile.read()
  vxfile.close()

  # Unpack that binary data to an array of floats, knowing that it is
  # in big-endian format. Why? God only knows.
  nvals = len(vx_raw_data)/4
  arr = np.zeros(nvals)
  for i in range(nvals):
    arr[i] = struct.unpack('>f',vx_raw_data[4*i:4*(i+1)])[0]
  vx = arr.reshape((ny,nx))
  ind = np.where(vx == -2e9)
  vx[ind] = nodatavalue

  # Go through the same rigmarole for the y-velocity file
  vyfile = open(filename+".vy","rb")
  vy_raw_data = vyfile.read()
  vyfile.close()
  arr = np.zeros(nvals)
  for i in range(nvals):
    arr[i] = struct.unpack('>f',vy_raw_data[4*i:4*(i+1)])[0]
  vy = arr.reshape((ny,nx))
  ind = np.where(vy == -2e9)
  vy[ind] = nodatavalue

  # Do the same thing for the files containing the errors
  exfile = open(filename+".ex","rb")
  ex_raw_data = exfile.read()
  exfile.close()
  arr = np.zeros(nvals)
  for i in range(nvals):
    arr[i] = struct.unpack('>f',ex_raw_data[4*i:4*(i+1)])[0]
  ex = arr.reshape((ny,nx))
  ind = np.where(ex == -2e9)
  ex[ind] = nodatavalue

  eyfile = open(filename+".ey","rb")
  ey_raw_data = eyfile.read()
  eyfile.close()
  arr = np.zeros(nvals)
  for i in range(nvals):
    arr[i] = struct.unpack('>f',ey_raw_data[4*i:4*(i+1)])[0]
  ey = arr.reshape((ny,nx))
  ind = np.where(ey == -2e9)
  ey[ind] = nodatavalue

  # Calculate velocity magnitude
  ind = np.where(vx != -2e9)
  v = np.zeros_like(vx)
  v[:,:] = nodatavalue
  v[ind] = (vx[ind]**2+vy[ind]**2)**(0.5)
    
  return (x,y,v,vx,vy,ex,ey,time,interval)

def readvelocity(DIR,track,file):
  
  '''
  x,y,v,vx,vy,ex,ey,time,interval = readvelocity(DIR,track,file)
  '''
  
   
  # Filename
  filename = DIR+track+"/"+file
  # Get time
  if os.path.exists(filename+".meta"):
    time,interval = readtime(filename)
    year,month,day = datelib.fracyear_to_date(time)
    date = "%04d%02d%02d" % (year,month,day)
  else:
    time = float('NaN')
    interval = float('NaN')
    date = float('NaN')
  
  if os.path.isfile(DIR+"TIF/"+track+"_v.tif"):
    x,y,v = geotifflib.read(DIR+"TIF/"+track+"_v.tif")
    x,y,vx = geotifflib.read(DIR+"TIF/"+track+"_vx.tif")
    x,y,vy = geotifflib.read(DIR+"TIF/"+track+"_vy.tif")
    x,y,ex = geotifflib.read(DIR+"TIF/"+track+"_ex.tif")
    x,y,ey = geotifflib.read(DIR+"TIF/"+track+"_ey.tif")
      # Check to see if there are geotiff files in the subdirectory. If not, read the binary data.
  elif os.path.isfile(DIR+"TIF/"+track+"_"+date+"_v.tif"):
    x,y,v = geotifflib.read(DIR+"TIF/"+track+"_"+date+"_v.tif")
    x,y,vx = geotifflib.read(DIR+"TIF/"+track+"_"+date+"_vx.tif")
    x,y,vy = geotifflib.read(DIR+"TIF/"+track+"_"+date+"_vy.tif")
    x,y,ex = geotifflib.read(DIR+"TIF/"+track+"_"+date+"_ex.tif")
    x,y,ey = geotifflib.read(DIR+"TIF/"+track+"_"+date+"_ey.tif")
  else:
    print "Unpacking binary velocity file ",track
    x,y,v,vx,vy,ex,ey,time,interval = readbinary(filename)
     
  return (x,y,v,vx,vy,ex,ey,time,interval)
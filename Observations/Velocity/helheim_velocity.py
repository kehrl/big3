#Various scripts for processing and pulling Helheim velocity data:

## velocity_at_points(xpt,ypt): get velocity timeseries at a particular point
## velocity_at_flowline(x,y): get velocity timeseries along flowline

#LMK, UW, 04/01/2014

import os
import math
import sys
import scipy.interpolate
import numpy as np
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Helheim/Observations/IceFronts"))
import geodat, helheim_icefronts

def velocity_at_points(xpt,ypt):
  # Finds velocity at nearest gridcell to xpt, ypt

  try:
    n = len(xpt)
  except:
    n = 0

  DIR_TSX = os.path.join(os.getenv("HOME"),"Data/Velocity/TSX/Helheim/")
  DIR_RADARSAT = os.path.join(os.getenv("HOME"),"Data/Velocity/RADARSAT/Helheim/")

  #################
  # LOAD TSX Data #
  #################

  DIRS=os.listdir(DIR_TSX)
  tpt=[]
  
  # Get number of velocity files
  m=0
  for DIR in DIRS:
    if DIR.startswith('track'):
      m = m+1

  vpt=np.zeros([m,n])
  tpt=np.zeros([m,1])
  count=0
  for j in range(0,len(DIRS)):
    DIR=DIRS[j]
    if DIR.startswith('track'):
      infile=DIR_TSX+DIR
      x,y,v,vx,vy,ex,ey,time=geodat.readvelocity(infile+"/mosaicOffsets")
      tpt[count]=time
      
      for i in range(0,n):
        try:
          xind = (abs(x-xpt[i])).argmin()
          yind = (abs(y-ypt[i])).argmin()
        except:
      	  xind = (abs(x-xpt)).argmin()
          yind = (abs(y-ypt)).argmin()

        vpt[count,i] = v[yind,xind]
        #print (abs(x-xpt)).min()
        
      count = count + 1
  
  # Sort arrays by time
  tpt_tsx = tpt
  vpt_tsx = vpt
  
  ######################
  # Load RADARSAT data #
  ######################
  DIRS=os.listdir(DIR_RADARSAT)
  
  # Get number of velocity files
  m=0
  for DIR in DIRS:
    if DIR.startswith('winter'):
      m = m+1

  vpt=np.zeros([m,n])
  tpt=np.zeros([m,1])
  count=0
  for j in range(0,len(DIRS)):
    DIR=DIRS[j]
    if DIR.startswith('winter'):
      infile=DIR_RADARSAT+DIR
      x,y,v,vx,vy,ex,ey,time=geodat.readvelocity(infile+"/mosaicOffsets")
      tpt[count]=time
      
      for i in range(0,n):
        try:
          xind = (abs(x-xpt[i])).argmin()
          yind = (abs(y-ypt[i])).argmin()
      	except:
      	  xind = (abs(x-xpt)).argmin()
          yind = (abs(y-ypt)).argmin()

        vpt[count,i] = v[yind,xind]
        #print (abs(x-xpt)).min()
        
      count = count + 1
  
  tpt_radarsat = tpt
  vpt_radarsat = vpt
  
  tpt_all = np.row_stack([tpt_radarsat,tpt_tsx])
  vpt_all = np.row_stack([vpt_radarsat,vpt_tsx])
  
  # Sort arrays by time  
  sortind=np.argsort(tpt_all,0)
  tpt_all = tpt_all[sortind[:,0]]
  vpt_all = vpt_all[sortind[:,0],:]
      	  
  return vpt_all,tpt_all

def velocity_along_flowline(xf,yf,dists):      
  # Uses a linear interpolation to flowline coordinates xf, yf
  
  DIR_TSX = os.path.join(os.getenv("HOME"),"Data/Velocity/TSX/Helheim/")
  DIR_RADARSAT = os.path.join(os.getenv("HOME"),"Data/Velocity/RADARSAT/Helheim/")
  
  ###############################################################################
  # Load terminus profiles so we can cutoff velocities in front of the terminus #
  ###############################################################################

  term_values, term_time = helheim_icefronts.distance_along_flowline(xf,yf,dists)

  #################
  # LOAD TSX Data #
  #################

  DIRS=os.listdir(DIR_TSX)
  tpt=[]
  
  # Get number of velocity files
  n = len(xf)
  m=0
  for DIR in DIRS:
    if DIR.startswith('track'):
      m = m+1

  velocities=np.zeros([n,m])
  velocities[:,:] = 'NaN'
  times=np.zeros([m,1])
  count=0
  for j in range(0,len(DIRS)):
    DIR=DIRS[j]
    if DIR.startswith('track'):
      errors = np.zeros([m])
      infile=DIR_TSX+DIR
      x,y,v,vx,vy,ex,ey,time=geodat.readvelocity(infile+"/mosaicOffsets")
      times[count]=time
      
      # Set up grid for interpolation
      fv = scipy.interpolate.RegularGridInterpolator([y,x],v)
      
      # Find flowline coordinates within grid and behind terminus
      terminus = 82000 #np.interp(time,term_time,term_values)
      minind = np.max([(abs(np.min(x) - xf)).argmin(),(abs(np.max(y) - yf)).argmin()])
      maxind = np.min([(abs(np.max(x) - xf)).argmin(),(abs(np.min(y) - yf)).argmin(),(abs(terminus-dists)).argmin()])
      if (xf[minind] < np.min(x)) or (yf[minind] > np.max(y)):
        minind = minind+1
      if (yf[maxind] < np.min(y)) or (xf[maxind] > np.max(x)) or (dists[maxind] > terminus):
        maxind = maxind-1
      coords = np.array([yf[minind:maxind],xf[minind:maxind]]).T
      
      velocities[minind:maxind,count] = fv(coords,'linear')
    
      if DIR.endswith('33047') or DIR.endswith('33214'):
        print "Still ignoring some directories"
        velocities[:,count] = 'NaN'
      
      count = count + 1
  
  # Sort arrays by time
  time_tsx = times
  velocities_tsx = velocities
  
  ######################
  # Load RADARSAT data #
  ######################
  DIRS=os.listdir(DIR_RADARSAT)
  
  # Get number of velocity files
  m=0
  for DIR in DIRS:
    if DIR.startswith('winter'):
      m = m+1

  velocites=np.zeros([n,m])
  times=np.zeros([m,1])
  count=0
  for j in range(0,len(DIRS)):
    DIR=DIRS[j]
    if DIR.startswith('winter'):
      infile=DIR_RADARSAT+DIR
      x,y,v,vx,vy,ex,ey,time=geodat.readvelocity(infile+"/mosaicOffsets")
      times[count]=time
      
      # Set up grid for interpolation
      f = scipy.interpolate.RegularGridInterpolator([y,x],v)
      
      # Find flowline coordinates within grid and behind terminus

      minind = np.max([(abs(np.min(x) - xf)).argmin(),(abs(np.max(y) - yf)).argmin()])
      maxind = np.min([(abs(np.max(x) - xf)).argmin(),(abs(np.min(y) - yf)).argmin(),(abs(dists-terminus)).argmin()])
      print minind,maxind
      if (xf[minind] < np.min(x)) or (yf[minind] > np.max(y)):
        minind = minind+1
      if (yf[maxind] < np.min(y)) or (xf[maxind] > np.max(x)) or (dists[maxind] > terminus):
        maxind = maxind-1
      coords = np.array([yf[minind:maxind],xf[minind:maxind]]).T
      
      velocities[minind:maxind,count] = f(coords)

      count = count + 1
  
  time_radarsat = times
  velocities_radarsat = velocities
      	  
  return velocities_tsx, time_tsx






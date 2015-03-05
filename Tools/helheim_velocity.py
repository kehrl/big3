#Various scripts for processing and pulling Helheim velocity data:

## velocity_at_eulpoints(xpt,ypt): get velocity timeseries at an Eulerian point
## velocity_at_lagpoints(x,y,dists,pts): get velocity timeseries at a Lagrangian point
## velocity_along_flowline(x,y): get velocity timeseries along flowline

#LMK, UW, 04/01/2014

import os
import math
import sys
import scipy.interpolate
import numpy as np
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Helheim/Tools"))
import geodat, helheim_icefronts

#########################################################################################
def velocity_at_eulpoints(xpt,ypt):
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
  ept=np.zeros([m,n])
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
        ept[count,i] = math.sqrt(ex[yind,xind]**2+ey[yind,xind]**2)
        #print (abs(x-xpt)).min()
        
      count = count + 1
  
  # Sort arrays by time
  tpt_tsx = tpt
  ept_tsx = ept
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
  ept=np.zeros([m,n])
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
        ept[count,i] = math.sqrt(ex[yind,xind]**2+ey[yind,xind]**2)
        #print (abs(x-xpt)).min()
        
      count = count + 1
  
  tpt_radarsat = tpt
  vpt_radarsat = vpt
  ept_radarsat = ept
  
  tpt_all = np.row_stack([tpt_radarsat,tpt_tsx])
  vpt_all = np.row_stack([vpt_radarsat,vpt_tsx])
  ept_all = np.row_stack([ept_radarsat,ept_tsx])
  
  # Sort arrays by time  
  sortind=np.argsort(tpt_all,0)
  tpt_all = tpt_all[sortind[:,0]]
  vpt_all = vpt_all[sortind[:,0],:]
  ept_all = ept_all[sortind[:,0],:]
      	  
  return vpt_all,tpt_all, ept_all
  
#########################################################################################
def velocity_along_flowline(xf,yf,dists):      
  # Uses a nearest interpolation to flowline coordinates xf, yf
  
  DIR_TSX = os.path.join(os.getenv("HOME"),"Data/Velocity/TSX/Helheim/")
  DIR_RADARSAT = os.path.join(os.getenv("HOME"),"Data/Velocity/RADARSAT/Helheim/")
  
  ###############################################################################
  # Load terminus profiles so we can cutoff velocities in front of the terminus #
  ###############################################################################

  term_values, term_time = helheim_icefronts.distance_along_flowline(xf,yf,dists,'icefront')

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

  velocities=np.zeros([m,n])
  velocities[:,:] = 'nan'
  times=np.zeros([m,1])
  termini=np.zeros([m,1])
  count=0
  for j in range(0,len(DIRS)):
    DIR=DIRS[j]
    if DIR.startswith('track'):
      errors = np.zeros([m])
      infile=DIR_TSX+DIR
      x,y,v,vx,vy,ex,ey,time=geodat.readvelocity(infile+"/mosaicOffsets")
      times[count]=time
      
      # Set up grid for interpolation
      fv = scipy.interpolate.RegularGridInterpolator([y,x],v,method='nearest',bounds_error=False)
      
      # Find flowline coordinates behind terminus
      terminus = np.interp(time,term_time,term_values)
      ind = np.where(dists > terminus)
        
      coords = np.array([yf,xf]).T
      velocities[count,:] = fv(coords)
      velocities[count,ind] = 'nan'
      
      termini[count]=terminus
      
      #if DIR.endswith('33047') or DIR.endswith('33214'):
      #  print "Still ignoring some directories"
      #  velocities[:,count] = 'nan'
      
      count = count + 1
  
  # Sort arrays by time
  time_tsx = times
  velocities_tsx = velocities
  term_tsx = termini
  del velocities, time, termini
  
  ######################
  # Load RADARSAT data #
  ######################
  DIRS=os.listdir(DIR_RADARSAT)
  
  # Get number of velocity files
  m=0
  for DIR in DIRS:
    if DIR.startswith('winter'):
      m = m+1  

  velocities=np.zeros([m,n])
  velocities[:,:] = 'nan'
  times=np.zeros([m,1])
  termini = np.zeros([m,1])
  count=0
  for j in range(0,len(DIRS)):
    DIR=DIRS[j]
    if DIR.startswith('winter'):
      infile=DIR_RADARSAT+DIR
      x,y,v,vx,vy,ex,ey,time=geodat.readvelocity(infile+"/mosaicOffsets")
      times[count]=time
      
      # Set up grid for interpolation
      fv = scipy.interpolate.RegularGridInterpolator([y,x],v,method='nearest',bounds_error=False)
      
      # Find flowline coordinates behind terminus
      terminus = np.interp(time,term_time,term_values)
      ind = np.where(dists > terminus)

      coords = np.array([yf,xf]).T
      
      velocities[count,:] = fv(coords)
      velocities[count,ind] = 'nan'
      
      termini[count]=terminus

      count = count + 1
  
  time_radarsat = times
  velocities_radarsat = velocities
  term_radarsat = termini
  
  print np.shape(velocities_tsx), np.shape(velocities_radarsat)    	  
  tpt_all = np.row_stack([time_radarsat,time_tsx])
  term_all = np.row_stack([term_radarsat,term_tsx])
  vpt_all = np.row_stack([velocities_radarsat,velocities_tsx]).T
  
  # Sort arrays by time  
  sortind=np.argsort(tpt_all,0)
  tpt_sort = tpt_all[sortind[:,0]]
  term_sort = term_all[sortind[:,0]]
  vpt_sort = vpt_all[:,sortind[:,0]]  	  
  
  #return velocities_tsx,time_tsx,term_tsx    	  
  return vpt_sort,tpt_sort,term_sort


##########################################################################################
def velocity_at_lagpoints(xf,yf,dists,pts):
  
  
  DIR_TSX = os.path.join(os.getenv("HOME"),"Data/Velocity/TSX/Helheim/")
  DIR_RADARSAT = os.path.join(os.getenv("HOME"),"Data/Velocity/RADARSAT/Helheim/")
  
  ###############################################################################
  # Load terminus profiles so we can calculate lagrangian points #
  ###############################################################################

  term_values, term_time = helheim_icefronts.distance_along_flowline(xf,yf,dists,'icefront')

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

  velocities=np.zeros([m,len(pts)])
  velocities[:,:] = 'nan'
  times=np.zeros([m,1])
  positions = np.zeros([m,len(pts)])
  positions[:,:] = 'nan'
  xpts_all = np.zeros([m,len(pts)])
  xpts_all[:,:] = 'NaN'
  ypts_all = np.zeros([m,len(pts)])
  ypts_all[:,:] = 'NaN'
  error=np.zeros([m,len(pts)])
  count=0
  for j in range(0,len(DIRS)):
    DIR=DIRS[j]
    if DIR.startswith('track'):
      errors = np.zeros([m])
      infile=DIR_TSX+DIR
      x,y,v,vx,vy,ex,ey,time=geodat.readvelocity(infile+"/mosaicOffsets")
      times[count]=time
      
      # Set up grid for interpolation
      fv = scipy.interpolate.RegularGridInterpolator([y,x],v,method='nearest',bounds_error=False)
      fex = scipy.interpolate.RegularGridInterpolator([y,x],ex,method='nearest',bounds_error=False)
      fey = scipy.interpolate.RegularGridInterpolator([y,x],ey,method='nearest',bounds_error=False)
          
      # Use terminus position to determine what points along the flowline we should be estimating velocities
      terminus = np.interp(time,term_time,term_values)
      flowdists = terminus-pts
      
      xpts = np.interp(flowdists,dists,xf)
      ypts = np.interp(flowdists,dists,yf)
      
      ind = np.where((xpts > np.min(x)) & (xpts < np.max(x)) & (ypts < np.max(y)) & (ypts > np.min(y)))
      velocities[count,ind[0]] = fv(np.column_stack([ypts[ind[0]],xpts[ind]]))
      positions[count,:] = flowdists
      xpts_all[count,:] = xpts
      ypts_all[count,:] = ypts
      for i in range(0,len(ind[0])):
        error[count,ind[0][i]]=math.sqrt(fex(np.column_stack([ypts[ind[0][i]],xpts[ind[0][i]]]))**2+fey(np.column_stack([ypts[ind[0][i]],xpts[ind[0][i]]]))**2)
      
      count = count + 1
  
  # Sort arrays by time
  time_tsx = times
  velocities_tsx = velocities
  error_tsx = error
  pos_tsx = positions
  xpt_tsx = xpts_all
  ypt_tsx = ypts_all
  del velocities, time, positions, xpts_all, error, ypts_all
  
  ######################
  # Load RADARSAT data #
  ######################
  DIRS=os.listdir(DIR_RADARSAT)
  
  # Get number of velocity files
  m=0
  for DIR in DIRS:
    if DIR.startswith('winter'):
      m = m+1

  velocities=np.zeros([m,len(pts)])
  times=np.zeros([m,1])
  error=np.zeros([m,len(pts)])
  positions = np.zeros([m,len(pts)])
  positions[:,:] = 'nan'
  xpts_all = np.zeros([m,len(pts)])
  xpts_all = 'nan'
  ypts_all = np.zeros([m,len(pts)])
  ypts_all = 'nan'
  count=0
  for j in range(0,len(DIRS)):
    DIR=DIRS[j]
    if DIR.startswith('winter'):
      infile=DIR_RADARSAT+DIR
      x,y,v,vx,vy,ex,ey,time=geodat.readvelocity(infile+"/mosaicOffsets")
      times[count]=time
      
      # Set up grid for interpolation
      fv = scipy.interpolate.RegularGridInterpolator([y,x],v,method='nearest',bounds_error=False)
      fex = scipy.interpolate.RegularGridInterpolator([y,x],ex,method='nearest',bounds_error=False)
      fey = scipy.interpolate.RegularGridInterpolator([y,x],ey,method='nearest',bounds_error=False)
      
      # Use terminus position to determine what points along the flowline we should be estimating velocities
      terminus = np.interp(time,term_time,term_values)
      flowdists = terminus - pts
      
      xpts = np.interp(flowdists,dists,xf)
      ypts = np.interp(flowdists,dists,yf)
      
      ind = np.where((xpts > np.min(x)) & (xpts < np.max(x)) & (ypts < np.max(y)) & (ypts > np.min(y)))
      velocities[count,ind[0]] = fv(np.column_stack([ypts[ind[0]],xpts[ind]]))
      positions[count,:] = flowdists
      xpts_all[count,:] = xpts
      ypts_all[count,:] = ypts
      
      for i in range(0,len(ind[0])):
        error[count,ind[0][i]]=math.sqrt(fex(np.column_stack([ypts[ind[0][i]],xpts[ind[0][i]]]))**2+fey(np.column_stack([ypts[ind[0][i]],xpts[ind[0][i]]]))**2)
      
      count = count + 1
  
  time_radarsat = times
  velocities_radarsat = velocities
  error_radarsat = error
  pos_radarsat = positions
  xpt_tsx = xpts_all
  ypt_tsx = ypts_all

  del positions,error,times,velocities,xpts_all,ypts_all

  tpt_all = np.row_stack([time_radarsat,time_tsx])
  vpt_all = np.row_stack([velocities_radarsat,velocities_tsx])
  ept_all = np.row_stack([error_radarsat,error_tsx])
  dists_all = np.row_stack([pos_radarsat,pos_tsx])
  ypt_all = np.row_stack([ypt_radarsat,ypt_tsx])
  xpt_all = np.row_stack([xpt_radarsat,xpt_tsx])
  
  # Sort arrays by time  
  sortind=np.argsort(tpt_all,0)
  tpt_all = tpt_all[sortind[:,0]]
  vpt_all = vpt_all[sortind[:,0],:]
  ept_all = ept_all[sortind[:,0],:]
  xpt_all = ypt_all[sortind[:,0],:]
  xpt_all = ypt_all[sortind[:,0],:]
  dists_all = dists_all[sortind[:,0],:]
      	  
  return vpt_all,tpt_all, ept_all, dists_all, xpts_all, ypts_all 


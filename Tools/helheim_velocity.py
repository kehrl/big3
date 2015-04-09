#Various scripts for processing and pulling Helheim velocity data:

## velocity_at_eulpoints(xpt,ypt): get velocity magnitude timeseries at an Eulerian point
## divergence_at_eulpoints(xpt,ypt): get divx,divy at Eulerian points through time
## velocity_at_lagpoints(x,y,dists,pts): get velocity timeseries at a Lagrangian point
## velocity_along_flowline(x,y): get velocity timeseries along flowline
## inversion_3D(file_velocity_in, dir_velocity_out): output velocity data for ElmerSolver

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

  # Finds velocity at nearest gridcell to xpt, ypt for all velocity maps. Output is 
  # velocity at xpt, ypt through time.

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
  
  # Find velocity along flowline with coordinates xf, yf. The variable "dists" (distance 
  # along flowline) is used to determine which points to throw out in front if the ice front.
  # The output is velocity along the flowline through time.
  
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
    	  
  return vpt_sort,tpt_sort,term_sort


##########################################################################################
def velocity_at_lagpoints(xf,yf,dists,pts):
  
  # Find velocity at lagrangian points with distance "pts" behind the glacier terminus.
  # Output is velocity through time.
  
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
  xpts_all[:,:] = 'nan'
  ypts_all = np.zeros([m,len(pts)])
  ypts_all[:,:] = 'nan'
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
  xpt_radarsat = xpts_all
  ypt_radarsat = ypts_all

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
  xpt_all = xpt_all[sortind[:,0],:]
  ypt_all = ypt_all[sortind[:,0],:]
  dists_all = dists_all[sortind[:,0],:]
      	  
  return vpt_all,tpt_all, ept_all, dists_all, xpt_all, ypt_all 

#########################################################################################
def divergence_at_eulpoints(xpt,ypt):
  
  # Finds divx, divy at nearest gridcell to xpt, ypt through time. The output is vx, vy, 
  # divx, and divy.

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

  divxpt=np.zeros([m,n])
  divypt=np.zeros([m,n])
  vpt=np.zeros([m,n])
  vxpt=np.zeros([m,n])
  vypt=np.zeros([m,n])
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
          vpt[count,i] = v[yind,xind]
          vxpt[count,i] = vx[yind,xind]
          vypt[count,i] = vy[yind,xind]
          if (x[xind] > xpt[i]):
            xind=xind-1
          if (y[yind] > ypt[i]):
            yind=yind-1
        except:
      	  xind = (abs(x-xpt)).argmin()
          yind = (abs(y-ypt)).argmin()
          vpt[count,i] = v[yind,xind]
          vxpt[count,i] = vx[yind,xind]
          vypt[count,i] = vy[yind,xind]
          if (x[xind] > xpt):
            xind=xind-1
          if (y[yind] > ypt):
            yind=yind-1
        
        divxpt[count,i] = (vx[yind,xind+1]-vx[yind,xind])/(x[xind+1]-x[xind])
        divypt[count,i] = (vy[yind+1,xind]-vy[yind,xind])/(y[yind+1]-y[yind])
        
      count = count + 1
  
  # Sort arrays by time
  tpt_tsx = tpt
  vpt_tsx = vpt
  vxpt_tsx = vxpt
  vypt_tsx = vypt
  divxpt_tsx = divxpt
  divypt_tsx = divypt
  
  ######################
  # Load RADARSAT data #
  ######################
  DIRS=os.listdir(DIR_RADARSAT)
  
  # Get number of velocity files
  m=0
  for DIR in DIRS:
    if DIR.startswith('winter'):
      m = m+1

  divxpt=np.zeros([m,n])
  divypt=np.zeros([m,n])
  vpt=np.zeros([m,n])
  vxpt=np.zeros([m,n])
  vypt=np.zeros([m,n])
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
          vpt[count,i] = v[yind,xind]
          vxpt[count,i] = vx[yind,xind]
          vypt[count,i] = vy[yind,xind]
          if (x[xind] > xpt[i]):
            xind=xind-1
          if (y[yind] > ypt[i]):
            yind=yind-1
        except:
      	  xind = (abs(x-xpt)).argmin()
          yind = (abs(y-ypt)).argmin()
          vpt[count,i] = v[yind,xind]
          vxpt[count,i] = vx[yind,xind]
          vypt[count,i] = vy[yind,xind]
          if (x[xind] > xpt):
            xind=xind-1
          if (y[yind] > ypt):
            yind=yind-1

        divxpt[count,i] = (vx[yind,xind+1]-vx[yind,xind])/(x[xind+1]-x[xind])
        divypt[count,i] = (vy[yind+1,xind]-vy[yind,xind])/(y[yind+1]-y[yind])
        
      count = count + 1
  
  tpt_radarsat = tpt
  vpt_radarsat = vpt
  vxpt_radarsat = vxpt
  vypt_radarsat = vypt
  divxpt_radarsat = divxpt
  divypt_radarsat = divypt
  
  tpt_all = np.row_stack([tpt_radarsat,tpt_tsx])
  vpt_all = np.row_stack([vpt_radarsat,vpt_tsx])
  vxpt_all = np.row_stack([vxpt_radarsat,vxpt_tsx])
  vypt_all = np.row_stack([vypt_radarsat,vypt_tsx])
  divxpt_all = np.row_stack([divxpt_radarsat,divxpt_tsx])
  divypt_all = np.row_stack([divypt_radarsat,divypt_tsx])
  
  # Sort arrays by time  
  sortind=np.argsort(tpt_all,0)
  tpt_all = tpt_all[sortind[:,0]]
  vpt_all = vpt_all[sortind[:,0],:]
  vxpt_all = vxpt_all[sortind[:,0],:]
  vypt_all = vypt_all[sortind[:,0],:]
  divxpt_all = divxpt_all[sortind[:,0],:]
  divypt_all = divypt_all[sortind[:,0],:]
      	  
  return vpt_all,vxpt_all, vypt_all, divxpt_all, divypt_all, tpt_all
  
#########################################################################################

def inversion_3D(file_velocity_in,dir_velocity_out):

  # This code outputs velocity measurements for Elmer to use during inversions. It combines 
  # "file_velocity_in" with a Greenland velocity map to create a continuous product.

  # Large velocity map
  file_velocity_in_global = os.path.join(os.getenv("HOME"),"Data/Velocity/Random/Greenland/track-07to10")
  x1,y1,v1,vx1,vy1,ex1,ey1,time=geodat.readvelocity(file_velocity_in_global+"/mosaicOffsets")

  # Individual velocity map
  x2,y2,v2,vx2,vy2,ex2,ey2,time=geodat.readvelocity(file_velocity_in+"/mosaicOffsets")

  # Chop larger velocity maps to desired size
  xmin = np.argmin(abs(x1-200000))
  xmax = np.argmin(abs(x1-320000))+1
  ymin = np.argmin(abs(y1--2620000))
  ymax = np.argmin(abs(y1--2510000))+1

  x1 = x1[xmin:xmax]
  y1 = y1[ymin:ymax]
  v1 = v1[ymin:ymax,xmin:xmax]
  vx1 = vx1[ymin:ymax,xmin:xmax]
  vy1 = vy1[ymin:ymax,xmin:xmax]
  ex1 = ex1[ymin:ymax,xmin:xmax]
  ey1 = ey1[ymin:ymax,xmin:xmax]
  
  xmin = np.argmin(abs(x2-200000))
  xmax = np.argmin(abs(x2-320000))+1
  ymin = np.argmin(abs(y2--2620000))
  ymax = np.argmin(abs(y2--2510000))+1

  x2 = x2[xmin:xmax]
  y2 = y2[ymin:ymax]
  v2 = v2[ymin:ymax,xmin:xmax]
  vx2 = vx2[ymin:ymax,xmin:xmax]
  vy2 = vy2[ymin:ymax,xmin:xmax]
  ex2 = ex2[ymin:ymax,xmin:xmax]
  ey2 = ey2[ymin:ymax,xmin:xmax]

  # Add individual velocity map to final product
  xfinal = np.arange(np.min(x1),np.max(x1)+x2[1]-x2[0],x2[1]-x2[0])
  yfinal = np.arange(np.min(y1),np.max(y1)+y2[1]-y2[0],y2[1]-y2[0])
  xgridfinal,ygridfinal = np.meshgrid(xfinal,yfinal)
  vx_final = np.zeros_like(xgridfinal)
  vy_final = np.zeros_like(xgridfinal)
  v_final = np.zeros_like(xgridfinal)
  vx_final[:,:] = 'NaN'
  vy_final[:,:] = 'NaN'
  v_final[:,:] = 'NaN'
  
  xind1 = np.argmin(abs(xfinal-np.min(x2)))
  xind2 = np.argmin(abs(xfinal-np.max(x2)))+1
  yind1 = np.argmin(abs(yfinal-np.min(y2)))
  yind2 = np.argmin(abs(yfinal-np.max(y2)))+1
  
  vx_final[yind1:yind2,xind1:xind2] = vx2
  vy_final[yind1:yind2,xind1:xind2] = vy2
  v_final[yind1:yind2,xind1:xind2] = v2
  
  
  # Fill in spots in velocity map that are not covered by the individual velocity map 
  nonnan=np.where(~(np.isnan(vx1)))
  xgrid1,ygrid1 = np.meshgrid(x1,y1)
  vx_fill = scipy.interpolate.griddata(np.column_stack([ygrid1[nonnan],xgrid1[nonnan]]),vx1[nonnan],(ygridfinal,xgridfinal),method='linear')
  vy_fill = scipy.interpolate.griddata(np.column_stack([ygrid1[nonnan],xgrid1[nonnan]]),vx1[nonnan],(ygridfinal,xgridfinal),method='linear')
  v_fill = scipy.interpolate.griddata(np.column_stack([ygrid1[nonnan],xgrid1[nonnan]]),v1[nonnan],(ygridfinal,xgridfinal),method='linear')

  nans = np.where(np.isnan(vx_final))
  vx_final[nans]=vx_fill[nans]
  vy_final[nans]=vy_fill[nans]
  v_final[nans]=v_fill[nans]

  # Write out velocities to files
  fid = open(dir_velocity_out+"/UDEM.xy","w")
  fid.write('{}\n{}\n'.format(len(xfinal),len(yfinal)))
  for i in range(0,len(xfinal)):
    for j in range(0,len(yfinal)):
      fid.write('{} {} {}\n'.format(xfinal[i],yfinal[j],vx_final[j,i]))
  fid.close()
  
  fid = open(dir_velocity_out+"/VDEM.xy","w")
  fid.write('{}\n{}\n'.format(len(xfinal),len(yfinal)))
  for i in range(0,len(xfinal)):
    for j in range(0,len(yfinal)):
      fid.write('{} {} {}\n'.format(xfinal[i],yfinal[j],vy_final[j,i]))
  fid.close()

  fid = open(dir_velocity_out+"/VMAG.xy","w")
  fid.write('{}\n{}\n'.format(len(xfinal),len(yfinal)))
  for i in range(0,len(xfinal)):
    for j in range(0,len(yfinal)):
      fid.write('{} {} {}\n'.format(xfinal[i],yfinal[j],v_final[j,i]))
  fid.close()

  return 1
'''
Various scripts for processing and pulling the Helheim and Kanger velocity data:

  velocity_at_eulpoints(xpt,ypt): get velocity magnitude timeseries 
  			at an Eulerian point
  divergence_at_eulpoints(xpt,ypt): get divx,divy at Eulerian points 
  			through time
  velocity_at_lagpoints(x,y,dists,pts): get velocity timeseries at 
  			a Lagrangian point
  velocity_along_flowline(x,y): get velocity timeseries along flowline
  inversion_3D(file_velocity_in, dir_velocity_out): output velocity 
  			data for ElmerSolver
  inversion_2D(x,y,dist,file_velocity_in, dir_velocity_out): output 
  			velocity along flowline for ElmerSolver
 
LMK, UW, 04/01/2014
'''
import os
import math
import shutil
import sys
import scipy.interpolate
import numpy as np
import netCDF4
import geodatlib, icefrontlib, geotifflib, datelib, masklib, jdcal
from scipy import stats

#########################################################################################
def convert_binary_to_geotiff(glacier):

  '''
I'm getting tired of unpacking Ian's binary velocity files every time 
I need to use them, so I've set up a script to convert all of them to 
geotifflib. Then the "geodat" module checks to see if there are geotiffs 
before unpacking the binary files.
  '''

  DIRTOP_TSX = os.path.join(os.getenv("DATA_HOME"),"Velocity/TSX/"+glacier+"/")
  DIRTOP_RADARSAT = os.path.join(os.getenv("DATA_HOME"),"Velocity/RADARSAT/Greenland/")

  
  # TSX files
  files = os.listdir(DIRTOP_TSX)
  for file in files:
    if file.startswith('track'):
      print file
      # Load binary data
      x,y,v,vx,vy,vz,ex,ey,time,interval = geodatlib.readbinary(DIRTOP_TSX+file+"/mosaicOffsets",nodatavalue=-2.0e9,read_vz=True)

      year,month,day = datelib.fracyear_to_date(time)
    
      # Set up date label for geotiff file
      date = "%04d%02d%02d" % (year,month,day)
    
      # Save as geotiff
      geotifflib.write_from_grid(x,y,np.flipud(v),-2.0e9,DIRTOP_TSX+"TIF/"+file+"_"+date+"_v.tif")
      geotifflib.write_from_grid(x,y,np.flipud(vx),-2.0e9,DIRTOP_TSX+"TIF/"+file+"_"+date+"_vx.tif")
      geotifflib.write_from_grid(x,y,np.flipud(vy),-2.0e9,DIRTOP_TSX+"TIF/"+file+"_"+date+"_vy.tif")
      geotifflib.write_from_grid(x,y,np.flipud(ex),-2.0e9,DIRTOP_TSX+"TIF/"+file+"_"+date+"_ex.tif")
      geotifflib.write_from_grid(x,y,np.flipud(ey),-2.0e9,DIRTOP_TSX+"TIF/"+file+"_"+date+"_ey.tif")
      geotifflib.write_from_grid(x,y,np.flipud(vz),-2.0e9,DIRTOP_TSX+"TIF/"+file+"_"+date+"_vz.tif")

  
  # RADARSAT files    
  files = os.listdir(DIRTOP_RADARSAT)
  for file in files:
    if file.startswith('winter'):
      print file
      # Load binary data
      x,y,v,vx,vy,ex,ey,time,interval = geodatlib.readbinary(DIRTOP_RADARSAT+file+"/mosaicOffsets")
    
      # Save as geotiff

      geotifflib.write_from_grid(x,y,np.flipud(v),-2.0e9,DIRTOP_RADARSAT+"TIF/"+file+"_v.tif")
      geotifflib.write_from_grid(x,y,np.flipud(vx),-2.0e9,DIRTOP_RADARSAT+"TIF/"+file+"_vx.tif")
      geotifflib.write_from_grid(x,y,np.flipud(vy),-2.0e9,DIRTOP_RADARSAT+"TIF/"+file+"_vy.tif")
      geotifflib.write_from_grid(x,y,np.flipud(ex),-2.0e9,DIRTOP_RADARSAT+"TIF/"+file+"_ex.tif")
      geotifflib.write_from_grid(x,y,np.flipud(ey),-2.0e9,DIRTOP_RADARSAT+"TIF/"+file+"_ey.tif")
      
  return 1

#########################################################################################
def velocity_grid(glacier,xmin=-np.Inf,xmax=np.Inf,ymin=-np.Inf,ymax=np.Inf,resolution=100):

  DIR_TSX = os.path.join(os.getenv("DATA_HOME"),"Velocity/TSX/"+glacier+"/")
  
  dx = dy = float(resolution)
  nx = int(np.ceil((xmax-xmin)/dx)+1)
  x = np.linspace(xmin,(nx-1)*dx+xmin,nx)
  ny = int(np.ceil((ymax-ymin)/dx)+1)
  y = np.linspace(ymin,(ny-1)*dy+ymin,ny)
  xgrid,ygrid = np.meshgrid(x,y) 
  coords = np.column_stack([ygrid.flatten(),xgrid.flatten()])
 
  #################
  # LOAD TSX Data #
  #################

  DIRs=os.listdir(DIR_TSX)
  
  # Get number of velocity files
  nt=0
  for DIR in DIRs:
    if DIR.startswith('track'):
      nt = nt+1

  # Set up variables
  velgrid = np.zeros([ny,nx,nt])
  time = np.zeros(nt)
  ergrid = np.zeros([ny,nx,nt])

  # Load velocity and mask
  count = 0
  for j in range(0,len(DIRs)):
    DIR=DIRs[j]
    if DIR.startswith('track'):
      # Load velocity
      x1,y1,v1,vx1,vy1,ex1,ey1,time_file,interval1 = geodatlib.readvelocity(DIR_TSX,DIR,"mosaicOffsets")
      
      time[count] = time_file
      year,month,day = datelib.fracyear_to_date(time_file)
      
      xind1 = np.argmin(abs(x1-xmin))
      xind2 = np.argmin(abs(x1-xmax))+1
      yind1 = np.argmin(abs(y1-ymin))
      yind2 = np.argmin(abs(y1-ymax))+1
      
      # Load velocity
      try:
        # If the input and output grids have the same dimensions...
        velgrid[:,:,count] = v1[yind1:yind2,xind1:xind2]
      except:
        # Otherwise interpolate onto output grid
        f_dem = scipy.interpolate.RegularGridInterpolator([y1,x1],v1,bounds_error = False,method='linear',fill_value=float('nan'))
        v_flatten = f_dem(coords)
    
        # Reshape to grid
        velgrid[:,:,count] = np.reshape(v_flatten,(ny,nx))
          
      count = count+1 
  
  # Sort velocities
  sortind = np.argsort(time)
  time = time[sortind]
  velgrid = velgrid[:,:,sortind]

  return x,y,velgrid,time
    
#########################################################################################
def velocity_at_eulpoints(xpt,ypt,glacier,data='all',xy_velocities='False'):

  '''
  Finds velocity at nearest gridcell to xpt, ypt for all 
  velocity maps. Output is velocity at xpt, ypt through time.

  vpt_sort,tpt_sort,ept_sort = velocity_at_eulpoints(xpt,ypt,glacier,data='all',xy_velocities='False')

  # Inputs:
  # xpt, ypt: coordinates of flowline
  # glacier: glacier name
  # data: TSX, Radarsat, or all data
  # xy_velocities: True or False; do you want the x,y velocities too?
  '''
  # Select data type
  if data == 'all':
    data = ['TSX','RADARSAT']
  elif data == 'RADARSAT':
    data = ['RADARSAT']
  elif data == 'TSX':
    data = ['TSX']
  else:
    print "Unknown data type"

  ###################
  # Load velocities #
  ###################
  
  # Find velocity files to be imported
  files = []
  dirs = []
  for type in data:
    if type == 'RADARSAT':
      DIRTOP = os.path.join(os.getenv("DATA_HOME"),"Velocity/RADARSAT/Greenland/")
    elif type == 'TSX':
      DIRTOP = os.path.join(os.getenv("DATA_HOME"),"Velocity/TSX/"+glacier+"/")
    

    DIRs=os.listdir(DIRTOP)
    for DIR in DIRs:
      if DIR.startswith('track') or DIR.startswith('winter'):
        files.append(DIR)
        dirs.append(DIRTOP)
  
  # Load velocities
  m = len(files)
  try:
    n = len(xpt)
  except:
    n = 1
  
  # Set up variables
  velocities = np.zeros([m,n])
  velocities[:,:] = 'nan'
  velocities_x = np.zeros([m,n])
  velocities_x[:,:] = 'nan'
  velocities_y = np.zeros([m,n])
  velocities_y[:,:] = 'nan'
  error = np.zeros([m,n])
  error[:,:] = 'nan'
  times=np.zeros([m])
  termini=np.zeros([m])
  count=0
  
  for i in range(0,len(files)):
    x,y,v,vx,vy,ex,ey,time,interval = geodatlib.readvelocity(''.join(dirs[i]),''.join(files[i]),"mosaicOffsets")
    if 'winter' in files[i]:
      time = float('20'+files[i][-2:])
    times[count]=time
      
    # Set up grid for interpolation
    fv = scipy.interpolate.RegularGridInterpolator([y,x],v,method='linear',bounds_error=False)
    fex = scipy.interpolate.RegularGridInterpolator([y,x],ex,method='linear',bounds_error=False)
    fey = scipy.interpolate.RegularGridInterpolator([y,x],ey,method='linear',bounds_error=False)
    fvx = scipy.interpolate.RegularGridInterpolator([y,x],vx,method='linear',bounds_error=False)
    fvy = scipy.interpolate.RegularGridInterpolator([y,x],vy,method='linear',bounds_error=False)
      
    # Find velocities 
    velocities[count,:] = fv(np.array([ypt,xpt]).T) 
    velocities_x[count,:] = fvx(np.array([ypt,xpt]).T)
    velocities_y[count,:] = fvy(np.array([ypt,xpt]).T)
    error[count,:] = velocities[count,:]*np.sqrt((fex(np.array([ypt,xpt]).T)/fvx(np.array([ypt,xpt]).T))**2+(fey(np.array([ypt,xpt]).T)/fvy(np.array([ypt,xpt]).T))**2)        
    
    count = count + 1
    
  # Sort arrays by time  
  sortind=np.argsort(times,0)
  tpt_sort = times[sortind]
  vpt_sort = velocities[sortind,:]
  ept_sort = error[sortind,:]
  vxpt_sort = velocities_x[sortind,:]
  vypt_sort = velocities_y[sortind,:]
  
  if xy_velocities == 'True':
    return vpt_sort,tpt_sort,ept_sort,vxpt_sort,vypt_sort
  else:  	  
    return vpt_sort,tpt_sort,ept_sort
  
#########################################################################################
def velocity_along_flowline(xf,yf,dists,glacier,cutoff='terminus',data='all'):      
  
  '''
  Find velocity along flowline with coordinates xf, yf. The variable "dists" (distance 
  along flowline) is used to determine which points to throw out in front if the ice front.
  The output is velocity along the flowline through time.
  
  Inputs:
  xf,yf: flowline positions
  dists: distances along flowline
  glacier: Kanger or Helheim
  cutoff: cutoff velocities in front of terminus if set to 'terminus'
  data: 'all' data, or just 'TSX' or 'RADARSAT' 
  
  Outputs:
  vpt_sort: velocities along flowline
  tpt_sort: time for velocities
  term_sort: interpolated terminus positions for velocities
  '''
  
  if data == 'all':
    data = ['TSX','RADARSAT']
  elif data == 'RADARSAT':
    data = ['RADARSAT']
  elif data == 'TSX':
    data = ['TSX']
  else:
    print "Unknown data type"
  
  ###############################################################################
  # Load terminus profiles so we can cutoff velocities in front of the terminus #
  ###############################################################################

  term_values, term_time = icefrontlib.distance_along_flowline(xf,yf,dists,glacier,type='icefront')
  
  ###################
  # Load velocities #
  ###################
  
  # Find velocity files to be imported
  files = []
  dirs = []
  for type in data:
    if type == 'RADARSAT':
      DIRTOP = os.path.join(os.getenv("DATA_HOME"),"Velocity/RADARSAT/Greenland/")
    elif type == 'TSX':
      DIRTOP = os.path.join(os.getenv("DATA_HOME"),"Velocity/TSX/"+glacier+"/")
       
    DIRs=os.listdir(DIRTOP)
    for DIR in DIRs:
      if DIR.startswith('track') or DIR.startswith('winter'):
        files.append(DIR)
        dirs.append(DIRTOP)
  
  # Load velocities
  m = len(files)
  n = len(xf)
  
  velocities=np.zeros([m,n])
  velocities[:,:] = 'nan'
  times=np.zeros(m)
  termini=np.zeros(m)
  count=0
  for i in range(0,len(files)):
    x,y,v,vx,vy,ex,ey,time,interval = geodatlib.readvelocity(''.join(dirs[i]),''.join(files[i]),"mosaicOffsets")
    if 'winter' in files[i]:
      time = float('20'+files[i][-2:])
    times[count]=time
      
    # Set up grid for interpolation
    fv = scipy.interpolate.RegularGridInterpolator([y,x],v,method='linear',bounds_error=False)
     
    # Find velocities
    velocities[count,:] = fv(np.array([yf,xf]).T) 
      
    # Find flowline coordinates behind terminus 
    terminus = np.interp(time,term_time,term_values)  
      
    if cutoff == 'terminus':
      ind = np.where(dists > terminus)
      velocities[count,ind] = 'nan'
      
    termini[count]=terminus
      
    count = count + 1
    
  # Sort arrays by time  
  sortind=np.argsort(times,0)
  tpt_sort = times[sortind]
  term_sort = termini[sortind]
  vpt_sort = velocities[sortind,:].T  	  
    	  
  # Print warning if removing points in front of ice front
  if cutoff == 'terminus':
    print "Cleaning up velocity points by removing points in front of ice front."
    print "You can change this setting by setting `cutoff = 'none'.'"    	  
    	  
  return vpt_sort,tpt_sort,term_sort


##########################################################################################
def velocity_at_lagpoints(xf,yf,dists,pts,glacier,data='all'):
  
  # Find velocity at lagrangian points with distance "pts" behind (or in front) of the 
  # glacier terminus.
  # Output is velocity through time.
  
  # Select data type
  if data == 'all':
    data = ['TSX','RADARSAT']
  elif data == 'RADARSAT':
    data = ['RADARSAT']
  elif data == 'TSX':
    data = ['TSX']
  else:
    print "Unknown data type"

  ###########################
  # Load terminus positions #
  ###########################
  
  term_values, term_time = icefrontlib.distance_along_flowline(xf,yf,dists,glacier,type='icefront')
  
  ###################
  # LOAD velocities #
  ###################
  
  # Find velocity files to be imported
  files = []
  dirs = []
  for type in data:
    if type == 'RADARSAT':
      DIRTOP = os.path.join(os.getenv("DATA_HOME"),"Velocity/RADARSAT/Greenland/")
    elif type == 'TSX':
      DIRTOP = os.path.join(os.getenv("DATA_HOME"),"Velocity/TSX/"+glacier+"/")
    
    DIRs=os.listdir(DIRTOP)
    for DIR in DIRs:
      if DIR.startswith('track') or DIR.startswith('winter'):
        files.append(DIR)
        dirs.append(DIRTOP)
  
  # Load velocities
  m = len(files)
  try:
    n = len(pts)
  except:
    n = 1
  
  velocities = np.zeros([m,n])
  velocities[:,:] = 'nan'
  positions = np.zeros([m,n])
  positions[:,:] = 'nan'
  xpts_all = np.zeros([m,n])
  ypts_all = np.zeros([m,n])
  error = np.zeros([m,n])
  error[:,:] = 'nan'
  times = np.zeros([m])
  intervals = np.zeros([m])
  termini = np.zeros([m])
  count=0
  for i in range(0,len(files)):
    x,y,v,vx,vy,ex,ey,time,interval = geodatlib.readvelocity(''.join(dirs[i]),''.join(files[i]),"mosaicOffsets")
    if 'winter' in files[i]:
      time = float('20'+files[i][-2:])
    times[count]=time
    intervals[count] = interval
      
    # Set up grid for interpolation
    fv = scipy.interpolate.RegularGridInterpolator([y,x],v,method='linear',bounds_error=False)
    fvx = scipy.interpolate.RegularGridInterpolator([y,x],vx,method='linear',bounds_error=False)
    fvy = scipy.interpolate.RegularGridInterpolator([y,x],vy,method='linear',bounds_error=False)
    fex = scipy.interpolate.RegularGridInterpolator([y,x],ex,method='linear',bounds_error=False)
    fey = scipy.interpolate.RegularGridInterpolator([y,x],ey,method='linear',bounds_error=False)
      
    # Find terminus position
    terminus = np.interp(time,term_time,term_values)
    flowdists = terminus+pts
      
    xpts = np.interp(flowdists,dists,xf)
    ypts = np.interp(flowdists,dists,yf)
    
    # Find velocities 
    velocities[count,:] = fv(np.array([ypts,xpts]).T)  
    error[count,:] = velocities[count,:]*np.sqrt((fex(np.array([ypts,xpts]).T)/fvx(np.array([ypts,xpts]).T))**2+(fey(np.array([ypts,xpts]).T)/fvy(np.array([ypts,xpts]).T))**2)        
    
     
    positions[count,:] = flowdists
    xpts_all[count,:] = xpts
    ypts_all[count,:] = ypts

    count = count + 1
  
  # Sort arrays by time  
  sortind=np.argsort(times,0)
  tpt_all = np.column_stack([times[sortind],intervals[sortind]])
  vpt_all = velocities[sortind,:]
  ept_all = error[sortind,:]
  xpt_all = xpts_all[sortind,:]
  ypt_all = ypts_all[sortind,:]
  dists_all = positions[sortind,:]
      	  
  return vpt_all, tpt_all, ept_all, dists_all, xpt_all, ypt_all 

#########################################################################################
def tsx_near_time(time,glacier,just_filename = False):

  '''
  
  x,y,vx,vy,v,time = tsx_near_time(time,glacier,just_filename = False)
  
  Find TSX data closest to "time".
  
  Inputs:
  time: time that you want data
  glacier: glacier name (Kanger or Helheim)
  just_filename: option to only return the filename

  Outputs:
  x,y: grid coordinates
  vx,vy: x and y velocities
  v: velocity magnitudes for gird
  time: time of transect

  '''

  DIR_TSX = os.path.join(os.getenv("DATA_HOME"),"Velocity/TSX/"+glacier+"/")

  DIRs=os.listdir(DIR_TSX)
  tpt=[]
  
  best_track = []
  min_diff = 1.0
  for DIR in DIRs:
    if DIR.startswith('track'):
      tsx_time,interval = geodatlib.readtime(DIR_TSX+DIR+"/mosaicOffsets")
      if abs(tsx_time-time) < min_diff:
        min_diff = abs(tsx_time-time)
        best_time = tsx_time
        best_track = DIR

  
  if just_filename:
    year,month,day = datelib.fracyear_to_date(best_time)
    
    return DIR_TSX+'TIF/'+best_track+'_'+"%04d%02d%02d" % (year,month,day),best_time
  
  else:
    # Return the closest velocity profile
    x,y,v,vx,vy,ex,ey,time,interval = geodatlib.readvelocity(DIR_TSX,best_track,"/mosaicOffsets")
    
    return x,y,vx,vy,v,time


#########################################################################################
def variability(glacier,time1,time2):

  ''
  ''

  DIR_TSX = os.path.join(os.getenv("DATA_HOME"),"Velocity/TSX/"+glacier+"/")

  if glacier == 'Helheim':
    xmin = 270000.0
    xmax = 354900.0
    ymin = -2601000.0
    ymax = -2541000.0
  elif glacier == 'Kanger':
    xmin = 457000.0
    xmax = 517000.0
    ymin = -2319100.0
    ymax = -2247100.0
  
  dx = dy = 100.
  nx = int(np.ceil((xmax-xmin)/dx)+1)
  x = np.linspace(xmin,(nx-1)*dx+xmin,nx)
  ny = int(np.ceil((ymax-ymin)/dx)+1)
  y = np.linspace(ymin,(ny-1)*dy+ymin,ny)
  xgrid,ygrid = np.meshgrid(x,y) 
  coords = np.column_stack([ygrid.flatten(),xgrid.flatten()])
 
  #################
  # LOAD TSX Data #
  #################

  DIRs=os.listdir(DIR_TSX)
  
  # Get number of velocity files
  nt=0
  for DIR in DIRs:
    if DIR.startswith('track'):
      nt = nt+1

  # Set up variables
  velgrid = np.zeros([ny,nx,nt])
  mask = np.zeros([ny,nx,nt])
  velgrid_mask = np.zeros([ny,nx,nt])
  time = np.zeros(nt)
  ergrid = np.zeros([ny,nx,nt])

  # Load velocity and mask
  count = 0
  for j in range(0,len(DIRs)):
    DIR=DIRs[j]
    if DIR.startswith('track'):
      # Load velocity
      x1,y1,v1,vx1,vy1,ex1,ey1,time_file,interval1 = geodatlib.readvelocity(DIR_TSX,DIR,"mosaicOffsets")
      
      time[count] = time_file
      year,month,day = datelib.fracyear_to_date(time_file)
      
      xind1 = np.argmin(abs(x1-xmin))
      xind2 = np.argmin(abs(x1-xmax))+1
      yind1 = np.argmin(abs(y1-ymin))
      yind2 = np.argmin(abs(y1-ymax))+1
      
      # Load velocity
      try:
        # If the input and output grids have the same dimensions...
        velgrid[:,:,count] = v1[yind1:yind2,xind1:xind2]
      except:
        # Otherwise interpolate onto output grid
        f_dem = scipy.interpolate.RegularGridInterpolator([y1,x1],v1,bounds_error = False,method='linear',fill_value=float('nan'))
        v_flatten = f_dem(coords)
    
        # Reshape to grid
        velgrid[:,:,count] = np.reshape(v_flatten,(ny,nx))
    
      # Load mask
      date = "%04d%02d%02d" % (year,month,day)
      maskfile = DIR_TSX+'TIF/'+DIR+'_'+date+'_'+'mask.tif'
      if os.path.isfile(maskfile):
        xmask,ymask,mask[:,:,count] = geotifflib.read(maskfile)
      else:
        xmask,ymask,mask[:,:,count] = masklib.load_grid(glacier,xmin,xmax,ymin,ymax,dx,icefront_time=time1)
        geotifflib.write_from_grid(xmask,ymask,np.flipud(mask[:,:,count]),float('nan'),maskfile)
      
      velgrid_mask[:,:,count] = np.array(velgrid[:,:,count])
      velgrid_mask[mask[:,:,count]==1,count] = float('nan')
      
      count = count+1 
  
  del count,maskfile,date,xind1,yind1,xind2,yind2,year,month,x1,y1,vx1,vy1,ex1,ey1,time_file,interval1
  
  # Throw out obvious outliers
  ind = np.where(velgrid > 16.0e3)
  velgrid[ind[0],ind[1],ind[2]] = float('nan')
  velgrid_mask[ind[0],ind[1],ind[2]] = float('nan')
  print "Throwing out velocities above 16 km/yr to deal with outliers in Kanger record"

  # Only keep data that falls between time1 and time2, and sort that data by time
  sortind = np.argsort(time)
  time = time[sortind]
  velgrid_mask = velgrid_mask[:,:,sortind]
  velgrid = velgrid[:,:,sortind]
  
  ind = np.where((time > time1) & (time < time2))[0]
  velgrid_mask = velgrid_mask[:,:,ind]
  time = time[ind]
  velgrid = velgrid[:,:,ind]

  # Get average and std values
  velmean = np.nanmean(velgrid_mask,axis=2)
  
  # Get linear trends
  veltrend = np.zeros_like(velmean)
  veltrend_time1 = np.zeros_like(velmean)
  veltrend_time2 = np.zeros_like(velmean)
  veltrend_count = np.zeros_like(velmean)
  veltrend_p = np.zeros_like(velmean)
  veltrend_error = np.zeros_like(velmean)
  veltrend_r = np.zeros_like(velmean)
  veltrend_intercept = np.zeros_like(velmean)
  veltrend_p[:,:] = float('nan')
  veltrend[:,:] = float('nan')
  veltrend_error[:,:] = float('nan')
  veltrend_r[:,:] = float('nan')
  veltrend_intercept[:,:] = float('nan')
  for j in range(0,len(y)):
    for i in range(0,len(x)):
      nonnan = np.where((~(np.isnan(velgrid_mask[j,i,:]))))[0]
      if len(nonnan) > 0.75*len(time):
        if (np.floor(np.min(time[nonnan]))==time1) and np.ceil(np.max(time[nonnan]))==time2:
          slope,intercept,r,p,std_err = stats.linregress(time[nonnan],velgrid_mask[j,i,nonnan])
          veltrend_count[j,i] = len(nonnan)
          veltrend[j,i] = slope
          veltrend_p[j,i] = p
          veltrend_error[j,i] = std_err
          veltrend_time1[j,i] = np.min(time[nonnan])
          veltrend_time2[j,i] = np.max(time[nonnan])
          veltrend_r[j,i] = r
          veltrend_intercept[j,i] = intercept
        
  # Detrend velocity timeseries      
  veldetrend = np.zeros_like(velgrid_mask)
  for i in range(0,len(time)):
    trend = veltrend_intercept+time[i]*veltrend
    veldetrend[:,:,i] = velgrid_mask[:,:,i]-trend
  
  # Calculate range of observed values
  velrange = np.zeros_like(velmean)
  velrange[:,:] = float('nan')
  for i in range(0,len(x)):
    for j in range(0,len(y)):
      nonnan = np.where(~(np.isnan(veldetrend[j,i,:])))[0]
      if len(nonnan) > 1:
        velrange[j,i] = np.max(veldetrend[j,i,nonnan])-np.min(veldetrend[j,i,nonnan])
  
  # Remove insignifcant trends      
  ind = np.where(veltrend_p > 0.05)
  veltrend[ind] = float('nan')
  veltrend_error[ind] = float('nan')

  # Get number of nonnan velocities for each pixel
  velcount = np.zeros([ny,nx])
  for j in range(0,ny):
    for i in range(0,nx):
      nonnan = len(np.where(~(np.isnan(velgrid_mask[j,i,:])))[0])
      velcount[j,i] = nonnan          
  
  sortind = np.argsort(time)
  velgrid_mask = velgrid_mask[:,:,sortind]
  time = time[sortind]
  
  return x,y,velgrid_mask,veltrend,veldetrend,velrange,velcount,veltrend_error,time

#########################################################################################
def divergence_at_eulpoints(xpt,ypt):
  
  # Finds divx, divy at nearest gridcell to xpt, ypt through time. The output is vx, vy, 
  # divx, and divy.

  try:
    n = len(xpt)
  except:
    n = 0

  DIR_TSX = os.path.join(os.getenv("DATA_HOME"),"Velocity/TSX/"+glacier+"/")
  DIR_RADARSAT = os.path.join(os.getenv("DATA_HOME"),"Velocity/RADARSAT/Greenland/")

  #################
  # LOAD TSX Data #
  #################

  DIRs=os.listdir(DIR_TSX)
  tpt=[]
  
  # Get number of velocity files
  m=0
  for DIR in DIRs:
    if DIR.startswith('track'):
      m = m+1

  divxpt=np.zeros([m,n])
  divypt=np.zeros([m,n])
  vpt=np.zeros([m,n])
  vxpt=np.zeros([m,n])
  vypt=np.zeros([m,n])
  tpt=np.zeros([m,1])
  count=0
  for j in range(0,len(DIRs)):
    DIR=DIRs[j]
    if DIR.startswith('track'):
      x,y,v,vx,vy,ex,ey,time,interval = geodatlib.readvelocity(DIR_TSX,DIR,"mosaicOffsets")
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
  DIRs=os.listdir(DIR_RADARSAT)
  
  # Get number of velocity files
  m=0
  for DIR in DIRs:
    if DIR.startswith('winter'):
      m = m+1

  divxpt=np.zeros([m,n])
  divypt=np.zeros([m,n])
  vpt=np.zeros([m,n])
  vxpt=np.zeros([m,n])
  vypt=np.zeros([m,n])
  tpt=np.zeros([m,1])
  count=0
  for j in range(0,len(DIRs)):
    DIR=DIRs[j]
    if DIR.startswith('winter'):
      print "Loading ",dir
      infile=DIR_RADARSAT+DIR
      x,y,v,vx,vy,ex,ey,time,interval = geodatlib.readvelocity(DIR_RADARSAT,DIR,"mosaicOffsets")
      tpt[count] = float('20'+DIR[9:])
      
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

def inversion_3D(glacier,x,y,time,dir_velocity_out='none',blur=False,dx='none'):

  '''
  Inputs:
  x : list of x coordinates for grid interpolation
  y : list of y coordinates for grid interpolation
  time : primary time for velocities, which will be filled in with other data
  file_velocity_in : velocity file for interpolation
  dir_velocity_out : directory for outputting the velocity
  
  Outputs:
  u : velocity in x-dir on grid defined by x,y
  y : velocity in y-dir on grid defined by x,y
  '''
  
  xmin = np.min(x)-5.0e3
  xmax = np.max(x)+5.0e3
  ymin = np.min(y)-5.0e3
  ymax = np.max(y)+5.0e3
  
  OUTDIR = os.path.join(os.getenv("DATA_HOME"),"Velocity/MosaicVelocities/"+glacier)
  
  # Large velocity map to fill in gaps in smaller velocity map
  file_velocity_all = os.path.join(os.getenv("DATA_HOME"),"Velocity/TSX/"+glacier+"/TIF/all-2008-2016")

  # If the region is bigger than what is covered by the TSX stripmaps, then we use Ian's big 
  # inSAR velocity map
  file_velocity_global = os.path.join(os.getenv("DATA_HOME"),"Velocity/Random/Greenland/AllGLVel/mosaicOffsets")

  year,month,day = datelib.fracyear_to_date(time)
  date = "%04d%02d%02d" % (year,month,day)

  # Individual velocity map for time step
  if time <= 2008:
    HOWATDIR = os.path.join(os.getenv("DATA_HOME"),"Velocity/Howat/"+glacier+"/")
    # Use Howat velocity maps
    print date
    if glacier == 'Kanger':
      if '200707' in date:
        filename1 = HOWATDIR+"OPT_E68.80N_2007-08/OPT_E68.80N_2007-08"
        filename2 = HOWATDIR+"OPT_E68.80N_2007-07/OPT_E68.80N_2007-07"
      elif '200107' in date:
        filename1 = HOWATDIR+"OPT_E68.80N_2001-07/OPT_E68.80N_2001-07"
        filename2 = filename1
      elif ('200308' in date) or ('200307' in date):
        filename1 = HOWATDIR+"OPT_E68.80N_2003-07/OPT_E68.80N_2003-07"
        filename2 = HOWATDIR+"OPT_E68.80N_2001-07/OPT_E68.80N_2001-07"
      elif '200506' in date:
        filename1 = HOWATDIR+"OPT_E68.80N_2005-06/OPT_E68.80N_2005-06"
        filename2 = filename1
      elif '200508' in date:
        filename1 = HOWATDIR+"OPT_E68.80N_2005-08/OPT_E68.80N_2005-08"
        filename2 = filename1
      elif '200605' in date:
        filename2 = HOWATDIR+"OPT_E68.80N_2006-05/OPT_E68.80N_2006-05"
        filename1 = HOWATDIR+"OPT_E68.80N_2006-04/OPT_E68.80N_2006-04"
      elif '200607' in date:
        filename2 = HOWATDIR+"OPT_E68.80N_2006-07/OPT_E68.80N_2006-07"
        filename1 = HOWATDIR+"OPT_E68.80N_2006-06/OPT_E68.80N_2006-06"
      elif '200609' in date:
        filename1 = HOWATDIR+"OPT_E68.80N_2006-09/OPT_E68.80N_2006-09"
        filename2 = filename1
    elif glacier == 'Helheim':
      if '200709' in date:
        filename1 = HOWATDIR+"OPT_E66.50N_2007-08/OPT_E66.50N_2007-08"
        filename2 = HOWATDIR+"OPT_E66.50N_2007-09/OPT_E66.50N_2007-09"
      elif '200408' in date:
        filename1 = HOWATDIR+"OPT_E66.50N_2004-08/OPT_E66.50N_2004-08"
        filename2 = HOWATDIR+"OPT_E66.50N_2004-07/OPT_E66.50N_2004-07"
      elif '200508' in date:
        filename1 = HOWATDIR+"OPT_E66.50N_2005-09/OPT_E66.50N_2005-09"
        filename2 = HOWATDIR+"OPT_E66.50N_2005-07/OPT_E66.50N_2005-07"
      elif '200608' in date:
        filename1 = HOWATDIR+"OPT_E66.50N_2006-09/OPT_E66.50N_2006-09"
        filename2 = HOWATDIR+"OPT_E66.50N_2006-07/OPT_E66.50N_2006-07"

    files_vx = ' '+filename1+'.vx.tif '+filename2+'.vx.tif '+\
            file_velocity_all+'_vx.tif'+' '+file_velocity_global+'_vx.tif'
    files_vy = ' '+filename1+'.vy.tif '+filename2+'.vy.tif '+\
            file_velocity_all+'_vy.tif'+' '+file_velocity_global+'_vy.tif'
  else:
    # Use TSX velocity maps  
    filename1,time1 = tsx_near_time(time,glacier,just_filename=True)
    filename2,time2 = tsx_near_time(time-11/365.,glacier,just_filename=True)
    filename3,time3 = tsx_near_time(time+11/365.,glacier,just_filename=True)
  
    if abs(time-time2) < abs(time-time3): 
      files_vx = ' '+filename1+'_vx.tif'+' '+filename2+'_vx.tif'+\
  	' '+filename3+'_vx.tif'+' '+file_velocity_all+'_vx.tif'+' '+file_velocity_global+'_vx.tif'
      files_vy = ' '+filename1+'_vy.tif'+' '+filename2+'_vy.tif'+\
  	' '+filename3+'_vy.tif'+' '+file_velocity_all+'_vy.tif'+' '+file_velocity_global+'_vy.tif'
    else:
      files_vx = ' '+filename1+'_vx.tif'+' '+filename3+'_vx.tif'+\
        ' '+filename2+'_vx.tif'+' '+file_velocity_all+'_vx.tif'+' '+file_velocity_global+'_vx.tif'
      files_vy = ' '+filename1+'_vy.tif'+' '+filename3+'_vy.tif'+\
        ' '+filename2+'_vy.tif'+' '+file_velocity_all+'_vy.tif'+' '+file_velocity_global+'_vy.tif'
  

  CURRENTDIR = os.getcwd()
  os.chdir(OUTDIR)
  filename_vx = 'mosaic-'+date+'-vx'
  filename_vy = 'mosaic-'+date+'-vy'
  if dx == 'none':
    os.system('dem_mosaic --hole-fill-length 5 --t_projwin '+str(xmin)+' '+str(ymin)+' '+str(xmax)+\
  		' '+str(ymax)+' --priority-blending-length 10 -o'+filename_vx+files_vx)
    os.system('dem_mosaic --hole-fill-length 5 --t_projwin '+str(xmin)+' '+str(ymin)+' '+str(xmax)+\
  		' '+str(ymax)+' --priority-blending-length 10 -o'+filename_vy+files_vy) 
  else:
    os.system('dem_mosaic --hole-fill-length 5 --t_projwin '+str(xmin)+' '+str(ymin)+' '+str(xmax)+\
  		' '+str(ymax)+' --tr '+str(dx)+' --priority-blending-length 10 -o'+filename_vx+files_vx)
    os.system('dem_mosaic --hole-fill-length 5 --t_projwin '+str(xmin)+' '+str(ymin)+' '+str(xmax)+\
  		' '+str(ymax)+' --tr '+str(dx)+' --priority-blending-length 10 -o'+filename_vy+files_vy) 
  
  xu,yu,uu = geotifflib.read(filename_vx+"-tile-0.tif")
  xv,yv,vv = geotifflib.read(filename_vy+"-tile-0.tif")

  if (blur == True) and (dx == 'none'):
    print "Blurring DEM over 17 pixels (roughly 1.5km in each direction)..."
    # 17 pixel gaussian blur
    vx_blur = scipy.ndimage.filters.gaussian_filter(uu,sigma=2,truncate=4)
    vy_blur = scipy.ndimage.filters.gaussian_filter(vv,sigma=2,truncate=4)
  else:
    vx_blur = uu
    vy_blur = vv
  
  os.chdir(CURRENTDIR)
  
  # Calculate velocity magnitude
  vmag = np.sqrt(vx_blur**2+vy_blur**2)
  
  ######################################################
  # Write out velocities to files for inversion solver #
  ######################################################

  # Interpolate to input grid
  xgrid,ygrid = np.meshgrid(x,y)
  fu = scipy.interpolate.RegularGridInterpolator((yu,xu),vx_blur,method='linear')
  vx = fu((ygrid,xgrid))
  fv = scipy.interpolate.RegularGridInterpolator((yv,xv),vy_blur,method='linear')
  vy = fv((ygrid,xgrid))
  
  if dir_velocity_out != 'none':
    #files = os.listdir(OUTDIR):
    #for file in files:
    #  if file.startswith('mosaic-'+date):
    #    shutil.copy(file,dir_velocity_out)
    
    # File for velocity in x-dir
    fidu = open(dir_velocity_out+"/udem.xy","w")
    fidu.write('{}\n{}\n'.format(len(x),len(y)))
  
    # File for velocity in y-dir
    fidv = open(dir_velocity_out+"/vdem.xy","w")
    fidv.write('{}\n{}\n'.format(len(x),len(y)))
  
    for i in range(0,len(x)):
      for j in range(0,len(y)):
        fidu.write('{} {} {}\n'.format(x[i],y[j],vx[j,i]))
        fidv.write('{} {} {}\n'.format(x[i],y[j],vy[j,i]))
  
    fidv.close()
    fidu.close()
  
  return vx,vy

def inversion_2D(x,y,d,glacier,time,dir_velocity_out,filt_len='none'):
  
  xmin = np.min(x)-2.0e3
  xmax = np.max(x)+2.0e3
  ymin = np.min(y)-2.0e3
  ymax = np.max(y)+2.0e3
  
  OUTDIR = os.path.join(os.getenv("DATA_HOME"),"Velocity/MosaicVelocities/"+glacier)
  
  # Large velocity map to fill in gaps in smaller velocity map
  file_velocity_global = os.path.join(os.getenv("DATA_HOME"),"Velocity/Random/Greenland/AllGLVel/mosaicOffsets")

  # Individual velocity map
  filename1,time1 = tsx_near_time(time,glacier,just_filename=True)
  filename2,time2 = tsx_near_time(time-0.1,glacier,just_filename=True)
  filename3,time3 = tsx_near_time(time+0.1,glacier,just_filename=True)
  year,month,day = datelib.fracyear_to_date(time1)
  date = "%04d%02d%02d" % (year,month,day)
  
  files_vx = ' '+filename1+'_vx.tif'+' '+filename2+'_vx.tif'+\
  		' '+filename3+'_vx.tif'+' '+file_velocity_global+'_vx.tif'
  files_vy = ' '+filename1+'_vy.tif'+' '+filename2+'_vy.tif'+\
  		' '+filename3+'_vy.tif'+' '+file_velocity_global+'_vy.tif'
  
  CURRENTDIR = os.getcwd()
  os.chdir(OUTDIR)
  filename_vx = 'mosaic-'+date+'-vx'
  if not(os.path.isfile(filename_vx+'-tile-0.tif')):
    os.system('dem_mosaic --t_projwin '+str(xmin)+' '+str(ymin)+' '+str(xmax)+\
  		' '+str(ymax)+' --priority-blending-length 10 -o'+filename_vx+files_vx)
  filename_vy = 'mosaic-'+date+'-vy'
  if not(os.path.isfile(filename_vy+'-tile-0.tif')):
    os.system('dem_mosaic --t_projwin '+str(xmin)+' '+str(ymin)+' '+str(xmax)+\
  		' '+str(ymax)+' --priority-blending-length 10 -o'+filename_vy+files_vy)
  
  xu,yu,uu = geotifflib.read(filename_vx+"-tile-0.tif")
  xv,yv,vv = geotifflib.read(filename_vy+"-tile-0.tif")  
  
  fu = scipy.interpolate.RegularGridInterpolator((yu,xu),uu)
  vx = fu((y,x))
  fv = scipy.interpolate.RegularGridInterpolator((yv,xv),vv)
  vy = fv((y,x))
  
  vnonnan = np.sqrt(vx**2+vy**2)
  
  # Filter velocities
  if filt_len != 'none':
    cutoff=(1/filt_len)/(1/(np.diff(d[1:3])*2))
    b,a=scipy.signal.butter(4,cutoff,btype='low')
    filtered=scipy.signal.filtfilt(b,a,vnonnan)
  else:
    filtered = np.array(vcomb)

  # Write out the velocity data
  fid = open(dir_velocity_out+"velocity.dat",'w')
  R=len(filtered)
  fid.write('{0}\n'.format(R))
  for j in range(0,R):
    fid.write('{} {}\n'.format(d[j],filtered[j]))
  fid.close()

  return filtered

def rosenau_landsat_at_pts(xpt,ypt,glacier,xy_velocities='False'):
  
  if glacier == 'Helheim':
    file = os.path.join(os.getenv("DATA_HOME"),"Velocity/Rosenau/Helheim/GRL_003_all.EPSG3413.vel_md.nc")
  elif glacier == 'Kanger':
    file = os.path.join(os.getenv("DATA_HOME"),"Velocity/Rosenau/Kanger/GRL_004_all.EPSG3413.vel_md.nc")

  data = netCDF4.Dataset(file)
  x = data.variables['x'][:]
  y = data.variables['y'][:]

  # Select data to load
  i1 = np.argmin(abs(np.min(xpt)-x-1e3))
  i2 = np.argmin(abs(np.max(xpt)-x+1e3))
  j1 = np.argmin(abs(np.min(ypt)-y-1e3))
  j2 = np.argmin(abs(np.max(ypt)-y+1e3))
  x = x[i1:i2]
  y = y[j1:j2]

  vx = data.variables['vx'][:,j1:j2,i1:i2]
  vy = data.variables['vy'][:,j1:j2,i1:i2]
  v = np.sqrt(vx**2+vy**2)
  time = datelib.date_to_fracyear(1970,1,1)+data.variables['time'][:]/365.25

  try:
    n = len(xpt)
  except:
    n = 1
  
  nt = len(time)
  tpt = time
  vxpt = np.zeros([nt,n])
  vypt = np.zeros([nt,n])
  vpt = np.zeros([nt,n])
  
  for k in range(0,n):
    i = np.argmin(abs(xpt[k]-x))
    j = np.argmin(abs(ypt[k]-y))
    vxpt[:,k] = vx[:,j,i]*365.25
    vypt[:,k] = vy[:,j,i]*365.25
    vpt[:,k] = v[:,j,i]*365.25
  
  vxpt[vpt==0] = float('nan')
  vxpt[vpt==0] = float('nan')
  vpt[vpt==0] = float('nan')
    
  if xy_velocities == 'True':
    return vpt,tpt,vxpt,vypt
  else:  	  
    return vpt,tpt

def howat_optical_at_pts(xpt,ypt,glacier,xy_velocities='False'):

  DIR = os.path.join(os.getenv("DATA_HOME"),"Velocity/Howat/"+glacier+"/")
  
  dirs = os.listdir(DIR)
  
  m=0
  for dir in dirs:
    if dir.startswith('OPT'):
      m = m+1
  try:
    n = len(xpt)
  except:
    n = 1
  
  # Set up variables
  velocities = np.zeros([m,n])
  velocities[:,:] = 'nan'
  velocities_x = np.zeros([m,n])
  velocities_x[:,:] = 'nan'
  velocities_y = np.zeros([m,n])
  velocities_y[:,:] = 'nan'
  error = np.zeros([m,n])
  error[:,:] = 'nan'
  times=np.zeros([m,2])
  
  count = 0
  for dir in dirs:
    if dir.startswith('OPT'):
      metafile = open(DIR+dir+'/'+dir+'.meta',"r")
      lines = metafile.readlines()
      metafile.close()
      jdates = []
      jdates.append(float(lines[0][36:48]))
      if len(lines[0]) > 50:
        jdates.append(float(lines[0][49:61]))
        if len(lines[0]) > 63:
          jdates.append(float(lines[0][62:74]))
          if len(lines[0]) > 75:
            jdates.append(float(lines[0][75:87]))
            if len(lines[0]) > 88:
              jdates.append(float(lines[0][88:100]))
              if len(lines[0]) > 101:
                jdates.append(float(lines[0][101:113]))
    
      # Get date
      times_all = []
      for jdate in jdates:
    	  year,month,day,fracday=jdcal.jd2gcal(jdate,0)
    	  times_all.append(datelib.date_to_fracyear(year,month,day+fracday))
    
      times[count,0] = np.mean(times_all)
      times[count,1] = np.max(times_all)-np.min(times_all)
    
      x,y,vx = geotifflib.read(DIR+dir+'/'+dir+'.vx.tif',no_data_value=-99999)
      x,y,vy = geotifflib.read(DIR+dir+'/'+dir+'.vy.tif',no_data_value=-99999)
      x,y,ex = geotifflib.read(DIR+dir+'/'+dir+'.ex.tif',no_data_value=-99999)
      x,y,ey = geotifflib.read(DIR+dir+'/'+dir+'.ey.tif',no_data_value=-99999)
      v = np.sqrt(vx**2+vy**2)
    
      fv = scipy.interpolate.RegularGridInterpolator([y,x],v,method='linear',bounds_error=False)
      fex = scipy.interpolate.RegularGridInterpolator([y,x],ex,method='linear',bounds_error=False)
      fey = scipy.interpolate.RegularGridInterpolator([y,x],ey,method='linear',bounds_error=False)
      fvx = scipy.interpolate.RegularGridInterpolator([y,x],vx,method='linear',bounds_error=False)
      fvy = scipy.interpolate.RegularGridInterpolator([y,x],vy,method='linear',bounds_error=False)
      
      # Find velocities 
      velocities[count,:] = fv(np.array([ypt,xpt]).T) 
      velocities_x[count,:] = fvx(np.array([ypt,xpt]).T)
      velocities_y[count,:] = fvy(np.array([ypt,xpt]).T)
      error[count,:] = velocities[count,:]*np.sqrt((fex(np.array([ypt,xpt]).T)/fvx(np.array([ypt,xpt]).T))**2+(fey(np.array([ypt,xpt]).T)/fvy(np.array([ypt,xpt]).T))**2)        
    
      count = count + 1
    
  # Sort arrays by time  
  sortind=np.argsort(times[:,0],0)
  tpt_sort = times[sortind]
  vpt_sort = velocities[sortind,:]
  ept_sort = error[sortind,:]
  vxpt_sort = velocities_x[sortind,:]
  vypt_sort = velocities_y[sortind,:]
  
  if xy_velocities == 'True':
    return vpt_sort,tpt_sort,ept_sort,vxpt_sort,vypt_sort
  else:  	  
    return vpt_sort,tpt_sort,ept_sort

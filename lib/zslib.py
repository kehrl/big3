'''
This script contains all elevation products for Helheim and Kanger, so that we are consistent
between scripts. 

LMK, UW, 11/12/2015

'''

import os
import sys
import numpy as np
import icefrontlib, coordlib, masklib, distlib,datelib, geotifflib, floatlib
import scipy.interpolate
import scipy.signal as signal
import shapely.geometry
import scipy.ndimage
from scipy import stats

def smith_db_query(xpts,ypts,glacier,verticaldatum='geoid',data='all',maxdist=200.):

  '''
  This code has not been fully developed. It will need to be improved 
  if we're going to use it.
  '''

  if glacier == 'Helheim':
    queryfile = os.path.join(os.getenv("DATA_HOME"),\
    	"Elevation/Smith_db_queries/helheim_ptdb_20150921.mat")
  elif glacier == 'Kanger':
    queryfile = os.path.join(os.getenv("DATA_HOME"),\
    	"Elevation/Smith_db_queries/helheim_ptdb_20150921.mat")

  import h5py
  f = h5py.File(queryfile, 'r')
  key = f.keys()[0]
  g = f[key]
  names = g.keys()
  #Limit to the ones we care about
  names = ['x', 'y', 'z', 'time', 'sensor']
  nc = len(names)
  nr = g.get(names[0]).size
  #Need float64 here for time ordinal
  a = np.zeros((nr, nc), dtype=np.float64)
  print "Importing h5 mat points to np array"
  for i,name in enumerate(names):
    a[:,i] = np.array(g.get(name))

  return zpt,zptstd,time

def gimp_at_pts(xpts,ypts,glacier,verticaldatum):
  
  '''
  elev = gimp_at_pts(xpts,ypts,glacier,verticaldatum)
  Find GIMP surface elevation at points xpts,ypts.
  
  Inputs:
  xpts,ypts: points where you want interpolation
  glacier: glacier name
  verticaldatum: ellipsoid or geoid
  
  Outputs:
  elev: 1-d numpy array of GIMP surface elevations at xpts, ypts
  '''
  

  # Read grid of gimp dem
  x,y,z = gimp_grid(np.min(xpts)-2.0e3,np.max(xpts)+2.0e3,np.min(ypts)-2.0e3,\
  			np.max(ypts)+2.0e3,glacier=glacier,verticaldatum=verticaldatum)
  
  # Interpolate DEM to points
  f = scipy.interpolate.RegularGridInterpolator((y,x),z,method="linear")
  zs = f(np.column_stack([ypts,xpts]))
  
  return zs
  
def gimp_grid(xmin,xmax,ymin,ymax,glacier,verticaldatum):
  
  '''
  x,y,elev = gimp_grid(xmin,xmax,ymin,ymax,glacier,verticaldatum)
  
  Find GIMP surface elevation in grid defined by xmin,xmax,ymin,ymax.
  
  Inputs:
  xmin,xmax,ymin,ymax: extent of grid (grid spacing will be the same as the original GIMP product)
  glacier: glacier name
  verticaldatum: ellipsoid or geoid
  
  Outputs:
  x,y: arrays of x and y coordinates for grid
  elev: grid of GIMP surface elevations 
  '''
  
  # Get correct tile of the GIMP DEM
  if glacier == 'Helheim':
    subset = 'gimpdem3_1'
  elif glacier == 'Kanger':
    subset = 'gimpdem4_2'
  else: 
    sys.exit("Unknown glacier.")
  
  # Load correct file for geoid or ellipsoid heights
  if verticaldatum == 'geoid':
    fileend = '-adj.tif'
  else:
    fileend = '.tif'
  
  # File name  
  file = os.path.join(os.getenv("DATA_HOME"),'Elevation/Gimp/'+subset+fileend)
    
  # Get GIMP DEM
  [x,y,zs]=geotifflib.read(file,xmin,xmax,ymin,ymax)
    
  return x,y,zs

def atm(years='all',verticaldatum='geoid'):

  '''
  atm = atm(years='all',verticaldatum='geoid')
  
  Get ATM data for particular years.
  
  Inputs:
  years: years that you want ATM data
  verticaldatum: ellipsoid or geoid
  
  Outputs:
  atm: dictionary of ATM surface elevations for years; column 1 is x coordinate, column 2 
  		is y coordinate, column 3 is surface elevation according to the chosen vertical datum 
  '''
  
  # Directory for ATM
  ATMDIR = os.path.join(os.getenv("DATA_HOME"),'Elevation/ATM/HelheimKanger/')
  DIRs = os.listdir(ATMDIR)
  
  atm = {}
  
  # Load only desired data
    
  for DIR in DIRs:
    if (DIR.startswith('2') or DIR.startswith('1')) and ((years == 'all') or (DIR[0:4] in years)):
      print "Loading ATM", DIR
      x=[]
      y=[]
      z=[]
      files = os.listdir(ATMDIR+DIR)
      for file in files:
        if 'nadir' in file and not(file.endswith('.xml')):
          try:
            data=np.loadtxt(ATMDIR+DIR+'/'+file,comments='#')
          except:
            data=np.loadtxt(ATMDIR+DIR+'/'+file,comments='#',delimiter=',')
          xfile = data[:,2]
          yfile = data[:,1]
          zfile = data[:,3]
          x = np.hstack([x,xfile])
          y = np.hstack([y,yfile])
          z = np.hstack([z,zfile])
      x2,y2 = coordlib.convert(x-360,y,4326,3413)
      if DIR[0:4] in atm.keys():
        print "Already data from that year, consider changing how you have labeled the directories"
      else:
        atm[DIR] = np.column_stack([x2,y2,z])
  
  # Choose vertical datum
  if verticaldatum=='geoid':
    for key in atm.keys():
      atm[key][:,2] = coordlib.geoidheight(atm[key][:,0],atm[key][:,1],atm[key][:,2])
  elif verticaldatum=='ellipsoid':
    atm = atm
  else: 
    sys.exit("Don't recognize that vertical datum.")
    
  return atm

def atm_along_flowline(xpts,ypts,glacier,years='all',cutoff='none',maxdist=200,method='closest',verticaldatum='geoid',filt_len='none'):
  
  '''
  pts = atm_along_flowline(xpts,ypts,glacier,years='all',cutoff='none',maxdist=200,verticaldatum='geoid')
  
  Get atm surface elevations along a flowline. Select points that are within "maxdist" of the flowline.
  
  Inputs:
  xpts,ypts: coordinates of flowline
  glacier: glacier name
  years: years that you want data
  cutoff: cutoff surface elevations in front of glacier terminus (options: 'none' or 'terminus')
  maxdist: maximum distance between ATM point and flowline point to include in output
  verticaldatum: ellipsoid or geoid
  
  Outputs:
  pts: dictionary of ATM points along flowline, maybe this should be changed? or you should 
  just use "atm_at_pts" which returns an array of surface elevations at the chosen points?
  
  '''
  
  # Get ATM data
  data = atm(years,verticaldatum)
  
  # Get terminus position if we want to throw data points out in front of terminus
  dists = distlib.transect(xpts,ypts)
  if cutoff == 'terminus':
    term_values, term_time = icefrontlib.distance_along_flowline(xpts,ypts,dists,glacier,'icefront')
  
  # Set up output as dictionary
  pts = {}
  
  dates = data.keys()
  R = len(xpts)
  for date in dates:
    z = np.zeros(R); z[:] = float('NaN')
    zstd = np.zeros(R); zstd[:] = float('NaN')
    d = np.zeros(R); d[:] = float('NaN')
    for i in range(0,R):
      if method == 'closest':
        ind = np.argmin((xpts[i]-data[date][:,0])**2 + (ypts[i]-data[date][:,1])**2)
        xatm = data[date][ind,0]
        yatm = data[date][ind,1]
        d[i] = distlib.between_pts(xpts[i],ypts[i],xatm,yatm)
        if d[i] < maxdist:
          z[i] = data[date][ind,2]
      elif method == 'average':
        ind = np.where(np.sqrt((xpts[i]-data[date][:,0])**2 + (ypts[i]-data[date][:,1])**2)<maxdist)[0]
        nonnan = np.where(~(np.isnan(data[date][ind,2])))[0]
        if len(nonnan) > 0.8*len(ind):
          z[i] = np.nanmean(data[date][ind,2])
          zstd[i] = np.nanstd(data[date][ind,2])
    # Filter if desired
    if filt_len != 'none':
      cutoff=(1/filt_len)/(1/((dists[1]-dists[0])*2))
      b,a=signal.butter(4,cutoff,btype='low')
      print "Filtering with a cutoff of ", filt_len
      for i in range(0,len(dates)):
        nonnan = np.where(~(np.isnan(z)))[0]
        if len(nonnan) > 20:
          z[nonnan]=signal.filtfilt(b,a,z[nonnan])
        
    if cutoff == 'terminus':
      # Get fractional year
      time = datelib.date_to_fracyear(float(date[0:4]),float(date[4:6]),float(date[6:8]))
    
      # Get terminus position when ATM was collected
      termpos = np.interp(time,term_time,term_values)
      ind = np.where(dists > termpos)[0]
      z[ind] = float('NaN') 
      zstd[ind] = float('NaN')
    
    # Save elevations if there are any
    if len(np.where(~(np.isnan(z)))[0]) > 0:
      pts[date] = np.column_stack([xpts,ypts,z,zstd])
  
  if cutoff == 'terminus':
    print "Cleaning up DEM points by removing points in front of ice front."
    print "You can change this setting by setting `cutoff = 'none'.'"
	    
  return pts
  
def atm_at_pts(xpts,ypts,glacier,years='all',maxdist=200,verticaldatum='geoid',method='average',cutoff='none'):

  '''
  zpts,time = atm_at_pts(xpts,ypts,glacier,years='all',maxdist=200,verticaldatum='geoid')
  
  Find ATM surface elevations at xpts,ypts.
  
  Inputs:
  xpts,ypts: desired coordinates for ATM surface elevations
  glacier: glacier name
  years: years that you want data
  maxdist: maximum allowable distance between point and ATM data
  verticaldatum: ellipsoid or geoid
  
  Outputs:
  zpts: array of surface elevations at xpts,ypts
  time: time for each array
  '''

  # Get ATM data
  atm = atm_along_flowline(xpts,ypts,glacier,years=years,maxdist=maxdist,cutoff='none',verticaldatum=verticaldatum,method=method)
  
  # Get dates
  dates = np.sort(atm.keys())
  time = np.zeros(len(dates))
  zpts = np.zeros([len(dates),len(xpts)])
  zptstd = np.zeros([len(dates),len(xpts)])
  
  # Load date and data
  for i in range(0,len(dates)):
    date = dates[i]
    time[i] = datelib.date_to_fracyear(float(date[0:4]),float(date[4:6]),float(date[6:8]))
    zpts[i,:] = atm[date][:,2]
    zptstd[i,:] = atm[date][:,3]

  return zpts,zptstd,time

def lvis(years='all',verticaldatum='geoid'):

  '''
  atm = atm(years='all',verticaldatum='geoid')
  
  Get ATM data for particular years.
  
  Inputs:
  years: years that you want ATM data
  verticaldatum: ellipsoid or geoid
  
  Outputs:
  atm: dictionary of lvis surface elevations for years; column 1 is x coordinate, column 2 
  		is y coordinate, column 3 is surface elevation according to the chosen vertical datum 
  '''
  
  # Directory for LVIS
  LVISDIR = os.path.join(os.getenv("DATA_HOME"),'Elevation/LVIS/HelheimKanger/')
  DIRs = os.listdir(LVISDIR)
  
  lvis = {}
  
  # Load only desired data
    
  for DIR in DIRs:
    if (DIR.startswith('2') or DIR.startswith('1')) and ((years == 'all') or (DIR[0:4] in years)):
      print "Loading LVIS", DIR
      x=[]
      y=[]
      z=[]
      files = os.listdir(LVISDIR+DIR)
      for file in files:
        if not(file.endswith('.xml')):
          try:
            data=np.loadtxt(LVISDIR+DIR+'/'+file,comments='#')
          except:
            data=np.loadtxt(LVISDIR+DIR+'/'+file,comments='#',delimiter=',')
          xfile = data[:,3]
          yfile = data[:,4]
          zfile = data[:,5]
          x = np.hstack([x,xfile])
          y = np.hstack([y,yfile])
          z = np.hstack([z,zfile])
      x2,y2 = coordlib.convert(x-360,y,4326,3413)
      if DIR[0:4] in lvis.keys():
        print "Already data from that year, consider changing how you have labeled the directories"
      else:
        lvis[DIR] = np.column_stack([x2,y2,z])
  
  # Choose vertical datum
  if verticaldatum=='geoid':
    for key in lvis.keys():
      lvis[key][:,2] = coordlib.geoidheight(lvis[key][:,0],lvis[key][:,1],lvis[key][:,2])
  elif verticaldatum=='ellipsoid':
    lvis = lvis
  else: 
    sys.exit("Don't recognize that vertical datum.")
    
  return lvis
  
def lvis_along_flowline(xpts,ypts,glacier,years='all',cutoff='none',maxdist=200,method='closest',verticaldatum='geoid',filt_len='none'):
  
  '''
  pts = lvis_along_flowline(xpts,ypts,glacier,years='all',cutoff='none',maxdist=200,verticaldatum='geoid')
  
  Get lvis surface elevations along a flowline. Select points that are within "maxdist" of the flowline.
  
  Inputs:
  xpts,ypts: coordinates of flowline
  glacier: glacier name
  years: years that you want data
  cutoff: cutoff surface elevations in front of glacier terminus (options: 'none' or 'terminus')
  maxdist: maximum distance between lvis point and flowline point to include in output
  verticaldatum: ellipsoid or geoid
  
  Outputs:
  pts: dictionary of lvis points along flowline, maybe this should be changed? or you should 
  just use "lvis_at_pts" which returns an array of surface elevations at the chosen points?
  
  '''
  
  # Get lvis data
  data = lvis(years,verticaldatum)
  
  # Get terminus position if we want to throw data points out in front of terminus
  dists = distlib.transect(xpts,ypts)
  if cutoff == 'terminus':
    term_values, term_time = icefrontlib.distance_along_flowline(xpts,ypts,dists,glacier,'icefront')
  
  # Set up output as dictionary
  pts = {}
  
  dates = data.keys()
  R = len(xpts)
  for date in dates:
    z = np.zeros(R); z[:] = float('NaN')
    zstd = np.zeros(R); zstd[:] = float('NaN')
    d = np.zeros(R); d[:] = float('NaN')
    for i in range(0,R):
      if method == 'closest':
        ind = np.argmin((xpts[i]-data[date][:,0])**2 + (ypts[i]-data[date][:,1])**2)
        xlvis = data[date][ind,0]
        ylvis = data[date][ind,1]
        d[i] = distlib.between_pts(xpts[i],ypts[i],xlvis,ylvis)
        if d[i] < maxdist:
          z[i] = data[date][ind,2]
      elif method == 'average':
        ind = np.where(np.sqrt((xpts[i]-data[date][:,0])**2 + (ypts[i]-data[date][:,1])**2)<maxdist)[0]
        nonnan = np.where(~(np.isnan(data[date][ind,2])))[0]
        if len(nonnan) > 0.8*len(ind):
          z[i] = np.nanmean(data[date][ind,2])
          zstd[i] = np.nanstd(data[date][ind,2])
    # Filter if desired
    if filt_len != 'none':
      cutoff=(1/filt_len)/(1/((dists[1]-dists[0])*2))
      b,a=signal.butter(4,cutoff,btype='low')
      print "Filtering with a cutoff of ", filt_len
      for i in range(0,len(dates)):
        nonnan = np.where(~(np.isnan(z)))[0]
        if len(nonnan) > 20:
          z[nonnan]=signal.filtfilt(b,a,z[nonnan])
        
    if cutoff == 'terminus':
      # Get fractional year
      time = datelib.date_to_fracyear(float(date[0:4]),float(date[4:6]),float(date[6:8]))
    
      # Get terminus position when LVIS was collected
      termpos = np.interp(time,term_time,term_values)
      ind = np.where(dists > termpos)[0]
      z[ind] = float('NaN') 
      zstd[ind] = float('NaN')
    
    # Save elevations if there are any
    if len(np.where(~(np.isnan(z)))[0]) > 0:
      pts[date] = np.column_stack([xpts,ypts,z,zstd])
  
  if cutoff == 'terminus':
    print "Cleaning up DEM points by removing points in front of ice front."
    print "You can change this setting by setting `cutoff = 'none'.'"
	    
  return pts
  
def lvis_at_pts(xpts,ypts,glacier,years='all',maxdist=200,verticaldatum='geoid',method='average',cutoff='none'):

  '''
  zpts,time = lvis_at_pts(xpts,ypts,glacier,years='all',maxdist=200,verticaldatum='geoid')
  
  Find LVIS surface elevations at xpts,ypts.
  
  Inputs:
  xpts,ypts: desired coordinates for LVIS surface elevations
  glacier: glacier name
  years: years that you want data
  maxdist: maximum allowable distance between point and LVIS data
  verticaldatum: ellipsoid or geoid
  
  Outputs:
  zpts: array of surface elevations at xpts,ypts
  time: time for each array
  '''

  # Get lvis data
  lvis = lvis_along_flowline(xpts,ypts,glacier,years=years,maxdist=maxdist,cutoff='none',verticaldatum=verticaldatum,method=method)
  
  # Get dates
  dates = np.sort(lvis.keys())
  time = np.zeros(len(dates))
  zpts = np.zeros([len(dates),len(xpts)])
  zptstd = np.zeros([len(dates),len(xpts)])
  
  # Load date and data
  for i in range(0,len(dates)):
    date = dates[i]
    time[i] = datelib.date_to_fracyear(float(date[0:4]),float(date[4:6]),float(date[6:8]))
    zpts[i,:] = lvis[date][:,2]
    zptstd[i,:] = lvis[date][:,3]

  return zpts,zptstd,time

def dem_error(glacier,filename,sensor):
  
  '''
  error = dem_error(glacier,date,sensor)
  
  Get absolute median error from the co-registration for the DEM on date "date".
  
  Inputs:
  glacier: glacier name
  date: date for error
  sensor: sensor name (Worldview, TDM, Spirit)
  
  Outputs:
  error: error value for chosen date and sensor combination
  '''
  
  if sensor == 'SPIRIT':
    error = 4.0
  else:
    if sensor == 'TDM':
      file = os.path.join(os.getenv("DATA_HOME"),'Elevation/TDM/'+glacier+'/'+glacier+'_TDM_robust_stats_after_fn.csv')
    elif sensor == 'WV':
      file = os.path.join(os.getenv("DATA_HOME"),'Elevation/Worldview/'+glacier+'/'+glacier+'_robust_stats_after_fn.csv')

    fid = open(file)
    lines = fid.readlines()
    for line in lines:
      p = line.split(",")
      if sensor == 'TDM':
        if p[0][0:25] == str(filename):
          error = float(p[-3])
      else:
        if p[0] == str(filename):
          error = float(p[-3])
         
  return error
   
def dem_grid(glacier,xmin=-np.Inf,xmax=np.Inf,ymin=-np.Inf,ymax=np.Inf,
		resolution=32,data='all',verticaldatum='geoid',years='all',method='linear',
		return_error=False):

  '''
  x,y,zs_nonnan,time_nonnan=dem_grid(glacier,xmin=0,xmax=0,ymin=0,ymax=0,resolution=32,data='all',verticaldatum='geoid')
  
  Output an array of worldview grids.
  
  Inputs
  xmin,xmax,ymin,ymax: extent for output grids
  data: type of data you want (WV, TDM, SPIRIT, or all)
  glacier: Kanger or Helheim
  verticaldatum: ellipsoid or geoid
  return_error: return NMAD co-registration errors if set to true
  
  Output: 
  x,y,zs_nonnan,time_nonnan: numpy array of surface elevations on the grid x,y for times_nonnan
  '''
  
  # Set up output grid
  dx = dy = float(resolution)
  nx = np.ceil((xmax-xmin)/dx)+1
  x = np.linspace(xmin,(nx-1)*dx+xmin,nx)
  ny = np.ceil((ymax-ymin)/dx)+1
  y = np.linspace(ymin,(ny-1)*dy+ymin,ny)
  xgrid,ygrid = np.meshgrid(x,y)
  
  # Data directories  
  WVDIR = os.path.join(os.getenv("DATA_HOME"),'Elevation/Worldview/'+glacier+'/')
  WVDIRs = os.listdir(WVDIR)
  
  TDMDIR = os.path.join(os.getenv("DATA_HOME"),'Elevation/TDM/'+glacier+'/')
  TDMDIRs = os.listdir(TDMDIR)
  
  SPIRITDIR = os.path.join(os.getenv("DATA_HOME"),'Elevation/SPIRIT/'+glacier+'/')
  SPIRITDIRs = os.listdir(SPIRITDIR)
  
  # Find file ending based on whether we want elevation relative to geoid or ellipsoid
  if verticaldatum == 'geoid':
    filestring = 'trans-adj.tif'
    spiritstring = 'hae-trans_source-DEM-adj.tif'
  else:
    filestring = 'trans.tif'
    spiritstring = 'hae-trans_source-DEM.tif'
  
  # Find dates where we have data in the desired region
  glacier_extent = shapely.geometry.Polygon([(xmin,ymin),(xmin,ymax),(xmax,ymax),(xmax,ymin)])
  
  # For Worldview data...
  dates=[]
  if ('WV' in data) or (data == 'all'):
    for DIR in WVDIRs:
      if (DIR[0:13] not in dates) and (DIR.startswith('2')) and (DIR.endswith(filestring)):
        wvxmin,wvxmax,wvymin,wvymax = geotifflib.extent(WVDIR+DIR)
        wv_extent = shapely.geometry.Polygon([(wvxmin,wvymin),(wvxmin,wvymax),(wvxmax,wvymax),(wvxmax,wvymin)])
        if glacier_extent.intersects(wv_extent):
          if not(years) or (years=='all'):
            dates.append(DIR[0:13])
          else:
            if len(years) == 4:
              if DIR[0:4] in years:
                dates.append(DIR[0:13])
            else:
              if DIR[0:8] in years:
                dates.append(DIR[0:13])
  # For TDM data...
  if ('TDM' in data) or (data == 'all'):
    for DIR in TDMDIRs:
      if DIR.endswith(filestring):
        tdxmin,tdxmax,tdymin,tdymax = geotifflib.extent(TDMDIR+DIR)
        td_extent = shapely.geometry.Polygon([(tdxmin,tdymin),(tdxmin,tdymax),(tdxmax,tdymax),(tdxmax,tdymin)])
        if glacier_extent.intersects(td_extent):
          if not(years) or (years=='all'):
            dates.append(DIR[0:13])
          else:
            if len(years) == 4:
              if DIR[0:4] in years:
                dates.append(DIR[0:13])
            else:
              if DIR[0:8] in years:
                dates.append(DIR[0:813])
  # For SPIRIT data...              
  if ('SPIRIT' in data) or (data == 'all'):
    for DIR in SPIRITDIRs:
      if DIR.endswith(spiritstring):
        sprxmin,sprxmax,sprymin,sprymax = geotifflib.extent(SPIRITDIR+DIR)
        spr_extent = shapely.geometry.Polygon([(sprxmin,sprymin),(sprxmin,sprymax),(sprxmax,sprymax),(sprxmax,sprymin)])
        if glacier_extent.intersects(spr_extent):
          if not(years) or (years=='all'):
            dates.append(DIR[0:13])
          else:
            if len(years) == 4:
              if DIR[0:4] in years:
                dates.append(DIR[0:13])
            else:
              if DIR[0:8] in years:
                dates.append(DIR[0:13])

  # Load data
  time = np.zeros(len(dates))
  zs = np.zeros([ny,nx,len(time)])
  zs[:,:] = float('NaN')
  error = np.zeros(len(dates))
  for i in range(0,len(dates)):
    date = dates[i]
    n = 1 # Count number of files for that date
    ngrid = np.ones([ny,nx])
    time[i] = datelib.date_to_fracyear(int(date[0:4]),int(date[4:6]),float(date[6:8]))
    for DIR in np.r_[TDMDIRs,WVDIRs,SPIRITDIRs]:
      if DIR.startswith(date) and (DIR.endswith(filestring) or DIR.endswith(spiritstring)):
        print "Loading ", DIR
        if os.path.isfile(WVDIR+DIR):
          xwv,ywv,zwv = geotifflib.read(WVDIR+DIR) 
          error[i] = dem_error(glacier,DIR[0:47],'WV')
        elif os.path.isfile(TDMDIR+DIR):
          xwv,ywv,zwv = geotifflib.read(TDMDIR+DIR) 
          error[i] = dem_error(glacier,DIR[0:25],'TDM')
        else:
          xwv,ywv,zwv = geotifflib.read(SPIRITDIR+DIR)
          error[i] = dem_error(glacier,date,'SPIRIT')
        
        # Interpolate onto output grid
        zwv_dem = scipy.interpolate.RegularGridInterpolator([ywv,xwv],zwv,bounds_error = False,method=method,fill_value=float('nan'))
        zwv_flattened = zwv_dem(np.column_stack([ygrid.flatten(),xgrid.flatten()]))
    
        # Reshape to grid
        zwv_on_grid = np.reshape(zwv_flattened,(len(y),len(x)))
        
        if n > 1: # there is more than one file for that date, so we need to treat it differently
          nonnan = np.array(np.where(~(np.isnan(zwv_on_grid))))
          for [k,j] in nonnan.T:
            if np.isnan(zs[k,j,i]):
              zs[k,j,i] = zwv_on_grid[k,j]
            elif ~(np.isnan(zwv_on_grid[k,j])):
              zs[k,j,i] = zwv_on_grid[k,j]#(zs[k,j,i]*(ngrid[k,j]-1)+zwv_on_grid[k,j])/ngrid[k,j]
              ngrid[k,j] = ngrid[k,j]+1
        else:
          zs[:,:,i] = zwv_on_grid
          n = n+1
      
  # It's possible we still pulled a DEM full of NaNs in the desired region. If so, let's chuck those.
  nonnan=[]
  for i in range(0,len(time)):
    if len(np.where(~(np.isnan(zs[:,:,i])))[0]) > 5:
      nonnan.append(i)
  zs_nonnan = zs[:,:,nonnan]
  time_nonnan = time[nonnan]
  error_nonnan = error[nonnan]    
  
  # Sort by time
  sortind = np.argsort(time_nonnan)
  time_nonnan = time_nonnan[sortind]
  zs_nonnan = zs_nonnan[:,:,sortind]
  error_nonnan = error_nonnan[sortind]
  
  if return_error == False:    
    return x,y,zs_nonnan,time_nonnan
  else:
    return x,y,zs_nonnan,time_nonnan,error_nonnan

def dem_along_flowline(xpts,ypts,glacier,years='all',cutoff='terminus',verticaldatum='geoid',filt_len='none',method='linear',data='all',return_error=False):

  '''
  zs,times = dem_along_flowline(xpts, ypts, glacier ,years='all', 
  		cutoff='terminus', verticaldatum='geoid', filt_len='none')
  
  Outputs worldview surface elevations along a flowline defined by xpts,ypts. Same as "atm_along_flowline" except
  for worldview data. See that function for more information.
  '''

  # Data directories
  WVDIR = os.path.join(os.getenv("DATA_HOME"),"Elevation/Worldview/"+glacier+"/")
  WVDIRs = os.listdir(WVDIR)
  
  TDMDIR = os.path.join(os.getenv("DATA_HOME"),"Elevation/TDM/"+glacier+"/")
  TDMDIRs = os.listdir(TDMDIR)
  
  SPIRITDIR = os.path.join(os.getenv("DATA_HOME"),"Elevation/SPIRIT/"+glacier+"/")
  SPIRITDIRs = os.listdir(SPIRITDIR)

  # Find file ending based on whether we want elevation relative to geoid or ellipsoid
  if verticaldatum == 'geoid':
    filestring = 'trans-adj.tif'
    spiritstring = 'hae-trans_source-DEM-adj.tif'
  elif verticaldatum == 'ellipsoid':
    spiritstring = 'hae-trans_source-DEM.tif'
    filestring = 'trans.tif'
  else:
    sys.exit("Unknown vertical datum.")
  
  # Load ice front positions so we can toss data in front of terminus
  dists = distlib.transect(xpts,ypts)
  if cutoff == 'terminus': 
    term_values, term_time = icefrontlib.distance_along_flowline(xpts,ypts,dists,glacier,type='icefront')
  
  
# Find dates where we have data in the desired region
  dates=[]
  #For Worldview...
  if ('WV' in data) or (data == 'all'):
    for DIR in WVDIRs:
      if (DIR[0:8] not in dates) and DIR.startswith('2') and DIR.endswith(filestring):
        xmin,xmax,ymin,ymax = geotifflib.extent(WVDIR+DIR)
        within = np.where((xpts > xmin) & (xpts < xmax) & (ypts > ymin) & (ypts < ymax))[0]
        if len(within) > 0:
          if not(years) or (years=='all'):
            dates.append(DIR[0:8])
          else: 
            if len(years) == 4:
              if DIR[0:4] in years:
                dates.append(DIR[0:8])
            else:
              if DIR[0:8] in years:
                dates.append(DIR[0:8])
  # For TDM...
  if ('TDM' in data) or (data == 'all'):
    for DIR in TDMDIRs:
      if DIR.endswith(filestring):
        xmin,xmax,ymin,ymax = geotifflib.extent(TDMDIR+DIR)
        within = np.where((xpts > xmin) & (xpts < xmax) & (ypts > ymin) & (ypts < ymax))[0]
        if len(within) > 0:
          if not(years) or (years=='all'):
            dates.append(DIR[0:8])
          else:
            if len(years) == 4:
              if DIR[0:4] in years:
                dates.append(DIR[0:8])
            else:
              if DIR[0:8] in years:
                dates.append(DIR[0:8])
  # For SPIRIT...
  if ('SPIRIT' in data) or (data == 'all'):
    for DIR in SPIRITDIRs:
      if DIR.endswith(spiritstring):
        xmin,xmax,ymin,ymax = geotifflib.extent(SPIRITDIR+DIR)
        within = np.where((xpts > xmin) & (xpts < xmax) & (ypts > ymin) & (ypts < ymax))[0]
        if len(within) > 0:
          if not(years) or (years=='all'):
            dates.append(DIR[0:8])
          else:
            if len(years) == 4:
              if DIR[0:4] in years:
                dates.append(DIR[0:8])
            else:
              if DIR[0:8] in years:
                dates.append(DIR[0:8])
  

  # Set up output 
  zs = np.zeros((len(dates),len(xpts)))
  zs[:,:] = float('nan')
  times = np.zeros((len(dates),2))
  error = np.zeros(len(dates))
  for i in range(0,len(dates)):
    date = dates[i]
    times[i,0] = datelib.date_to_fracyear(float(date[0:4]),float(date[4:6]),float(date[6:]))
    times[i,1] = date
    for DIR in np.r_[TDMDIRs,WVDIRs,SPIRITDIRs]:
      if DIR.startswith(date) and (DIR.endswith(filestring) or DIR.endswith(spiritstring)): 
        #print "Loading data from "+DIR+"\n"
        if os.path.isfile(WVDIR+DIR):
          x,y,z = geotifflib.read(WVDIR+DIR) 
          error[i] = dem_error(glacier,DIR[0:47],'WV')
        elif os.path.isfile(TDMDIR+DIR):
          x,y,z = geotifflib.read(TDMDIR+DIR) 
          error[i] = dem_error(glacier,DIR[0:25],'TDM')
        else:
          x,y,z = geotifflib.read(SPIRITDIR+DIR)
          error[i] = dem_error(glacier,date,'SPIRIT')

        dem = scipy.interpolate.RegularGridInterpolator([y,x],z,method=method)
    
        # Find points that fall within the DEM
        ind = np.where((xpts > np.min(x)) & (xpts < np.max(x)) & (ypts > np.min(y)) & (ypts < np.max(y)))
        if ind:
          # Get elevation at coordinates
          zpts = np.zeros(len(xpts)); zpts[:] = 'NaN'
          zpts[ind] = dem(np.column_stack([ypts[ind],xpts[ind]]))
          
          # Get fractional year
          time = datelib.date_to_fracyear(float(date[0:4]),float(date[4:6]),float(date[6:8]))

          if cutoff == 'terminus':
            # Get terminus position at time of worldview image
            termpos = np.interp(time,term_time,term_values)
            ind = np.where(dists > termpos)[0]
            zpts[ind] = float('NaN')
          
          # Save points to output dictionary
          if len(np.where(~(np.isnan(zpts)))[0]) < 1:
            print "Not loading ",DIR
          elif date in times[:,1]:
            nans = np.where((np.isnan(pts[date][:,2])))
            zs[i,nans] = zpts[nans] 
          else:
            zs[i,:] = zpts
  
  # Sort by time
  sortind = np.argsort(times[:,0])
  times = times[sortind,:]
  zs = zs[sortind,:]
  error = error[sortind]
  
  # Print warning if removing points in front of ice front
  if cutoff == 'terminus':
    print "Cleaning up DEM points by removing points in front of ice front."
    print "You can change this setting by setting `cutoff = 'none'.'"
  
  # Filter if desired
  if filt_len != 'none':
    cutoff=(1/filt_len)/(1/((dists[1]-dists[0])*2))
    b,a=signal.butter(4,cutoff,btype='low')
    print "Filtering with a cutoff of ", filt_len
    for i in range(0,len(dates)):
      nonnan = np.where(~(np.isnan(zs[i,:])))[0]
      if len(nonnan) > 20:

        zs[i,nonnan]=signal.filtfilt(b,a,zs[i,nonnan])
  
  if return_error == False:    
    return zs,times
  else:
    return zs,times,error
  
def dem_at_pts(xpts,ypts,glacier,years='all',verticaldatum='geoid',cutoff='none',method='linear',radius=200,data='all'):

  '''
  zpts,zpts_error,time = dem_at_pts(xpts,ypts,glacier,years='all',
  	verticaldatum='geoid',cutoff='none',method='linear',radius=200,data='all'):

  
  Interpolates worldview surface elevations to points xpts,ypts.
  
  Inputs:
  xpts,ypts: points where we want surface elevations
  glacier: glacier name
  years: years that we want data
  verticaldatum: geoid or ellipsoid
  cutoff: cutoff elevations in front of terminus ('terminus' or 'none')
  method: 'linear' extrapolation or 'average' value for a region defined by radius
  radius: radius for average value (only necessary if method is 'average')
  data: all, WV, TDM, or SPIRIT
  
  Outputs:
  zpts: array of surface elevations for points xpts,ypts
  time: time of the arrays
  '''
  
  if method == 'linear':
    zpts,times,zpts_error = dem_along_flowline(xpts,ypts,glacier,years=years,cutoff=cutoff,verticaldatum=verticaldatum,method=method,data=data,return_error=True)
    time = times[:,0]
  elif method == 'average':
    xmin = np.min(xpts)-radius*3
    ymin = np.min(ypts)-radius*3
    xmax = np.max(xpts)+radius*3
    ymax = np.max(ypts)+radius*3
    
    xwv,ywv,zwv,timewv,zpts_error = dem_grid(glacier,xmin,xmax,ymin,ymax,years='all',verticaldatum='geoid',return_error=True,data=data)
    N = len(timewv)
    
    time = timewv
    zpts = np.zeros([N,len(xpts)])
    zpts[:,:] = float('nan')
    
    xwv_grid,ywv_grid = np.meshgrid(xwv,ywv) 
    xwv_flattened = xwv_grid.flatten()
    ywv_flattened = ywv_grid.flatten()
    
    for j in range(0,len(xpts)):
      ind = np.where(np.sqrt((xpts[j] - xwv_flattened)**2+(ypts[j] - ywv_flattened)**2) < radius)[0]
      for i in range(0,N):
        nonnan = np.where(~(np.isnan(zwv[:,:,i].flatten()[ind])))[0]
        if len(nonnan) > 0.9*len(ind):
          zpts[i,j] = np.nanmean(zwv[:,:,i].flatten()[ind])

    
  return zpts,zpts_error,time

def dem_at_lagpts(xf,yf,dists,zb,pts,glacier,years='all',verticaldatum='geoid',cutoff='none',method='linear',radius=200,data='all',rho_i=917.,rho_sw=1020.):

  '''
  zpts,zpts_haf,zpts_error,time = dem_at_pts(xpts,ypts,glacier,years='all',
  	verticaldatum='geoid',cutoff='none',method='linear',radius=200,data='all'):

  
  Interpolates worldview surface elevations to lagrangian points behind flowline.
  
  Inputs:
  xf,yf: x,y coordinates for flowline
  dists: distances along flowline
  zb: bed elevations (given relative to the geoid) to calculate height above flotation
  pts: lagrangian distances where we want surface elevation
  glacier: glacier name
  years: years that we want data
  verticaldatum: geoid or ellipsoid
  cutoff: cutoff elevations in front of terminus ('terminus' or 'none')
  method: 'linear' extrapolation or 'average' value for a region defined by radius
  radius: radius for average value (only necessary if method is 'average')
  data: all, WV, TDM, or SPIRIT
  
  Outputs:
  zpts: array of surface elevations for points xpts,ypts
  zpts_hab: array of height above flotation for distances "pts" behind terminus
  zpts_error
  time: time of the arrays
  '''
  
  # Get surface elevation along flowline for each timestep
  if method == 'linear':
    zf,times,zf_error = dem_along_flowline(xf,yf,glacier,years=years,cutoff=cutoff,verticaldatum=verticaldatum,method=method,data=data,return_error=True)
    time = times[:,0]
  elif method == 'average':
    xmin = np.min(xf)-radius*3
    ymin = np.min(yf)-radius*3
    xmax = np.max(xf)+radius*3
    ymax = np.max(yf)+radius*3
    
    xwv,ywv,zwv,timewv,zf_error = dem_grid(glacier,xmin,xmax,ymin,ymax,years='all',verticaldatum='geoid',return_error=True,data=data)
    N = len(timewv)
    
    time = timewv
    zf = np.zeros([N,len(xf)])
    zf[:,:] = float('nan')
    
    xwv_grid,ywv_grid = np.meshgrid(xwv,ywv) 
    xwv_flattened = xwv_grid.flatten()
    ywv_flattened = ywv_grid.flatten()
    
    for j in range(0,len(xf)):
      ind = np.where(np.sqrt((xf[j] - xwv_flattened)**2+(yf[j] - ywv_flattened)**2) < radius)[0]
      for i in range(0,N):
        nonnan = np.where(~(np.isnan(zwv[:,:,i].flatten()[ind])))[0]
        if len(nonnan) > 0.9*len(ind):
          zf[i,j] = np.nanmean(zwv[:,:,i].flatten()[ind])
    del ind

  # Get terminus positions
  terminus_val,terminus_time = icefrontlib.distance_along_flowline(xf,yf,dists,glacier,type='icefront')
  
  # Get surface elevations at desired distances from terminus
  N = len(pts)
  zpts = np.zeros([len(time),N])
  zpts[:,:] = float('nan')
  zpts_haf = np.zeros([len(time),N])
  zpts_haf[:,:] = float('nan')
  termpts = np.zeros(len(time))
  zpts_error = np.zeros(len(time))
  for j in range(0,N):
    for i in range(0,len(time)):
      termpts[i] = np.interp(time[i],terminus_time,terminus_val)
      minind = np.argmin(abs(dists-(termpts[i]+pts[j])))
      zpts[i,j] = zf[i,minind]
      zpts_haf[i,j] = zf[i,minind]-floatlib.height(zb[minind],rho_i=rho_i,rho_sw=rho_sw)
      zpts_error[i] = zf_error[i]
  
    
  return zpts,zpts_haf,zpts_error,termpts,time


def grid_near_time(time,glacier,verticaldatum='geoid'):

  '''
  Find the elevation grid closest in time to the input time.
  
  Inputs:
  time: time that you want the grid
  glacier: glacier name
  verticaldatum: geoid or ellipsoid
  
  Outputs:
  x,y: grid coordinates
  zs: surface elevations
  besttime: time of grid that is closest to input time
  '''

  # Worldview data
  WVDIR = os.path.join(os.getenv("DATA_HOME"),"Elevation/Worldview/"+glacier+"/")
  WVDIRs = os.listdir(WVDIR)
  TDMDIR = os.path.join(os.getenv("DATA_HOME"),"Elevation/TDM/"+glacier+"/")
  TDMDIRs = os.listdir(TDMDIR)


  # Find file ending based on whether we want elevation relative to geoid or ellipsoid
  if verticaldatum == 'geoid':
    filestring = 'trans-adj.tif'
  else:
    filestring = 'trans.tif'
  
# Find dates where we have data in the desired region
  besttime = 0
  for DIR in WVDIRs:
    if DIR.startswith('2') and DIR.endswith(filestring):
      demtime = datelib.date_to_fracyear(int(DIR[0:4]),int(DIR[4:6]),int(DIR[6:8]))
      if abs(demtime-time) < abs(besttime-time):
        besttime = demtime
        bestfile = WVDIR+DIR
  for DIR in TDMDIRs:
    if DIR.startswith('2') and DIR.endswith(filestring):
      demtime = datelib.date_to_fracyear(int(DIR[0:4]),int(DIR[4:6]),int(DIR[6:8]))
      if abs(demtime-time) < abs(besttime-time):
        besttime = demtime
        bestfile = TDMDIR+DIR
       
  x,y,zs = geotifflib.read(bestfile)

  return x,y,zs,besttime

def dem_continuous(glacier,xmin,xmax,ymin,ymax,date,verticaldatum='geoid',blur=False):

  '''
  xg,yg,zs = dem_continuous(glacier, date, xmin, xmax, ymin, ymax,
  				verticaldatum='geoid', fillin='none')
  
  The worldview DEMs do not cover the entire glacier extent, so 
  this function combines individual DEMs to produce a continuous DEM 
  for different time periods.
  
  Inputs:
  glacier: glacier name
  date: date for terminus DEM (see function for options)
  xmin,xmax,ymin,ymax: desired output grid size
  verticaldatum: geoid or elliipsoid
  
  Outputs:
  xg,yg: x and y coordinates for output grid
  zg: continuous surface DEM 
  '''
  
  # Get correct tile of the GIMP DEM
  if glacier == 'Helheim':
    subset = 'gimpdem3_1'
  elif glacier == 'Kanger':
    subset = 'gimpdem4_2'
  else: 
    sys.exit("Unknown glacier.")
  
  # Load correct file for geoid or ellipsoid heights
  if verticaldatum == 'geoid':
    fileend = 'trans-adj.tif'
  elif verticaldatum == 'ellipsoid':
    fileend = 'trans.tif'
  else:
    sys.exit("Unknown vertical datum.")
  
  # Gimp filename  
  gimpfile = os.path.join(os.getenv("DATA_HOME"),'Elevation/Gimp/'+subset+fileend[5:])
  
  # Directories for high res DEMs
  OUTDIR = os.path.join(os.getenv("DATA_HOME"),"Elevation/MosaicDEMs/"+glacier+"/")
  WVDIR = os.path.join(os.getenv("DATA_HOME"),"Elevation/Worldview/"+glacier+"/")
  TDMDIR = os.path.join(os.getenv("DATA_HOME"),"Elevation/TDM/"+glacier+"/")
  
  # Select dates for worldview images 

  
  if glacier == 'Helheim':
    if date == '20110319':
      dates = ['20110319','20110615']
    elif date == '20110828':
      dates = ['20110828','20110824','20110615']
    elif date == '20120624':
      dates = ['20120624','20120629','20120513']
    elif date == '20140509':
      dates = ['20140509','2014611','20140408','20130508']
    elif date == '20140731':
      dates = ['20140731','20130508']
    elif date == ['20141016']:
      dates = ['20141016','20130508']
    else:
      dates = np.array([date])
     
  elif glacier == 'Kanger':
    if date == '20110712':
      dates = ['20110712','20110808','20110823']
      dates_backup = ['20110824']
  
  
  files = ''
  dirs_wv = os.listdir(WVDIR)
  dirs_tdm = os.listdir(TDMDIR)
  for d in dates:
    for dir in dirs_tdm:
      if (d in dir) and (dir.endswith(fileend)):
        files = files+' '+TDMDIR+dir
    for dir in dirs_wv:
      if (d in dir) and (dir.endswith(fileend)):
        files = files+' '+WVDIR+dir

  #files = files+' '+gimpfile  
  
  CURRENTDIR = os.getcwd()
  os.chdir(OUTDIR)
  filename = 'mosaic-'+date+'-'+verticaldatum
  #if not(os.path.isfile(filename+'-tile-0.tif')):
  os.system('dem_mosaic --hole-fill-length 5 --t_projwin '+str(xmin)+' '+str(ymin)+' '+str(xmax)+\
  		' '+str(ymax)+' --priority-blending-length 40 -o'+filename+files)
  
  xg,yg,zs = geotifflib.read(filename+"-tile-0.tif")
  
  if blur == True:
    print "Blurring DEM over 17 pixels (roughly 500 m in each direction)..."
    # 17 pixel gaussian blur
    zs_blur = scipy.ndimage.filters.gaussian_filter(zs,sigma=2,truncate=4)
  else:
    zs_blur = zs
  
  os.chdir(CURRENTDIR)
  
  return xg,yg,zs_blur
  
def dem_continuous_flowline(xf,yf,dists,glacier,date,verticaldatum='geoid',filt_len='none'):

  '''
  zflow = dem_continuous_flowline(xf,yf,glacier,date,verticaldatum='geoid',fillin='none')

  Same as "dem_continuous", except output the coordinates along a flowline rather than as a 
  grid. See "dem_continuous" for more information.

  '''
  
  # Get grid for that date 
  xg,yg,zgrid = dem_continuous(glacier,np.min(x)-1e3,np.max(x)+1e3,np.min(y)-1e3,np.max(y)+1e3,date,verticaldatum=verticaldatum,blur=False)
  
  # Interpolate grid onto flowline
  dem = scipy.interpolate.RegularGridInterpolator([yg,xg],zgrid)
   
  zflow = dem((yf,xf)) 
  
  # Filter flowline if filt_len is not set to none
  if filt_len != 'none':
    cutoff=(1/filt_len)/(1/((dists[1]-dists[0])*2))
    b,a=signal.butter(4,cutoff,btype='low')
    zflow_filtered = signal.filtfilt(b,a,zflow)
  else:
    zflow_filtered = zflow
   
  return zflow_filtered
  
def variability(glacier,time1,time2,verticaldatum='geoid',resolution=32.,data='all'):

  '''
  '''
  
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

  # Load DEMs
  # Set up output grid
  dx = dy = float(resolution)
  nx = np.ceil((xmax-xmin)/dx)+1
  x = np.linspace(xmin,(nx-1)*dx+xmin,nx)
  ny = np.ceil((ymax-ymin)/dx)+1
  y = np.linspace(ymin,(ny-1)*dy+ymin,ny)
  xgrid,ygrid = np.meshgrid(x,y)
    
  WVDIR = os.path.join(os.getenv("DATA_HOME"),'Elevation/Worldview/'+glacier+'/')
  WVDIRs = os.listdir(WVDIR)
  TDMDIR = os.path.join(os.getenv("DATA_HOME"),'Elevation/TDM/'+glacier+'/')
  TDMDIRs = os.listdir(TDMDIR)
  SPIRITDIR = os.path.join(os.getenv("DATA_HOME"),'Elevation/SPIRIT/'+glacier+'/')
  SPIRITDIRs = os.listdir(SPIRITDIR)
  
  # Find file ending based on whether we want elevation relative to geoid or ellipsoid
  if verticaldatum == 'geoid':
    filestring = 'trans-adj.tif'
  else:
    filestring = 'trans.tif'
  
  # Find dates where we have data in the desired region
  dates=[]
  glacier_extent = shapely.geometry.Polygon([(xmin,ymin),(xmin,ymax),(xmax,ymax),(xmax,ymin)])
  if 'WV' in data or data=='all':
    for DIR in WVDIRs:
      if (DIR[0:8] not in dates) and DIR.startswith('2') and (DIR.endswith(filestring)):
        wvxmin,wvxmax,wvymin,wvymax = geotifflib.extent(WVDIR+DIR)
        wv_extent = shapely.geometry.Polygon([(wvxmin,wvymin),(wvxmin,wvymax),(wvxmax,wvymax),(wvxmax,wvymin)])
        time_file = datelib.date_to_fracyear(int(DIR[0:4]),int(DIR[4:6]),float(DIR[6:8]))
        if glacier_extent.intersects(wv_extent):
          if (time_file > time1) and (time_file < time2):
            dates.append(DIR[0:8])
        del wvxmin,wvxmax,wvymin,wvymax,wv_extent
  if 'TDM' in data or data=='all':
    for DIR in TDMDIRs:
      if (DIR.endswith(filestring)) and (DIR[0:8] not in dates):
        TDMmin,TDMmax,tdymin,tdymax = geotifflib.extent(TDMDIR+DIR)
        td_extent = shapely.geometry.Polygon([(TDMmin,tdymin),(TDMmin,tdymax),(TDMmax,tdymax),(TDMmax,tdymin)])
        time_file = datelib.date_to_fracyear(int(DIR[0:4]),int(DIR[4:6]),float(DIR[6:8]))
        if glacier_extent.intersects(td_extent):
          if (time_file > time1) and (time_file < time2):
            dates.append(DIR[0:8])
        del td_extent,TDMmin,TDMmax,tdymin,tdymax
  if 'SPIRIT' in data or data=='all':
    for DIR in SPIRITDIRs:
      if DIR.endswith(filestring) and (DIR[0:8] not in dates):
        sprxmin,sprxmax,sprymin,sprymax = geotifflib.extent(SPIRITDIR+DIR)
        spr_extent = shapely.geometry.Polygon([(sprxmin,sprymin),(sprxmin,sprymax),(sprxmax,sprymax),(sprxmax,sprymin)])
        time_file = datelib.date_to_fracyear(int(DIR[0:4]),int(DIR[4:6]),float(DIR[6:8]))
        if glacier_extent.intersects(spr_extent):
          if (time_file > time1) and (time_file < time2):
            dates.append(DIR[0:8])
        del spr_extent, sprxmin,sprxmax,sprymin,sprymax
  del time_file
   
  # Load data
  time = np.zeros(len(dates))
  zs = np.zeros([ny,nx,len(time)])
  mask = np.zeros([ny,nx,len(time)])
  for i in range(0,len(dates)):
    date = dates[i]
    n = 1 # Count number of files for that date
    ngrid = np.ones([ny,nx])
    time[i] = datelib.date_to_fracyear(int(date[0:4]),int(date[4:6]),float(date[6:8]))
    for DIR in np.r_[TDMDIRs,WVDIRs]:
      if DIR.startswith(date) and DIR.endswith(filestring): 
        # Read file
        if  os.path.isfile(WVDIR+DIR):
          maskfile = WVDIR+date+'_mask.tif'
          xwv,ywv,zwv = geotifflib.read(WVDIR+DIR) 
          zwv[zwv==0] = float('NaN')
        elif os.path.isfile(TDMDIR+DIR):
          maskfile = TDMDIR+date+'_mask.tif'
          xwv,ywv,zwv = geotifflib.read(TDMDIR+DIR) 
          zwv[zwv==0] = float('NaN')
        else:
          maskfile = SPIRITDIR+date+'_mask.tif'
          xwv,ywv,zwv = geotifflib.read(SPIRITDIR+DIR) 
          zwv[zwv==0] = float('NaN')          
        
        # Interpolate onto output grid
        zwv_dem = scipy.interpolate.RegularGridInterpolator([ywv,xwv],zwv,bounds_error = False,method='linear',fill_value=float('nan'))
        zwv_flattened = zwv_dem(np.column_stack([ygrid.flatten(),xgrid.flatten()]))
    
        # Reshape to grid
        zwv_on_grid = np.reshape(zwv_flattened,(len(y),len(x)))
        
        if n > 1: # there is more than one file for that date, so we need to treat it differently
          nonnan = np.array(np.where(~(np.isnan(zwv_on_grid))))
          for [k,j] in nonnan.T:
            if np.isnan(zs[k,j,i]):
              zs[k,j,i] = zwv_on_grid[k,j]
            else:
              zs[k,j,i] = (zs[k,j,i]*(ngrid[k,j]-1)+zwv_on_grid[k,j])/ngrid[k,j]
        else:
          zs[:,:,i] = zwv_on_grid
          n = n+1
          
          # Set up mask file
          if os.path.isfile(maskfile):
            xmask,ymask,mask[:,:,i] = geotifflib.read(maskfile)
          else:
            xmask,ymask,mask[:,:,i] = masklib.load_grid(glacier,xmin,xmax,ymin,ymax,dx,icefront_time=time[i])
            geotifflib.write_from_grid(xmask,ymask,np.flipud(mask[:,:,i]),float('nan'),maskfile)   
    
    zs[zs==0] = float('NaN')
  
  del xwv,ywv,zwv,maskfile,n,ngrid
      
  # It's possible we still pulled a DEM full of NaNs in the desired region. If so, let's chuck those.
  nonnan=[]
  for i in range(0,len(time)):
    if len(np.where(~(np.isnan(zs[:,:,i])))[0]) > 5:
      nonnan.append(i)
  zs = zs[:,:,nonnan]
  time = time[nonnan]    

  # Mask out non-ice surface elevations
  #allmask = np.nanmax(mask,axis=2)
  #zs[allmask==1,:] = float('nan')
  for i in range(0,len(time)):
    zs[mask[:,:,i]==1,i] = float('nan') 
  
  # Sort data by time
  sortind = np.argsort(time)
  time = time[sortind]
  zs = zs[:,:,sortind]
  del sortind
  
  # Get standard deviation and mean of all elevations
  zsstd = np.nanstd(zs,axis=2)
  zsmean = np.nanmean(zs,axis=2)
  
  # Get trends
  zstrend = np.zeros_like(zsmean)
  zstrend_time1 = np.zeros_like(zsmean)
  zstrend_time2 = np.zeros_like(zsmean)
  zstrend_p = np.zeros_like(zsmean)
  zstrend_count = np.zeros_like(zsmean)
  zstrend_error = np.zeros_like(zsmean)
  zstrend_r = np.zeros_like(zsmean)
  zstrend_intercept = np.zeros_like(zsmean)
  zstrend[:,:] = float('nan')
  zstrend_p[:,:] = float('nan')
  zstrend_error[:,:] = float('nan')
  zstrend_r[:,:] = float('nan')
  zstrend_intercept[:,:] = float('nan')
  for j in range(0,len(y)):
    for i in range(0,len(x)):
      nonnan = np.where(~(np.isnan(zs[j,i,:])))[0]
      zstrend_count[j,i] = len(nonnan)
      if len(nonnan) > 1:
        if (np.floor(np.min(time[nonnan]))==time1) and np.ceil(np.max(time[nonnan]))==time2: 
        # don't want to calculate a trend if measurements don't cover the entire time period
          zstrend_time1[j,i] = np.min(time[nonnan])
          zstrend_time2[j,i] = np.max(time[nonnan])
          slope,intercept,r,p,std_err = stats.linregress(time[nonnan],zs[j,i,nonnan])
          zstrend[j,i] = slope
          zstrend_intercept[j,i] = intercept
          zstrend_error[j,i] = std_err
          zstrend_p[j,i] = p
          zstrend_r[j,i] = r

  # Remove trend from all observations
  zsdetrend = np.zeros_like(zs)
  zsdetrend[:,:] = float('nan')
  nonnan = np.where(~(np.isnan(zstrend)))
  for i in range(0,len(time)):
    trend = zstrend_intercept+zstrend*time[i]
    #zsdetrend[:,:,i] = zs[:,:,i]
    zsdetrend[:,:,i] = zs[:,:,i]-trend

  # Calculate range of observed values
  zsrange = np.zeros_like(zsmean)
  zsrange[:,:] = float('nan')
  for i in range(0,len(x)):
    for j in range(0,len(y)):
      nonnan = np.where(~(np.isnan(zsdetrend[j,i,:])))[0]
      if len(nonnan) > 1:
        zsrange[j,i] = np.max(zsdetrend[j,i,nonnan])-np.min(zsdetrend[j,i,nonnan])
  
  # Remove insignificant trends from trend matrix
  ind = np.where(zstrend_p > 0.05)
  zstrend[ind] = float('nan')
  
  # Get number of nonnan velocities for each pixel
  zscount = np.zeros([len(y),len(x)])
  for j in range(0,len(y)):
    for i in range(0,len(x)):
      nonnan = np.where(~(np.isnan(zs[j,i,:])))[0]
      zscount[j,i] = len(nonnan)   
  
  return x,y,zs,zstrend,zsdetrend,zsrange,zscount,time


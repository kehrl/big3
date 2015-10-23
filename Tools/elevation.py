# This script contains all elevation products for Helheim and Kanger, so that we are consistent
# between scripts. 

# Functions:
# gimp_pts(xpts,ypts,verticaldatum): gimp elevations at xpts,ypts
# gimp_grid(xmin,xmax,ymin,ymax,glacier,verticaldatum): gimp elevations in the grid specified by
#	xmin,xmax,ymin,ymax
# atm(year,verticaldatum): ATM data for chosen year or "all" years
# atm(xpts,ypts,years,maxdist,verticaldatum): get ATM data at specific points given by xpts, ypts

import os
import sys
import numpy as np
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"Util/Modules"))
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"BigThreeGlaciers/Tools"))
import icefronts, glacier_flowline, coords
import coords, geotiff
import scipy.interpolate, dist, fracyear
import scipy.signal as signal
import shapely.geometry
import scipy.ndimage

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
  [x,y,zs]=geotiff.read(file,xmin,xmax,ymin,ymax)
    
  return x,y,zs

def atm(years,verticaldatum):

  '''
  atm = atm(years,verticaldatum)
  
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
      x2,y2 = coords.convert(x-360,y,4326,3413)
      if DIR[0:4] in atm.keys():
        print "Already data from that year, consider changing how you have labeled the directories"
      else:
        atm[DIR] = np.column_stack([x2,y2,z])
  
  # Choose vertical datum
  if verticaldatum=='geoid':
    for key in atm.keys():
      atm[key][:,2] = coords.geoidheight(atm[key][:,0],atm[key][:,1],atm[key][:,2])
  elif verticaldatum=='ellipsoid':
    atm = atm
  else: 
    print "Don't recognize that vertical datum, defaulting to geoid"
    
  return atm

def atm_along_flowline(xpts,ypts,glacier,years='all',cutoff='none',maxdist=200,method='closest',verticaldatum='geoid'):
  
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
  if cutoff == 'terminus':
    dists = dist.transect(xpts,ypts)
    term_values, term_time = icefronts.distance_along_flowline(xpts,ypts,dists,glacier,'icefront')
  
  # Set up output as dictionary
  pts = {}
  
  dates = data.keys()
  R = len(xpts)
  for date in dates:
    z = np.zeros(R); z[:] = float('NaN')
    d = np.zeros(R); d[:] = float('NaN')
    for i in range(0,R):
      if method == 'closest':
        ind = np.argmin((xpts[i]-data[date][:,0])**2 + (ypts[i]-data[date][:,1])**2)
        xatm = data[date][ind,0]
        yatm = data[date][ind,1]
        d[i] = dist.between_pts(xpts[i],ypts[i],xatm,yatm)
        if d[i] < maxdist:
          z[i] = data[date][ind,2]
      elif method == 'average':
        ind = np.where(np.sqrt((xpts[i]-data[date][:,0])**2 + (ypts[i]-data[date][:,1])**2)<maxdist)[0]
        nonnan = np.where(~(np.isnan(data[date][ind,2])))[0]
        if len(nonnan) > 0.8*len(ind):
          z[i] = np.nanmean(data[date][ind,2])
        
    if cutoff == 'terminus':
      # Get fractional year
      time = fracyear.date_to_fracyear(float(date[0:4]),float(date[4:6]),float(date[6:8]))
    
      # Get terminus position when ATM was collected
      termpos = np.interp(time,term_time,term_values)
      ind = np.where(dists > termpos)[0]
      z[ind] = float('NaN') 
    
    # Save elevations if there are any
    if len(np.where(~(np.isnan(z)))[0]) > 0:
      pts[date] = np.column_stack([xpts,ypts,z,d])
  
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
  
  # Load date and data
  for i in range(0,len(dates)):
    date = dates[i]
    time[i] = fracyear.date_to_fracyear(float(date[0:4]),float(date[4:6]),float(date[6:8]))
    zpts[i,:] = atm[date][:,2]

  return zpts,time
   
def worldview_grid(glacier,xmin=-np.Inf,xmax=np.Inf,ymin=-np.Inf,ymax=np.Inf,resolution=32,years='all',verticaldatum='geoid'):

  '''
  x,y,zs_nonnan,time_nonnan=worldview_grid(glacier,xmin=0,xmax=0,ymin=0,ymax=0,resolution=32,years='all',verticaldatum='geoid')
  
  Output an array of worldview grids.
  
  Inputs
  xmin,xmax,ymin,ymax: extent for output grids
  years: year that you want data or "all" for all years 
  glacier: Kanger or Helheim
  verticaldatum: ellipsoid or geoid
  
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
    
  WVDIR = os.path.join(os.getenv("DATA_HOME"),'Elevation/Worldview/'+glacier+'/')
  TDXDIR = os.path.join(os.getenv("DATA_HOME"),'Elevation/TDX/'+glacier+'/')
  WVDIRs = os.listdir(WVDIR)
  TDXDIRs = os.listdir(TDXDIR)
  
  # Find file ending based on whether we want elevation relative to geoid or ellipsoid
  if verticaldatum == 'geoid':
    filestring = 'trans-adj.tif'
  else:
    filestring = 'trans.tif'
  
  # Find dates where we have data in the desired region
  dates=[]
  glacier = shapely.geometry.Polygon([(xmin,ymin),(xmin,ymax),(xmax,ymax),(xmax,ymin)])
  for DIR in WVDIRs:
    if (DIR[0:8] not in dates) and DIR.startswith('2') and (DIR.endswith(filestring)):
      wvxmin,wvxmax,wvymin,wvymax = geotiff.extent(WVDIR+DIR)
      wv_extent = shapely.geometry.Polygon([(wvxmin,wvymin),(wvxmin,wvymax),(wvxmax,wvymax),(wvxmax,wvymin)])
      if glacier.intersects(wv_extent):
        if not(years) or (years=='all'):
          dates.append(DIR[0:8])
        else:
          if len(years) == 4:
            if DIR[0:4] in years:
              dates.append(DIR[0:8])
          else:
            if DIR[0:8] in years:
              dates.append(DIR[0:8])
  for DIR in TDXDIRs:
    if DIR.endswith(filestring):
      tdxmin,tdxmax,tdymin,tdymax = geotiff.extent(TDXDIR+DIR)
      td_extent = shapely.geometry.Polygon([(tdxmin,wvymin),(tdxmin,tdymax),(tdxmax,wvymax),(tdxmax,tdymin)])
      if glacier.intersects(td_extent):
        if not(years) or (years=='all'):
          dates.append(DIR[0:8])
        else:
          if len(years) == 4:
            if DIR[0:4] in years:
              dates.append(DIR[0:8])
          else:
            if DIR[0:8] in years:
              dates.append(DIR[0:8])
 
  # Load data
  time = np.zeros(len(dates))
  zs = np.zeros([ny,nx,len(time)])
  for i in range(0,len(dates)):
    date = dates[i]
    n = 1 # Count number of files for that date
    time[i] = fracyear.date_to_fracyear(int(date[0:4]),int(date[4:6]),float(date[6:8]))
    for DIR in np.r_[TDXDIRs,WVDIRs]:
      if DIR.startswith(date) and DIR.endswith(filestring): 
        # Read file
        if  os.path.isfile(WVDIR+DIR):
          xwv,ywv,zwv = geotiff.read(WVDIR+DIR) 
          zwv[zwv==0] = float('NaN')
        else:
          xwv,ywv,zwv = geotiff.read(TDXDIR+DIR) 
          zwv[zwv==0] = float('NaN')
        
        # Interpolate onto output grid
        zwv_dem = scipy.interpolate.RegularGridInterpolator([ywv,xwv],zwv,bounds_error = False,method='linear',fill_value=float('nan'))
        zwv_flattened = zwv_dem(np.column_stack([ygrid.flatten(),xgrid.flatten()]))
    
        # Reshape to grid
        zwv_on_grid = np.reshape(zwv_flattened,(len(y),len(x)))
        zwv_on_grid[zwv_on_grid=='NaN']=0
        
        if n > 1: # there is more than one file for that date, so we need to treat it differently
          zs[:,:,i] = (zs[:,:,i]*(n-1)+zwv_on_grid)/n
        else:
          zs[:,:,i] = zwv_on_grid
          n = n+1
          
    zs[zs==0] = 'NaN' 
      
  # It's possible we still pulled a DEM full of NaNs in the desired region. If so, let's chuck those.
  nonnan=[]
  for i in range(0,len(time)):
    if len(np.where(~(np.isnan(zs[:,:,i])))[0]) > 5:
      nonnan.append(i)
  zs_nonnan = zs[:,:,nonnan]
  time_nonnan = time[nonnan]    
  
  # Sort by time
  sortind = np.argsort(time_nonnan)
  time_nonnan = time_nonnan[sortind]
  zs_nonnan = zs_nonnan[:,:,sortind]
      
  return x,y,zs_nonnan,time_nonnan

def worldview_along_flowline(xpts,ypts,glacier,years='all',cutoff='terminus',verticaldatum='geoid',filt_len='none',method='linear'):

  '''
  pts = worldview_along_flowline(xpts, ypts, glacier ,years='all', 
  		cutoff='terminus', verticaldatum='geoid', filt_len='none')
  
  Outputs worldview surface elevations along a flowline defined by xpts,ypts. Same as "atm_along_flowline" except
  for worldview data. See that function for more information.
  '''

  # Worldview data
  WVDIR = os.path.join(os.getenv("DATA_HOME"),"Elevation/Worldview/"+glacier+"/")
  WVDIRs = os.listdir(WVDIR)
  TDXDIR = os.path.join(os.getenv("DATA_HOME"),"Elevation/TDX/"+glacier+"/")
  TDXDIRs = os.listdir(TDXDIR)


  # Find file ending based on whether we want elevation relative to geoid or ellipsoid
  if verticaldatum == 'geoid':
    filestring = 'trans-adj.tif'
  else:
    filestring = 'trans.tif'
  
  # Set up output dictionary
  pts = {}
  
  # Load ice front positions so we can toss data in front of terminus
  if cutoff == 'terminus':
    dists = dist.transect(xpts,ypts)
    term_values, term_time = icefronts.distance_along_flowline(xpts,ypts,dists,glacier,type='icefront')
  
  
# Find dates where we have data in the desired region
  dates=[]
  for DIR in WVDIRs:
    if (DIR[0:8] not in dates) and DIR.startswith('2') and DIR.endswith(filestring):
      xmin,xmax,ymin,ymax = geotiff.extent(WVDIR+DIR)
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
  for DIR in TDXDIRs:
    if DIR.endswith(filestring):
      xmin,xmax,ymin,ymax = geotiff.extent(TDXDIR+DIR)
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


  for date in dates:
    for DIR in np.r_[TDXDIRs,WVDIRs]:
      if DIR.startswith(date) and (DIR.endswith(filestring)): 
        #print "Loading data from "+DIR+"\n"
        if os.path.isfile(WVDIR+DIR):
          x,y,z = geotiff.read(WVDIR+DIR)
          z[z == 0] = float('NaN')
        else:
          x,y,z = geotiff.read(TDXDIR+DIR)

        dem = scipy.interpolate.RegularGridInterpolator([y,x],z,method=method)
    
        # Find points that fall within the DEM
        ind = np.where((xpts > np.min(x)) & (xpts < np.max(x)) & (ypts > np.min(y)) & (ypts < np.max(y)))
        if ind:
          # Get elevation at coordinates
          zpts = np.zeros(len(xpts)); zpts[:] = 'NaN'
          zpts[ind] = dem(np.column_stack([ypts[ind],xpts[ind]]))
          
          # Get fractional year
          time = fracyear.date_to_fracyear(float(date[0:4]),float(date[4:6]),float(date[6:8]))

          if cutoff == 'terminus':
            # Get terminus position at time of worldview image
            termpos = np.interp(time,term_time,term_values)
            ind = np.where(dists > termpos)[0]
            zpts[ind] = 'NaN' 
          
          # Save points to output dictionary
          if len(np.where(~(np.isnan(zpts)))[0]) < 1:
            print "Not loading ",DIR
          elif date in pts.keys():
            nans = np.where((np.isnan(pts[date][:,2])))
            pts[date][nans,2] = zpts[nans] 
          else:
            pts[date] = np.column_stack([xpts,ypts,zpts])
  
  # Print warning if removing points in front of ice front
  if cutoff == 'terminus':
    print "Cleaning up DEM points by removing points in front of ice front."
    print "You can change this setting by setting `cutoff = 'none'.'"
  
  # Filter if desired
  if filt_len != 'none':
    cutoff=(1/filt_len)/(1/((dists[1]-dists[0])*2))
    b,a=signal.butter(4,cutoff,btype='low')
    for key in pts.keys():
      nonnan = np.where(~(np.isnan(pts[key][:,2])))[0]
      if len(nonnan) > 20:
        pts[key][nonnan,2]=signal.filtfilt(b,a,pts[key][nonnan,2])
  
  return pts

def worldview_at_pts(xpts,ypts,glacier,years='all',verticaldatum='geoid',cutoff='none',method='linear',radius=500):

  '''
  zpts,time = worldview_at_pts(xpts,ypts,glacier,
  				years='all',verticaldatum='geoid')
  
  Interpolates worldview surface elevations to points xpts,ypts.
  
  Inputs:
  xpts,ypts: points where we want surface elevations
  glacier: glacier name
  years: years that we want data
  verticaldatum: geoid or ellipsoid
  cutoff: cutoff elevations in front of terminus ('terminus' or 'none')
  method: 'linear' extrapolation or 'average' value for a region defined by radius
  radius: radius for average value (only necessary if method is 'average')
  
  Outputs:
  zpts: array of surface elevations for points xpts,ypts
  time: time of the arrays
  '''
  
  if method == 'linear':
    wv = worldview_along_flowline(xpts,ypts,glacier,years=years,cutoff=cutoff,verticaldatum=verticaldatum,method=method)
    time = np.zeros(len(dates))
    zpts = np.zeros([len(dates),len(xpts)])
    dates = np.sort(wv.keys())
    for i in range(0,len(dates)):
      date = dates[i]
      time[i] = fracyear.date_to_fracyear(float(date[0:4]),float(date[4:6]),float(date[6:8]))
      zpts[i,:] = wv[date][:,2]
  elif method == 'average':
    xmin = np.min(xpts)-radius*3
    ymin = np.min(ypts)-radius*3
    xmax = np.max(xpts)+radius*3
    ymax = np.max(ypts)+radius*3
    
    xwv,ywv,zwv,timewv = worldview_grid(glacier,xmin,xmax,ymin,ymax,years='all',verticaldatum='geoid')
    N = len(timewv)
    
    time = timewv
    zpts = np.zeros([N,len(xpts)])
    zpts_std = np.zeros([N,len(xpts)])
    zpts[:,:] = float('nan')
    zpts_std[:,:] = float('nan')
    
    xwv_grid,ywv_grid = np.meshgrid(xwv,ywv) 
    xwv_flattened = xwv_grid.flatten()
    ywv_flattened = ywv_grid.flatten()
    
    for j in range(0,len(xpts)):
      ind = np.where(np.sqrt((xpts[j] - xwv_flattened)**2+(ypts[j] - ywv_flattened)**2) < radius)[0]
      for i in range(0,N):
        nonnan = np.where(~(np.isnan(zwv[:,:,i].flatten()[ind])))[0]
        if len(nonnan) > 0.8*len(ind):
          zpts[i,j] = np.nanmean(zwv[:,:,i].flatten()[ind])
          zpts_std[i,j] = np.nanvar(zwv[:,:,i].flatten()[ind])
      
        
    
  return zpts,zpts_std,time

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
  TDXDIR = os.path.join(os.getenv("DATA_HOME"),"Elevation/TDX/"+glacier+"/")
  TDXDIRs = os.listdir(TDXDIR)


  # Find file ending based on whether we want elevation relative to geoid or ellipsoid
  if verticaldatum == 'geoid':
    filestring = 'trans-adj.tif'
  else:
    filestring = 'trans.tif'
  
# Find dates where we have data in the desired region
  besttime = 0
  for DIR in WVDIRs:
    if DIR.startswith('2') and DIR.endswith(filestring):
      demtime = fracyear.date_to_fracyear(int(DIR[0:4]),int(DIR[4:6]),int(DIR[6:8]))
      if abs(demtime-time) < abs(besttime-time):
        besttime = demtime
        bestfile = WVDIR+DIR
  for DIR in TDXDIRs:
    if DIR.startswith('2') and DIR.endswith(filestring):
      demtime = fracyear.date_to_fracyear(int(DIR[0:4]),int(DIR[4:6]),int(DIR[6:8]))
      if abs(demtime-time) < abs(besttime-time):
        besttime = demtime
        bestfile = TDXDIR+DIR
       
  x,y,zs = geotiff.read(bestfile)

  return x,y,zs,besttime

def dem_continuous(glacier,date,verticaldatum='geoid',fillin=False,blur=False):

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
  fillin: fill in no data locations
  
  Outputs:
  xg,yg: x and y coordinates for output grid
  zg: continuous surface DEM 
  '''
  
  # Get correct tile of the GIMP DEM
  if glacier == 'Helheim':
    subset = 'gimpdem3_1'
    xmin = 233000.0
    xmax = 318000.0
    ymin = -2601000.0
    ymax = -2515000.0
  elif glacier == 'Kanger':
    subset = 'gimpdem4_2'
    xmin = 420000.0
    xmax = 500000.0
    ymin = -2320000.0
    ymax = -2210000.0
  else: 
    sys.exit("Unknown glacier.")
  
  # Load correct file for geoid or ellipsoid heights
  if verticaldatum == 'geoid':
    fileend = '-adj.tif'
  else:
    fileend = '.tif'
  
  # Gimp filename  
  gimpfile = os.path.join(os.getenv("DATA_HOME"),'Elevation/Gimp/'+subset+fileend)
  
  # Directories for high res DEMs
  OUTDIR = os.path.join(os.getenv("DATA_HOME"),"Elevation/MosaicDEMs/"+glacier+"/")
  WVDIR = os.path.join(os.getenv("DATA_HOME"),"Elevation/Worldview/"+glacier+"/")
  TDXDIR = os.path.join(os.getenv("DATA_HOME"),"Elevation/Worldview/"+glacier+"/")
  
  # Select dates for worldview images 
  if glacier == 'Helheim':
    if date == '20120624':
      dates = ['20120624','20120513','20120520']
    elif date == '20110319':
      dates = ['20110319','20110615']
      dates_backup = ['20110628']
    elif date == '20110628':
      dates = ['20110628','20110615']
      dates_backup = ['20110319']
  elif glacier == 'Kanger':
    if date == '20110712':
      dates = ['20110712','20110808','20110823']
      dates_backup = ['20110824']
  
  
  files = ''
  dirs_wv = os.listdir(WVDIR)
  for d in dates:
    for dir in dirs_wv:
      if (d in dir) and (dir.endswith(fileend)):
        files = files+' '+WVDIR+dir
  files = files+' '+gimpfile  
  
  CURRENTDIR = os.getcwd()
  os.chdir(OUTDIR)
  filename = 'mosaic-'+date+'-'+verticaldatum
  if not(os.path.isfile(filename+'-tile-0.tif')):
    os.system('dem_mosaic --t_projwin '+str(xmin)+' '+str(ymin)+' '+str(xmax)+\
  		' '+str(ymax)+' --priority-blending-length 40 -o'+filename+files)
  
  xg,yg,zs = geotiff.read(filename+"-tile-0.tif")
  
  if blur == True:
    print "Blurring DEM over 17 pixels (roughly 500 m in each direction)..."
    # 17 pixel gaussian blur
    zs_blur = scipy.ndimage.filters.gaussian_filter(zs,sigma=2,truncate=4)
  
  os.chdir(CURRENTDIR)
  
  return xg,yg,zs
  
def dem_continuous_flowline(xf,yf,glacier,date,verticaldatum='geoid',fillin='none',filt_len='none'):

  '''
  zflow = dem_continuous_flowline(xf,yf,glacier,date,verticaldatum='geoid',fillin='none')

  Same as "dem_continuous", except output the coordinates along a flowline rather than as a 
  grid. See "dem_continuous" for more information.

  '''
  
  # Get grid for that date 
  xg,yg,zgrid = dem_continuous(glacier,date,verticaldatum,fillin,blur=False)
  
  # Interpolate grid onto flowline
  dem = scipy.interpolate.RegularGridInterpolator([yg,xg],zgrid)
   
  zflow = dem((yf,xf)) 
  
  # Filter flowline if filt_len is not set to none
  if filt_len != 'none':
    print "fix this"
    cutoff=(1/filt_len)/(1/((dists[1]-dists[0])*2))
    b,a=signal.butter(4,cutoff,btype='low')
    zflow_filtered = signal.filtfilt(b,a,zflow)
  else:
    zflow_filtered = zflow
   
  return zflow_filtered


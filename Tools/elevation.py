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
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
sys.path.append(os.path.join(os.getenv("HOME"),"Code/BigThreeGlaciers/Tools"))
import icefronts, glacier_flowline
import coords, geotiff
import scipy.interpolate, dist, fracyear
import scipy.signal as signal
import shapely.geometry

def gimp_at_pts(xpts,ypts,glacier,verticaldatum):
  
  if glacier == 'Helheim':
    file = os.path.join(os.getenv("HOME"),'Data/Elevation/Gimp/gimpdem3_1.tif')
  elif glacier == 'Kanger':
    file = os.path.join(os.getenv("HOME"),'Data/Elevation/Gimp/gimpdem4_2.tif')
  else: 
    sys.exit("Unknown glacier.")
  
  # Read GIMP DEM
  [x,y,z]=geotiff.read(file)
  
  # Interpolate DEM to points
  f = scipy.interpolate.RegularGridInterpolator((y,x),z,method="linear")
  zs = f(np.column_stack([ypts,xpts]))
  
  # Choose vertical datum
  if verticaldatum == "ellipsoid":
    elev = zs
  elif verticaldatum == "geoid":
    geoidheight = coords.geoidheight(xpts,ypts)
    elev = zs - geoidheight
  else:
    print "Unknown datum, defaulting to height above ellipsoid"
  
  return elev
  
def gimp_grid(xmin,xmax,ymin,ymax,glacier,verticaldatum):
  
  if glacier == 'Helheim':
    file = os.path.join(os.getenv("HOME"),'Data/Elevation/Gimp/gimpdem3_1.tif')
  elif glacier == 'Kanger':
    file = os.path.join(os.getenv("HOME"),'Data/Elevation/Gimp/gimpdem4_2.tif')
  else: 
    sys.exit("Unknown glacier.")
    
  # Get GIMP DEM
  [x,y,zs]=geotiff.read(file,xmin,xmax,ymin,ymax)
  
  if verticaldatum == "ellipsoid":
    elev = np.array(zs)
  elif verticaldatum == "geoid":
    (xgrid,ygrid)=np.meshgrid(x,y)
    geoidheight = coords.geoidheight(xgrid.flatten(),ygrid.flatten())
    zg = np.reshape(geoidheight,(len(y),len(x)))
    elev = zs - zg
  else: 
    print "Unrecognized vertical datum."
    
  return x,y,elev

def atm(years,verticaldatum):
  import coords
  
  # Directory for ATM
  ATMDIR = os.path.join(os.getenv("HOME"),'Data/Elevation/ATM/HelheimKanger/')
  DIRs = os.listdir(ATMDIR)
  
  atm = {}
  
  # Load only desired data
    
  for DIR in DIRs:
    if (DIR.startswith('2') or DIR.startswith('1')) and ((years == 'all') or (DIR[0:4] in years)):
      print "Loading",DIR,"..."
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
      geoidheight = coords.geoidheight(atm[key][:,0],atm[key][:,1])
      atm[key][:,2] = atm[key][:,2] - geoidheight
  elif verticaldatum=='ellipsoid':
    atm = atm
  else: 
    print "Don't recognize that vertical datum, defaulting to geoid"
    
  return atm

def atm_along_flowline(xpts,ypts,glacier,years='all',cutoff='none',maxdist=200,verticaldatum='geoid'):
  
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
    z = np.zeros(R); z[:] = 'NaN'
    d = np.zeros(R); d[:] = 'NaN'
    for i in range(0,R):
      ind = np.argmin((xpts[i]-data[date][:,0])**2 + (ypts[i]-data[date][:,1])**2)
      xatm = data[date][ind,0]
      yatm = data[date][ind,1]
      d[i] = dist.between_pts(xpts[i],ypts[i],xatm,yatm)
      if d[i] < maxdist:
        z[i] = data[date][ind,2]
        
    if cutoff == 'terminus':
      # Get fractional year
      time = fracyear.date_to_fracyear(float(date[0:4]),float(date[4:6]),float(date[6:8]))
    
      # Get terminus position when ATM was collected
      termpos = np.interp(time,term_time,term_values)
      ind = np.where(dists > termpos)[0]
      z[ind] = 'NaN' 
    
    # Save elevations if there are any
    if len(np.where(~(np.isnan(z)))[0]) > 0:
      pts[date] = np.column_stack([xpts,ypts,z,d])
  
  if cutoff == 'terminus':
    print "Cleaning up DEM points by removing points in front of ice front."
    print "You can change this setting by setting `cutoff = 'none'.'"
	    
  return pts
  
def atm_at_pts(xpts,ypts,glacier,years='all',maxdist=200,verticaldatum='geoid'):

  # Get ATM data
  atm = atm_along_flowline(xpts,ypts,glacier,years,'none',maxdist,verticaldatum)
  
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
   
def worldview_grid(glacier,xmin=0,xmax=0,ymin=0,ymax=0,resolution=32,years='all',verticaldatum='geoid'):

  # Inputs
  # xmin,xmax,ymin,ymax: extent for output grids
  # years: year that you want data or "all" for all years 
  # glacier: Kanger or Helheim
  # verticaldatum: ellipsoid or geoid
  
  # Output: 
  # x,y,zs_nonnan,time_nonnan: numpy array of surface elevations on the grid x,y for times_nonnan

  # Set up output grid
  dx = dy = float(resolution)
  nx = np.ceil((xmax-xmin)/dx)+1
  x = np.linspace(xmin,(nx-1)*dx+xmin,nx)
  ny = np.ceil((ymax-ymin)/dx)+1
  y = np.linspace(ymin,(ny-1)*dy+ymin,ny)
  xgrid,ygrid = np.meshgrid(x,y)
    
  WVDIR = os.path.join(os.getenv("HOME"),'/Users/kehrl/Data/Elevation/Worldview/'+glacier+'/')
  DIRs = os.listdir(WVDIR)
  
  # Find dates where we have data in the desired region
  dates=[]
  glacier = shapely.geometry.Polygon([(xmin,ymin),(xmin,ymax),(xmax,ymax),(xmax,ymin)])
  for DIR in DIRs:
    if (DIR[0:8] not in dates) and DIR.startswith('2') and (DIR.endswith('_32m.tif') or DIR.endswith('32m_trans.tif')):
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
 
  # Load data
  time = np.zeros(len(dates))
  zs = np.zeros([ny,nx,len(time)])
  for i in range(0,len(dates)):
    date = dates[i]
    n = 1 # Count number of files for that date
    time[i] = fracyear.date_to_fracyear(int(date[0:4]),int(date[4:6]),float(date[6:8]))
    for DIR in DIRs:
      if DIR.startswith(date) and (DIR.endswith('_32m.tif') or DIR.endswith('32m_trans.tif')): 
        # Read file
        xwv,ywv,zwv = geotiff.read(WVDIR+DIR) 
        zwv[zwv==0] = 'NaN'
        
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
  
  # Set elevation relative to vertical datum
  geoidheight = coords.geoidheight(xgrid.flatten(),ygrid.flatten())
  geoidheight = np.reshape(geoidheight,(ny,nx))
  if verticaldatum == 'geoid':
    for i in range(0,len(time)):
      zs[:,:,i] = zs[:,:,i] - geoidheight
      
  # It's possible we still pulled a DEM full of NaNs in the desired region. If so, let's chuck those.
  nonnan=[]
  for i in range(0,len(time)):
    if len(np.where(~(np.isnan(zs[:,:,i])))[0]) > 5:
      nonnan.append(i)
  zs_nonnan = zs[:,:,nonnan]
  time_nonnan = time[nonnan]    
      
  return x,y,zs_nonnan,time_nonnan

def worldview_along_flowline(xpts,ypts,glacier,years='all',cutoff='terminus',verticaldatum='geoid',filt_len='none'):

  # Worldview data
  WVDIR = os.path.join(os.getenv("HOME"),"/Users/kehrl/Data/Elevation/Worldview/"+glacier+"/")
  DIRs = os.listdir(WVDIR)
  
  # Set up output dictionary
  pts = {}
  
  # Load ice front positions so we can toss data in front of terminus
  if cutoff == 'terminus':
    dists = dist.transect(xpts,ypts)
    term_values, term_time = icefronts.distance_along_flowline(xpts,ypts,dists,glacier,'icefront')
  
  
# Find dates where we have data in the desired region
  dates=[]
  for DIR in DIRs:
    if (DIR[0:8] not in dates) and DIR.startswith('2') and DIR.endswith('.tif'):
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

  for date in dates:
    for DIR in DIRs:
      if DIR.startswith(date) and (DIR.endswith('_32m.tif') or DIR.endswith('32m_trans.tif')): 
        #print "Loading data from "+DIR+"\n"
        x,y,z = geotiff.read(WVDIR+DIR)
        z[z == 0] ='NaN'

        dem = scipy.interpolate.RegularGridInterpolator([y,x],z)
    
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
	    
  # Choose if elevation should be relative to geoid or ellipsoid
  if verticaldatum == 'geoid':
    geoidheight = coords.geoidheight(xpts,ypts)
    for key in pts.keys():
      pts[key][:,2] = pts[key][:,2] - geoidheight
  elif verticaldatum == 'ellipsoid':
     pts = pts
  else:
    print "Unknown datum, defaulting to height above ellipsoid"
  
  # Filter if desired
  if filt_len != 'none':
    cutoff=(1/filt_len)/(1/((dists[1]-dists[0])*2))
    b,a=signal.butter(4,cutoff,btype='low')
    for key in pts.keys():
      nonnan = np.where(~(np.isnan(pts[key][:,2])))[0]
      if len(nonnan) > 20:
        pts[key][nonnan,2]=signal.filtfilt(b,a,pts[key][nonnan,2])
  
  return pts

def worldview_at_pts(xpts,ypts,glacier,years='all',verticaldatum='geoid'):

  wv = worldview_along_flowline(xpts,ypts,glacier,years,'none',verticaldatum)
  
  dates = np.sort(wv.keys())
  time = np.zeros(len(dates))
  zpts = np.zeros([len(dates),len(xpts)])
  
  for i in range(0,len(dates)):
    date = dates[i]
    time[i] = fracyear.date_to_fracyear(float(date[0:4]),float(date[4:6]),float(date[6:8]))
    zpts[i,:] = wv[date][:,2]

  return zpts,time

def dem_continuous(xmin,xmax,ymin,ymax,glacier,date,verticaldatum='geoid',fillin='none'):

  # Load gimp DEM
  xg,yg,zg = gimp_grid(xmin,xmax,ymin,ymax,glacier,verticaldatum)
  
  # Create output
  zs = np.zeros_like(zg)
  zs = np.array(zg)
  
  # Select dates for worldview images 
  if glacier == 'Helheim':
    if date == '20120624':
      dates = ['20120624','20120513','20120520']
      dates_backup = ['20110319','20110615']
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
  
  # Load worldview images and update dates file in case there are multiple 
  # DEMs for a single date
  xwv,ywv,zwv,timewv = worldview_grid(glacier,np.min(xg),np.max(xg),np.min(yg),np.max(yg),resolution=xg[1]-xg[0],years=dates,verticaldatum=verticaldatum)
  xwvb,ywvb,zwvb,timewvb = worldview_grid(glacier,np.min(xg),np.max(xg),np.min(yg),np.max(yg),resolution=xg[1]-xg[0],years=dates_backup,verticaldatum=verticaldatum)  
    
  # Use worldview DEM
  for i in range(0,len(xg)):
    for j in range(0,len(yg)):
      nonnan = np.where(~(np.isnan(zwv[j,i,:])))[0]
      if len(nonnan) > 0:
        zs[j,i] = np.mean(zwv[j,i,nonnan])
      elif (len(np.where(~(np.isnan(zwvb[j,i,:])))[0])) > 0 and (fillin != 'none'):
        nonnan = np.where(~(np.isnan(zwvb[j,i,:])))[0]
        zs[j,i] = np.mean(zwvb[j,i,nonnan])
  
  return xg,yg,zs
  
def dem_continuous_flowline(xf,yf,glacier,date,verticaldatum='geoid',fillin='none'):

  # Get grid dimensions for dem_continuous
  xmin = np.min(xf)-10.0e3
  xmax = np.max(xf)+10.0e3
  ymin = np.min(yf)-10.0e3
  ymax = np.max(yf)+10.0e3
  
  # Get grid for that date 
  xg,yg,zgrid = dem_continuous(xmin,xmax,ymin,ymax,glacier,date,verticaldatum,fillin)
  
  # Interpolate grid onto flowline
  dem = scipy.interpolate.RegularGridInterpolator([yg,xg],zgrid)
   
  zflow = dem((yf,xf)) 
   
  return zflow
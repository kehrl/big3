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
import icefronts
import coords, geotiff
import scipy.interpolate, jdcal, dist

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
      os.chdir(ATMDIR+DIR)
      files = os.listdir(ATMDIR+DIR)
      for file in files:
        print file
        if 'nadir' in file and not(file.endswith('.xml')):
          try:
            data=np.loadtxt(file,comments='#')
          except:
            data=np.loadtxt(file,comments='#',delimiter=',')
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
    print "Don't recognize that vertical datum, defaulting to ellipsoid"
    
  return atm

def atm_along_flowline(xpts,ypts,glacier,years='all',cutoff='none',maxdist=500,verticaldatum='geoid'):
  
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
      year = float(date[0:4])
      day = jdcal.gcal2jd(year,float(date[4:6]),float(date[6:8]))
      day2 = jdcal.gcal2jd(year,12,31)
      day1 = jdcal.gcal2jd(year-1,12,31)
      doy = day[1]+day[0]-day1[1]-day1[0]
      time = ( year + doy/(day2[1]+day2[0]-day1[0]-day1[1])) 
    
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
   
def worldview_grid(years,glacier,verticaldatum):

  # Inputs
  # years: year that you want data or "all" for all years 
  # resolution: options are 2 or 32
  # verticaldatum: ellipsoid or geoid
  
  # Output: dictionary including all available data for the chosen years, with the keys 
  #  indicates the dates of each grid

  worldview = {}
  
  WVDIR = os.path.join(os.getenv("HOME"),'/Users/kehrl/Data/Elevation/Worldview/HelheimKanger/')
  DIRs = os.listdir(WVDIR)
  
  os.chdir(WVDIR)
  
  # Find dates where we have data in the desired region
  dates=[]
  for DIR in DIRs:
    if (DIR[0:8] not in dates) and DIR.startswith('2'):
      xmin,xmax,ymin,ymax = geotiff.extent(DIR)
      if (np.max(x) > xmin) and (np.min(x) < xmax) and (np.max(y) > ymin) and (np.min(y) > ymax):
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
      if DIR.startswith(date) and DIR.endswith('_32m_trans.tif'): 
        print "Loading data from "+DIR+"\n"
        x,y,zs = geotiff.read(DIR) 
        zs[zs==0] = 'NaN'
        worldview[date]=x,y,zs
  
  if verticaldatum == 'geoid':
     for key in worldview.keys():
       nx = len(worldview[key][0])
       ny = len(worldview[key][1])
       x_grid,y_grid = np.meshgrid(worldview[key][0],worldview[key][1])
       geoidheight = coords.geoidheight(x_grid.flatten(),y_grid.flatten())
       worldview[key][2][:,:] = worldview[key][2] - np.reshape(geoidheight,(ny,nx))
  elif verticaldatum == 'ellipsoid':
    worldview = worldview
  else:
    print "Unknown datum, defaulting to ellipsoid"
  
    
  return worldview

def worldview_along_flowline(xpts,ypts,glacier,years='all',cutoff='terminus',verticaldatum='geoid'):

  # Worldview data
  WVDIR = os.path.join(os.getenv("HOME"),"/Users/kehrl/Data/Elevation/Worldview/HelheimKanger/")
  DIRs = os.listdir(WVDIR)
  
  os.chdir(WVDIR)
  
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
      xmin,xmax,ymin,ymax = geotiff.extent(DIR)
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
      if DIR.startswith(date) and DIR.endswith('_32m.tif'): 
        #print "Loading data from "+DIR+"\n"
        x,y,z = geotiff.read(DIR)
        z[z == 0] ='NaN'

        dem = scipy.interpolate.RegularGridInterpolator([y,x],z)
    
        # Find points that fall within the DEM
        ind = np.where((xpts > np.min(x)) & (xpts < np.max(x)) & (ypts > np.min(y)) & (ypts < np.max(y)))
        if ind:
          # Get elevation at coordinates
          zpts = np.zeros(len(xpts)); zpts[:] = 'NaN'
          zpts[ind] = dem(np.column_stack([ypts[ind],xpts[ind]]))
          
          # Get fractional year
          year = float(date[0:4])
          day = jdcal.gcal2jd(year,float(date[4:6]),float(date[6:8]))
          day2 = jdcal.gcal2jd(year,12,31)
          day1 = jdcal.gcal2jd(year-1,12,31)
          doy = day[1]+day[0]-day1[1]-day1[0]
          time = ( year + doy/(day2[1]+day2[0]-day1[0]-day1[1])) 
          
          if cutoff == 'terminus':
            # Get terminus position at time of worldview image
            termpos = np.interp(time,term_time,term_values[:,0])
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
  
  return pts

def worldview_at_pts(xpts,ypts,glacier,years='all',verticaldatum='geoid'):

  wv = worldview_along_flowline(xpts,ypts,glacier,years,'none',verticaldatum)
  
  dates = np.sort(wv.keys())
  time = np.zeros(len(dates))
  zpts = np.zeros([len(dates),len(xpts)])
  
  for i in range(0,len(dates)):
    date = dates[i]
    year = float(date[0:4])
    day = jdcal.gcal2jd(year,float(date[4:6]),float(date[6:8]))
    day2 = jdcal.gcal2jd(year,12,31)
    day1 = jdcal.gcal2jd(year-1,12,31)
    doy = day[1]+day[0]-day1[1]-day1[0]
    time[i] = ( year + doy/(day2[1]+day2[0]-day1[0]-day1[1]))  
    zpts[i,:] = wv[date][:,2]

  return zpts,time

def dem_continuous(xmin,xmax,ymin,ymax,glacier,date,verticaldatum):

  # Load gimp DEM
  xg,yg,zg = gimp_grid(xmin,xmax,ymin,ymax,glacier,verticaldatum)
  
  # Create output
  zs = np.zeros_like(zg)
  zs = np.array(zg)
  
  # Load worldview DEMs 
  if date == '20120624' and glacier == 'Helheim':
    dates = ['20120624','20120513','20120520']
    dates_backup = ['20110319','20110615']
    wv = worldview_grid(dates,glacier,verticaldatum)
    wv_backup = worldview_grid(dates_backup,glacier,verticaldatum)
  
  # Interpolate worldview DEMs onto gimp DEM
  xgrid,ygrid = np.meshgrid(xg,yg)
  
  zwv_on_grid = np.zeros([len(yg),len(xg),len(dates)])
  n = 0
  for date in dates:
    # Get data
    xwv = wv[date][0]
    ywv = wv[date][1]
    zwv = wv[date][2]
    zwv[zwv == 0] = 'NaN'
    
    # Interpolate
    zwv_dem = scipy.interpolate.RegularGridInterpolator([ywv,xwv],zwv,bounds_error = False)
    zwv_flattened = zwv_dem(np.column_stack([ygrid.flatten(),xgrid.flatten()]))
    
    # Reshape to grid
    zwv_on_grid[:,:,n] = np.reshape(zwv_flattened,(len(yg),len(xg)))
    n = n+1
  
  zwv_on_grid_backup = np.zeros([len(yg),len(xg),len(dates_backup)])
  n = 0  
  for date in dates_backup:
    # Get data
    xwv = wv_backup[date][0]
    ywv = wv_backup[date][1]
    zwv = wv_backup[date][2]
    zwv[zwv == 0] = 'NaN'
    
    # Interpolate
    zwv_dem = scipy.interpolate.RegularGridInterpolator([ywv,xwv],zwv,bounds_error = False)
    zwv_flattened = zwv_dem(np.column_stack([ygrid.flatten(),xgrid.flatten()]))
    
    # Reshape to grid
    zwv_on_grid_backup[:,:,n] = np.reshape(zwv_flattened,(len(yg),len(xg)))
    n = n+1  
    
  # Find extent of primary worldview images  
  xmin_wv = np.min([np.min(wv[date][0]) for date in dates])
  xmax_wv = np.max([np.max(wv[date][0]) for date in dates])
  ymin_wv = np.max([np.min(wv[date][1]) for date in dates])
  ymax_wv = np.max([np.max(wv[date][1]) for date in dates])
  
  xind1 = np.argmin(abs(xg-xmin_wv))
  xind2 = np.argmin(abs(xg-xmax_wv))
  yind1 = np.argmin(abs(yg-ymin_wv))
  yind2 = np.argmin(abs(yg-ymax_wv))
  
  # Use worldview DEM
  for i in range(xind1,xind2):
    for j in range(yind1,yind2):
      nonnan = np.where(~(np.isnan(zwv_on_grid[j,i,:])))[0]
      if len(nonnan) > 0:
        zs[j,i] = np.mean(zwv_on_grid[j,i,nonnan])
      elif len(np.where(~(np.isnan(zwv_on_grid_backup[j,i,:])))[0]) > 0:
        nonnan = np.where(~(np.isnan(zwv_on_grid_backup[j,i,:])))[0]
        zs[j,i] = zwv_on_grid_backup[j,i,nonnan]
        print zg[j,i], zs[j,i]
  
  return xg,yg,zs
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

def gimp_pts(xpts,ypts,glacier,verticaldatum):
  
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
    elev = zs
  elif verticaldatum == "geoid":
    (xgrid,ygrid)=np.meshgrid(x,y)
    geoidheight = coords.geoidheight(xgrid.flatten(),ygrid.flatten())
    zg = np.reshape(geoidheight,(len(y),len(x)))
    elev = zs - zg
    
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

def atm_at_pts(xpts,ypts,years,maxdist,verticaldatum):
  
  # Get ATM data
  data = atm(years,verticaldatum)
  
  # Set up output as dictionary
  pts = {}
  
  keys = data.keys()
  R = len(xpts)
  for key in keys:
    z = np.zeros(R); z[:] = 'NaN'
    d = np.zeros(R); d[:] = 'NaN'
    for i in range(0,R):
      ind = np.argmin((xpts[i]-data[key][:,0])**2 + (ypts[i]-data[key][:,1])**2)
      xatm = data[key][ind,0]
      yatm = data[key][ind,1]
      d[i] = dist.between_pts(xpts[i],ypts[i],xatm,yatm)
      if d[i] < maxdist:
        z[i] = data[key][ind,2]
    pts[key] = np.column_stack([xpts,ypts,z,d])
  
  return pts
   
def worldview_grid(years,resolution,verticaldatum):

  # Inputs
  # years: year that you want data or "all" for all years 
  # resolution: options are 2 or 32
  # verticaldatum: ellipsoid or geoid
  
  # Output: dictionary including all available data for the chosen years, with the keys 
  #  indicates the dates of each grid

  import coords

  worldview = {}
  
  WVDIR = os.path.join(os.getenv("HOME"),'/Users/kehrl/Data/Elevation/Worldview/Helheim/')
  DIRs = os.listdir(WVDIR)
  os.chdir(WVDIR)
  
  dates=[]
  for DIR in DIRs:
    if (DIR[0:8] not in dates) and DIR.startswith('2'):
      if not(years) or (years=='all'):
        dates.append(DIR[0:8])
      else: 
        if DIR[0:4] in years:
          dates.append(DIR[0:8])
 
  for date in dates:
    if resolution == 32:
      for DIR in DIRs:
        if DIR.startswith(date) and DIR.endswith('_32m_trans.tif'): 
          print "Loading data from "+DIR+"\n"
          worldview[date] = geotiff.read(DIR) 
    elif resolution == 2:
      for DIR in DIRs:
        if DIR.startswith(date) and DIR.endswith('_tr4x_align'):
          print "Loading data from "+DIR+"\n"
          worldview[date] = geotiff.read(DIR+"/"+DIR[0:56]+"-trans_reference-DEM.tif")
  
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

def worldview_pts(xpts,ypts,glacier,resolution,years,verticaldatum):

  # Worldview data
  WVDIR = os.path.join(os.getenv("HOME"),'/Users/kehrl/Data/Elevation/Worldview/Helheim/')
  DIRs = os.listdir(WVDIR)
  os.chdir(WVDIR)
  
  # Load ice front positions so we can toss data in front of terminus
  dists = dist.transect(xpts,ypts)
  term_values, term_time = icefronts.distance_along_flowline(xpts,ypts,dists,'icefront',glacier)
  
  dates=[]
  for DIR in DIRs:
    if (DIR[0:8] not in dates) and DIR.startswith('2'):
      if not(years) or (years=='all'):
        dates.append(DIR[0:8])
      else: 
        if DIR[0:4] in years:
          dates.append(DIR[0:8])
 
  n = 0
  time = np.zeros(len(dates))
  zpts = np.zeros([len(xpts),len(dates)])
  zpts[:,:] = 'NaN'
  for date in dates:
    if resolution == 32:
      for DIR in DIRs:
        if DIR.startswith(date) and DIR.endswith('_32m_trans.tif'): 
          print "Loading data from "+DIR+"\n"
          data = geotiff.read(DIR)
          x = data[0]
          y = data[1]
          z = data[2]
          z[z == 0] ='NaN'

          dem = scipy.interpolate.RegularGridInterpolator([y,x],z)
    
          # Find points that fall within the DEM
          ind = np.where((xpts > np.min(x)) & (xpts < np.max(x)) & (ypts > np.min(y)) & (ypts < np.max(y)))
          if ind:
            # Get fractional year
            year = float(date[0:4])
            day = jdcal.gcal2jd(year,float(date[4:6]),float(date[6:8]))
            day2 = jdcal.gcal2jd(year,12,31)
            day1 = jdcal.gcal2jd(year-1,12,31)
            doy = day[1]+day[0]-day1[1]-day1[0]
            time[n] = ( year + doy/(day2[1]+day2[0]-day1[0]-day1[1])) 
            
            # Get terminus position at time of worldview image
            terminus = np.interp(time[n],term_time,term_values[:,0])
            
            # Get elevation at coordinates
            zpts[ind,n] = dem(np.column_stack([ypts[ind],xpts[ind]]))
            ind = np.where(dists > terminus)
            zpts[ind,n] = 'NaN' 
            n = n+1
    	    
    elif resolution == 2:
      for DIR in DIRs:
        if DIR.startswith(date) and DIR.endswith('_tr4x_align'):
          print "Loading data from "+DIR+"\n"
          data = geotiff.read(DIR+"/"+DIR[0:56]+"-trans_reference-DEM.tif")
          x = data[0]
          y = data[1]
          z = data[2]
          z[z == 0] ='NaN'

          dem = scipy.interpolate.RegularGridInterpolator([y,x],z)
    
          # Find points that fall within the DEM
          ind = np.where((xpts > np.min(x)) & (xpts < np.max(x)) & (ypts > np.min(y)) & (ypts < np.max(y)))
          if ind:
            # Get fractional year
            year = float(date[0:4])
            day = jdcal.gcal2jd(year,float(date[4:6]),float(date[6:8]))
            day2 = jdcal.gcal2jd(year,12,31)
            day1 = jdcal.gcal2jd(year-1,12,31)
            doy = day[1]+day[0]-day1[1]-day1[0]
            time[n] = ( year + doy/(day2[1]+day2[0]-day1[0]-day1[1])) 
            
            # Get terminus position at time of worldview image
            termind = np.argmin(abs(term_time-time[n]))
            if term_time[termind] < time[n]:
              terminus = np.max(term_values[termind],term_values[termind+1])
            else:
              terminus = np.max(term_values[termind-1],term_values[termind])
            
            # Get elevation at coordinates
            zpts[ind,n] = dem(np.column_stack([ypts[ind],xpts[ind]]))
            ind = np.where(dists > terminus)
            zpts[ind,n] = 'NaN'
    	    n = n+1

  # Choose if elevation should be relative to geoid or ellipsoid
  if verticaldatum == 'geoid':
    geoidheight = coords.geoidheight(xpts,ypts)
    for i in range(0,len(time)):
      zpts[:,i] = zpts[:,i] - geoidheight
  elif verticaldatum == 'ellipsoid':
     zpts = zpts
  else:
    print "Unknown datum, defaulting to height above ellipsoid"
  return zpts,time
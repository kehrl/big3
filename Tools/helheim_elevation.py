# This script contains all elevation products for Helheim, so that we are consistent
# between scripts. 

# Functions:
# gimp_pts(xpts,ypts,verticaldatum): gimp elevations at xpts,ypts
# gimp_grid(xmin,xmax,ymin,ymax,verticaldatum): gimp elevations in the grid specified by
#	xmin,xmax,ymin,ymax
# atm(year,verticaldatum): atm data for chosen year or "all" years

import os
import sys
import numpy as np
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Helheim/Tools"))
import helheim_icefronts
import coords, geotiff
import scipy.interpolate, jdcal, dist

def gimp_pts(xpts,ypts,verticaldatum):
  file = os.path.join(os.getenv("HOME"),'Data/Elevation/Gimp/gimpdem3_1.tif')
  [x,y,z]=geotiff.read(file)
  
  f = scipy.interpolate.RegularGridInterpolator((y,x),z,method="linear")
  zs = f(np.column_stack([ypts,xpts]))
  
  if verticaldatum == "ellipsoid":
    elev = zs
  elif verticaldatum == "geoid":
    geoidheight = coords.geoidheight(xpts,ypts)
    elev = zs - geoidheight
  else:
    print "Unknown datum, defaulting to height above ellipsoid"
  
  return elev
  
def gimp_grid(xmin,xmax,ymin,ymax,verticaldatum):
  
  file = os.path.join(os.getenv("HOME"),'Data/Elevation/Gimp/gimpdem3_1.tif')
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
  
  atm = {}
  
  ATMDIR = os.path.join(os.getenv("HOME"),'Data/Elevation/ATM/Helheim/')
  DIRs = os.listdir(ATMDIR)
  
  if (not years) or (years=='all'):
    for DIR in DIRs:
      if DIR.startswith('2') or DIR.startswith('1'):
        x=[]
        y=[]
        z=[]
        os.chdir(ATMDIR+DIR)
        files = os.listdir(ATMDIR+DIR)
        for file in files:
          if (file.endswith('nadir5seg')) or (file.endswith('nadir3seg')):
            data = np.loadtxt(file)
          elif (file.endswith('nadir3seg_50pt.csv')):
            data = np.loadtxt(file,skiprows=10)
          xfile = data[:,2]
          yfile = data[:,1]
          zfile = data[:,3]
          x=np.hstack([x,xfile])
          y=np.hstack([y,yfile])
          z=np.hstack([z,zfile])
        print "Loading ", DIR, "..."
        x2,y2 = coords.convert(x-360,y,4326,3413)
        if DIR[0:4] in atm.keys():
          print "Already data from that year, consider changing how you have labeled the directories"
        else:
          atm[DIR] = np.column_stack([x2,y2,z])
  else:
    for DIR in DIRs:
      if DIR.startswith(years):
        x=[]
        y=[]
        z=[]
        os.chdir(ATMDIR+DIR)
        files = os.listdir(ATMDIR+DIR)
        for file in files:
          if 'nadir' in file:
            try:
              data=np.loadtxt(file,comments='#')
            except:
              data=np.loadtxt(file,comments='#',delimiter=',')
              
            xfile = data[:,2]
            yfile = data[:,1]
            zfile = data[:,3]
            x=np.hstack([x,xfile])
            y=np.hstack([y,yfile])
            z=np.hstack([z,zfile])
        x2,y2 = coords.convert(x-360,y,4326,3413)
        if DIR[0:4] in atm.keys():
          print "Already data from that year, consider changing how you have labeled the directories"
        else:
          atm[DIR] = np.column_stack([x2,y2,z])
  
  if verticaldatum=='geoid':
    for key in atm.keys():
      geoidheight = coords.geoidheight(atm[key][:,0],atm[key][:,1])
      atm[key][:,2] = atm[key][:,2] - geoidheight
  elif verticaldatum=='ellipsoid':
    atm = atm
  else: 
    print "Don't recognize that vertical datum, defaulting to ellipsoid"
    
  return atm

def atm_at_pts(x,y,years):

  atm = atm(years)
  
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

def worldview_pts(xpts,ypts,resolution,years,verticaldatum):

  # Worldview data
  WVDIR = os.path.join(os.getenv("HOME"),'/Users/kehrl/Data/Elevation/Worldview/Helheim/')
  DIRs = os.listdir(WVDIR)
  os.chdir(WVDIR)
  
  # Load ice front positions so we can toss data in front of terminus
  dists = dist.transect(xpts,ypts)
  term_values, term_time = helheim_icefronts.distance_along_flowline(xpts,ypts,dists,'icefront')
  
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
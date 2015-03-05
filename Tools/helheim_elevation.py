import os
import sys
import numpy as np
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
import coords, geotiff
import scipy.interpolate, jdcal

def atm(years):

  atm = {}
  
  ATMDIR = os.path.join(os.getenv("HOME"),'/Users/kehrl/Data/Elevation/ATM/Helheim/')
  DIRs = os.listdir(ATMDIR)
  
  if (not years) or (years=='all'):
    print "not working"
    for DIR in DIRs:
      if DIR.startswith('2') or DIR.startswith('1'):
        x=[]
        y=[]
        z=[]
        os.chdir(ATMDIR+DIR)
        files = os.listdir(ATMDIR+DIR)
        for file in files:
          if (file.endswith('nadir5seg')) or (file.endswith('nadir3seg')):
            data=np.loadtxt(file)
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
    
  return atm

def atm_at_pts(x,y,years):

  atm = atm(years)
  
  return pts
 
  
def worldview(years,resolution):

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
            worldview[date] = geotiff.read(DIR) 
      elif resolution == 2:
        for DIR in DIRs:
          if DIR.startswith(date) and DIR.endswith('_tr4x_align'):
            worldview[date] = geotiff.read(DIR+"/"+DIR[0:56]+"-trans_reference-DEM.tif")
    
  return worldview


def worldview_at_pts(xpts,ypts,resolution,years):

  data = worldview(years,resolution)

  dates = data.keys()
  n = 0
  time = np.zeros(len(dates))
  zpts = np.zeros([len(xpts),len(dates)])
  zpts[:,:] = 'NaN'
  for date in dates:
    # Get DEM measurements for that date
    x = data[date][0]
    y = data[date][1]
    z = data[date][2]
    z[z == 0] ='NaN'
    
    # Create DEM interpolation product
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
      
      # Get elevation at coordinates
      zpts[ind,n] = dem(np.column_stack([ypts[ind],xpts[ind]]))  
      n = n+1
      
  return time,zpts
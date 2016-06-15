# This module deals with the climate products.

import os
import numpy as np
import netCDF4, jdcal
import coordlib, datelib
import scipy.signal as signal
import scipy
from scipy.spatial import cKDTree

def racmo_grid(xmin,xmax,ymin,ymax,variable,epsg=3413,maskvalues='ice'):

  ''' 
  Pull all values for RACMO smb, t2m, zs, or runoff values for the region defined 
  by xmin,xmax,ymin,ymax.
  
  xrac_subset,yrac_subset,var_subset,time = racmo_grid(xmin,xmax,ymin,ymax,
  		variable,epsg=3413,mask='ice')
  
  Inputs:
  xmin,xmax,ymin,ymax : region where you want racmo data
  variable : what variable you want (runoff, t2m, zs, smb)
  maskvalues : if you want only 'ice' or 'notice' or 'both' values
  
  Outputs:
  xrac_subet,yrac_subset : output x,y
  var_subset: value of chosen variable at these points
  time : time
  '''

  # RACMO data
  if variable != 'zs':
    files = [(os.path.join(os.getenv("DATA_HOME"),"Climate/RACMO/2015_09_Laura_Kehrl/RACMO2.3_GRN11_"+variable+"_daily_2001_2010.nc")), \
           (os.path.join(os.getenv("DATA_HOME"),"Climate/RACMO/2015_09_Laura_Kehrl/RACMO2.3_GRN11_"+variable+"_daily_2011_2014.nc")), \
           (os.path.join(os.getenv("DATA_HOME"),"Climate/RACMO/2015_09_Laura_Kehrl/RACMO2.3_GRN11_"+variable+"_daily_2015.nc"))]
  else:
    files = [(os.path.join(os.getenv("DATA_HOME"),"Climate/RACMO/ZS_ZGRN_V5_1960-2014_detrended_2day.nc"))]

  rec1 = netCDF4.Dataset(files[0])
  if variable != 'zs':
    rec2 = netCDF4.Dataset(files[1])
    rec3 = netCDF4.Dataset(files[2])
  mask = netCDF4.Dataset(os.path.join(os.getenv("DATA_HOME"),"Climate/RACMO/2015_09_Laura_Kehrl/RACMO23_masks_ZGRN11.nc")).variables['icemask'][:]

  # Load RACMO data
  lat = np.array(rec1.variables['lat'][:])
  lon = np.array(rec1.variables['lon'][:])
  var1 = np.array(rec1.variables[variable][:])
  daysfrom1950_1 = np.array(rec1.variables['time'][:])
  if variable != 'zs':
    var2 = np.array(rec2.variables[variable][:])
    daysfrom1950_2 = np.array(rec2.variables['time'][:])
    var3 = np.array(rec3.variables[variable][:])
    if variable != 't2m':
      var3 = np.array(var3)/(60*60*24.0)
    days2015 = np.array(rec3.variables['time'][:])

  # Convert date to fractional year
  startday1950 = jdcal.gcal2jd(1950,1,1)
  Nt1 = len(daysfrom1950_1)
  if variable != 'zs':
    Nt2 = len(daysfrom1950_2)
    Nt3 = len(days2015)
    time = np.zeros(Nt1+Nt2+Nt3)
    for i in range(0,Nt1):
      year,month,day,fracday = jdcal.jd2gcal(startday1950[0],startday1950[1]+daysfrom1950_1[i])
      time[i] = datelib.date_to_fracyear(year,month,day) 
    for i in range(0,Nt2):
      year,month,day,fracday = jdcal.jd2gcal(startday1950[0],startday1950[1]+daysfrom1950_2[i])
      time[i+Nt1] = datelib.date_to_fracyear(year,month,day)
    for i in range(0,Nt3):
      time[i+Nt1+Nt2] = datelib.doy_to_fracyear(2015,1+days2015[i])
  else:
    time = daysfrom1950_1 
    time = time[0:-71]
    var1 = var1[0:-71,:,:]
    
  
  # Convert lat,lon to epsg 3413
  xrac,yrac = coordlib.convert(lon,lat,4326,epsg)

  # Find x,y indices that fall within the desired grid and check to make sure that the chosen
  # indices fall on the ice mask (mask == 1) 
  if maskvalues == 'ice':
    xind = np.where((xrac >= xmin) & (xrac <= xmax) & (mask == 1))
  elif maskvalues == 'notice':
    xind = np.where((xrac >= xmin) & (xrac <= xmax) & (mask == 0))
  elif maskvalues == 'both':
    xind = np.where((xrac >= xmin) & (xrac <= xmax))
  xrac_subset = xrac[xind]
  yrac_subset = yrac[xind]
  if variable != 'zs':
    var1_subset = var1[:,:,xind[0],xind[1]]  
    var2_subset = var2[:,:,xind[0],xind[1]]
    var3_subset = var3[:,xind[0],xind[1]]
  else:
    var1_subset = var1[:,xind[0],xind[1]] 
  mask_subset = mask[xind[0],xind[1]]
  
  if maskvalues == 'ice':
    yind = np.where((yrac_subset >= ymin) & (yrac_subset <= ymax) & (mask_subset == 1))
  elif maskvalues == 'notice':
    yind = np.where((yrac_subset >= ymin) & (yrac_subset <= ymax) & (mask_subset == 0))
  elif maskvalues == 'both':
    yind = np.where((yrac_subset >= ymin) & (yrac_subset <= ymax))
  xrac_subset = xrac_subset[yind]
  yrac_subset = yrac_subset[yind]
  if variable != 'zs':  
    var1_subset = var1_subset[:,:,yind]
    var2_subset = var2_subset[:,:,yind]
    var3_subset = var3_subset[:,yind]
    var_subset = np.row_stack([var1_subset[:,0,0,:],var2_subset[:,0,0,:],var3_subset[:,0,:]])
  else:
    var1_subset = var1_subset[:,yind]
    var_subset = var1_subset[:,0,:]
  
  if variable == 't2m':
    # Convert Kelvin to Celsius
    var_subset=var_subset
  elif variable == 'smb' or variable == 'precip' or variable == 'runoff':
    # If variable is smb, convert kg m-2 s-1 to kg m-2 d-1
    var_subset=var_subset*(60*60*24.0)
  
  return xrac_subset,yrac_subset,var_subset,time

def racmo_interpolate_to_cartesiangrid(x,y,variable,epsg=3413,maskvalues='ice',timing='mean'):


  ''' 
  Pull all values for RACMO smb, t2m, zs, or runoff values for the region defined 
  by arrays x,y.
  
  xrac,yrac,var,time = racmo_grid(xmin,xmax,ymin,ymax,
  		variable,epsg=3413,mask='ice')
  
  Inputs:
  x,y : grid to interpolate RACMO values onto
  variable : what variable you want (runoff, t2m, zs, smb)
  maskvalues : if you want only 'ice' or 'notice' or 'both' values
  
  Outputs:
  var: value of chosen variable at these points
  time : time
  '''

  # RACMO data
  if variable != 'zs':
    files = [(os.path.join(os.getenv("DATA_HOME"),"Climate/RACMO/2015_09_Laura_Kehrl/RACMO2.3_GRN11_"+variable+"_daily_2001_2010.nc")), \
           (os.path.join(os.getenv("DATA_HOME"),"Climate/RACMO/2015_09_Laura_Kehrl/RACMO2.3_GRN11_"+variable+"_daily_2011_2014.nc")), \
           (os.path.join(os.getenv("DATA_HOME"),"Climate/RACMO/2015_09_Laura_Kehrl/RACMO2.3_GRN11_"+variable+"_daily_2015.nc"))]
  else:
    files = [(os.path.join(os.getenv("DATA_HOME"),"Climate/RACMO/ZS_ZGRN_V5_1960-2014_detrended_2day.nc"))]

  rec1 = netCDF4.Dataset(files[0])
  if variable != 'zs':
    rec2 = netCDF4.Dataset(files[1])
    rec3 = netCDF4.Dataset(files[2])
  mask = netCDF4.Dataset(os.path.join(os.getenv("DATA_HOME"),"Climate/RACMO/2015_09_Laura_Kehrl/RACMO23_masks_ZGRN11.nc")).variables['icemask'][:]

  # Load RACMO data
  lat = np.array(rec1.variables['lat'][:])
  lon = np.array(rec1.variables['lon'][:])
  var1 = np.array(rec1.variables[variable][:])
  daysfrom1950_1 = np.array(rec1.variables['time'][:])
  if variable != 'zs':
    var2 = np.array(rec2.variables[variable][:])
    daysfrom1950_2 = np.array(rec2.variables['time'][:])
    var3 = np.array(rec3.variables[variable][:])
    if variable != 't2m':
      var3 = np.array(var3)/(60*60*24.0)
    days2015 = np.array(rec3.variables['time'][:])

  # Convert date to fractional year
  startday1950 = jdcal.gcal2jd(1950,1,1)
  Nt1 = len(daysfrom1950_1)
  if variable != 'zs':
    Nt2 = len(daysfrom1950_2)
    Nt3 = len(days2015)
    time = np.zeros(Nt1+Nt2+Nt3)
    for i in range(0,Nt1):
      year,month,day,fracday = jdcal.jd2gcal(startday1950[0],startday1950[1]+daysfrom1950_1[i])
      time[i] = datelib.date_to_fracyear(year,month,day) 
    for i in range(0,Nt2):
      year,month,day,fracday = jdcal.jd2gcal(startday1950[0],startday1950[1]+daysfrom1950_2[i])
      time[i+Nt1] = datelib.date_to_fracyear(year,month,day)
    for i in range(0,Nt3):
      time[i+Nt1+Nt2] = datelib.doy_to_fracyear(2015,1+days2015[i])
  else:
    time = daysfrom1950_1 
    time = time[0:-71]
    var1 = var1[0:-71,:,:]
  
  # Combine variables into same file  
  var = np.row_stack([var1[:,0,:,:],var2[:,0,:,:],var3])
  
  # Convert lat,lon to epsg 3413
  xrac,yrac = coordlib.convert(lon,lat,4326,epsg)

  # Get dimensions of output grid
  nx = len(x)
  ny = len(y)
  
  if maskvalues == 'ice':
    ind = np.where((mask == 1))
  elif maskvalues == 'notice':
    ind = np.where((mask == 0))
  elif maskvalues == 'both':
    ind = np.where((mask == 1) | (mask == 0))
  
  # Make a KD-tree so we can do range searches fast
  xracflat = xrac[ind]
  yracflat = yrac[ind]
  tree = cKDTree(np.column_stack([xracflat,yracflat]))
  
  # Make a gridded data set from the model output
  if timing == 'mean':
    vargrid = np.zeros([len(y),len(x)])
    varflat = np.mean(var,axis=0)[ind]
    time = np.mean(time)
    
  # For each point in the grid,
  for i in range(ny):
    for j in range(nx):
      L = tree.query_ball_point( (x[j], y[i]), 10.e3 )
      
      # Initialize the weights to 0
      weights = 0.0
      
      # For all the nearby model points,
      for l in L:
        xp = xracflat[l]
        yp = yracflat[l]
        
        # find the distance to the current point and the
        # appropriate weight
        r = np.sqrt( (x[j] - xp)**2 + (y[i] - yp)**2 )
        w = (10.e3/(r))**3
        weights += w
        
        vargrid[i, j] += w * varflat[l]
        
      vargrid[i,j] /= weights

  if variable == 'smb' or variable == 'precip' or variable == 'runoff':
    # If variable is smb, convert kg m-2 s-1 to kg m-2 d-1
    vargrid=vargrid*(60*60*24.0)
  
  return time,vargrid

  
def racmo_at_pts(xpt,ypt,variable,filt_len='none',epsg=3413,method='nearest',maskvalues='ice'):

  '''
  Calculate RACMO timeseries at points xpts,ypts. You can either use the nearest 
  RACMO grid point or linearly interpolate from the nearest points.
  
  xrac_subset,yrac_subset,filtered,time = racmo_at_pts(xpt,ypt,variable,
  		filt_len='none',epsg=3413,method='nearest',mask='ice')
  		
  Inputs:
  xpt,ypt : points where you want the RACMO timeseries
  variable : what RACMO variable you want (t2m, zs, runoff, smb)
  filt_len : whether or not the timeseries should be filtered and if it should be filtered over how many days
  method : method for calculating the RACMO timeseries at xpt, ypt (linear or nearest)
  mask : should we only grab points for 'ice' or 'notice' or 'both'
  
  Outputs :
  xrac_subset, yrac_subset : coordinates for nearest points
  filtered : filtered (or not filtered) timeseries
  time : time
  '''

  xrac,yrac,var,time = racmo_grid(np.min(xpt)-20e3,np.max(xpt)+20e3,np.min(ypt)-20e3,np.max(ypt)+20e3,variable,epsg=3413,maskvalues=maskvalues)

  # Number of timesteps
  Nt = len(time)
 
  # Number of points
  if isinstance(xpt,(float,int)):
    var_subset = np.zeros(Nt)
  else:
    Npt = len(xpt)
    var_subset = np.zeros([Nt,Npt])
    # Set up output variables
    xrac_subset = np.zeros(Npt)
    yrac_subset = np.zeros(Npt)
    
  if method == 'nearest':
    if isinstance(xpt,(float,int)):
      ind = np.argmin((xrac-xpt)**2.0 + (yrac-ypt)**2.0)
      xrac_subset = xrac[ind]
      yrac_subset = yrac[ind]
      var_subset = var[:,ind]
    else:
      for j in range(0,Npt):
        ind = np.argmin((xrac-xpt[j])**2.0 + (yrac-ypt[j])**2.0)
        xrac_subset[j] = xrac[ind]
        yrac_subset[j] = yrac[ind]
        var_subset[:,j] = var[:,ind]
  elif method == 'linear':
    for j in range(0,Nt):
      var_subset[j,:] = scipy.interpolate.griddata((xrac,yrac),var[j,:],(xpt,ypt),method='linear')
    
  # Filter the timeseries
  if filt_len != 'none':
    print "Filtered timeseries for variable",variable
    cutoff=(1/(filt_len/365.25))/(1/(np.diff(time[1:3])*2))
    b,a=signal.butter(3,cutoff,btype='low')
    if isinstance(xpt,(float,int)):
      filtered = signal.filtfilt(b,a,var_subset)
    else:
      filtered = np.zeros([Nt,Npt])
      for j in range(0,Npt):
        filtered[:,j] = signal.filtfilt(b,a,var_subset[:,j])
  else:
    filtered = var_subset
    

  
  return xrac_subset,yrac_subset,filtered,time


def seasonlength(time,values,variable):

  # Get years
  year = np.unique(np.floor(time))
    
  # Calculate length of melt season according to different variables
  if variable == 'runoff':
    threshold = 1.0e-3
  elif variable == 'smb':
    threshold = 0.0
    values = -values
  elif variable == 't2m':
    threshold = 0.0
  
  # Set up output variables
  firstday = np.zeros(len(year)) # first day of melt season in fractional years
  lastday = np.zeros(len(year)) # last day of melt season in fractional years
  meltlength = np.zeros(len(year)) # length of melt season in days
  total = np.zeros(len(year))
    
  for i in range(0,len(year)):
    # Get indices for that year
    ind = np.where(np.floor(time) == year[i])[0]
    day1 = 0
    day2 = 0
    yearlen = datelib.date_to_doy(year[i],12,31)[1]
    # Get runoff values > 0 (which is when melt is happening)
    for j in range(0,len(ind)):
      if (values[ind[j]]) > threshold:
        if (day1 == 0):
          day1 = time[ind[j]-1]
          day1ind = ind[j]
          lastind = j
          day2 = time[ind[j]]
          day2ind = ind[j]
        elif (day1 != 0) and (j-lastind)==1.0:
          day2 = time[ind[j]]
          day2ind = ind[j]
          lastind = j
        elif j-lastind > 1.0:
          if (day2-day1)*yearlen > meltlength[i]:
            meltlength[i] = (day2-day1)*yearlen
            firstday[i] = day1
            lastday[i] = day2
            day1ind_final = day1ind
            day2ind_final = day2ind
            day1 = 0
            day2 = 0
          else:
            day1 = 0
            day2 = 0
      elif j == len(ind)-1:
        if (day2-day1)*yearlen > meltlength[i]:
          meltlength[i] = (day2-day1)*yearlen
          firstday[i] = day1
          lastday[i] = day2  
          day1ind_final = day1ind
          day2ind_final = day2ind
    total[i] = np.sum(values[day1ind_final:day2ind_final])
                        
  return year,firstday,lastday,meltlength,total
  
def calculate_PDD(time,T):
  
  PPD = np.zeros(len(T))
  
  for i in range(0,len(time)):
    if i == 0:
      year = np.floor(time[i])
    if np.floor(time[i])==year:
      if T[i] > 0:
        PPD[i] = PPD[i-1]+T[i]
      elif i > 0:
        PPD[i] = PPD[i-1]
    else:
      year = np.floor(time[i]) 
      if T[i] > 0:
        PPD[i] = T[i]
    

  return time,PDD
  
def cumsmb(xpt,ypt,epsg=3413,rho_i=900.0,method='nearest',maskvalues='ice'):

  '''
  Cumulative surface mass balance.
  '''

  xrac,yrac,smb,time = racmo_at_pts(xpt,ypt,'smb',\
  		filt_len='none',epsg=epsg,method=method,maskvalues=maskvalues)


  if isinstance(xpt,(float,int)):
    zs_smb = np.zeros(len(time))
    for i in range(1,len(time)):
      zs_smb[i] = zs_smb[i-1]+smb[i]
  else:
    zs_smb = np.zeros([len(time),len(xpt)])
    for i in range(1,len(time)):
      zs_smb[i,:] = zs_smb[i-1,:]+smb[i,:]
  
  zs_smb = zs_smb/rho_i
  
  return xrac,yrac,zs_smb,time

def SIF_at_pts(xpts,ypts,epsg=3413,filt_len='none',variable='sif'):

  #file = os.path.join(os.getenv("DATA_HOME"),"Climate/SeaIce/METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2_1447191850310.nc")
  #file = os.path.join(os.getenv("DATA_HOME"),"Climate/SeaIce/METOFFICE-GLO-SST-REANALYSIS_1.nc")
  #file = os.path.join(os.getenv("DATA_HOME"),"Climate/SeaIce/METOFFICE-GLO-SST-REANALYSIS_2.nc")
  file = os.path.join(os.getenv("DATA_HOME"),"Climate/SeaIce/METOFFICE-GLO-SST-V2.nc")


  rec = netCDF4.Dataset(file)
  
  lat = rec.variables['lat'][:]
  lon = rec.variables['lon'][:]
  time = rec.variables['time'][:]
  mask = rec.variables['mask'][:]
  if variable == 'sif':
    sif = rec.variables['sea_ice_fraction'][:]
  elif variable == 'sst':
    sif = rec.variables['analysed_sst'][:]
  
  lon_grid,lat_grid = np.meshgrid(lon,lat)
  
  xsif,ysif = coordlib.convert(lon_grid,lat_grid,4326,epsg)
  
  try:
    bestind = np.zeros(len(xpts))
    print "need to fix code for multiple points"
  except:
    mindist = 500.e3
    for j in range(0,len(lat)):
      for i in range(0,len(lon)):
        if (mask[100,j,i] != 2) and (xsif[j,i] > xpts):
          d = np.sqrt((xsif[j,i]-xpts)**2+(ysif[j,i]-ypts)**2)
          if d < mindist:
            mindist = d
            bestind = [j,i]
    
  # Convert time to fractional years
  # Why they have time listed as seconds since Jan. 1, 1981 is beyond me. 
  fractime = 1981. + time/(365.25*24*60*60)

  # Filter the timeseries
  if filt_len != 'none':
    print "Filtered timeseries for "+variable
    cutoff=(1/(filt_len/365.25))/(1/(np.diff(fractime[1:3])*2))
    b,a=signal.butter(3,cutoff,btype='low')
    filtered = signal.filtfilt(b,a,sif[:,bestind[0],bestind[1]])
  else:
    filtered = sif[:,bestind[0],bestind[1]]

  return xsif[bestind[0],bestind[1]],ysif[bestind[0],bestind[1]],filtered,fractime
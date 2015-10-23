# This is a module to help load the RACMO timeseries. 

import os
import sys
import numpy as np
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"Util/Modules"))
import netCDF4, jdcal
import coords, fracyear
import scipy.signal as signal

def racmo_grid(xmin,xmax,ymin,ymax,variable,epsg=3413):

  # RACMO data
  files = [(os.path.join(os.getenv("DATA_HOME"),"Climate/RACMO/2015_09_Laura_Kehrl/RACMO2.3_GRN11_"+variable+"_daily_2001_2010.nc")), \
           (os.path.join(os.getenv("DATA_HOME"),"Climate/RACMO/2015_09_Laura_Kehrl/RACMO2.3_GRN11_"+variable+"_daily_2011_2014.nc"))]

  rec1 = netCDF4.Dataset(files[0])
  rec2 = netCDF4.Dataset(files[1])
  mask = netCDF4.Dataset(os.path.join(os.getenv("DATA_HOME"),"Climate/RACMO/2015_09_Laura_Kehrl/RACMO23_masks_ZGRN11.nc")).variables['icemask'][:]


  # Load RACMO data
  lat = np.array(rec1.variables['lat'][:])
  lon = np.array(rec1.variables['lon'][:])
  var1 = np.array(rec1.variables[variable][:])
  var2 = np.array(rec2.variables[variable][:])
  daysfrom1950_1 = np.array(rec1.variables['time'][:])
  daysfrom1950_2 = np.array(rec2.variables['time'][:])

  # Convert date to fractional year
  startday1950 = jdcal.gcal2jd(1950,1,1)
  Nt1 = len(daysfrom1950_1)
  Nt2 = len(daysfrom1950_2)
  time = np.zeros(Nt1+Nt2)
  for i in range(0,Nt1):
    year,month,day,fracday = jdcal.jd2gcal(startday1950[0],startday1950[1]+daysfrom1950_1[i])
    time[i] = fracyear.date_to_fracyear(year,month,day)
  for i in range(0,Nt2):
    year,month,day,fracday = jdcal.jd2gcal(startday1950[0],startday1950[1]+daysfrom1950_2[i])
    time[i+Nt1] = fracyear.date_to_fracyear(year,month,day)
  
  # Convert lat,lon to epsg 3413
  xrac,yrac = coords.convert(lon,lat,4326,epsg)

  # Find x,y indices that fall within the desired grid and check to make sure that the chosen
  # indices fall on the ice mask (mask == 1) 
  xind = np.where((xrac >= xmin) & (xrac <= xmax) & (mask == 1))
  xrac_subset = xrac[xind]
  yrac_subset = yrac[xind]
  var1_subset = var1[:,:,xind[0],xind[1]]
  var2_subset = var2[:,:,xind[0],xind[1]]
  mask_subset = mask[xind[0],xind[1]]

  yind = np.where((yrac_subset >= ymin) & (yrac_subset <= ymax) & (mask_subset == 1))
  xrac_subset = xrac_subset[yind]
  yrac_subset = yrac_subset[yind]
  var1_subset = var1_subset[:,:,yind]
  var2_subset = var2_subset[:,:,yind]
  
  var_subset = np.row_stack([var1_subset[:,0,0,:],var2_subset[:,0,0,:]])
  
  if variable == 't2m':
    # Convert Kelvin to Celsius
    var_subset=var_subset-273.15
  elif variable == 'smb' or variable == 'precip' or variable == 'runoff':
    # If variable is smb, convert kg m-2 s-1 to kg m-2 d-1
    var_subset=var_subset*(60*60*24.0)

  
  return xrac_subset,yrac_subset,var_subset,time

  
def racmo_at_pts(xpt,ypt,variable,filt_len='none',epsg=3413,method='nearest'):

  xrac,yrac,var,time = racmo_grid(np.min(xpt)-20e3,np.max(xpt)+20e3,np.min(ypt)-20e3,np.max(ypt)+20e3,variable,epsg=3413)
 
  # Number of points
  Npt = xpt.shape
  Nt = len(time)
  
  # Set up output variables
  xrac_subset = np.zeros(Npt)
  yrac_subset = np.zeros(Npt)
  
  if not Npt:
    var_subset = np.zeros([Nt,0])
  else:
    var_subset = np.zeros([Nt,Npt])

  
  if method == 'nearest':
    if not Npt:
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
    if not Npt:
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
    yearlen = fracyear.date_to_doy(year[i],12,31)[1]
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
  
#def cumsmb(xpt,ypt,epsg=3413,method='nearest'):
#
#  '''
#  Cumulative surface mass balance.
#  '''
#
#  xrac_subset,yrac_subset,filtered,time = racmo_at_pts(xpt,ypt,variable,\
#  		filt_len='none',epsg=epsg,method=method)
#
#
#return
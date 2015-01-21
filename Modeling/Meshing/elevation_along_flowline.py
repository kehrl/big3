# This file takes the Helheim IceBridge data and plots it along a flowline.

import os
import sys
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
import numpy as np
import matplotlib.pyplot as plt
import icebridge
import geotiff
from scipy import interpolate

# Moving average
def movingaverage(interval, window_size):
  window = np.ones(int(window_size))/float(window_size)
  return np.convolve(interval, window, 'same')

##########
# Inputs #
##########

Glacier = "Helheim"

# ATM data
DIRATM=os.path.join(os.getenv("HOME"),"Data/Elevation/IceBridge/"+Glacier)
DIRATMs = [os.path.join(DIRATM,f) for f in os.listdir(DIRATM) if os.path.isdir(os.path.join(DIRATM,f)) ]

# Worldview Data
DIRWV=os.path.join(os.getenv("HOME"),"Data/DEMs/Worldview/"+Glacier+"/")

# Flowline
flowline=np.genfromtxt(os.path.join(os.getenv("HOME"),"Models/Helheim/Meshes/Flowline/Worldview_Advance/Inputs/flowline.dat"),skiprows=1)

################
# Get ATM data #
################

# Load all ATM data
files=[]
atm={}
for DIR in DIRATMs:
  allfiles=os.listdir(DIR)
  for file in allfiles:
    if file.endswith('nadir3seg'):
      files.append(os.path.join(DIR,file))

  atm[DIR[-8:-1]] = icebridge.read(files)


# Get ATM elevations along flowline
atm_flowline={}
xmin=np.min(flowline[:,1])
xmax=np.max(flowline[:,1])
ymin=np.min(flowline[:,2])
ymax=np.max(flowline[:,2])
for key in atm.keys():
  # Setup output variable
  atm_flowline[key] = np.zeros_like(flowline[:,0])
  atm_flowline[key][:] = 'NaN'
  
  # Find ATM elevation measurements within 500 m of flowline
  data=atm[key][(atm[key][:,0] >= xmin) & (atm[key][:,0] <= xmax) & (atm[key][:,1] >= ymin) &(atm[key][:,1] <= ymax)]
  for i in range(0,len(data)):
    mindist=(np.sqrt((data[i,0]-flowline[:,1])**2+(data[i,1]-flowline[:,2])**2)).min()
    if mindist < 500:
      minind=(np.sqrt((data[i,0]-flowline[:,1])**2+(data[i,1]-flowline[:,2])**2)).argmin()
      atm_flowline[key][minind]=data[i,2]
  
  atm_flowline[key]=movingaverage(atm_flowline[key],5)
  # Filter results
  #nonnan=np.where(~np.isnan(atm_flowline[key]))
  #satisfied=0
  #while not(satisfied):
  #  
  #  satisfied=1
  #surf_filt_len=float(1000)
  #cutoff=(1/surf_filt_len)/(1/(np.diff(flowline[1:3,0])*2))
  #b,a=signal.butter(4,cutoff,btype='low')
  #atm_flowline[key]=

######################
# Get Worldview Data #
###################### 

files=[f for f in os.listdir(DIRWV) if f.endswith('.tif') ]

wv_flowline={}
for file in files:
  [xs,ys,zs] = geotiff.read(DIRWV+file)
  zs_dem = interpolate.RectBivariateSpline(ys,xs,zs)
  zs_interp = zs_dem.ev(flowline[:,2],flowline[:,1])
  zs_interp[zs_interp < 40] = 'NaN'
  surf_filt_len=float(1000)
  cutoff=(1/surf_filt_len)/(1/(np.diff(flowline[1:3,0])*2))
  #b,a=signal.butter(4,cutoff,btype='low')
  wv_flowline[file[0:8]]=movingaverage(zs_interp,5)
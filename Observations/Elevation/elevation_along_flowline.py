# This file takes the Helheim IceBridge data and plots it along a flowline.

import os
import sys
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules/"))
sys.path.append(os.path.join(os.getenv("HOME"),"Code/BigThreeGlaciers/Tools/"))
import glacier
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

Glacier = 'Helheim'

x,y,dists = 

################
# Get ATM data #
################

 = helheim_elevation.atm(years)


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
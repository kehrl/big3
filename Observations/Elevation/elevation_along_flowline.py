# This file takes the Helheim IceBridge data and plots it along a flowline.

import os
import sys
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules/"))
sys.path.append(os.path.join(os.getenv("HOME"),"Code/BigThreeGlaciers/Tools/"))
import glacier_flowline, elevation
import numpy as np
import matplotlib.pyplot as plt


##########
# Inputs #
##########

glacier = 'Kanger'

x,y,dists = glacier_flowline.load(glacier)

dists_eul = [0.0,5.0,10.0,20.0,30.0] # kilometers

time1 = 2008 #start time for plot
time2 = 2015.5 # end time for plot

################
# Get terminus #
################

terminus_val, terminus_time = icefronts.distance_along_flowline(x,y,dists,glacier,type='icefront')

# Chop to desired time interval
indt = np.where((terminus_time > time1) & (terminus_time < time2))
terminus_time = terminus_time[indt[0]]
terminus_val = terminus_val[indt[0]]
del indt

# Reference for terminus position
terminus = np.mean(terminus_val)

ind=[]
for i in range(0,len(dists_eul)):
  ind.append( (abs(dists - (terminus - dists_eul[i]*1e3))).argmin() )

################
# Get ATM data #
################

atm_data = elevation.atm_along_flowline(x,y,glacier,'all',cutoff='terminus',maxdist=500,verticaldatum='geoid')
zpt_atm,time_atm = elevation.atm_at_pts(x[ind],y[ind],glacier,'all',

######################
# Get Worldview Data #
###################### 

wv_data = elevation.worldview_along_flowline(x,y,glacier,years='all',cutoff='terminus',verticaldatum='geoid')
zpt_wv,time_wv = elevation.worldview_at_pts(x[ind],y[ind],glacier,'all','geoid')

#########################################
# Make plot of elevation along flowline #
#########################################

dates = np.sort(atm_data.keys())
for date in dates:
  plt.plot(dists,atm_data[date][:,2],label=date,linewidth=1.5)
plt.legend()

dates = np.sort(wv_data.keys())
for date in dates:
  plt.plot(dists,wv_data[date][:,2],label=date,linewidth=1.5)
plt.legend()
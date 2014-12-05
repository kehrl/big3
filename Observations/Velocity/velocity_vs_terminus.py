# This file compares Helheim's terminus position to its velocity.

import os
import sys
import numpy as np
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
import helheim_velocity, helheim_icefronts
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import AutoMinorLocator
import shapefile
import jdcal
from shapely.geometry import LineString

##########
# Inputs #
##########

# Flowline
MESHNAME = "Worldview_Advance"
DIRM = os.path.join(os.getenv("HOME"),"Models/Helheim/Meshes/Flowline/"+MESHNAME+"/")
flowline = np.loadtxt(DIRM+"Inputs/flowline.dat",skiprows=1)

# Locations for velocities
dists = [0,5,10,20,28.9] # kilometers
terminus = 82 # location of terminus in flowline (km)

##################################
# Get velocities at those points # 
##################################
ind=[]
for i in range(0,len(dists)):
  ind.append( (abs(flowline[:,0] - 1e3*(terminus - dists[i]))).argmin() )
  
velocity_val,velocity_time = helheim_velocity.velocity_at_points(flowline[ind,1],flowline[ind,2])
velocity_val = velocity_val[0:-2,:] # Ignore bad velocities
velocity_time = velocity_time[0:-2,:] # Ignore bad velocities

##################
# Get ice fronts #
##################
terminus_val, terminus_time = helheim_icefronts.distance_along_flowline(flowline[:,1],flowline[:,2],flowline[:,0])

##############
# Make plots #
##############

# Plot terminus, velocity through time
plt.figure(figsize=(32,12))
gs = matplotlib.gridspec.GridSpec(3,1)

# Get mean velocities
plt.subplot(gs[:2, :]) 
#plt.plot([2000,2014],[0,0],'k')
colors=['b','c','g','y','r']
for i in (range(0,len(velocity_val[0,:]))):
  meanvel = np.mean(velocity_val[~(np.isnan(velocity_val[:,i])),i])
  plt.plot(velocity_time[~(np.isnan(velocity_val[:,i]))],(velocity_val[~(np.isnan(velocity_val[:,i])),i]-velocity_val[0,i])/1e3,'o-',color=colors[i],label=str(dists[i])+' km',linewidth=2,markersize=10)
ax = plt.gca()
plt.xticks(range(2000,2014),fontsize=2)
x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
ax.xaxis.set_major_formatter(x_formatter)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.set_xticklabels([])
plt.xlim([2012,2013])
plt.ylim([-0.1,4.2])
plt.yticks([0,1,2,3,4],fontsize=30) 
#plt.legend(loc=2,fontsize=30)
plt.ylabel('Relative glacier velocity \n (km/yr)',fontsize=28)
ax.tick_params('both', length=20, width=2, which='major')
ax.tick_params('both', length=10, width=1, which='minor')

plt.subplot(gs[2, :])
plt.plot(terminus_time,(terminus_val-82000)/1e3,'ko-',linewidth=2,markersize=10)
#plt.plot([2008,2014],[(terminus_val[0]-82500)/1e3,(terminus_val[0]-82500)/1e3])
ax = plt.gca()
ax.xaxis.set_major_formatter(x_formatter)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
plt.xticks(range(2000,2014),fontsize=30)
plt.xlim([2012,2013])  
plt.ylim([-1.5,7])
plt.yticks([0,2,4,6],fontsize=28)
plt.ylabel('Terminus position \n (km)',fontsize=30)
ax.tick_params('both', length=20, width=2, which='major')
ax.tick_params('both', length=10, width=1, which='minor')
plt.tight_layout()
plt.subplots_adjust(hspace=0,wspace=0)

# Plot velocity vs. terminus position
plt.figure()
terminus_interp = np.interp(velocity_time,terminus_time,terminus_val)
for i in [0]:
  meanvel = np.mean(velocity_val[~(np.isnan(velocity_val[:,i])),i])
  plt.plot(terminus_interp-np.mean(terminus_interp),(velocity_val[:,i]-meanvel)/1e3,'.',color=colors[i],label=str(dists[i])+' km')
plt.xlabel('Terminus position (km)')
plt.ylabel('Glacier velocity (km/yr)')

# Plot terminus position along flowline
plt.figure()
floatheight = -(flowline[:,3])*(1.024/0.917-1)
ind = np.where((terminus_time > 2008.5))
terminus_interp = terminus_val[ind]
velocity_interp = np.interp(terminus_time[ind],velocity_time[:,0],velocity_val[:,0])
surf_interp = np.interp(terminus_interp,flowline[:,0],floatheight)
bed_interp = np.interp(terminus_interp,flowline[:,0],flowline[:,3])
plt.plot(flowline[0:-15,0],flowline[0:-15,4],'k')
plt.plot([flowline[-16,0],flowline[-16,0]],[flowline[-16,4],60],'k')
plt.plot(flowline[:,0],floatheight,'k--')
plt.plot(flowline[:,0],flowline[:,3],'k')
for i in range(0,len(surf_interp)):
  velcolor = 'b'
  #if velocity_interp[i] < 8000:
  #  velcolor='b'
  #elif velocity_interp[i] < 9000:
  #  velcolor='m'
  #elif velocity_interp[i] > 9000:
  #  velcolor='r'
  plt.plot([terminus_interp[i],terminus_interp[i]],[surf_interp[i],bed_interp[i]],color=velcolor,linewidth=1.5)
plt.xlim([65000,87000])
plt.ylim([-1000,800])
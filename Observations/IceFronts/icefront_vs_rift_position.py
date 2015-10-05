# This file compares the location of rift and terminus positions

import os
import sys
import numpy as np
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Helheim/Observations/IceFronts"))
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Helheim/Observations/Bed"))
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Helheim/Observations/Velocity"))
import helheim_velocity, helheim_icefronts, helheim_bed
import matplotlib.pyplot as plt
import geotiff
import math

###########
# Options #
###########

MESHNAME = "High"
DIRM = os.path.join(os.getenv("MODEL_HOME"),"Helheim/Meshes/Flowline/"+MESHNAME+"/")
flowline = np.loadtxt(DIRM+"Inputs/flowline.dat",skiprows=1)

################
# Plot options #
################

cresisbeds = 1 # plot cresis bed profiles
icefronts = 0 # plot ice fronts
rifts = 0 # plot rifts
xlim1=78 #km along flowline
xlim2=90 #km along flowline

image = geotiff.readrgb(os.path.join(os.getenv("DATA_HOME"),"Imagery/Landsat/Helheim/TIF/20140611135917_LC82320132014162LGN00.tif"))

# Bed profiles
years=['2013','2012','2009b','2008','2001']
colors=['c','b','r','m','g']

############
# Get data #
############

bed2001 = helheim_bed.cresis('2001')
bed2001dist = np.interp(bed2001[:,0],flowline[:,1],flowline[:,0])

terminus_val, terminus_time = helheim_icefronts.distance_along_flowline(flowline[:,1],flowline[:,2],flowline[:,0],'icefront')
rift_val, rift_time = helheim_icefronts.distance_along_flowline(flowline[:,1],flowline[:,2],flowline[:,0],'rift')

#terminus_val, terminus_time = helheim_icefronts.distance_along_flowline(flowline2001bed[:,0],flowline2001bed[:,1],flowline2001bed_dist,'icefront')
#rift_val, rift_time = helheim_icefronts.distance_along_flowline(flowline2001bed[:,0],flowline2001bed[:,1],flowline2001bed_dist,'rift')

terminus_bed = np.interp(terminus_val,bed2001dist,bed2001[:,2])
rift_bed = np.interp(rift_val,bed2001dist,bed2001[:,2])
rift_surf = np.interp(rift_val,flowline[:,0],flowline[:,4])
terminus_surf = np.interp(terminus_val,flowline[:,0],flowline[:,4])

###########################
# Make cross section plot #
###########################

plt.figure()
if icefronts == 1:
  for i in range(0,len(terminus_time)):
    plt.plot([terminus_val[i]/1e3,terminus_val[i]/1e3],[terminus_bed[i],terminus_surf[i]],'0.75',linewidth=1.5)
if rifts == 1:
  for i in range(0,len(rift_time)):
    plt.plot([rift_val[i]/1e3,rift_val[i]/1e3],[rift_bed[i],rift_surf[i]],'r',linewidth=1.5)
plt.plot(flowline[:,0]/1e3,flowline[:,4],'k',linewidth=1.5)
plt.plot(flowline[:,0]/1e3,flowline[:,3],'k',linewidth=1.5,label='Mesh bed')

if cresisbeds == 1: 
  for i in range(0,len(years)):
    bed = helheim_bed.cresis(years[i])
    nonzero = np.where(bed[:,2] < -20)
    dist = np.interp(bed[:,0],flowline[:,1],flowline[:,0])
    zb = bed[nonzero,2]
    plt.plot(dist[nonzero]/1e3,bed[nonzero,2].T,color=colors[i],label=years[i],linewidth=1.5)
  plt.legend()
else:
  plt.plot(bed2001dist/1e3,bed2001[:,2],'g',label='2001')
  plt.legend()
  
plt.xlim([xlim1,xlim2])
plt.ylim([-800,400])
plt.ylabel('Elevation (masl)')
plt.xlabel('Distance along flowline (km)')
plt.savefig('/Users/kehrl/Bigtmp/Helheim_Jan26/fig_bed_profile.pdf',FORMAT='PDF')
plt.close()

######################
# Make overview plot #
######################
plt.figure()
image_extent = [np.min(image[0]),np.max(image[0]),np.max(image[1][np.nonzero(image[1])]),np.min(image[1])]
plt.imshow(image[2],extent=image_extent)
plt.gca().invert_yaxis()
plt.axis('equal')
plt.xlim([300000,314000])
plt.ylim([-2582000,-2572000])
ind = np.where((flowline[:,0] > xlim1) & (flowline[:,0] < xlim2))
plt.plot(flowline[ind[0],1],flowline[ind[0],2],'k',linewidth=3)
if cresisbeds == 1:
  for i in range(0,len(years)):
    bed = helheim_bed.cresis(years[i])
    dist = np.interp(bed[:,0],flowline[:,1],flowline[:,0])/1e3
    ind = np.where((dist > xlim1) & (dist < xlim2))
    plt.plot(bed[ind[0],0],bed[ind[0],1],color=colors[i],label=years[i],linewidth=3)
  plt.legend()

plt.plot(flowline[:,1],flowline[:,2],'k',linewidth=1.5)
plt.title('CreSIS bed profiles')
plt.savefig('/Users/kehrl/Bigtmp/Helheim_Jan26/fig_bed_map.pdf',FORMAT='PDF')
plt.close()
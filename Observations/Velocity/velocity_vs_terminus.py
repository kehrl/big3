# This file compares Helheim's terminus position to its velocity.

import os
import sys
import numpy as np
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Helheim/Observations/IceFronts"))
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Helheim/Observations/Velocity"))
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Helheim/Observations/Bed"))
import helheim_velocity, helheim_icefronts, helheim_bed
import matplotlib.pyplot as plt
import matplotlib, geotiff
from matplotlib.ticker import AutoMinorLocator
import shapefile
import jdcal
import pylab
from shapely.geometry import LineString
import matplotlib.cm as cmx
import matplotlib.colors as colors

##########
# Inputs #
##########

# Flowline
MESHNAME = "High"
DIRM = os.path.join(os.getenv("HOME"),"Models/Helheim/Meshes/Flowline/"+MESHNAME+"/")
flowline = np.loadtxt(DIRM+"Inputs/flowline.dat",skiprows=1)

# Locations for velocities
dists = [0,5,10,20,28.9] # kilometers
terminus = 82 # location of terminus in flowline (km)

# Image for plotting
image = geotiff.readrgb(os.path.join(os.getenv("HOME"),"Data/Imagery/Landsat/Helheim/TIF/20140611135917_LC82320132014162LGN00.tif"))

################
# Plot options #
################

time1 = 2009 #start time for plot
time2 = 2015 # end time for plot
seasonal = 1 # plot seasonal bars, to highlight seasonal trends
normalizedvelocity = 0 # divide by average velocity
terminusbed = 1

##################################
# Get velocities at those points # 
##################################

ind=[]
for i in range(0,len(dists)):
  ind.append( (abs(flowline[:,0] - 1e3*(terminus - dists[i]))).argmin() )

velocity_val,velocity_time, velocity_error = helheim_velocity.velocity_at_points(flowline[ind,1],flowline[ind,2])
velocity_val[81,0] = 'NaN' # Ignore bad velocities
velocity_val[80,1] = 'NaN' # Ignore bad velocities
velocity_val[80:82,2] = 'NaN' # Ignore bad velocities
velocitypoints = np.column_stack([flowline[ind,1],flowline[ind,2]])

indt = np.where((velocity_time > time1) & (velocity_time < time2))
velocity_time = velocity_time[indt[0]]
velocity_val = velocity_val[indt[0],:]
del indt

##################
# Get ice fronts #
##################

terminus_val, terminus_time = helheim_icefronts.distance_along_flowline(flowline[:,1],flowline[:,2],flowline[:,0],'icefront')
rift_val, rift_time = helheim_icefronts.distance_along_flowline(flowline[:,1],flowline[:,2],flowline[:,0],'rift')
indt = np.where((terminus_time > time1) & (terminus_time < time2))
terminus_time = terminus_time[indt[0]]
terminus_val = terminus_val[indt[0]]
indt = np.where((rift_time > time1) & (rift_time < time2))
rift_val = rift_val[indt[0]]
rift_time = rift_time[indt[0]]
del indt

###########
# Get bed #
###########

bed2001 = helheim_bed.cresis('2001')
bed2001dist = np.interp(bed2001[:,0],flowline[:,1],flowline[:,2])

############################################
# Plot: velocity vs. terminus through time #
############################################

# Plot terminus, velocity through time
plt.figure(figsize=(12,7))
if terminusbed == 1:
  gs = matplotlib.gridspec.GridSpec(4,1)
else:
  gs = matplotlib.gridspec.GridSpec(3,1)

# Get mean velocities
plt.subplot(gs[:2, :]) 
#plt.plot([2000,2014],[0,0],'k')
coloptions=['b','c','g','y','r']
if normalizedvelocity:
  for i in (range(0,len(velocity_val[0,:]))):
    meanvel = np.mean(velocity_val[~(np.isnan(velocity_val[:,i])),i])
    plt.plot(velocity_time[~(np.isnan(velocity_val[:,i]))],(velocity_val[~(np.isnan(velocity_val[:,i])),i])/meanvel,'o-',color=coloptions[i],label=str(-dists[i])+' km',linewidth=2,markersize=5)
  #plt.ylim([0.7,1.4])
  #plt.yticks([0.8,1.0,1.2,1.4])
  plt.ylabel('Change from mean velocity',fontsize=12)
else:
  for i in (range(0,len(velocity_val[0,:]))):
    meanvel = np.mean(velocity_val[~(np.isnan(velocity_val[:,i])),i])
    plt.plot(velocity_time[~(np.isnan(velocity_val[:,i]))],(velocity_val[~(np.isnan(velocity_val[:,i])),i]-velocity_val[0,i])/1e3,'o-',color=coloptions[i],label=str(-dists[i])+' km',linewidth=2,markersize=5)
  #plt.ylim([-0.1,4.2])
  #plt.yticks([0,1,2,3,4],fontsize=10) 
  plt.ylabel('Relative velocity \n (km/yr)',fontsize=12)
ax = plt.gca()
plt.xticks(range(2000,2016),fontsize=12)
x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
ax.xaxis.set_major_formatter(x_formatter)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.set_xticklabels([])
plt.xlim([time1,time2])
plt.legend(loc=2,fontsize=10)
ax.tick_params('both', length=20, width=2, which='major')
ax.tick_params('both', length=10, width=1, which='minor')
if seasonal:
  xTickPos = np.linspace(time1,time2,(time2-time1)*2+1)
  xTickPos = xTickPos[:-1]
  ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.9','w'],linewidth=0)

plt.subplot(gs[2, :])
plt.plot(terminus_time,(terminus_val-82000)/1e3,'ko-',linewidth=2,markersize=5)
plt.plot(rift_time,(rift_val-82000)/1e3,'ro',linewidth=2,markersize=5)
#plt.plot([2008,2014],[(terminus_val[0]-82500)/1e3,(terminus_val[0]-82500)/1e3])
ax = plt.gca()
ax.xaxis.set_major_formatter(x_formatter)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
plt.xticks(range(2000,2016),fontsize=12)
plt.xlim([time1,time2]) 
plt.yticks(range(0,4), fontsize=12)
#plt.ylim([-1.5,7])
#plt.yticks([0,2,4,6],fontsize=10)
plt.ylabel('Terminus position \n (km)',fontsize=12)
#ax.arrow(2001,2.5,0,1.5,head_width=0.05, head_length=0.1, fc='k', ec='k',linewidth=3)
#ax.arrow(2001,1.5,0,-1.5,head_width=0.05, head_length=0.1, fc='k', ec='k',linewidth=3)
#ax.text(2001.1,0.5,'Retreat',fontsize=36,fontweight='bold')
#ax.text(2001.1,3.0,'Advance',fontsize=36,fontweight='bold')
ax.tick_params('both', length=20, width=2, which='major')
ax.tick_params('both', length=10, width=1, which='minor')
if seasonal == 1:
  xTickPos = np.linspace(time1,time2,(time2-time1)*2+1)
  xTickPos = xTickPos[:-1]
  ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]),bottom=min(plt.ylim()), color=['0.9','w'],linewidth=0)
if terminusbed == 1:
  ax.set_xticklabels([])

if terminusbed == 1:
  terminus_bed = np.interp(terminus_val,bed2001dist,bed2001[:,2])
  plt.subplot(gs[3, :])
  ax = plt.gca()
  plt.xticks(range(2000,2016),fontsize=12)
  ax.set_xticklabels([])
  ax.xaxis.set_major_formatter(x_formatter)
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  ax.yaxis.set_minor_locator(AutoMinorLocator(2))
  plt.xlim([time1,time2]) 
  ax.tick_params('both', length=20, width=2, which='major')
  ax.tick_params('both', length=10, width=1, which='minor')
  plt.ylabel('Bed at terminus \n (m)',fontsize=12)
  plt.plot(terminus_time,terminus_bed,'ko-',linewidth=2,markersize=5)
  plt.yticks([-700,-600,-500])
  if seasonal:
    xTickPos = np.linspace(time1,time2,(time2-time1)*2+1)
    xTickPos = xTickPos[:-1]
    ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.9','w'],linewidth=0)
 
plt.tight_layout()
plt.subplots_adjust(hspace=0.1,wspace=0)
plt.savefig('fig_velocity_time_'+str(time1)+'to'+str(time2)+'.pdf',FORMAT='PDF')
plt.close()

#########################################
# Plot overview map for previous figure #
#########################################
plt.figure()
image_extent = [np.min(image[0]),np.max(image[0]),np.max(image[1][np.nonzero(image[1])]),np.min(image[1])]
plt.imshow(image[2],extent=image_extent)
plt.gca().invert_yaxis()
plt.axis('equal')
plt.xlim([np.min(velocitypoints[:,0])-2e3,317000])
plt.ylim([-2583000,np.max(velocitypoints[:,1])+2e3])
termx,termy,termt = helheim_icefronts.load_all(time1,time2,'icefront')
riftx,rifty,riftt = helheim_icefronts.load_all(time1,time2,'rift')
plt.plot(termx,termy,color='0.75')
plt.plot(riftx,rifty,color='r')
plt.plot(flowline[:,1],flowline[:,2],'k',linewidth=1.5)
plt.title('Terminus positions from '+str(time1)+' to '+str(time2))
for i in range(0,len(coloptions)):
  plt.plot(velocitypoints[i,0],velocitypoints[i,1],'o',color=coloptions[i],markersize=10)
plt.savefig('fig_velocity_map_'+str(time1)+'to'+str(time2)+'.pdf',FORMAT='PDF')
plt.close()

# Plot velocity vs. terminus position
#plt.figure(figsize=(5,5))
#gs = matplotlib.gridspec.GridSpec(3,1)

#plt.subplot(gs[:2, :]) 
#ind = np.where(terminus_time > 2008)
#ind2 = np.where(velocity_time > 2008)
#terminus_interp = np.interp(velocity_time[ind2],terminus_time[ind],terminus_val[ind])
#plt.plot((terminus_interp)/1e3-82,(velocity_val[ind2[0],0])/1e3,'bo',label='Measured',markeredgecolor='b')
#plt.plot(term_rkm/1e3-82,velocities_rkm[:,0]/1e3,'ro',label='Simulated',markeredgecolor='r')
#ax = plt.gca()
#plt.xlim([0,4])
#plt.xticks([0,1,2,3,4],fontsize=24)
#plt.yticks([7,8,9,10],fontsize=24)
#plt.ylim([6.5,10.5])
#ax.xaxis.set_minor_locator(AutoMinorLocator(2))
#ax.yaxis.set_minor_locator(AutoMinorLocator(2))
#ax.tick_params('both', length=10, width=1, which='major')
#ax.tick_params('both', length=5, width=1, which='minor')
#plt.ylabel('Velocity (km/yr)',fontsize=28)
#plt.legend(fontsize=20,numpoints=1,handlelength=0)
#plt.xlabel('Terminus (km)',fontsize=28)
#plt.tight_layout()

#plt.subplot(gs[2, :]) 
#ax = plt.gca()
#plt.plot((flowline[:,0]-82000)/1e3,flowline[:,3],'k',linewidth=2)
#plt.plot((flowline[:,0]-82000)/1e3,flowline[:,4],'k',linewidth=2)
#plt.xlim([0,4])
#plt.yticks([-700,-600,-500],fontsize=24)
#plt.ylim([-700,-500])
#plt.ylabel('Bed (m asl)',fontsize=30)
#plt.xlabel('Terminus (km)',fontsize=30)
#ax.xaxis.set_minor_locator(AutoMinorLocator(2))
#ax.yaxis.set_minor_locator(AutoMinorLocator(2))
#ax.tick_params('both', length=10, width=1, which='major')
#ax.tick_params('both', length=5, width=1, which='minor')
#plt.xticks([0,1,2,3,4],fontsize=24)
#plt.tight_layout()
#plt.subplots_adjust(hspace=0.1,wspace=0)


# Plot terminus position along flowline
#plt.figure(figsize=(8,4))
#pylab.rcParams['xtick.major.pad']='8'
#pylab.rcParams['ytick.major.pad']='8'
#jet = cm = plt.get_cmap('RdYlBu_r')
#cNorm  = colors.Normalize(vmin=6, vmax=10)
#scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
#CS3 = plt.contourf([[0,0],[0,0]],np.arange(6,10.1,.1), cmap=jet)

#floatheight = -(flowline[:,3])*(1.024/0.917-1)
#ind = np.where((terminus_time > 2008))
#terminus_interp = terminus_val[ind]
#velocity_interp = np.interp(terminus_time[ind],velocity_time[:,0],velocity_val[:,0])
#surf_interp = np.interp(terminus_interp,flowline[:,0],floatheight)
#bed_interp = np.interp(terminus_interp,flowline[:,0],flowline[:,3])

#plt.plot(flowline[0:-74,0]/1e3-82,flowline[0:-74,4],'k',linewidth=3)

#for i in range(0,len(surf_interp)):
  #velcolor = scalarMap.to_rgba(velocity_interp[i]/1e3)
#  plt.plot([terminus_interp[i]/1e3-82,terminus_interp[i]/1e3-82],[surf_interp[i],bed_interp[i]],color='0.8',linewidth=2)
#plt.plot(flowline[:,0]/1e3-82,floatheight,'k--',linewidth=3,label='Flotation height')
#pylab.fill_between(flowline[:,0]/1e3-82,-10000,flowline[:,3],facecolor=[139/255.0,90/255.0,43/255.0],alpha=0.5)
#plt.plot(flowline[:,0]/1e3-82,flowline[:,3],'k',linewidth=3)

#plt.plot([flowline[-16,0]/1e3-82,flowline[-16,0]/1e3-82],np.array([flowline[-16,4]-10,flowline[-16,3]]),'k',linewidth=2)
#plt.plot(terminus1_rkm['coord1']/1e3-82,terminus1_rkm['coord2'],'k',linewidth=3)
#plt.plot([terminus1_rkm['coord1'][0]/1e3-82-1,terminus1_rkm['coord1'][0]/1e3-82-1],[np.interp(terminus1_rkm['coord1'][0]-1000,flowline[:,0],flowline[:,4]),np.interp(terminus1_rkm['coord1'][0]-1000,flowline[:,0],flowline[:,3])],'r',linewidth=3)
#plt.plot([terminus1_rkm['coord1'][0]/1e3-82-0.25,terminus1_rkm['coord1'][0]/1e3-82-0.25],[np.interp(terminus1_rkm['coord1'][0]-250,flowline[:,0],flowline[:,4]),np.interp(terminus1_rkm['coord1'][0]-250,flowline[:,0],flowline[:,3])],'b',linewidth=3)

#ax = plt.gca()
#ax.tick_params('both', length=10, width=1, which='major')
#plt.xlim([-12,4.7])
#plt.ylim([-1000,800])
#plt.xticks([-12,-8,-4,0,4],fontsize=24)
#plt.yticks([-1000,-500,0,500],fontsize=24)
#cbar = plt.colorbar(CS3,ticks=range(6,12))
#cbar.set_label('Velocity near terminus (km/yr)',fontsize=20)
#cbar.ax.set_ylim([cbar.norm(6), cbar.norm(11)])
#plt.legend(fontsize=20,handlelength=0)
#plt.xlabel('Distance along flowline (km)',fontsize=26)
#plt.ylabel('Elevation (m asl)',fontsize=26)
#plt.tight_layout()
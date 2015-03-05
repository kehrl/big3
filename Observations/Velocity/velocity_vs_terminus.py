# This file compares Helheim's terminus position to its velocity.

import os
import sys
import numpy as np
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Helheim/Tools"))
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
Lagrangian = 1
dists_eul = [0,5,10,20,30] # kilometers
dists_lag = np.array([1.0,5.0,10.0,20.0,30.0])
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

############################
# Get velocities at points # 
############################

# Lagrangian points

vellag_val,vellag_time,vellag_error,pos = helheim_velocity.velocity_at_lagpoints(flowline[:,1],flowline[:,2],flowline[:,0],np.array(dists_lag)*1e3)
vellag_val[80,2:4]='NaN'
vellag_val[81,:]='NaN'

indt = np.where((vellag_time > time1) & (vellag_time < time2))
vellag_time = vellag_time[indt[0]]
vellag_val = vellag_val[indt[0],:]
pos=pos[indt[0],:]

del indt

# Eulerian points

ind=[]
for i in range(0,len(dists_eul)):
  ind.append( (abs(flowline[:,0] - 1e3*(terminus - dists_eul[i]))).argmin() )

veleul_val,veleul_time, veleul_error = helheim_velocity.velocity_at_eulpoints(flowline[ind,1],flowline[ind,2])
veleul_val[81,0] = 'NaN' # Ignore bad velocities
veleul_val[80,1] = 'NaN' # Ignore bad velocities
veleul_val[80:82,2] = 'NaN' # Ignore bad velocities
velocitypoints = np.column_stack([flowline[ind,1],flowline[ind,2]])

indt = np.where((veleul_time > time1) & (veleul_time < time2))
veleul_time = veleul_time[indt[0]]
veleul_val = veleul_val[indt[0],:]

del indt

if Lagrangian == 1:
  velocity_val = vellag_val
  velocity_time = vellag_time
  dists = dists_lag
else:
  velocity_val = veleul_val
  velocity_time = veleul_time
  dists_eul

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
bed2001dist = np.interp(bed2001[:,0],flowline[:,1],flowline[:,0])

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
  xTickPos = np.linspace(time1-0.25,time2-0.25,(time2-time1)*2+1)
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
  xTickPos = np.linspace(time1-0.25,time2-0.25,(time2-time1)*2+1)
  xTickPos = xTickPos[:-1]
  ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]),bottom=min(plt.ylim()), color=['0.9','w'],linewidth=0)
if terminusbed == 1:
  ax.set_xticklabels([])

if terminusbed == 1:
  terminus_bed = np.interp(terminus_val,bed2001dist,bed2001[:,2])
  plt.subplot(gs[3, :])
  ax = plt.gca()
  ax.tick_params('both', length=20, width=2, which='major')
  ax.tick_params('both', length=10, width=1, which='minor')
  plt.ylabel('Bed at terminus \n (m)',fontsize=12)
  plt.plot(terminus_time,terminus_bed,'ko-',linewidth=2,markersize=5)
  plt.xticks(range(2000,2016),fontsize=12)
  labels=[]
  for i in range(2000,2017):
    labels.append('Jan \n'+str(i))
  ax.set_xticklabels(labels)
  plt.xlim([time1,time2]) 
  #ax.xaxis.set_major_formatter(x_formatter)
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  ax.yaxis.set_minor_locator(AutoMinorLocator(2))
  plt.yticks([-700,-600,-500])
  if seasonal:
    xTickPos = np.linspace(time1-0.25,time2-0.25,(time2-time1)*2+1)
    xTickPos = xTickPos[:-1]
    ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.9','w'],linewidth=0)

plt.subplots_adjust(hspace=0.1,wspace=0) 
plt.tight_layout()
#plt.savefig(os.path.join(os.getenv("HOME"),'Bigtmp/fig_velocity_time_'+str(time1)+'to'+str(time2)+'.pdf'),FORMAT='PDF')
#plt.close()

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
plt.savefig(os.path.join(os.getenv("HOME"),'Bigtmp/fig_velocity_map_'+str(time1)+'to'+str(time2)+'.pdf'),FORMAT='PDF')
plt.close()

###################################################
# Plot terminus position, bed depth, vs. velocity #
###################################################

plt.figure(figsize=(10,5))

# Terminus position vs. velocity
terminus_interped = np.interp(velocity_time,terminus_time,terminus_val)
bed_interped = np.interp(terminus_interped,bed2001dist,bed2001[:,2])


plt.subplot(121) 
ind=np.where(~np.isnan(velocity_val[:,0]))
R = np.corrcoef(terminus_interped[ind,0],velocity_val[ind,0])
plt.plot(terminus_interped/1e3-82,velocity_val[:,0]/1e3,'k.')
plt.text(3,10,"R$^2$ = %.2f" % R[0,1]**2)
plt.xlabel('Terminus position (km)')
plt.ylabel('Glacier velocity (km/yr)')
plt.xticks([0,1,2,3,4])

plt.subplot(122)
plt.plot(bed_interped,velocity_val[:,0]/1e3,'k.')
R = np.corrcoef(bed_interped[ind,0],velocity_val[ind,0])
plt.xlabel('Bed at terminus (m)')
plt.ylabel('Glacier velocity (km/yr)')
plt.xticks([-700,-600,-500])
plt.text(-550,10,"R$^2$ = %.2f" % R[0,1]**2)

###########################
# Plot velocity radargram #
###########################

flowline_v,flowline_t,termini = helheim_velocity.velocity_along_flowline(flowline[:,1],flowline[:,2],flowline[:,0])
flowline_tint=np.linspace(2009,2015,1001)
term_int=np.zeros([len(flowline_tint),1])
flowline_vint = np.zeros([len(flowline[:,0]),len(flowline_tint)])

for i in range(0,len(flowline_tint)):
  ind = np.argmin(abs(flowline_tint[i]-flowline_t))
  term_int[i] = termini[ind]
  flowline_vint[:,i] = np.log10(flowline_v[:,ind])

plt.figure(figsize=(6,6))
ax=plt.gca()
plt.imshow(flowline_vint,extent=[flowline_tint[0],flowline_tint[-1],flowline[-1,0]/1e3,flowline[0,0]],interpolation='nearest',aspect='auto',norm=matplotlib.colors.LogNorm())
plt.fill_between(terminus_time,terminus_val/1e3,88,color='w')
plt.plot(terminus_time,terminus_val/1e3,'k.-')
plt.xlim([2009,2014.5])
ax.get_xaxis().get_major_formatter().set_useOffset(False)
plt.ylim([87,57])
plt.xticks([2009,2010,2011,2012,2013,2014,2015])
plt.clim(np.log10(4000),np.log10(9000))
cb=plt.colorbar()
cb.set_ticks([np.log10(4000),np.log10(5000),np.log10(6000),np.log10(7000),np.log10(8000),np.log10(9000)])
cb.set_ticklabels(['4','5','6','7','8','9'])
cb.set_label('Velocity (km/yr)')
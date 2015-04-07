# This code considers thinning near the terminus of Helheim Glacier due to 
# longitudinal/transverse strain. Do we expect the terminus to be grounded or floating?

# LMK, UW, 3/2/15

import os
import sys
import numpy as np
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Helheim/Tools"))
import helheim_velocity, helheim_icefronts, helheim_bed, helheim_elevation
import matplotlib.pyplot as plt
import geotiff, dist

##########
# Inputs #
##########
time1=2009
time2=2015

#############
# Load data #
#############

# Load Helheim bed
bed = helheim_bed.cresis('2001')
bed = bed[17:,:] #throw out points that don't follow a flowline
dists = dist.transect(bed[:,0],bed[:,1])-5165 # Last value to offset distances so similar zero point as that in "velocity_vs_terminus"

# Terminus positions
term, termt = helheim_icefronts.distance_along_flowline(bed[:,0],bed[:,1],dists,'icefront')
rift, riftt = helheim_icefronts.distance_along_flowline(bed[:,0],bed[:,1],dists,'rift')

# Load elevations
elev,elevt = helheim_elevation.worldview_at_pts(bed[:,0],bed[:,1],32,[])

# Calculate divergence
v,vx,vy,divx,divy,vt = helheim_velocity.divergence_at_eulpoints(bed[:,0],bed[:,1])
indices = np.where((vt > time1) & (vt < time2))
v=v[indices[0],:]
vx=vx[indices[0],:]
vy=vy[indices[0],:]
divx=divx[indices[0],:]
divy=divy[indices[0],:]
vt=vt[indices[0]]


# Bad transects (81,82,87) should maybe be removed.

# Divergence is noisy, so let's remove spikes. Right now I'm just removing spikes by using
# an ad hoc cutoff of > 1.5 in divy. This should be fixed if we want to develop this further.
ind = np.where(abs(divy) > 1.5)
divx[ind[0]] = 'NaN'
divy[ind[0]] = 'NaN'
v[ind[0]] = 'NaN'
vx[ind[0]] = 'NaN'
vy[ind[0]] = 'NaN'

# Delete points in front of the ice front
for i in range(0,len(vt)):
  terminus = np.interp(vt[i],termt,term[:,0]) 
  ind = np.where(dists > terminus)
  divx[i,ind[0]]='NaN'
  divy[i,ind[0]]='NaN'
  v[i,ind[0]] = 'NaN'
  vx[i,ind[0]] = 'NaN'
  vy[i,ind[0]] = 'NaN'
  
# Let's get average divergence over time interval from time1 to time2
avedivx = np.zeros([len(dists)])
avedivy = np.zeros([len(dists)])
avevx = np.zeros([len(dists)])
avevy = np.zeros([len(dists)])
avev = np.zeros([len(dists)])

for i in range(0,len(dists)):
  ind = np.where(~(np.isnan(divx[:,i])))
  avedivx[i] = np.mean(divx[ind[0],i])
  ind = np.where(~(np.isnan(divy[:,i])))
  avedivy[i] = np.mean(divy[ind[0],i])
  ind = np.where(~(np.isnan(v[:,i])))
  avevx[i] = np.mean(vx[ind[0],i])
  avevy[i] = np.mean(vy[ind[0],i])
  avev[i] = np.mean(v[ind[0],i])
#del divx, divy, vy, v

# Now let's start looking at shortening along the flowline due to horizontal stretching
timestep = 1/365.25 # in years
divz=-(divx+divy)
avedivz = -(avedivx+avedivy)

startind = 2 # Distance along flowline to set as original height for stretching (we should play around with this value)
startelev = np.mean(elev[startind,~(np.isnan(elev[startind,:]))])

# For averages
avedists=np.zeros([500])
avedists[:]='NaN'
aveheights=np.zeros([500])
aveheights[:] = 'NaN'
satisfied=0
n=0
pos=dists[startind]
avedists[0]=pos
Ho = startelev - bed[startind,2]
aveheights[0]=Ho
while not(satisfied):
  ind = np.argmin(abs(pos-dists))
  H=avedivz[ind]*Ho*timestep+Ho
  pos=pos+avev[ind]*timestep
  if (pos > 8000-5165) or (np.isnan(pos)):
    satisfied = 1
  n = n+1
  Ho=H
  avedists[n]=pos
  aveheights[n]=H 
  
aveelevs = np.zeros([500])
avebeds = np.zeros([500])
avebeds[:]='NaN'
aveelevs[:]='NaN'
meanelev = np.zeros([len(dists)])
for i in range(0,len(dists)):
  meanelev[i]=np.mean(elev[i,~(np.isnan(elev[i,:]))])
nonnan = np.where(~(np.isnan(avedists)))
bedint = np.interp(avedists[nonnan[0]],dists[:,0],bed[:,2])
elevint = np.interp(avedists[nonnan[0]],dists[:,0],meanelev)
aveelevs[nonnan[0]]=bedint+aveheights[nonnan[0]]
avebeds[nonnan[0]]=elevint-aveheights[nonnan[0]]
del bedint, elevint

# For all velocity transects
estdists=np.zeros([500,len(vt)])
estdists[:,:]='NaN'
estheights=np.zeros([500,len(vt)])
estheights[:,:]='NaN'
for i in range(0,len(vt)):
  satisfied=0
  n=0
  pos=dists[startind]
  Ho = startelev - bed[startind,2]
  while not(satisfied):
    ind = np.argmin(abs(pos-dists))
    H=divz[i,ind]*Ho*timestep+Ho
    pos=pos+v[i,ind]*timestep
    if (pos > 8000-5165) or (np.isnan(pos)):
      satisfied = 1
    n = n+1
    Ho=H
    estdists[n,i]=pos
    estheights[n,i]=H 

estelevs = np.zeros([500,len(vt)])
estbeds = np.zeros([500,len(vt)])
meanelev = np.zeros([len(dists)])
for i in range(0,len(dists)):
  meanelev[i]=np.mean(elev[i,~(np.isnan(elev[i,:]))])
for i in range(0,len(vt)):
  nonnan = np.where(~(np.isnan(estdists[:,i])))
  bedint = np.interp(estdists[nonnan[0],i],dists[:,0],bed[:,2])
  elevint = np.interp(estdists[nonnan[0],i],dists[:,0],meanelev)
  estelevs[nonnan[0],i]=bedint+estheights[nonnan[0],i]
  estbeds[nonnan[0],i]=elevint-estheights[nonnan[0],i]
  del bedint, elevint

###################  
# Make some plots #
###################

floatH=(1.02/0.9)*abs(bed[:,2])+bed[:,2]
morlighem = helheim_bed.morlighem(bed[:,0],bed[:,1])
morlighem[-41:]='NaN'

# Plot estimated surface elevations
plt.figure(figsize=(8,5))
ax=plt.gca()
plt.plot(estdists[:,0]/1e3,estelevs[:,0],'0.6',label='All divu=0')
plt.plot(estdists/1e3,estelevs,'0.6')
plt.plot(dists/1e3,elev[:,0],color=[1,0,0],label='Worldview elevation',linewidth=1.5)
plt.plot(dists/1e3,elev,color=[1,0,0],linewidth=1.5)
plt.plot(avedists/1e3,aveelevs,'k',label='Ave divu=0',linewidth=1.5)
plt.plot(dists/1e3,floatH,'k--',label='Flotation height')
plt.plot(dists/1e3,bed[:,2],'b',label='2001 bed',linewidth=1.5)
#plt.plot(dists/1e3,morlighem,'g',label='Morlighem bed',linewidth=1.5)
plt.plot(dists/1e3,gimp,'y',label='Gimp dem',linewidth=1.5)
plt.xlim([-5.3,4.2])
plt.ylim([-800,400])
plt.xlabel('Distance along flowline (km)')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#indices=np.where((termt > time1) & (termt < time2))
plt.plot(term/1e3,np.zeros(len(term))-200,'k.',label='Terminus positions')
#indices=np.where((riftt > time1) & (riftt < time2))
plt.plot(rift/1e3,np.zeros(len(rift))-100,'r.',label='Rift positions')
ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5),fontsize=10)
plt.ylabel('Elevation (m asl)')
plt.savefig(os.path.join(os.getenv("HOME"),'Bigtmp/fig_divelev_term_'+str(time1)+'to'+str(time2)+'.pdf'),FORMAT='PDF')

# Plot estimated bed elevations
plt.figure(figsize=(8,5))
ax=plt.gca()
plt.plot(estdists[:,0]/1e3,estbeds[:,0],'0.6',label='All divu=0')
plt.plot(estdists/1e3,estbeds,'0.6')
plt.plot(dists/1e3,elev[:,0],color=[1,0,0],label='Worldview elevation',linewidth=1.5)
plt.plot(dists/1e3,elev,color=[1,0,0],linewidth=1.5)
plt.plot(avedists/1e3,avebeds,'k',label='Ave divu=0',linewidth=1.5)
plt.plot(dists/1e3,floatH,'k--',label='Flotation height')
plt.plot(dists/1e3,bed[:,2],'b',label='2001 bed',linewidth=1.5)
plt.plot(dists/1e3,morlighem,'g',label='Morlighem bed',linewidth=1.5)
plt.xlim([-5.3,4])
plt.ylim([-900,300])
plt.xlabel('Distance along flowline (km)')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5),fontsize=10)
plt.ylabel('Elevation (m asl)')
plt.savefig(os.path.join(os.getenv("HOME"),'Bigtmp/fig_divbed_'+str(time1)+'to'+str(time2)+'.pdf'),FORMAT='PDF')


# Plot velocity and divergence terms
plt.figure(figsize=(8,5))
plt.subplot(311)
ax=plt.gca()
plt.plot(dists/1e3,v.T/1e3,'0.6')
plt.plot(dists/1e3,avev/1e3,'k',linewidth=1.5)
plt.ylabel('Glacier velocity (km/yr)')
plt.xlim([-5.3,4])
ax.set_xticklabels([])

plt.subplot(312)
ax=plt.gca()
plt.plot(dists/1e3,divx.T,'0.6')
plt.plot(dists/1e3,avedivx,'k',linewidth=1.5)
plt.ylabel('du/dx')
plt.xlim([-5.3,4])
plt.ylim([-9,9])
plt.yticks([-5,0,5])
ax.set_xticklabels([])

plt.subplot(313)
plt.plot(dists/1e3,divy.T,'0.6')
plt.plot(dists/1e3,avedivy,'k',linewidth=1.5)
plt.xlabel('Distance along flowline (km)')
plt.xlim([-5.3,4])
plt.yticks([-1,0,1])
plt.ylabel('dv/dy')

plt.savefig(os.path.join(os.getenv("HOME"),'Bigtmp/fig_veldiv_'+str(time1)+'to'+str(time2)+'.pdf'),FORMAT='PDF')


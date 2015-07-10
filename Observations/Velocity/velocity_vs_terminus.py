# This file compares Helheim's terminus position to its velocity.

import os
import sys
import numpy as np
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
sys.path.append(os.path.join(os.getenv("HOME"),"Code/BigThreeGlaciers/Tools"))
import velocity, icefronts, bed
import matplotlib.pyplot as plt
import matplotlib, geotiff
from matplotlib.ticker import AutoMinorLocator
import shapefile
import scipy.signal as signal
import jdcal
import pylab
from shapely.geometry import LineString
import matplotlib.cm as cmx
import matplotlib.colors as colors

##########
# Inputs #
##########

# Glacier
glacier = 'Helheim' # Options: Helheim, Kanger

# Locations for velocities
Lagrangian = 0
dists_eul = [0.0,5.0,10.0,20.0,30.0] # kilometers
dists_lag = np.array([1.0,5.0,10.0,20.0,30.0])

# Image for plotting
if glacier == "Helheim":
  image = geotiff.read(os.path.join(os.getenv("HOME"),"Data/Mosaics/Helheim/mosaicHelheim.2014-159.148.38725_1-20mgeo.tif"))
  image = geotiff.readrgb(os.path.join(os.getenv("HOME"),"Data/Imagery/ASTER/Helheim/20130812141113_AST2125988263_2.tif"))
elif glacier == "Kanger":
  image = geotiff.read(os.path.join(os.getenv("HOME"),"Data/Mosaics/Kanger/mosaicKang.2014-160.163.38740_1-20mgeo.tif"))
    
################
# Plot options #
################

time1 = 2008 #start time for plot
time2 = 2015.5 # end time for plot
seasonal = 1 # plot seasonal bars, to highlight seasonal trends
terminusbed = 0

############ 
# Flowline #
############

x,y,dists = glacier_flowline.load(glacier)

##################
# Get ice fronts #
##################

terminus_val, terminus_time = icefronts.distance_along_flowline(x,y,dists,glacier,type='icefront')
rift_val, rift_time = icefronts.distance_along_flowline(x,y,dists,glacier,type='rift')

# Chop to desired time interval
indt = np.where((terminus_time > time1) & (terminus_time < time2))
terminus_time = terminus_time[indt[0]]
terminus_val = terminus_val[indt[0]]
indt = np.where((rift_time > time1) & (rift_time < time2))
rift_val = rift_val[indt[0]]
rift_time = rift_time[indt[0]]
del indt

# Reference for terminus position
terminus = np.mean(terminus_val)

############################
# Get velocities at points # 
############################

# Eulerian points

ind=[]
for i in range(0,len(dists_eul)):
  ind.append( (abs(dists - (terminus - dists_eul[i]*1e3))).argmin() )

veleul_val,veleul_time, veleul_error = velocity.velocity_at_eulpoints(x[ind],y[ind],glacier)
velocitypoints = np.column_stack([x[ind],y[ind]])

# Chop to desired time interval
indt = np.where((veleul_time > time1) & (veleul_time < time2))
veleul_time = veleul_time[indt[0]]
veleul_val = veleul_val[indt[0],:]
del indt

# Lagrangian points

#vellag_val,vellag_time,vellag_error,vellag_dist,vellag_x,vellag_y = velocity.velocity_at_lagpoints(flowline[:,1],flowline[:,2],dists,np.array(dists_lag)*1e3,Helheim)
#vellag_val[80,2:4]='NaN'
#vellag_val[81,:]='NaN'

#indt = np.where((vellag_time > time1) & (vellag_time < time2))
#vellag_time = vellag_time[indt[0]]
#vellag_val = vellag_val[indt[0],:]
#vellag_dist=vellag_dist[indt[0],:]
#vellag_x=vellag_x[indt[0],:]
#vellag_y=vellag_y[indt[0],:]

#del indt

############################################
# Plot: velocity vs. terminus through time #
############################################

#plt.figure(figsize=(4.527,3.74))
plt.figure(figsize=(6.5,4))
gs = matplotlib.gridspec.GridSpec(4,1)

# Plot terminus
plt.subplot(gs[0, :])
nonnan = np.where(~(np.isnan(terminus_val)))[0]
plt.plot(terminus_time[nonnan],terminus_val[nonnan]/1e3-terminus/1e3,'ko--',linewidth=1,markersize=3)
plt.plot(rift_time,(rift_val-terminus)/1e3,'ro',linewidth=2,markersize=3)
ax = plt.gca()
x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
ax.xaxis.set_major_formatter(x_formatter)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
plt.xticks(range(2000,2016),fontsize=8,fontname="Arial")
plt.xlim([time1,time2])
plt.yticks(np.arange(-6,8,2),fontsize=8,fontname="Arial")
#plt.yticks([0,2,4,6],fontsize=10)
plt.ylabel('Terminus position \n (km)',fontsize=8,fontname="Arial")
ax.tick_params('both', length=10, width=1.5, which='major')
ax.tick_params('both', length=5, width=1, which='minor')
if seasonal:
  xTickPos = np.linspace(np.floor(time1)-0.25,np.floor(time2)-0.25,(np.floor(time2)-np.floor(time1))*2+1)
  ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.9','w'],linewidth=0)
ax.set_xticklabels([])
plt.ylim(np.floor((np.min(terminus_val)-terminus)/1e3),np.ceil((np.max(terminus_val)-terminus)/1e3))


# Plot velocities
plt.subplot(gs[1:, :]) 
#plt.plot([2000,2014],[0,0],'k')
coloptions=['b','c','g','y','r']
for i in range(0,len(dists_eul)):
  nonnan = np.where(~(np.isnan(veleul_val[:,i])))[0]
  plt.plot(veleul_time[nonnan],(veleul_val[nonnan,i])/1e3,'o--',color=coloptions[i],label=str(-dists_eul[i])+' km',linewidth=1,markersize=3)
plt.ylabel('Glacier velocity \n (km/yr)',fontsize=8,fontname="Arial")
ax = plt.gca()
nonnan = np.where(~np.isnan(veleul_val[:,0]))[0]
if glacier == 'Kanger':
  plt.yticks(range(0,12),fontsize=8,fontname="Arial")
elif glacier == 'Helheim':
  plt.yticks(range(0,12),fontsize=8,fontname="Arial")
x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
ax.xaxis.set_major_formatter(x_formatter)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
plt.legend(loc=2,fontsize=8,numpoints=1,handlelength=0.5)
matplotlib.rc('font',family="Arial",)
ax.tick_params('both', length=10, width=1.5, which='major')
ax.tick_params('both', length=5, width=1, which='minor')
labels=[]
plt.xticks(range(2000,2017),fontsize=8,fontname="Arial")
for i in range(2000,2017):
  labels.append('Jan \n'+str(i))
ax.set_xticklabels(labels)
plt.ylim([0,12])
if seasonal:
  xTickPos = np.linspace(np.floor(time1)-0.25,np.floor(time2)-0.25,(np.floor(time2)-np.floor(time1))*2+1)
  ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.9','w'],linewidth=0)
plt.xlim([time1,time2])


plt.subplots_adjust(hspace=0.05,wspace=0) 
plt.savefig(os.path.join(os.getenv("HOME"),'Bigtmp/'+glacier+'_veleul_time_'+str(time1)+'to'+str(time2)+'.pdf'),FORMAT='PDF')
plt.close()

#########################################
# Plot overview map for previous figure #
#########################################
plt.figure(figsize=(4,4))

jet = cm = plt.get_cmap('rainbow') 
cNorm  = colors.Normalize(vmin=0.0, vmax=1.0)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
Z = [[0,0],[0,0]]
levels = np.arange(0,365,1)
CS3 = plt.contourf(Z, levels, cmap=jet)
plt.clf()

image_extent = [np.min(image[0]),np.max(image[0]),np.max(image[1][np.nonzero(image[1])]),np.min(image[1])]
plt.imshow(image[2],extent=image_extent,cmap='Greys_r')
plt.gca().invert_yaxis()
plt.axis('equal')
plt.xlim([min(velocitypoints[:,0])-5.0e3,max(velocitypoints[:,0])+5.0e3])
plt.ylim([min(velocitypoints[:,1])-5.0e3,max(velocitypoints[:,1])+5.0e3])
ax = plt.gca()

termx,termy,termt = icefronts.load_all(time1,time2,'icefront',glacier)
riftx,rifty,riftt = icefronts.load_all(time1,time2,'rift',glacier)
for i in range(0,len(termx[0,:])):
  value = scalarMap.to_rgba(termt[i]-np.floor(termt[i]))
  plt.plot(termx[:,i],termy[:,i],color=value)
plt.plot(x,y,'k',linewidth=1.5)
for i in range(0,len(velocitypoints)):
  plt.plot(velocitypoints[i,0],velocitypoints[i,1],'o',color=coloptions[i],markersize=7)
ax.set_xticklabels([])
ax.set_yticklabels([])
plt.colorbar(CS3)

plt.savefig(os.path.join(os.getenv("HOME"),'Bigtmp/'+glacier+'_velocity_map_'+str(time1)+'to'+str(time2)+'.pdf'),FORMAT='PDF')
plt.close()

###########################
# Plot velocity radargram #
###########################

flowline_v,flowline_t,termini = velocity.velocity_along_flowline(x,y,dists,glacier)
flowline_tint=np.linspace(2008.5,2015,1051)
term_int=np.zeros([len(flowline_tint),1])
flowline_vint = np.zeros([len(dists),len(flowline_tint)])
flowline_vintlog = np.zeros([len(dists),len(flowline_tint)])
flowline_vintper = np.zeros([len(dists),len(flowline_tint)])

# Filter velocities
filt_len=1000.0 # meters
cutoff=(1/filt_len)/(1/(np.diff(dists[1:3])*2))
b,a=signal.butter(4,cutoff,btype='low')
flowline_vfilt=np.zeros_like(flowline_v)
for i in range(0,len(flowline_t)):
  nonnan = np.where(~(np.isnan(flowline_v[:,i])))[0]
  if len(nonnan) > 1:
    interped = np.interp(dists[nonnan[0]:nonnan[-1]],dists[nonnan],flowline_v[nonnan,i])
    flowline_vfilt[nonnan[0]:nonnan[-1],i]=signal.filtfilt(b,a,interped)
    flowline_vfilt[np.where(np.isnan(flowline_v))]='NaN'
    flowline_vfilt[flowline_vfilt==0]='NaN'

for i in range(0,len(flowline_tint)):
  ind = np.argmin(abs(flowline_tint[i]-flowline_t))
  term_int[i] = termini[ind]
  flowline_vint[:,i]=flowline_vfilt[:,ind]
  flowline_vintlog[:,i] = np.log10(flowline_vfilt[:,ind])


meanvel=np.zeros(len(dists))
for i in range(0,len(dists)):
  nonnan = np.where(~np.isnan(flowline_vint[i,:]))
  meanvel[i] = np.mean(flowline_vint[i,nonnan])  
  
for i in range(0,len(flowline_tint)):
  ind = np.argmin(abs(flowline_tint[i]-flowline_t))
  flowline_vintper[:,i]=((flowline_vfilt[:,ind]-meanvel))

# Plot linear scale
plt.figure(figsize=(6.5,5))
plt.subplot(121)
ax=plt.gca()
plt.imshow(flowline_vint/1e3,extent=[flowline_tint[0],flowline_tint[-1],(dists[-1]-terminus)/1e3,(dists[0]-terminus)/1e3],interpolation='nearest',aspect='auto')
plt.fill_between(terminus_time,(terminus_val-terminus)/1e3,5,color='w')
plt.plot(terminus_time,(terminus_val-terminus)/1e3,'ko--',markersize=3,linewidth=1)
ax.get_xaxis().get_major_formatter().set_useOffset(False)
plt.xticks(range(2000,2017),fontsize=8,fontname="Arial")
labels=[]
for i in range(2000,2017):
  labels.append('Jan \n'+str(i))
ax.set_xticklabels(labels)
#for i in np.array([4,6,8]):
#  junk = np.zeros(len(flowline_tint))
#  for j in range(0,len(flowline_tint)):
#    nonnan = np.where(~(np.isnan(flowline_vint[:,j])))
#    ind = np.argmin(abs(flowline_vint[nonnan,j]/1e3-i))
#    junk[j] = flowline[nonnan[0][ind],0]/1e3
#  plt.plot(flowline_tint,junk,'k',linewidth=1.5)
plt.xlim([2008.5,2014.5])
plt.clim(3,9)
cb=plt.colorbar(orientation="horizontal",fraction=0.07)
cb.set_ticks([3,4,5,6,7,8,9])
cb.set_ticklabels(['3','4','5','6','7','8','9'])
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
cb.set_label('Glacier velocity (km/yr)',fontsize=10,fontname="Arial")
plt.ylabel('Distance along flowline (km)',fontsize=10,fontname="Arial")
ax.tick_params('both', length=10, width=1.5, which='major')
ax.tick_params('both', length=5, width=1, which='minor')
plt.ylim([5,-40])
plt.yticks(fontsize=8,fontname="Arial")
#plt.text(2009,50,"a",color='w',fontsize=22,fontweight="bold")

# Plot percentage
plt.subplot(122)
ax=plt.gca()
plt.imshow(flowline_vintper/1e3,extent=[flowline_tint[0],flowline_tint[-1],(dists[-1]-terminus)/1e3,(dists[0]-terminus)/1e3],interpolation='nearest',aspect='auto',cmap='bwr')
plt.fill_between(terminus_time,(terminus_val-terminus)/1e3,5,color='w')
plt.plot(terminus_time,(terminus_val-terminus)/1e3,'ko--',markersize=3,linewidth=1)
ax.get_xaxis().get_major_formatter().set_useOffset(False)
plt.xticks(range(2000,2017),fontsize=8,fontname="Arial")
plt.yticks(fontsize=8,fontname="Arial")
labels=[]
for i in range(2000,2017):
  labels.append('Jan \n'+str(i))
ax.set_xticklabels(labels)
plt.xlim([2008.5,2014.5])
plt.clim(-0.5,0.5)
cb=plt.colorbar(orientation="horizontal",fraction=0.07)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
cb.set_ticks([-0.4,-0.2,0,0.2,0.4])
cb.set_label('Change from average (km/yr)',fontsize=10,fontname="Arial")
ax.set_yticklabels([])
plt.ylim([5,-40])
ax.tick_params('both', length=10, width=1.5, which='major')
ax.tick_params('both', length=5, width=1, which='minor')
#plt.text(2009,50,"b",fontsize=22,fontweight="bold")

plt.subplots_adjust(hspace=0,wspace=-0.1) 
plt.tight_layout()

plt.savefig(os.path.join(os.getenv("HOME"),'Bigtmp/'+glacier+'_velocity_radargram_lines_'+str(time1)+'to'+str(time2)+'.pdf'),FORMAT='PDF')
plt.close()

######################
# Plot bed, terminus #
######################

# Plot bed near the terminus
fig = plt.figure(figsize=(5,3))
ax1 = fig.add_subplot(111)
ax2 = ax1.twinx()
ax1.set_zorder(ax2.get_zorder()+1)
ax1.patch.set_visible(False)
if glacier == 'Helheim':
  #years = ['2001']
  years = ['2001','2006b','2008b','2009','2012','2013','2014']
elif glacier == 'Kanger':
  years = ['2008','2009a','2009b','2012','2013','2014']
for year in years:
  bedpts = bed.cresis(year,glacier,'geoid')
  minind = np.argmin(abs(x-np.min(bedpts[:,0])))
  maxind = np.argmin(abs(x-np.max(bedpts[:,0])))
  ind = range(minind,maxind)
  bed_interp = np.zeros(len(ind))
  for i in range(0,len(ind)):
    j = ind[i]
    k = np.argmin((x[j]-bedpts[:,0])**2+(y[j]-bedpts[:,1]))
    bed_interp[i] = bedpts[k,2]
  ax1.plot((dists[minind:maxind]-terminus)/1e3,bed_interp,label=year)

ax1.set_xlabel('Distance from terminus (km)',fontsize=10)
ax1.set_ylabel('Elevation (m asl)',fontsize=10)
ax1.tick_params(axis='x', labelsize=8)
ax1.tick_params(axis='y', labelsize=8)
ax1.set_xlim([-10,4])
plt.legend(loc=2,fontsize=8,ncol=2)
 
ax2.plot((terminus_val-terminus)/1e3,terminus_time,'.-',color='0.7')
ax2.set_ylabel('Time',fontsize=10,color='0.6')
for tl in ax2.get_yticklabels():
  tl.set_color('0.6')
ax2.get_yaxis().get_major_formatter().set_useOffset(False)
ax2.tick_params(axis='y', labelsize=8)
ax2.set_xlim([-10,4])
plt.tight_layout()

plt.savefig(os.path.join(os.getenv("HOME"),'Bigtmp/'+glacier+'_bed.pdf'),FORMAT='PDF')
plt.close()


# Plot map to show location of bed profiles
plt.figure(figsize=(4,4))
image_extent = [np.min(image[0]),np.max(image[0]),np.max(image[1][np.nonzero(image[1])]),np.min(image[1])]
plt.imshow(image[2],extent=image_extent,cmap='Greys_r')
plt.gca().invert_yaxis()
plt.axis('equal')
plt.xlim([min(velocitypoints[:,0])-5.0e3,max(velocitypoints[:,0])+5.0e3])
plt.ylim([min(velocitypoints[:,1])-5.0e3,max(velocitypoints[:,1])+5.0e3])
ax = plt.gca()

for year in years:
  bedpts = bed.cresis(year,glacier,'geoid')
  plt.plot(bedpts[:,0],bedpts[:,1],label=year)
plt.legend(fontsize=8,ncol=2)
ax.set_xticklabels([])
ax.set_yticklabels([])

plt.savefig(os.path.join(os.getenv("HOME"),'Bigtmp/'+glacier+'_bed_map.pdf'),FORMAT='PDF')
plt.close()


###################################################
# Plot terminus position, bed depth, vs. velocity #
###################################################

#plt.figure(figsize=(7,7))
#lagind=0
#eulind=0

# Terminus position vs. velocity
#terminus_interped = np.interp(veleul_time,terminus_time,terminus_val)
#bed_interped = np.interp(terminus_interped,dists,flowline[:,3])

#plt.subplot(221) 
#ax = plt.gca()
#ind=np.where(~np.isnan(veleul_val[:,eulind]))
#R = np.corrcoef(terminus_interped[ind,0],veleul_val[ind,0])
#plt.plot(terminus_interped/1e3-terminus,veleul_val[:,0]/1e3,'ko',markerfacecolor='k')
#plt.text(2.8,10.2,"R$^2$ = %.2f" % R[0,1]**2)
#plt.xlabel('Terminus position (km)')
#plt.ylabel('Eulerian velocity \n at '+str(int(np.round(dists_eul[eulind])))+' km (km/yr)')
#plt.xticks([0,1,2,3,4])
#plt.yticks(range(0,12))
#plt.ylim([np.min(np.floor(veleul_val[ind,eulind]/1e3)),np.max(np.ceil(veleul_val[ind,eulind]/1e3))])
#ax.set_xticklabels([])

#plt.subplot(222)
#ax=plt.gca()
#plt.plot(bed_interped,veleul_val[:,0]/1e3,'ko',markerfacecolor='k')
#R = np.corrcoef(bed_interped[ind,0],veleul_val[ind,0])
#plt.xlabel('Bed at terminus (m)')
#plt.ylabel('Eulerian velocity (km/yr)')
#plt.xticks([-700,-600,-500])
#plt.text(-580,10.2,"R$^2$ = %.2f" % R[0,1]**2)
#plt.yticks(range(0,12))
#plt.ylim([np.min(np.floor(veleul_val[ind,eulind]/1e3)),np.max(np.ceil(veleul_val[ind,eulind]/1e3))])
#ax.set_xticklabels([])
#ax.set_yticklabels([])

#plt.subplot(223)
#ax=plt.gca()
#ind=np.where(~np.isnan(vellag_val[:,0]))
#R = np.corrcoef(terminus_interped[ind,0],vellag_val[ind,1])
#plt.plot(terminus_interped/1e3-terminus,vellag_val[:,0]/1e3,'ko',markerfacecolor='k')
#plt.text(2.8,10.2,"R$^2$ = %.2f" % R[0,1]**2)
#plt.xlabel('Terminus position (km)')
#plt.ylabel('Lagrangian velocity \n at '+str(int(np.round(dists_lag[lagind])))+' km (km/yr)')
#plt.xticks([0,1,2,3,4])
#plt.yticks(range(0,12))
#plt.ylim([np.min(np.floor(veleul_val[ind,lagind]/1e3)),np.max(np.ceil(veleul_val[ind,lagind]/1e3))])


#plt.subplot(224)
#ax=plt.gca()
#plt.plot(bed_interped,vellag_val[:,0]/1e3,'ko',markerfacecolor='k')
#R = np.corrcoef(bed_interped[ind,0],vellag_val[ind,1])
#plt.xlabel('Bed at terminus (m)')
#plt.ylabel('Eulerian velocity (km/yr)')
#plt.xticks([-700,-600,-500])
#plt.text(-580,10.2,"R$^2$ = %.2f" % R[0,1]**2)
#ax.set_yticklabels([])
#plt.yticks(range(0,12))
#plt.ylim([np.min(np.floor(veleul_val[ind,lagind]/1e3)),np.max(np.ceil(veleul_val[ind,lagind]/1e3))])

#plt.subplots_adjust(hspace=0.1,wspace=0.1) 
#plt.tight_layout()
#plt.savefig(os.path.join(os.getenv("HOME"),'Bigtmp/fig_velocity_terminus_'+str(time1)+'to'+str(time2)+'.pdf'),FORMAT='PDF')
#plt.close()
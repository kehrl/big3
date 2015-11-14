# This file takes the Helheim IceBridge data and plots it along a flowline.

import os
import sys
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"Util/Modules/"))
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"BigThreeGlaciers/Tools/"))
import glacier_flowline, elevation, icefronts, flotation
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib


##########
# Inputs #
##########

args = sys.argv
glacier = args[1][:] # Options: Kanger, Helheim

if glacier == 'Helheim':
  x,y,zb,dists = glacier_flowline.load(glacier,verticaldatum='geoid',shapefilename='flowline_flightline',bedmodel='aniso',bedsmoothing=4)
else:
  x,y,zb,dists = glacier_flowline.load(glacier)

dists_eul = -1.0*np.array([2.0,5.0,10.0,15.0,20.0]) # kilometers

xref,yref,zref = elevation.gimp_grid(np.min(x)-10.0e3,np.max(x)+10.0e3,np.min(y)-10.0e3,np.max(y)+10.0e3,glacier,verticaldatum='geoid')

################
# Get terminus #
################

terminus_val, terminus_time = icefronts.distance_along_flowline(x,y,dists,glacier,type='icefront')

# Chop to desired time interval
#indt = np.where((terminus_time > time1) & (terminus_time < time2))
#terminus_time = terminus_time[indt[0]]
#terminus_val = terminus_val[indt[0]]
#del indt

# Reference for terminus position
#terminus = np.mean(terminus_val)

ind=[]
for i in range(0,len(dists_eul)):
  ind.append( (abs(dists - dists_eul[i]*1e3)).argmin() )

################
# Get ATM data #
################

atm_data = elevation.atm_along_flowline(x,y,glacier,'all',cutoff='none',maxdist=500,verticaldatum='geoid',filt_len='none')
zpt_atm,time_atm = elevation.atm_at_pts(x[ind],y[ind],glacier,years='all',maxdist=250,verticaldatum='geoid',method='average')

##############################
# Get Worldview and TDX Data #
############################## 

zs_dem,time_dem = elevation.dem_along_flowline(x,y,glacier,years='all',cutoff='none',verticaldatum='geoid',filt_len=0.5e3)
zpt_wv,zpt_std,time_wv = elevation.dem_at_pts(x[ind],y[ind],glacier,years='all',verticaldatum='geoid',cutoff='terminus',method='average',radius=250)

############################
# GIMP elevation at points #
############################

gimp = elevation.gimp_at_pts(x[ind],y[ind],glacier,'geoid')

#################################
# Plot elevation along flowline #
#################################

if glacier == 'Helheim':
  ind1 = np.where(time_dem[:,1] == 20130804)[0]
  ind2 = np.where(time_dem[:,1] == 20140127)[0]
elif glacier == 'Kanger':
  ind1 = np.where(time_dem[:,1] == 20140628)[0]
  ind2 = np.where(time_dem[:,1] == 20141023)[0]

fig = plt.figure(figsize=(4,3))
gs = matplotlib.gridspec.GridSpec(3,1)
matplotlib.rc('font',family='Arial')

coloptions=['r','b','g','limegreen','gold']

plt.subplot(gs[0:2])
for i in range(0,len(dists_eul)):
  plt.plot([dists_eul[i],dists_eul[i]],[0,750],c=coloptions[i],lw=0.75)
ax = plt.gca()
plt.plot(dists/1e3,flotation.height(zb),'k:',linewidth=1.5,label='Flotation')
#ax.fill_between(dists/1e3,flotation.height(zb-50), flotation.height(zb+50),
#    alpha=0.2,facecolor='0.7',edgecolor='0.7',antialiased=True)
if glacier == 'Helheim':
  plt.plot(dists/1e3,atm_data['20010521'][:,2],'k',label='2001-05-21',lw=1.5)
elif glacier == 'Kanger':
  plt.plot(dists/1e3,atm_data['20010520'][:,2],'k',label='2001-05-20',lw=1.5)
plt.plot(dists/1e3,atm_data['20050518'][:,2],'0.8',label='2005-05-18',lw=1.5)
date = str(int(time_dem[ind1,1]))
plt.plot(dists/1e3,zs_dem[ind1,:].T,color='r',linewidth=1.5,label=date[0:4]+'-'+date[4:6]+'-'+date[6:])
date = str(int(time_dem[ind2,1]))
plt.plot(dists/1e3,zs_dem[ind2,:].T,color='b',linewidth=1.5,label=date[0:4]+'-'+date[4:6]+'-'+date[6:])
plt.xticks(np.arange(-30,10,5),fontsize=8)
ax.set_xticklabels([])
plt.yticks(np.arange(-1000,1000,250),fontsize=8)
plt.xlim([-21,6])
if glacier == 'Helheim':
  plt.ylim([0,670])
  plt.text(-19.5,500,'b',fontsize=9,fontweight='bold')
elif glacier == 'Kanger':
  plt.ylim([0,750])
  plt.text(-19.5,560,'b',fontsize=9,fontweight='bold')
plt.ylabel('Elevation (m asl)',fontsize=8)
plt.legend(loc=1,fontsize=8,numpoints=1,handlelength=0.6,labelspacing=0.2)

plt.subplot(gs[-1])
for i in range(0,len(dists_eul)):
  plt.plot([dists_eul[i],dists_eul[i]],[-1200,-200],c=coloptions[i],lw=0.75)
plt.plot(dists/1e3,zb,color='k',linewidth=1.5,label='Bed')
plt.xlabel('Distance from mean terminus (km)',fontsize=8)
plt.xticks(np.arange(-30,10,5),fontsize=8)
plt.xlim([-21,6])
plt.text(-19.5,-500,'c',fontsize=9,fontweight='bold')
plt.yticks(np.arange(-1000,-250,250),fontsize=8)
plt.ylim([-1150,-300])

# Save figure
plt.tight_layout()
plt.subplots_adjust(wspace=0.04,hspace=0.04)
plt.savefig(os.path.join(os.getenv("HOME"),'Bigtmp/'+glacier+'_zs_flowline.pdf'),FORMAT='PDF')
plt.close()


##########################################
# Plot elevation through time since 2000 #
##########################################

plt.figure(figsize=(4,3))
ax = plt.gca()
coloptions=['b','c','g','y','r']
# Make symbols for legend
for i in range(0,len(ind)):
  plt.plot(0,0,'s',color=coloptions[i],label=str(int(dists_eul[i]))+' km')
plt.plot(0,0,'o',color='w',label='Worldview')
plt.plot(0,0,'^',color='w',label='ATM')
for i in range(1,len(ind)):
  plt.plot(time_wv,zpt_wv[:,i]-gimp[i],'o',markersize=5,color=coloptions[i])
  plt.plot(time_atm,zpt_atm[:,i]-gimp[i],'^',markersize=5,color=coloptions[i])
x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
plt.xlim([2001,2015])
plt.grid()
plt.ylim([-10,200])
plt.xticks(np.arange(2002,2016,2),fontsize=8)
plt.yticks(np.arange(-100,200,50),fontsize=8)
plt.ylabel('Elevation relative to GIMP DEM (m)',fontsize=8)
plt.legend(loc=1,labelspacing=0.3,numpoints=1,fontsize=8)
pos = ax.get_position() 
ax.set_position([pos.x0+0.05, pos.y0,  pos.width, pos.height*1.1])
ax.xaxis.set_major_formatter(x_formatter)
xTickPos = np.linspace(2000-0.25,2016-0.25,(2016-2000)*2+1)
ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.9','w'],linewidth=0)

# Save figure
plt.savefig(os.path.join(os.getenv("HOME"),'Bigtmp/'+glacier+'_zs_time_longterm.pdf'),FORMAT='PDF')
plt.close()

#####################
# Plot overview map # 
#####################

#plt.figure(figsize=(4,4))
#plt.imshow(zref,extent=[np.min(xref),np.max(xref),np.min(yref),np.max(yref)],origin='lower',clim=[0,1200],cmap='rainbow')
#plt.colorbar()
#CS=plt.contour(xref,yref,zref,[200,400,600,800,1000,1200],colors='k')
#plt.clabel(CS, inline=1, fontsize=8)
#plt.xlim([min(x[ind])-5.0e3,max(x[ind])+5.0e3])
#plt.ylim([min(y[ind])-10.0e3,max(y[ind])+5.0e3])
#plt.plot(x,y,'k--',linewidth=1.5)
#for i in range(0,len(ind)):
#  plt.plot(x[ind[i]],y[ind[i]],'o',color=coloptions[i],markersize=15)

#plt.savefig(os.path.join(os.getenv("HOME"),'Bigtmp/'+glacier+'_zs_overview.pdf'),FORMAT='PDF')
#plt.close()
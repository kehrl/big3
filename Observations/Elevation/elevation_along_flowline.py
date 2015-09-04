# This file takes the Helheim IceBridge data and plots it along a flowline.

import os
import sys
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules/"))
sys.path.append(os.path.join(os.getenv("HOME"),"Code/BigThreeGlaciers/Tools/"))
import glacier_flowline, elevation, icefronts, flotation
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib


##########
# Inputs #
##########

#glacier = 'Kanger'
glacier = 'Helheim'

x,y,zb,dists = glacier_flowline.load(glacier)
float = flotation.height(zb,rho_i=917.0,rho_sw=1020.0)

dists_eul = -1.0*np.array([2.0,5.0,10.0,20.0,30.0]) # kilometers

time1 = 2008 #start time for plot
time2 = 2015.5 # end time for plot

xref,yref,zref = elevation.gimp_grid(np.min(x)-10.0e3,np.max(x)+10.0e3,np.min(y)-10.0e3,np.max(y)+10.0e3,glacier,verticaldatum='geoid')

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
  ind.append( (abs(dists - dists_eul[i]*1e3)).argmin() )

################
# Get ATM data #
################

atm_data = elevation.atm_along_flowline(x,y,glacier,'all',cutoff='terminus',maxdist=200,verticaldatum='geoid')
zpt_atm,time_atm = elevation.atm_at_pts(x[ind],y[ind],glacier,years='all',maxdist=200,verticaldatum='geoid')

######################
# Get Worldview Data #
###################### 

wv_data = elevation.worldview_along_flowline(x,y,glacier,years='all',cutoff='terminus',verticaldatum='geoid',filt_len='none')
zpt_wv,time_wv = elevation.worldview_at_pts(x[ind],y[ind],glacier,years='all',verticaldatum='geoid')

############################
# GIMP elevation at points #
############################

gimp = elevation.gimp_at_pts(x[ind],y[ind],glacier,'geoid')

#################################
# Plot elevation along flowline #
#################################

dates = np.sort(np.concatenate([atm_data.keys(),wv_data.keys()]))

fig = plt.figure(figsize=(4.8,4))
ax = fig.add_subplot(111)
jet = cm = plt.get_cmap('jet') 
cNorm  = colors.Normalize(vmin=0, vmax=len(dates))
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

for i in range(0,len(dates)):
    colorVal = scalarMap.to_rgba(i)
    try:
      ax.plot((dists-terminus)/1e3,atm_data[dates[i]][:,2],':',linewidth=1.5,color=colorVal,label=dates[i][0:4]+'-'+dates[i][4:6]+'-'+dates[i][6:])
    except:
      ax.plot((dists-terminus)/1e3,wv_data[dates[i]][:,2],color=colorVal,linewidth=1.2,label=dates[i][0:4]+'-'+dates[i][4:6]+'-'+dates[i][6:])
ax.plot(dists/1e3,float,'k--',linewidth=1.5,label='Flotation')
handles,labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc=2,ncol=4,bbox_to_anchor=(-0.1, -0.15),labelspacing=0.1,borderpad=0.2,columnspacing=0.3,fontsize=9)
plt.ylabel('Elevation (m asl)',fontsize=9)
plt.xlabel('Distance from mean terminus (km)',fontsize=9)
plt.xticks(np.arange(-30,10,2),fontsize=9)
plt.xlim([-10,4.5])
plt.ylim([0,400])
plt.yticks(np.arange(0,500,100),fontsize=9)
pos = ax.get_position() 
ax.set_position([pos.x0, 0.4,  pos.width*1.1, pos.height*0.73])

# Save figure
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
plt.xticks(np.arange(2002,2016,2),fontsize=9)
plt.yticks(np.arange(-100,200,50),fontsize=9)
plt.ylabel('Elevation relative to GIMP DEM (m)',fontsize=9)
plt.legend(loc=1,labelspacing=0.3,numpoints=1,fontsize=9)
pos = ax.get_position() 
ax.set_position([pos.x0+0.05, pos.y0,  pos.width, pos.height*1.1])
ax.xaxis.set_major_formatter(x_formatter)
xTickPos = np.linspace(2000-0.25,2016-0.25,(2016-2000)*2+1)
ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.9','w'],linewidth=0)

# Save figure
plt.savefig(os.path.join(os.getenv("HOME"),'Bigtmp/'+glacier+'_zs_time_longterm.pdf'),FORMAT='PDF')
plt.close()

##########################################
# Plot elevation through time since 2011 #
##########################################

plt.figure(figsize=(5,3))
ax = plt.gca()
coloptions=['b','c','g','y','r']
# Make symbols for legend
startind = np.argmin(abs(time_wv-2015.0))
for i in range(0,len(ind)):
  plt.plot(0,0,'s',color=coloptions[i],label=str(int(dists_eul[i]))+' km')
plt.plot(0,0,'o',color='w',label='Worldview')
plt.plot(0,0,'^',color='w',label='ATM')
for i in range(0,len(ind)):
  plt.plot(time_wv,zpt_wv[:,i]-gimp[i],'o',markersize=5,color=coloptions[i])
  plt.plot(time_atm,zpt_atm[:,i]-gimp[i],'^',markersize=5,color=coloptions[i])
x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
ax.xaxis.set_major_formatter(x_formatter)
labels=[]
for i in range(2000,2017):
  labels.append('Jan \n'+str(i))
plt.xticks(range(2000,2017))
ax.set_xticklabels(labels,fontsize=9,fontname='Arial')
plt.xlim([2011,2015])
plt.grid()
#plt.yticks(np.arange(-60,60,20),fontsize=9)
plt.ylabel('Elevation relative to GIMP DEM (m)',fontsize=9)
plt.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0,labelspacing=0.3,numpoints=1,fontsize=9)
pos = ax.get_position() 
ax.set_position([pos.x0, pos.y0,  pos.width-0.15, pos.height*1.05])
plt.ylim([-20,20])
plt.yticks(fontsize=9)
xTickPos = np.linspace(2000-0.25,2016-0.25,(2016-2000)*2+1)
ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.9','w'],linewidth=0)

# Save figure
plt.savefig(os.path.join(os.getenv("HOME"),'Bigtmp/'+glacier+'_zs_time_shortterm.pdf'),FORMAT='PDF')
plt.close()

#####################
# Plot overview map # 
#####################

plt.figure(figsize=(4,4))
plt.imshow(zref,extent=[np.min(xref),np.max(xref),np.min(yref),np.max(yref)],origin='lower',clim=[0,1200],cmap='rainbow')
plt.colorbar()
CS=plt.contour(xref,yref,zref,[200,400,600,800,1000,1200],colors='k')
plt.clabel(CS, inline=1, fontsize=8)
plt.xlim([min(x[ind])-5.0e3,max(x[ind])+5.0e3])
plt.ylim([min(y[ind])-10.0e3,max(y[ind])+5.0e3])
plt.plot(x,y,'k--',linewidth=1.5)
for i in range(0,len(ind)):
  plt.plot(x[ind[i]],y[ind[i]],'o',color=coloptions[i],markersize=15)

plt.savefig(os.path.join(os.getenv("HOME"),'Bigtmp/'+glacier+'_zs_overview.pdf'),FORMAT='PDF')
plt.close()
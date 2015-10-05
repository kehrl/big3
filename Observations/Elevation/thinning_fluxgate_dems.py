# This script compares elevation change from the DEMs to anticipated thinning rates from 
# the fluxgate method.

# LMK, UW, 8/31/2015

import os
import sys
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"Util/Modules/"))
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"BigThreeGlaciers/Tools/"))
import glacier_flowline, elevation, icefronts, fluxgate
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib
import geotiff

glacier = 'Helheim' 

if glacier == 'Helheim':
  xmin = 285000.0
  xmax = 320000.0
  ymin = -2588000.0
  ymax = -2566000.0
  ximage,yimage,zimage = geotiff.read(os.path.join(os.getenv("DATA_HOME"),"Mosaics/Helheim/mosaicHelheim.2014-159.148.38725_1-20mgeo.tif"))
elif glacier == 'Kanger':
  xmin = 449800.0
  xmax = 503000.0
  ymin = -2302000.0
  ymax = -2266000.0
  ximage,yimage,zimage = geotiff.read(os.path.join(os.getenv("DATA_HOME"),"Mosaics/Kanger/mosaicKang.2014-160.163.38740_1-20mgeo.tif"))
    
###########################
# dH/dt from the fluxgate #
###########################

flux_time_gate1_smith,flux_dH_gate1_smith = fluxgate.fluxgate_thinning(glacier,"fluxgate1",bedsource='smith')
flux_time_gate1_cresis,flux_dH_gate1_cresis = fluxgate.fluxgate_thinning(glacier,"fluxgate1",bedsource='cresis')
x_gate1,y_gate1 = fluxgate.fluxbox_geometry(glacier,"fluxgate1")

flux_time_gate2_smith,flux_dH_gate2_smith = fluxgate.fluxgate_thinning(glacier,"fluxgate2",bedsource='smith')
flux_time_gate2_cresis,flux_dH_gate2_cresis = fluxgate.fluxgate_thinning(glacier,"fluxgate2",bedsource='cresis')
x_gate2,y_gate2 = fluxgate.fluxbox_geometry(glacier,"fluxgate2")

flux_time_gate3_smith,flux_dH_gate3_smith = fluxgate.fluxgate_thinning(glacier,"fluxgate3",bedsource='smith')
flux_time_gate3_cresis,flux_dH_gate3_cresis = fluxgate.fluxgate_thinning(glacier,"fluxgate3",bedsource='cresis')
x_gate3,y_gate3 = fluxgate.fluxbox_geometry(glacier,"fluxgate3")

######################
# dH/dt from WV DEMs #
######################

# Load WV DEMs	
xwv,ywv,zwv,timewv = elevation.worldview_grid(glacier,xmin,xmax,ymin,ymax,years='all',verticaldatum='ellipsoid')

wv_time_gate1,wv_dH_gate1 = fluxgate.dem_thinning(glacier,xwv,ywv,zwv,timewv,"fluxgate1")
wv_time_gate2,wv_dH_gate2 = fluxgate.dem_thinning(glacier,xwv,ywv,zwv,timewv,"fluxgate2")
wv_time_gate3,wv_dH_gate3 = fluxgate.dem_thinning(glacier,xwv,ywv,zwv,timewv,"fluxgate3")

###########################
# Plot dH/dt through time #
###########################

time1 = 2008.5
time2 = 2015.5
plt.figure(figsize=(6.5,6.5))
gs = matplotlib.gridspec.GridSpec(3,1)

plt.subplot(gs[0,:])
ax = plt.gca()
nonnan = np.where(~(np.isnan(flux_dH_gate3_smith[:,0])))[0]
plt.errorbar(flux_time_gate3_smith[nonnan],flux_dH_gate3_smith[nonnan,0],flux_dH_gate3_smith[nonnan,1],linestyle='dashed',marker='o',color='b',markersize=3)
plt.plot(flux_time_gate3_cresis,flux_dH_gate3_cresis[:,0],'k+',markersize=3)
plt.plot(wv_time_gate3[:,0],wv_dH_gate3[:,0],'ko',markersize=3)
plt.xticks(range(2000,2017),fontsize=8,fontname="Arial")
ax.set_xticklabels([])

plt.subplot(gs[1,:])
ax = plt.gca()
nonnan = np.where(~(np.isnan(flux_dH_gate1_smith[:,0])))[0]
plt.errorbar(flux_time_gate1_smith[nonnan],flux_dH_gate1_smith[nonnan,0],flux_dH_gate1_smith[nonnan,1],linestyle='dashed',marker='o',color='g',markersize=3)
plt.plot(flux_time_gate1_cresis,flux_dH_gate1_cresis[:,0],'k+',markersize=3)
plt.plot(wv_time_gate1[:,0],wv_dH_gate1[:,0],'ko',markersize=3)
plt.xticks(range(2000,2017),fontsize=8,fontname="Arial")
plt.ylabel("dH/dt (m/yr)",fontsize=9,fontname='Arial')
ax.set_xticklabels([])

plt.subplot(gs[2,:])
ax = plt.gca()
nonnan = np.where(~(np.isnan(flux_dH_gate2_smith[:,0])))[0]
plt.errorbar(flux_time_gate2_smith[nonnan],flux_dH_gate2_smith[nonnan,0],flux_dH_gate2_smith[nonnan,1],linestyle='dashed',marker='o',color='r',markersize=3)
plt.plot(flux_time_gate2_cresis,flux_dH_gate2_cresis[:,0],'k+',markersize=3)
plt.plot(wv_time_gate2[:,0],wv_dH_gate2[:,0],'ko',markersize=3)
x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
ax.xaxis.set_major_formatter(x_formatter)
labels=[]
for i in range(2000,2017):
  labels.append('Jan \n'+str(i))
plt.xticks(range(2000,2017))
ax.set_xticklabels(labels,fontsize=8,fontname='Arial')

for i in range(0,3):
  plt.subplot(gs[i,:])
  plt.xlim([time1,time2])
  plt.ylim([-100,100])
  plt.yticks(np.arange(-75,100,25),fontsize=8)
  ax = plt.gca()
  xTickPos = np.linspace(np.floor(time1)-0.25,np.floor(time2)-0.25,(np.floor(time2)-np.floor(time1))*2+1)
  ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.8','w'],linewidth=0)

plt.subplots_adjust(hspace=0.05,wspace=0) 
plt.savefig(os.path.join(os.getenv("HOME"),'Bigtmp/'+glacier+'_thinning_timeline.pdf'),FORMAT='PDF')
plt.close()

#####################
# Plot overview map #
#####################

plt.figure(figsize=(4,3))
plt.imshow(zimage,extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],origin='lower',cmap='Greys_r')
ax = plt.gca()
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
path = matplotlib.path.Path(np.column_stack([x_gate1,y_gate1]))
patch = matplotlib.patches.PathPatch(path,facecolor='g',edgecolor='k',alpha=0.4,lw=0)
ax.add_patch(patch)
plt.plot(x_gate1,y_gate1,'k',linewidth=1)
path = matplotlib.path.Path(np.column_stack([x_gate2,y_gate2]))
patch = matplotlib.patches.PathPatch(path,facecolor='r',edgecolor='k',alpha=0.4,lw=0)
ax.add_patch(patch)
plt.plot(x_gate2,y_gate2,'k',linewidth=1)
path = matplotlib.path.Path(np.column_stack([x_gate3,y_gate3]))
patch = matplotlib.patches.PathPatch(path,facecolor='b',edgecolor='k',alpha=0.4,lw=0)
ax.add_patch(patch)
plt.plot(x_gate3,y_gate3,'k',linewidth=1)
plt.xticks([])
plt.yticks([])

plt.savefig(os.path.join(os.getenv("HOME"),'Bigtmp/'+glacier+'_thinning_overview.pdf'),FORMAT='PDF',dpi=600)
plt.close()
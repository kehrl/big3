# Make some Helheim/Kanger plots for PARCA 2017

import datelib, geotifflib, glaclib, zslib, fluxlib, vellib, masklib, icefrontlib, climlib, floatlib, zslib, bedlib
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib
import cubehelix
from mpl_toolkits.basemap import Basemap
from matplotlib.ticker import AutoMinorLocator
from matplotlib.path import Path
from matplotlib.patches import PathPatch

#############
# Load data #
#############

# Flowlines
x_H,y_H,zb_H,dists_H = glaclib.load_flowline('Helheim',shapefilename='flowline_flightline',filt_len=2.0e3,bedsource='cresis')
x_K,y_K,zb_K,dists_K = glaclib.load_flowline('Kanger',shapefilename='flowline_flightline',filt_len=2.0e3,bedsource='cresis')

# Load DEMs	
xdem_H,ydem_H,zdem_H,timedem_H,errordem_H = zslib.dem_grid('Helheim',285000.0,320000.0,-2588000.0,-2566000.0,years='all',verticaldatum='ellipsoid',method='nearest',return_error=True)
xdem_K,ydem_K,zdem_K,timedem_K,errordem_K = zslib.dem_grid('Kanger',449800.0,503000.0,-2302000.0,-2266000.0,years='all',verticaldatum='ellipsoid',method='nearest',return_error=True)

# Ice-front positions
terminus_val_H, terminus_time_H = icefrontlib.distance_along_flowline(x_H,y_H,dists_H,'Helheim',type='icefront',time1=2000.,time2=2016.5)
terminus_val_K, terminus_time_K = icefrontlib.distance_along_flowline(x_K,y_K,dists_K,'Kanger',type='icefront',time1=2000.,time2=2016.5)

# Locations where we want velocities and surface elevations
dists_eul = -1*np.array([2.,5.,10.,15.,20.])
ind_eul_H=[]
for i in range(0,len(dists_eul)):
  ind_eul_H.append( (abs(dists_H - dists_eul[i]*1e3)).argmin() )
ind_eul_K=[]
for i in range(0,len(dists_eul)):
  ind_eul_K.append( (abs(dists_K - dists_eul[i]*1e3)).argmin() )

# Load velocities
vel_val_H,vel_time_H,vel_error_H = vellib.velocity_at_eulpoints(x_H[ind_eul_H],y_H[ind_eul_H],'Helheim',data='all')
vel_val_K,vel_time_K,vel_error_K = vellib.velocity_at_eulpoints(x_K[ind_eul_K],y_K[ind_eul_K],'Kanger',data='all')

# Load elevations
zpt_atm_H,zptstd_atm_H,time_atm_H = zslib.atm_at_pts(x_H[ind_eul_H],y_H[ind_eul_H],'Helheim',years='all',maxdist=200.,verticaldatum='geoid',method='average',cutoff='terminus')
zpt_dem_H,zpterror_dem_H,time_dem_H = zslib.dem_at_pts(x_H[ind_eul_H],y_H[ind_eul_H],'Helheim',years='all',verticaldatum='geoid',cutoff='terminus',method='average',radius=200.)
zpt_atm_K,zptstd_atm_K,time_atm_K = zslib.atm_at_pts(x_K[ind_eul_K],y_K[ind_eul_K],'Kanger',years='all',maxdist=200.,verticaldatum='geoid',method='average',cutoff='terminus')
zpt_dem_K,zpterror_dem_K,time_dem_K = zslib.dem_at_pts(x_K[ind_eul_K],y_K[ind_eul_K],'Kanger',years='all',verticaldatum='geoid',cutoff='terminus',method='average',radius=200.)

# Color options for plotting
coloptions=['r','b','g','limegreen','gold']


###############
# Main figure #
###############

for glacier in ['Kanger','Helheim']:

  vbars = 'none'
  plot_calving = 0

  if glacier == 'Helheim':
    x = x_H; y = y_H; zb = zb_H; dists = dists_H
    vel_val = vel_val_H; vel_time = vel_time_H
    terminus_val = terminus_val_H; terminus_time = terminus_time_H
    time_dem = time_dem_H; zpt_dem = zpt_dem_H; zpterror_dem = zpterror_dem_H
    time_atm = time_atm_H; zpt_atm = zpt_atm_H
    ind_eul = ind_eul_H
  elif glacier == 'Kanger':  
    x = x_K; y = y_K; zb = zb_K; dists = dists_K
    vel_val = vel_val_K; vel_time = vel_time_K
    terminus_val = terminus_val_K; terminus_time = terminus_time_K
    time_dem = time_dem_K; zpt_dem = zpt_dem_K; zpterror_dem = zpterror_dem_K
    time_atm = time_atm_K; zpt_atm = zpt_atm_K
    ind_eul = ind_eul_K
    
  plt.figure(figsize=(8,5))
  time1 = 2000.0; time2 = 2016.75; 
  years = np.arange(np.floor(time1),np.ceil(time2)+1)
 
  gs = matplotlib.gridspec.GridSpec(6,1)
  matplotlib.rc('font',family='Arial')

  # Plot terminus
  plt.subplot(gs[0, :])
  ax = plt.gca()

  for i in range(2000,2017):
    plt.plot([i,i],[-5,120],color='0.5',lw=0.5)    
  if vbars == 'runoff':
    for i in range(0,len(year_runoff)):
      path = matplotlib.path.Path([[day1_runoff[i],-5],[day1_runoff[i],5],[day2_runoff[i],5],[day2_runoff[i],-5]])
      patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
      ax.add_patch(patch)
  elif vbars == 'seasonal':
    for i in range(0,len(years)):
      path = matplotlib.path.Path([[years[i]+0.25,-6],[years[i]+0.25,12],[years[i]+0.75,12],[years[i]+0.75,-6]])
      patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
      ax.add_patch(patch)  

  nonnan = np.where(~(np.isnan(terminus_val)))[0]
  plt.plot(terminus_time[nonnan],terminus_val[nonnan]/1e3,'ko',linewidth=1,markersize=3.5)
  x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
  ax.xaxis.set_major_formatter(x_formatter)
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  plt.xticks(range(2000,2017),fontsize=12,fontname="Arial")
  #ax.xaxis.tick_top()
  labels=[]
  for i in range(2000,2017):
    labels.append('Jan \n'+str(i))
  plt.xticks(np.arange(2000,2017,2))
  ax.set_xticklabels([])
  #ax.set_xticklabels(labels,fontsize=11,fontname='Arial')
  plt.xlim([time1,time2])
  plt.yticks(np.arange(-6,8,2),fontsize=11,fontname="Arial")
  plt.ylabel('Terminus \n (km)',fontsize=12,fontname="Arial")
  ax.tick_params('both', length=6, width=1.25, which='major')
  ax.tick_params('x', length=6, width=1.25, which='minor')
  plt.ylim([-5,6])
  #plt.text(2008.25,-2.3,'a',fontsize=12,fontname='Arial',fontweight='bold')

  # Plot velocities
  ax = plt.subplot(gs[1:4, :]) 
  coloptions=['r','b','g','limegreen','gold']
  markoptions=['o','o','o','o','o','o']

  for i in range(2000,2017):
    plt.plot([i,i],[-5,120],color='0.5',lw=0.5)
  x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
  ax.xaxis.set_major_formatter(x_formatter)
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  ax.yaxis.set_minor_locator(AutoMinorLocator(2))
  matplotlib.rc('font',family="Arial",)
  ax.tick_params('both', length=6, width=1.25, which='major')
  ax.tick_params('both', length=6, width=1.25, which='minor')
  for i in range(0,len(dists_eul)):
    nonnan = np.where(~(np.isnan(vel_val[:,i])))[0]
    plt.plot(vel_time[nonnan],(vel_val[nonnan,i])/1e3,markoptions[i],color=coloptions[i],label=glacier[0]+'{0:02d}'.format(int(abs(dists_eul[i]))),markersize=5)
  plt.legend(loc=2,borderpad=0.3,fontsize=10,numpoints=1,handlelength=0.2,handletextpad=0.5,labelspacing=0.1,ncol=2,columnspacing=0.8)
  plt.yticks(np.arange(2,14,2),fontsize=12,fontname="Arial")
  plt.ylabel('Glacier speed \n (km/yr)',fontsize=12,fontname="Arial")
  labels=[]
  plt.xticks(np.arange(2000,2017,2),fontsize=12,fontname="Arial")
  ax.set_xticklabels([])
  if vbars == 'runoff':
    for i in range(0,len(year_runoff)):
      path = matplotlib.path.Path([[day1_runoff[i],-3],[day1_runoff[i],12],[day2_runoff[i],12],[day2_runoff[i],-3]])
      patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
      ax.add_patch(patch)  
  elif vbars == 'seasonal':
    for i in range(0,len(years)):
      path = matplotlib.path.Path([[years[i]+0.25,-3],[years[i]+0.25,12],[years[i]+0.75,12],[years[i]+0.75,-3]])
      patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
      ax.add_patch(patch)  
  plt.xlim([time1,time2])
  if glacier == 'Helheim':
    plt.ylim([3,11])
  elif glacier == 'Kanger':
    plt.ylim([1.5,13])
  #plt.text(2008.25,3.8,'b',fontsize=12,fontname='Arial',fontweight='bold')

  # Plot surface elevations
  plt.subplot(gs[4:, :])
  ax = plt.gca()
  for i in range(2000,2017):
    plt.plot([i,i],[-30,400],color='0.5',lw=0.5)
  # Set up points for legend
  plt.errorbar(0,0,capsize=1,yerr=0.5,fmt='o',color='k',markersize=3.5,label='DEM')
  plt.plot(0,0,'+',color='k',markersize=5,markeredgewidth=2,label='ATM')
  if glacier == 'Helheim':
    plt.plot(0,0,'rs',label='H02',markersize=3.5)
    plt.plot(0,0,'bs',label='H05',markersize=3.5)
  elif glacier == 'Kanger':
    plt.plot(0,0,'bs',label='K05',markersize=3.5)
    plt.plot(0,0,'gs',label='H10',markersize=3.5)
  if vbars == 'runoff':
    for i in range(0,len(year_runoff)):
      path = matplotlib.path.Path([[day1_runoff[i],-50],[day1_runoff[i],240],[day2_runoff[i],240],[day2_runoff[i],-50]])
      patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
      ax.add_patch(patch)
  elif vbars == 'seasonal':
    for i in range(0,len(years)):
      path = matplotlib.path.Path([[years[i]+0.25,-50],[years[i]+0.25,240],[years[i]+0.75,240],[years[i]+0.75,-50]])
      patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
      ax.add_patch(patch)  
  if glacier == 'Helheim':
    plt.fill_between([time1,time2],np.ones(2)*floatlib.height(zb[ind_eul[0]]-50)-floatlib.height(zb[ind_eul[0]]),np.ones(2)*floatlib.height(zb[ind_eul[0]]+50)-floatlib.height(zb[ind_eul[0]]),alpha=0.1,facecolor='r',edgecolor='r',antialiased=True,zorder=2)
    plt.plot([time1,time2],[0,0],'r:',linewidth=1.0)
    nonnan = np.where(~(np.isnan(zpt_atm[:,0])))[0]
    plt.plot(time_atm[nonnan],zpt_atm[nonnan,0]-floatlib.height(zb[ind_eul[0]]),'+',color=coloptions[0],markersize=6,markeredgewidth=2)    
    nonnan = np.where(~(np.isnan(zpt_dem[:,0])))[0]
    plt.plot(time_dem[nonnan],zpt_dem[nonnan,0]-floatlib.height(zb[ind_eul[0]]),'o',color=coloptions[0],markersize=5)
  elif glacier == 'Kanger':
  	plt.fill_between([time1,time2],np.ones(2)*floatlib.height(zb[ind_eul[1]]-50)-floatlib.height(zb[ind_eul[1]]),np.ones(2)*floatlib.height(zb[ind_eul[1]]+50)-floatlib.height(zb[ind_eul[1]]),alpha=0.1,facecolor='b',edgecolor='b',antialiased=True,zorder=2)
  	plt.plot([time1,time2],[0,0],'b:',linewidth=1.0)
  	nonnan = np.where(~(np.isnan(zpt_dem[:,1])))[0]  	
  	plt.plot(time_dem[nonnan],zpt_dem[nonnan,1]-floatlib.height(zb[ind_eul[1]]),'o',color=coloptions[1],markersize=5)
  	nonnan = np.where(~(np.isnan(zpt_atm[:,1])))[0]  
  	plt.plot(time_atm[nonnan],zpt_atm[nonnan,1]-floatlib.height(zb[ind_eul[1]]),'+',color=coloptions[1],markersize=6,markeredgewidth=2)
  x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
  ax.xaxis.set_major_formatter(x_formatter)
  plt.xticks(np.arange(2000,2017,2))
  labels=[]
  for i in np.arange(2000,2017,2):
    labels.append('Jan \n'+str(i))
  x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
  ax.xaxis.set_major_formatter(x_formatter)
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  ax.tick_params('both', length=6, width=1.25, which='major')
  ax.tick_params('x', length=6, width=1.25, which='minor')
  ax.set_xticklabels(labels,fontsize=12,fontname='Arial')
  plt.xlim([time1,time2])
  if glacier == 'Helheim':
    plt.ylabel('  Height above flotation \n at H02 (m)',fontsize=12,fontname="Arial")
    plt.yticks(np.arange(-20,200,20),fontsize=12,fontname="Arial")
    plt.ylim([-5,120])
  elif glacier == 'Kanger':
    plt.ylabel('  Height above flotation \n at K05 (m)',fontsize=12,fontname="Arial")
    plt.yticks(np.arange(-20,200,50),fontsize=12,fontname="Arial")
    plt.ylim([-25,205])
  ind = [0,3,1,2]
  handles, labels = ax.get_legend_handles_labels()
  labels = [labels[i] for i in ind]
  handles = [handles[i] for i in ind]
  plt.legend(handles,labels,loc=3,borderpad=0.3,handleheight=0.1,fontsize=10,numpoints=1,handlelength=0.4,labelspacing=0.05,ncol=2,columnspacing=0.7,handletextpad=0.5)
  #plt.text(2008.25,-2,'c',fontsize=12,fontname='Arial',fontweight='bold')

  ax2 = ax.twinx()
  if glacier == 'Helheim':
    nonnan = np.where(~(np.isnan(zpt_atm[:,1])))[0]
    ax2.plot(time_atm[nonnan],zpt_atm[nonnan,1]-floatlib.height(zb[ind_eul[1]]),'+',color=coloptions[1],markersize=6,markeredgewidth=2,label='ATM')
    nonnan = np.where(~(np.isnan(zpt_dem[:,1])))[0]
    ax2.plot(time_dem[nonnan],zpt_dem[nonnan,1]-floatlib.height(zb[ind_eul[1]]),'o',color=coloptions[1],markersize=5)
  elif glacier == 'Kanger':
    nonnan = np.where(~(np.isnan(zpt_dem[:,2])))[0]
    ax2.plot(time_dem[nonnan],zpt_dem[nonnan,2]-floatlib.height(zb[ind_eul[2]]),'o',color=coloptions[2],markersize=5,label='WV')
    nonnan = np.where(~(np.isnan(zpt_atm[:,2])))[0]
    ax2.plot(time_atm[nonnan],zpt_atm[nonnan,2]-floatlib.height(zb[ind_eul[2]]),'+',color=coloptions[2],markersize=6,markeredgewidth=2,label='ATM')
  ax2.set_xlim([time1,time2])
  ax2.tick_params('both', length=6, width=1.25, which='major')
  ax2.tick_params('x', length=3, width=1, which='minor')

  if glacier == 'Helheim':
    ax2.set_ylabel('    Height above flotation at H05 (m)',fontsize=12,fontname='Arial')
    ax2.set_yticks(np.arange(60,180,20))
    ax2.set_ylim([55,180])
    ax2.set_yticklabels([60,80,100,120,140,160],fontsize=12,fontname='arial')
  elif glacier == 'Kanger':
    ax2.set_yticks(np.arange(150,400,50))
    ax2.set_ylim([145,375])
    ax2.set_yticklabels(np.arange(150,400,50),fontsize=12,fontname='arial')
    ax2.set_ylabel('    Height above flotation at K10 (m)',fontsize=12,fontname='Arial')

  plt.tight_layout()
  plt.subplots_adjust(hspace=0.03,wspace=0,top=0.98,right=0.93,left=0.095,bottom=0.09) 
  plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/PARCA_"+glacier+"_all.pdf"),FORMAT='PDF',dpi=600)
  plt.close()

del k,i, path, patch, ax, labels, handles, time1, time2


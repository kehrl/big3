# Make some plots for my AGU 2016 presentation

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

# Flux gates
xgate_H,ygate_H = fluxlib.fluxbox_geometry('Helheim',"fluxgate1")

# Load DEMs	
xdem_H,ydem_H,zdem_H,timedem_H,errordem_H = zslib.dem_grid('Helheim',285000.0,320000.0,-2588000.0,-2566000.0,years='all',verticaldatum='ellipsoid',method='nearest',return_error=True)

# And flotation condition for Helheim
xmin = 283000.
xmax = 313000.
ymin = -2587000.
ymax = -2552000.
xwv_H,ywv_H,zwv_H,timewv_H = zslib.dem_grid('Helheim',xmin,xmax,ymin,ymax,years='all',resolution=32,verticaldatum='geoid')
xf_H,yf_H,zabovefloat_H = floatlib.extent(xwv_H,ywv_H,zwv_H,timewv_H,'Helheim',rho_i=917.0,rho_sw=1025.0,bedsource='cresis',verticaldatum='geoid')

# Thinning rates
time_H,dH_H,Q_H,hbar_H,ubar_H,smb_H,width_H,area_H = fluxlib.fluxgate_thinning('Helheim','fluxgate1','cresis',10.,timing='velocity')
dem_time_H,dem_dH_H = fluxlib.dem_thinning('Helheim',xdem_H,ydem_H,zdem_H,timedem_H,errordem_H,"fluxgate1")

# Ice-front positions
terminus_val_H, terminus_time_H = icefrontlib.distance_along_flowline(x_H,y_H,dists_H,'Helheim',type='icefront',time1=2000.,time2=2016.5)

# Locations where we want velocities and surface elevations
dists_eul = -1*np.array([2.,5.,10.,15.,20.])
ind_eul_H=[]
for i in range(0,len(dists_eul)):
  ind_eul_H.append( (abs(dists_H - dists_eul[i]*1e3)).argmin() )

# Load velocities
vel_val_H,vel_time_H,vel_error_H = vellib.velocity_at_eulpoints(x_H[ind_eul_H],y_H[ind_eul_H],'Helheim',data='all')

# Load elevations
zpt_atm_H,zptstd_atm_H,time_atm_H = zslib.atm_at_pts(x_H[ind_eul_H],y_H[ind_eul_H],'Helheim',years='all',maxdist=200.,verticaldatum='geoid',method='average',cutoff='terminus')
zpt_dem_H,zpterror_dem_H,time_dem_H = zslib.dem_at_pts(x_H[ind_eul_H],y_H[ind_eul_H],'Helheim',years='all',verticaldatum='geoid',cutoff='terminus',method='average',radius=200.)

# Filter length for 
filt_len = 11.

# Load smb from RACMO
xrac_H,yrac_H,smbrac_H,timerac_H = climlib.racmo_at_pts(np.mean(xgate_H),np.mean(ygate_H),'smb',filt_len=filt_len)

xrac_H,yrac_H,runrac_H,timerac_H = climlib.racmo_at_pts(np.mean(xgate_H),np.mean(ygate_H),'runoff',filt_len=filt_len)

ind = np.argmin(abs(dists_H-np.nanmax(terminus_val_H)))
xsst_H,ysst_H,sst_H,timesst_H = climlib.SIF_at_pts(np.mean(x_H[ind]),np.mean(y_H[ind]),filt_len=filt_len,variable='sst')

ind = np.argmin(abs(dists_H-np.nanmax(terminus_val_H)))
xsif_H,ysif_H,sif_H,timesif_H = climlib.SIF_at_pts(np.mean(x_H[ind]),np.mean(y_H[ind]),filt_len=filt_len)

yH_runoff,day1H_runoff,day2H_runoff,meltlengthH_runoff,totalH_runoff = climlib.seasonlength(timerac_H,runrac_H,'runoff')
yH_smb,day1H_smb,day2_smb,meltlengthH_smb,totalH_smb = climlib.seasonlength(timerac_H,smbrac_H,'smb')

# Color options for plotting
coloptions=['r','b','g','limegreen','gold']

############################
# Plot behavior since 2000 #
############################

# Helheim 

plt.figure(figsize=(4.5,5))
gs = matplotlib.gridspec.GridSpec(3,1)

plt.subplot(gs[0])
ax = plt.gca()
for child in ax.get_children():
  if isinstance(child, matplotlib.spines.Spine):
    child.set_color('k')
ax.tick_params(axis='both',colors='k',which='both')
ax.arrow(2007,3.3,0,0.6,head_width=0.1,fc='k',ec='k',lw=1.5,head_length=0.1)
ax.text(2007.3,3.3,'Advance',fontsize=16)
ax.arrow(2007,-2.7,0,-0.6,head_width=0.1,fc='k',ec='k',lw=1.5,head_length=0.1)
ax.text(2007.3,-3.3,'Retreat',fontsize=16)
ax.tick_params(axis='both',colors='k',which='both')
plt.yticks(np.arange(-4,6,2),fontsize=12)
plt.ylim([-4.5,5])
plt.xlim([2000,2016])
plt.xticks(np.arange(2000,2020,4))
plt.plot(terminus_time_H,terminus_val_H/1e3,'ro',markersize=3.5)
plt.ylabel('Terminus \n (km)',color='k',fontsize=14)
x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
ax.xaxis.set_major_formatter(x_formatter)
ax.set_xticklabels([])

plt.subplot(gs[1])
ax=plt.gca()
for child in ax.get_children():
  if isinstance(child, matplotlib.spines.Spine):
    child.set_color('k')
ax.tick_params(axis='both',colors='k',which='both')
plt.xlim([2000,2016])
plt.plot(vel_time_H,vel_val_H[:,1]/1e3,'ro',markersize=3.5)
plt.ylabel('Glacier speed \n (km/yr)',color='k',fontsize=14)
plt.yticks(np.arange(5,10,1),fontsize=12)
plt.ylim([4.5,8.5])
plt.xticks(np.arange(2000,2020,4))
x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
ax.xaxis.set_major_formatter(x_formatter)
ax.set_xticklabels([])

plt.subplot(gs[2])
ax=plt.gca()
for child in ax.get_children():
  if isinstance(child, matplotlib.spines.Spine):
    child.set_color('k')
ax.tick_params(axis='both',colors='k',which='both')
plt.xlim([2000,2016])
plt.plot(time_atm_H,zpt_atm_H[:,1],'ro',markersize=3.5)
plt.plot(time_dem_H,zpt_dem_H[:,1],'ro',markersize=3.5)
plt.ylabel('Surface elevation \n (m asl)',color='k',fontsize=14)
#plt.ylim([4.5,12.5])
plt.ylim([120,300])
plt.yticks(np.arange(120,300,40),fontsize=12)
plt.xticks(np.arange(2000,2020,4))
x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
ax.xaxis.set_major_formatter(x_formatter)

plt.tight_layout()
plt.subplots_adjust(wspace=0.04,hspace=0.04,top=0.98,right=0.95,left=0.2,bottom=0.05)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/AGU_Helheim_2000.pdf"),format="PDF")
plt.close()


#######################
# Plot bed elevations #
#######################

for glacier in ['Helheim']:
  if glacier == 'Helheim':
    atm_data = zslib.atm_along_flowline(x_H,y_H,glacier,years='20010521',verticaldatum='geoid')
    zsdem,junk = zslib.dem_along_flowline(x_H,y_H,glacier,years='20130804',verticaldatum='geoid')
    x = x_H; y = y_H; dists = dists_H
  elif glacier == 'Kanger':
    atm_data = zslib.atm_along_flowline(x_K,y_K,glacier,years='20010520',verticaldatum='geoid')
    zsdem,junk = zslib.dem_along_flowline(x_K,y_K,glacier,years='20130714',verticaldatum='geoid')
    x = x_K; y = y_K; dists = dists_K
  # Get radar thicknesses close to flightline
  cresis = bedlib.cresis('all',glacier)
  if glacier == 'Helheim':
    cresis2001 = bedlib.cresis('2001',glacier)
    cresis = np.row_stack([cresis,cresis2001])
  cutoff = 200.
  dcresis = []
  zcresis = []
  tcresis = []
  for i in range(0,len(cresis[:,0])):
    mindist = np.min(np.sqrt((cresis[i,0]-x)**2+(cresis[i,1]-(y))**2))
    if mindist < cutoff:
      minind = np.argmin(np.sqrt((cresis[i,0]-x)**2+(cresis[i,1]-(y))**2))
      dcresis.append(dists[minind])
      zcresis.append(cresis[i,2])
      tcresis.append(cresis[i,4])
  dcresis = np.array(dcresis)
  zcresis = np.array(zcresis)
  

  fig = plt.figure(figsize=(4,2.7))
  gs = matplotlib.gridspec.GridSpec(2,1)
  matplotlib.rc('font',family='Arial')
 
  plt.subplot(gs[0])
  ax = plt.gca()
  ind = np.argmin(abs(dists--5e3))
  if glacier == 'Helheim':
    plt.plot(dists/1e3,atm_data['20010521'][:,2],'b',label='21 May 2001',lw=1.2)
    ind = np.where((terminus_time_H > 2008.) & (terminus_time_H < 2016.))
    plt.plot([np.min(terminus_val_H[ind])/1e3,np.min(terminus_val_H[ind])/1e3],[-1500,1500],'0.5',lw=1.5)
    plt.plot([np.max(terminus_val_H[ind])/1e3,np.max(terminus_val_H[ind])/1e3],[-1500,1500],'0.5',lw=1.5)
    plt.plot(dists_H/1e3,zsdem[0,:].T,color='r',linewidth=1.2,label='04 Aug 2013')
    plt.plot(dists_H/1e3,floatlib.height(zb_H),'k',linewidth=1.2,label='Flotation',dashes=[2,2,2,2])
  plt.xticks(np.arange(-30,10,5),fontsize=10)
  ax.set_xticklabels([])
  plt.yticks(np.arange(-1000,1000,250),fontsize=10)
  plt.xlim([-21,6])
  if glacier == 'Helheim':
    plt.ylim([0,670])
  elif glacier == 'Kanger':
    plt.ylim([0,800])
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  ax.yaxis.set_minor_locator(AutoMinorLocator(2))
  ax.tick_params('both', length=6, width=1.25, which='major')
  ax.tick_params('both', length=3, width=1, which='minor')
  plt.ylabel('Elevation (m asl)                                    ',fontsize=10)
  plt.legend(loc=1,fontsize=10,borderpad=0.2,numpoints=1,handlelength=0.6,labelspacing=0.1,handletextpad=0.3,markerscale=2)

  plt.subplot(gs[-1])
  ax = plt.gca()
  plt.plot(dcresis/1e3,zcresis,'.',c='0.7',markersize=2.5,label='Measured')
  plt.yticks(np.arange(-1250,-250,250),fontsize=10)
  if glacier == 'Helheim':
    ind = np.where((terminus_time_H > 2008.) & (terminus_time_H < 2016.))[0]
    plt.plot([np.min(terminus_val_H[ind])/1e3,np.min(terminus_val_H[ind])/1e3],[-1500,1500],'0.5',lw=1.5)
    plt.plot([np.max(terminus_val_H[ind])/1e3,np.max(terminus_val_H[ind])/1e3],[-1500,1500],'0.5',lw=1.5)
    plt.plot(dists_H/1e3,zb_H,color='k',linewidth=1.2,label='Smoothed')
    plt.ylim([-1150,-300])
    plt.text(np.min(terminus_val_H[ind])/1e3+0.2,-1000,'2008-\n2016',fontsize=10,fontname='Arial')
  elif glacier == 'Kanger':
    ind = np.where((terminus_time_K > 2008.) & (terminus_time_K < 2016.))[0]
    plt.plot([np.min(terminus_val_K[ind])/1e3,np.min(terminus_val_K[ind])/1e3],[-1500,1500],'0.5',lw=1.5)
    plt.plot([np.max(terminus_val_K[ind])/1e3,np.max(terminus_val_K[ind])/1e3],[-1500,1500],'0.5',lw=1.5)
    plt.text(np.min(terminus_val_K[ind])/1e3+0.5,-1150,'2008-2016',fontsize=10,fontname='Arial')
    ind = np.argmin(abs(dists_K--5e3))
    #plt.plot(dists/1e3,zb,':',color='k',linewidth=1.5)
    plt.plot(dists_K[0:ind]/1e3,zb_K[0:ind],color='k',linewidth=1.5,label='Smoothed')
    plt.text(-3.2,-950,'??',fontsize=10,fontname='arial')
    plt.text(1.5,-500,'??',fontsize=10,fontname='arial')
    plt.ylim([-1300,-300])
  plt.xlabel('Distance along flowline (km)',fontsize=10)
  plt.xticks(np.arange(-30,10,5),fontsize=10)
  plt.xlim([-21,6])
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  ax.yaxis.set_minor_locator(AutoMinorLocator(2))
  ax.tick_params('both', length=6, width=1.25, which='major')
  ax.tick_params('both', length=3, width=1, which='minor')
  #plt.legend(loc=4,borderpad=0.2,fontsize=10,numpoints=1,handletextpad=0.3,handlelength=0.5,labelspacing=0.1,markerscale=2)
  plt.tight_layout()
  plt.subplots_adjust(wspace=0.04,hspace=0.04,top=0.98,right=0.98,left=0.135,bottom=0.15)
  plt.savefig(os.path.join(os.getenv("HOME"),'Bigtmp/AGU_'+glacier+'_zs_flowline.pdf'),FORMAT='PDF',dpi=600)
  plt.close()

del zsdem,atm_data,x,y,dists,ax,cresis,zcresis,dcresis,fig,gs,junk,ind,tcresis,cutoff


###############
# Main figure #
###############

for glacier in ['Helheim']:

  vbars = 'seasonal'
  plot_calving = 0

  if glacier == 'Helheim':
    x = x_H; y = y_H; zb = zb_H; dists = dists_H
    vel_val = vel_val_H; vel_time = vel_time_H
    terminus_val = terminus_val_H; terminus_time = terminus_time_H
    time_dem = time_dem_H; zpt_dem = zpt_dem_H; zpterror_dem = zpterror_dem_H
    time_atm = time_atm_H; zpt_atm = zpt_atm_H
    day1_runoff = day1H_runoff; day2_runoff = day2H_runoff; year_runoff = yH_runoff
    ind_eul = ind_eul_H
    smb = smb_H; dem_dH = dem_dH_H; dem_time = dem_time_H;
    time = time_H; dH = dH_H
    timesif = timesif_H; sif = sif_H
    timerac = timerac_H; runrac = runrac_H
    total_runoff = totalH_runoff

    demtimes = ['20130508','20130804','20131031','20140127']
      
  calvingstyle = icefrontlib.calving(glacier)

  for k in range(0,1):
    
    plt.figure(figsize=(7.45,6))
    time1 = 2008.0; time2 = 2016.75; 
    years = np.arange(np.floor(time1),np.ceil(time2)+1)
 
    gs = matplotlib.gridspec.GridSpec(7,1)
    matplotlib.rc('font',family='Arial')

    # Plot terminus
    plt.subplot(gs[0, :])
    ax = plt.gca()
    
    plt.ylim([-6,6])

    if vbars == 'runoff':
      for i in range(0,len(year_runoff)):
        path = matplotlib.path.Path([[day1_runoff[i],-5],[day1_runoff[i],5],[day2_runoff[i],5],[day2_runoff[i],-5]])
        patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
        ax.add_patch(patch)
    elif vbars == 'seasonal':
      for i in range(0,len(years)):
        path = matplotlib.path.Path([[years[i]+0.25,-3],[years[i]+0.25,12],[years[i]+0.75,12],[years[i]+0.75,-3]])
        patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
        ax.add_patch(patch)  
    if (plot_calving):
      ind = np.where((calvingstyle[:,1] == 'Domino'))[0]
      for i in ind:
        if i == ind[1]:
          plt.bar(float(calvingstyle[i,0])-0.025,10.0,width=0.035,color=[1,0.6,0.6],edgecolor='none',label='Non-tabular',bottom=-5)
        elif i != 0:
          plt.bar(float(calvingstyle[i,0])-0.025,10.0,width=0.035,color=[1,0.6,0.6],edgecolor='none',bottom=-5)
      ind = np.where((calvingstyle[:,1] == 'Mixed'))[0]
      for i in ind:
        if i == ind[1]:
          plt.bar(float(calvingstyle[i,0])-0.025,10.0,width=0.035,color=[0.4,0.8,0.6],edgecolor='none',label='Mixed',bottom=-5)
        elif i != 0:
          plt.bar(float(calvingstyle[i,0])-0.025,10.0,width=0.035,color=[0.4,0.8,0.6],edgecolor='none',bottom=-5)

      ind = np.where(calvingstyle[:,1] == 'Tabular')[0]
      for i in ind:
        if i == ind[1]:
          plt.bar(float(calvingstyle[i,0])-0.025,10.0,width=0.035,color=[0.6,0.6,1],edgecolor='none',label='Tabular',bottom=-5)
        elif i !=0:
          plt.bar(float(calvingstyle[i,0])-0.025,10.0,width=0.035,color=[0.6,0.6,1],edgecolor='none',bottom=-5)
      plt.legend(loc=2,fontsize=10,numpoints=1,handlelength=0.3,labelspacing=0.05,ncol=3,columnspacing=0.7,handletextpad=0.2,borderpad=0.25)
    nonnan = np.where(~(np.isnan(terminus_val)))[0]
    plt.plot(terminus_time[nonnan],terminus_val[nonnan]/1e3,'ko',linewidth=1,markersize=3.5)
    x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax.xaxis.set_major_formatter(x_formatter)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    plt.xticks(range(2000,2017),fontsize=12,fontname="Arial")
    ax.xaxis.tick_top()
    labels=[]
    for i in range(2000,2017):
      labels.append('Jan \n'+str(i))
    plt.xticks(range(2000,2017))
    ax.set_xticklabels(labels,fontsize=11,fontname='Arial')
    plt.xlim([time1,time2])
    plt.yticks(np.arange(-6,8,2),fontsize=11,fontname="Arial")
    plt.ylabel('Terminus (km)',fontsize=12,fontname="Arial")
    ax.tick_params('both', length=6, width=1.25, which='major')
    ax.tick_params('both', length=3, width=1, which='minor')
    plt.ylim([-3,3])
    plt.text(2008.25,-2.3,'a',fontsize=12,fontname='Arial',fontweight='bold')

    # Plot velocities
    ax = plt.subplot(gs[1:3, :]) 
    coloptions=['r','b','g','limegreen','gold']
    markoptions=['o','o','o','o','o','o']


    for i in range(0,len(dists_eul)):
      nonnan = np.where(~(np.isnan(vel_val[:,i])))[0]
      plt.plot(vel_time[nonnan],(vel_val[nonnan,i])/1e3,markoptions[i],color=coloptions[i],label=glacier[0]+'{0:02d}'.format(int(abs(dists_eul[i]))),markersize=4)
    plt.legend(loc=2,borderpad=0.3,fontsize=10,numpoints=1,handlelength=0.2,handletextpad=0.5,labelspacing=0.1,ncol=5,columnspacing=0.8)
    plt.yticks(range(2,12),fontsize=12,fontname="Arial")
    plt.ylabel('Glacier speed \n (km/yr)',fontsize=12,fontname="Arial")
    x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax.xaxis.set_major_formatter(x_formatter)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    matplotlib.rc('font',family="Arial",)
    ax.tick_params('both', length=6, width=1.25, which='major')
    ax.tick_params('both', length=3, width=1, which='minor')
    labels=[]
    plt.xticks(range(2000,2017),fontsize=12,fontname="Arial")
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
    plt.ylim([3.5,11])
    plt.text(2008.25,3.8,'b',fontsize=12,fontname='Arial',fontweight='bold')

    # Plot surface elevations
    plt.subplot(gs[-4:-2, :])
    ax = plt.gca()
    # Set up points for legend
    plt.errorbar(0,0,capsize=1,yerr=0.5,fmt='o',color='k',markersize=3.5,label='DEM')
    plt.plot(0,0,'+',color='k',markersize=5,markeredgewidth=2,label='ATM')

    plt.plot(0,0,'rs',label='H02',markersize=3.5)
    plt.plot(0,0,'bs',label='H05',markersize=3.5)
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
    #plt.fill_between([time1,time2],np.ones(2)*floatlib.height(zb[ind_eul[0]]-50)-floatlib.height(zb[ind_eul[0]]),np.ones(2)*floatlib.height(zb[ind_eul[0]]+50)-floatlib.height(zb[ind_eul[0]]),alpha=0.1,facecolor='r',edgecolor='r',antialiased=True,zorder=2)
    #plt.plot([time1,time2],[0,0],'r:',linewidth=1.0)
    plt.plot(time_atm,zpt_atm[:,0]-floatlib.height(zb[ind_eul[0]]),'+',color=coloptions[0],markersize=6,markeredgewidth=2)    
    plt.errorbar(time_dem,zpt_dem[:,0]-floatlib.height(zb[ind_eul[0]]),capsize=1,yerr=zpterror_dem,fmt='o',color=coloptions[0],markersize=4)
    x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax.xaxis.set_major_formatter(x_formatter)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.set_xticklabels([])
    plt.xticks(range(2000,2016),fontsize=12,fontname="Arial")
    plt.xlim([time1,time2])
    ax.tick_params('both', length=6, width=1.25, which='major')
    ax.tick_params('both', length=3, width=1, which='minor')
    plt.ylabel('  Height above flotation \n at H02 (m)',fontsize=12,fontname="Arial")
    plt.yticks(np.arange(-10,50,10),fontsize=12,fontname="Arial")
    plt.ylim([-5,50])
    ind = [0,3,1,2]
    handles, labels = ax.get_legend_handles_labels()
    labels = [labels[i] for i in ind]
    handles = [handles[i] for i in ind]
    plt.legend(handles,labels,loc=1,borderpad=0.3,handleheight=0.1,fontsize=10,numpoints=1,handlelength=0.4,labelspacing=0.05,ncol=4,columnspacing=0.7,handletextpad=0.5)
    plt.text(2008.25,-2,'c',fontsize=12,fontname='Arial',fontweight='bold')

    ax2 = ax.twinx()
    ax2.plot(time_atm,zpt_atm[:,1]-floatlib.height(zb[ind_eul[1]]),'+',color=coloptions[1],markersize=6,markeredgewidth=2,label='ATM')
    ax2.errorbar(time_dem,zpt_dem[:,1]-floatlib.height(zb[ind_eul[1]]),capsize=1,yerr=zpterror_dem,fmt='o',color=coloptions[1],markersize=4,label='WV')
    ax2.plot(time_dem,zpt_dem[:,1]-floatlib.height(zb[ind_eul[1]]),'o',color=coloptions[1],markersize=4)
    ax2.set_xlim([time1,time2])
    ax2.tick_params('both', length=6, width=1.25, which='major')
    ax2.tick_params('both', length=3, width=1, which='minor')
    ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.set_yticks(np.arange(60,120,10))
    ax2.set_ylim([55,110])
    ax2.set_yticklabels([60,70,80,90,100],fontsize=12,fontname='arial')
    ax2.set_ylabel('    Height above flotation at H05 (m)',fontsize=12,fontname='Arial')

    plt.subplot(gs[-2:,:])
    ax = plt.gca()
    plt.plot([time1,time2],[0,0],'k')
    ind = np.where(np.isnan(smb))
    ind2 = np.where(np.floor(time)==2015.)
    for i in ind:
      smb[i] = np.interp(time[i]-2016,time[ind2]-2015.,smb[ind2])
    plt.plot(timerac_H,smbrac_H/900.*100,'0.4',lw=2,label='SMB')
    plt.errorbar(time,(dH[:,0]+smb)/365.25*100,yerr=3,fmt='rs',markersize=4,capsize=1,lw=1.5,label='Flux-gate')
    plt.errorbar(dem_time[:,0],dem_dH[:,0]/365.25*100,xerr=dem_time[:,1],yerr=dem_dH[:,1]/365.25*100,fmt='bo',markersize=4,capsize=1,lw=0.5,label='DEM')
    plt.legend(loc=2,numpoints=1,ncol=3,handletextpad=0.5,fontsize=10,columnspacing=0.7,markerscale=1,handlelength=1.0,borderpad=0.3)
    plt.ylabel(r'dh/dt (cm/d)',fontname='Arial',fontsize=12)
    plt.yticks(np.arange(-40,30,10),fontsize=12,fontname='Arial')
    plt.ylim([-45,25])
    if vbars == 'runoff':
      for i in range(0,len(yH_runoff)):
        path = matplotlib.path.Path([[day1_runoff[i],-165],[day1_runoff[i],160],[day2_runoff[i],160],[day2_runoff[i],-165]])
        patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
        ax.add_patch(patch) 
    elif vbars == 'seasonal':
      for i in range(0,len(years)):
        path = matplotlib.path.Path([[years[i]+0.25,-165],[years[i]+0.25,160],[years[i]+0.75,160],[years[i]+0.75,-165]])
        patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
        ax.add_patch(patch)  
    plt.xticks(range(2000,2017))
    labels=[]
    for i in range(2000,2017):
      labels.append('Jan \n'+str(i))
    x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax.xaxis.set_major_formatter(x_formatter)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params('both', length=6, width=1.25, which='major')
    ax.tick_params('both', length=3, width=1, which='minor')
    ax.set_xticklabels(labels,fontsize=12,fontname='Arial')
    plt.xlim([time1,time2])
    plt.text(2008.25,-40,'d',fontsize=12,fontname='Arial',fontweight='bold')

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.03,wspace=0,top=0.93,right=0.91,left=0.1,bottom=0.06) 
    plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/AGU_"+glacier+"_all.pdf"),FORMAT='PDF',dpi=600)
    plt.close()

del k,i, path, patch, ax, labels, handles, time1, time2

########################
# Plot individual DEMs #
########################

for k in range(0,1):
    
  ind = np.where(np.isnan(zb))[0]
  ind2 = np.where(~(np.isnan(zb)))[0]
  zb[ind] = np.interp(dists[ind],dists[ind2],zb[ind2])
    
  # Get radar thicknesses close to flightline
  cresis = bedlib.cresis('all',glacier)
  if glacier == 'Helheim':
    cresis2001 = bedlib.cresis('2001',glacier)
    cresis = np.row_stack([cresis,cresis2001])
  cutoff = 200.
  dcresis = []
  zcresis = []
  tcresis = []
  for i in range(0,len(cresis[:,0])):
    mindist = np.min(np.sqrt((cresis[i,0]-x)**2+(cresis[i,1]-(y))**2))
    if mindist < cutoff:
      minind = np.argmin(np.sqrt((cresis[i,0]-x)**2+(cresis[i,1]-(y))**2))
      dcresis.append(dists[minind])
      zcresis.append(cresis[i,2])
      tcresis.append(cresis[i,4])
  dcresis = np.array(dcresis)
  zcresis = np.array(zcresis)

  # Plot glacier cross section
  labels=['a','b','c','d','e']
  sortind = [0,2,1,3]
  fig = plt.figure(figsize=(6,4.6))
  gs = matplotlib.gridspec.GridSpec(2,2)
  for k in range(0,len(demtimes)):
    plt.subplot(gs[sortind[k]])
    ax = plt.gca()
    plt.plot(dcresis/1e3,zcresis,'.',c='0.7',markersize=4)
    plt.plot(dists/1e3,floatlib.height(zb),'k:',lw=2)
    if k > 0:
      for j in range(0,k):
        zs_demtime,demtime = zslib.dem_along_flowline(x,y,glacier,years=demtimes[j],verticaldatum='geoid',filt_len=100.)
        ind = np.where(~(np.isnan(zs_demtime[0,:])))[0][-1]
        plt.plot(dists/1e3,zs_demtime[0,:],c='lightblue',lw=1.5)
        plt.plot(dists/1e3,floatlib.icebottom(zb,zs_demtime[0,:]),c='lightblue',lw=1.5)
        plt.plot([dists[ind]/1e3,dists[ind]/1e3],[floatlib.icebottom(zb,zs_demtime[0,:])[ind],zs_demtime[0,ind]],'lightblue',lw=1.5)
    zs_demtime,demtime = zslib.dem_along_flowline(x,y,glacier,years=demtimes[k],verticaldatum='geoid',filt_len=100.)
    plt.plot(dists/1e3,zs_demtime[0,:],'b',lw=1.5)
    plt.plot(dists/1e3,floatlib.icebottom(zb,zs_demtime[0,:]),'b',lw=1.5)
    ind = np.where(~(np.isnan(zs_demtime[0,:])))[0][-1]
    plt.plot([dists[ind]/1e3,dists[ind]/1e3],[floatlib.icebottom(zb,zs_demtime[0,:])[ind],zs_demtime[0,ind]],'b',lw=1.5)
    plt.plot(dists/1e3,zb,'k',lw=1.5)
    if (k == 1):
      plt.xlabel('                                            Distance along flowline (km)',fontsize=18,fontname='Arial')
    elif (k != 3):
      ax.set_xticklabels([])
    if (k == 1):
      plt.ylabel('                                Elevation (m asl)',fontsize=18,fontname='Arial')
    elif (k != 0):
      ax.set_yticklabels([])
    plt.xticks(np.arange(-20,10,5),fontsize=14,fontname='Arial')
    plt.xlim([-10,3])
    plt.yticks(np.arange(-1000,1000,250),fontsize=14,fontname='Arial')
    if glacier == 'Helheim':
      plt.ylim([-800,380])
      #plt.text(-2.5,280,demtimes[k][0:4]+' '+datelib.month_to_string(int(demtimes[k][4:6]))+' '+demtimes[k][6:8],fontsize=11,fontname='Arial')
    plt.text(-5,220,demtimes[k][0:4]+' '+datelib.month_to_string(int(demtimes[k][4:6]))+' '+demtimes[k][6:8],fontsize=18,fontname='Arial')
    plt.text(-9.7,-750,labels[k],fontsize=16,fontname='Arial',fontweight='bold')
  plt.tight_layout()
  plt.subplots_adjust(hspace=0.03,wspace=0.03,top=0.97,right=0.97,left=0.14,bottom=0.13) 
  plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/AGU_"+glacier+"_dem.pdf"),FORMAT='PDF',dpi=600)
  plt.close()

  del cresis,dcresis,zcresis 
  
#######################################
# Terminus position vs. thinning rate #
#######################################

term_dhdt = np.interp(time,terminus_time,terminus_val)
depth_dhdt = np.interp(term_dhdt,dists,zb)
Qbalance = 31.9e9

# Get indices for advance and retreat dh/dt
ind_retreat =[]
ind_advance = []
for i in range(1,len(term_dhdt-1)):
  if (term_dhdt[i] > term_dhdt[i-1]):
    ind_advance.append(i)
  else:
    ind_retreat.append(i)

# Schoof
C=1./1200.
C=((6.3e-25*(917.*9.8)**(4)*(1-917./1024.)**3.)/(4**3*9e5))**(1./(1./3.+1))
q_schoof_terminus = C*(floatlib.height(zb)-zb)**(19./4.)*yearinsec
q_schoof_time = np.interp(term_dhdt,dists,q_schoof_terminus)

gs = matplotlib.gridspec.GridSpec(1,2)
fig = plt.figure(figsize=(6,3.2))

plt.subplot(gs[0])
plt.plot([0,0],[-120,60],c=0.7,lw=1.5)
plt.plot([-2,2.5],[0,0],'k',lw=1.5)
plt.errorbar(term_dhdt[ind_advance]/1e3,dH[ind_advance,0]/365.25*100,yerr=3,fmt='o',color='r',markersize=4,label='advance')
plt.errorbar(term_dhdt[ind_retreat]/1e3,dH[ind_retreat,0]/365.25*100,yerr=3,fmt='o',color='b',markersize=4,label='retreat')
plt.plot(dists/1e3,(Qbalance/5700-q_schoof_terminus)/(365.25*1e4)*100,'k--',lw=2,label='Schoof 2007')
plt.xlabel('Terminus position\n (km)',fontsize=10,fontname='Arial')
plt.ylabel('Dynamic dh/dt (cm/d)',fontsize=10,fontname='Arial')
plt.xticks([-2,-1,0,1,2],fontsize=10,fontname='Arial')
plt.yticks(np.arange(-40,30,10),fontsize=10,fontname='Arial')
plt.ylim([-45,15])
plt.xlim([-2,2.2])
plt.text(-1.8,10,'a',fontname='Arial',fontsize=10,fontweight='bold')
plt.legend(loc=4,numpoints=1,handletextpad=0.2,fontsize=10,columnspacing=0.05,markerscale=1,handlelength=1.2,borderpad=0.25)


plt.subplot(gs[1])
ax = plt.gca()
plt.plot([np.interp(0,dists,zb),np.interp(0,dists,zb)],[-120,60],'k',lw=1.5)
plt.plot([-710,-550],[0,0],'k',lw=1.5)
plt.errorbar(depth_dhdt[ind_advance],dH[ind_advance,0]/365.25*100,yerr=3,fmt='o',color='r',markersize=4)
plt.errorbar(depth_dhdt[ind_retreat],dH[ind_retreat,0]/365.25*100,yerr=3,fmt='o',color='b',markersize=4)
#plt.errorbar(depth_dhdt,dH[:,0]/365.25*100,yerr=3,fmt='o',color='b',markersize=5)
ind = np.argsort(depth_dhdt)
plt.plot(depth_dhdt[ind],(Qbalance/5700-q_schoof_time[ind])/(365.25*1e4)*100,'k--',lw=2,label='Schoof 2007')
plt.xlabel('Bed elevation at terminus\n (m)',fontsize=10,fontname='Arial')
plt.xticks([-700,-650,-600,-550],fontsize=10,fontname='Arial')
plt.yticks(np.arange(-40,30,10),fontsize=10,fontname='Arial')
plt.ylim([-45,15])
plt.xlim([-710,-550])
ax.set_yticklabels([])
plt.text(-702,10,'b',fontname='Arial',fontsize=10,fontweight='bold')

plt.tight_layout()
plt.subplots_adjust(hspace=0.03,wspace=0.04,top=0.97,right=0.96,left=0.13,bottom=0.26) 
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/AGU_"+glacier+"_thinning.pdf"),FORMAT='PDF',dpi=600)
plt.close()

#######################
# Helheim ungrounding #
#######################

for i in range(0,2):
  if i == 0:
    ximage,yimage,image = geotifflib.readrgb(os.path.join(os.getenv("DATA_HOME"),"Mosaics/Helheim/mosaicHelheim.2013-128.148.32713_1-20mgeo.tif"))
    ind = np.argmin(abs(timewv_H-datelib.date_to_fracyear(2013,5,8)))
  else:
    ind = np.argmin(abs(timewv_H-datelib.date_to_fracyear(2014,1,27)))
    ximage,yimage,image = geotifflib.readrgb(os.path.join(os.getenv("DATA_HOME"),"Mosaics/Helheim/mosaicHelheim.2014-027.148.36721_1-20mgeo.tif"))
  data = zabovefloat_H[:,ind]
  

  fig = plt.figure(figsize=(2.5,2.5))
  matplotlib.rc('font',family='Arial')
  ax = plt.gca()
  ax.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower')
  ax.axis('equal')
  p = plt.scatter(xf_H,yf_H,c=data,lw=0,cmap='RdBu_r',vmin=-10,vmax=10)


  xmin = 304000.0
  xmax = 313600.0
  ymin = -2582000.0
  ymax = -2571000.0

  ax.set_xticklabels([])
  ax.set_yticklabels([])
  ax.set_yticks([])
  ax.set_xticks([])
  ax.set_ylim([ymin,ymax])
  ax.set_xlim([xmin,xmax])
  plt.tight_layout()
  
  xmin,xmax = plt.xlim()
  ymin,ymax = plt.ylim()
  path = matplotlib.path.Path([[0.5*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
  			[0.99*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
  			[0.99*(xmax-xmin)+xmin,0.69*(ymax-ymin)+ymin],
  			[0.5*(xmax-xmin)+xmin,0.69*(ymax-ymin)+ymin],
  			[0.5*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin]])
  patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1)
  ax.add_patch(patch)
  cbaxes = fig.add_axes([0.57, 0.91, 0.35, 0.02]) 
  cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=[-5,0,5]) 
  ax.text(xmin+0.56*(xmax-xmin),ymin+0.78*(ymax-ymin),'Hfloat (m)',fontsize=11,fontname='Arial')
  cb.ax.tick_params(labelsize=11)
  ax.plot([xmin+0.54*(xmax-xmin),xmin+0.54*(xmax-xmin)+1e3],[ymin+0.74*(ymax-ymin),ymin+0.74*(ymax-ymin)],'k',linewidth=1.5)
  ax.plot([xmin+0.54*(xmax-xmin),xmin+0.54*(xmax-xmin)],[ymin+0.74*(ymax-ymin),ymin+0.72*(ymax-ymin)],'k',linewidth=1.5)
  ax.plot([xmin+0.54*(xmax-xmin)+1e3,xmin+0.54*(xmax-xmin)+1e3],[ymin+0.74*(ymax-ymin),ymin+0.72*(ymax-ymin)],'k',linewidth=1.5)
  ax.text(xmin+0.58*(xmax-xmin)+1e3,ymin+0.71*(ymax-ymin),'1 km',fontsize=11,fontname='Arial')
  #ax.plot(x_H[ind_eul_H[1]],y_H[ind_eul_H[1]],'ro',markersize=5)

  plt.subplots_adjust(wspace=0,hspace=0,top=0.98,right=0.98,left=0.02,bottom=0.02)

  plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/AGU_Helheim_unground"+str(i)+".pdf"),FORMAT='PDF',dpi=600,transparent=True)
  plt.close()

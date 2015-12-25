# This file compares Helheim's terminus position to its vellib.

import os
import sys
import numpy as np
import vellib, icefrontlib, bedlib, glaclib, zslib, fluxlib, floatlib, climlib, masklib, demshadelib, geotifflib, datelib
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import AutoMinorLocator
import scipy.signal as signal
import cubehelix
from mpl_toolkits.basemap import Basemap

##########
# Inputs #
##########

# Get arguments
args = sys.argv
glacier = args[1][:] # Options: Kanger, Helheim

# Locations for velocities
dists_eul = -1.0*np.array([2,5.0,10.0,15.0,20.0]) # kilometers
dists_mel = np.array([-0.5])
#dists_lag = np.array([1.0,5.0,10.0,20.0,30.0])

# What bed to use for thinning gates through fluxgate method
if glacier == 'Helheim':
  bedsource = 'smith'
elif glacier == 'Kanger':
  bedsource = 'cresis'
    
################
# Plot options #
################

time1 = 2008.0 #start time for plot
time2 = 2016.0 # end time for plot
seasonal = 1 # plot seasonal bars, to highlight seasonal trends
normalized = 0
lagrangian = 0

# Which plots?
plot_overview = 1
plot_radargram = 0
plot_bed = 0
plot_images = 0
plot_years = 1

############ 
# Flowline #
############

x,y,zb,dists = glaclib.load_flowline(glacier,shapefilename='flowline_flightline',filt_len=2.0e3)

##################
# Get ice fronts #
##################

terminus_val, terminus_time = icefrontlib.distance_along_flowline(x,y,dists,glacier,type='icefront',time1=time1,time2=time2)
rift_val, rift_time = icefrontlib.distance_along_flowline(x,y,dists,glacier,type='rift',time1=time1,time2=time2)

########################
# Get calving behavior #
########################

calvingstyle = icefrontlib.calving(glacier)

############################
# Get velocities at points # 
############################

# Find where points are located along flowline
ind_eul=[]
for i in range(0,len(dists_eul)):
  ind_eul.append( (abs(dists - dists_eul[i]*1e3)).argmin() )

# Load velocities for glacier and ice melange
if lagrangian == 1:
  vel_val,vel_time,vel_error,vel_dists,vel_x,vel_y = vellib.velocity_at_lagpoints(x,y,dists,dists_eul*1e3,glacier)
else:
  vel_val,vel_time,vel_error = vellib.velocity_at_eulpoints(x[ind_eul],y[ind_eul],glacier)

velmel_val,velmel_time,velmel_error,velmel_dists,velmel_x,velmel_y = vellib.velocity_at_lagpoints(x,y,dists,dists_mel*1e3,glacier)
velocitypoints = np.column_stack([x[ind_eul],y[ind_eul]])

# Chop to desired time interval
indt = np.where((velmel_time[:,0] > time1) & (velmel_time[:,0] < time2))[0]
vel_time = vel_time[indt]
vel_val = vel_val[indt,:]
velmel_val = velmel_val[indt,:]
velmel_time = velmel_time[indt,:]
del indt

# Get rid of velocities in vel_val that are in front of the ice front 
# (we don't want the melange speed to be used accidentally).
interped = np.interp(vel_time,terminus_time,terminus_val)
vel_val[interped < dists_eul[0]*1.0e3,-1:] = float('NaN')
del interped

###############################################
# Get thinning rates inferred from flux gates #
###############################################

if glacier == 'Helheim':
  xmin = 287000.0
  xmax = 320000.0
  ymin = -2588000.0
  ymax = -2560000.0
elif glacier == 'Kanger':
  xmin = 449800.0
  xmax = 503000.0
  ymin = -2302000.0
  ymax = -2266000.0

xdem,ydem,zdem,timedem,errordem = zslib.dem_grid(glacier,xmin,xmax,ymin,ymax,years='all',verticaldatum='ellipsoid',return_error=True)
dem_time,dem_dH = fluxlib.dem_thinning(glacier,xdem,ydem,zdem,timedem,errordem,"fluxgate3",type='rate')
flux_time,flux_dH = fluxlib.fluxgate_thinning(glacier,"fluxgate3",bedsource=bedsource)
xflux,yflux = fluxlib.fluxbox_geometry(glacier,"fluxgate3")

xrac,yrac,smbrac,timerac = climlib.racmo_at_pts(np.mean(xflux),np.mean(yflux),'smb',filt_len=14.0)
xrac,yrac,runrac,timerac = climlib.racmo_at_pts(np.mean(xflux),np.mean(yflux),'runoff',filt_len=14.0)
xrac,yrac,zsrac,timeraczs = climlib.racmo_at_pts(np.mean(xflux),np.mean(yflux),'zs',filt_len=14.0)
xsif,ysif,sif,timesif = climlib.SIF_at_pts(np.mean(xflux),np.mean(yflux),filt_len=14.)

dH_time,dH_flux,dH_dem,dH_smb = fluxlib.compare_thinning_rates(dem_time,dem_dH,flux_time,flux_dH,timerac,smbrac,rho_i=900.0)

################################
# Get elevations near terminus #
################################

zpt_atm,zptstd_atm,time_atm = zslib.atm_at_pts(x[ind_eul],y[ind_eul],glacier,years='all',maxdist=250.,verticaldatum='geoid',method='average',cutoff='terminus')

zpt_dem,zpterror_dem,time_dem = zslib.dem_at_pts(x[ind_eul],y[ind_eul],glacier,years='all',verticaldatum='geoid',cutoff='terminus',method='average',radius=250.)
zpt_wv,zpterror_wv,time_wv = zslib.dem_at_pts(x[ind_eul],y[ind_eul],glacier,years='all',verticaldatum='geoid',cutoff='terminus',method='average',radius=250.,data='WV')
zpt_tdm,zpterror_tdm,time_tdm = zslib.dem_at_pts(x[ind_eul],y[ind_eul],glacier,years='all',verticaldatum='geoid',cutoff='terminus',method='average',radius=250.,data='TDM')
zpt_spirit,zpterror_spirit,time_spirit = zslib.dem_at_pts(x[ind_eul],y[ind_eul],glacier,years='all',verticaldatum='geoid',cutoff='terminus',method='average',radius=250.,data='SPIRIT')


# Get rid of elevations that are in front of the ice front 
interped = np.interp(time_dem,terminus_time,terminus_val)
zpt_dem[interped < dists_eul[0]*1.0e3+250,0] = float('NaN')
interped = np.interp(time_wv,terminus_time,terminus_val)
zpt_wv[interped < dists_eul[0]*1.0e3+250,0] = float('NaN')
interped = np.interp(time_tdm,terminus_time,terminus_val)
zpt_tdm[interped < dists_eul[0]*1.0e3+250,0] = float('NaN')
interped = np.interp(time_spirit,terminus_time,terminus_val)
zpt_spirit[interped < dists_eul[0]*1.0e3+250,0] = float('NaN')
interped = np.interp(time_atm,terminus_time,terminus_val)
zpt_atm[interped < dists_eul[0]*1.0e3+250,0] = float('NaN')
del interped

###############
# Plot images #
###############

if plot_images == 1:
  
  DIRLANDSAT = os.path.join(os.getenv("DATA_HOME"),"Imagery/Landsat/"+glacier+"/TIF/")
  DIRWV = os.path.join(os.getenv("DATA_HOME"),"Imagery/Worldview/"+glacier+"/")
  DIRTSX = os.path.join(os.getenv("DATA_HOME"),"Mosaics/"+glacier+"/")
  
  if glacier == 'Helheim':
    xmin = 304000.0
    xmax = 314000.0
    ymin = -2582500.0
    ymax = -2572500.0
    imagefiles = ['mosaicHelheim.2013-216.148.34049_1-20mgeo.tif',
             'mosaicHelheim.2013-304.148.35385_1-20mgeo.tif',
             'mosaicHelheim.2014-027.148.36721_1-20mgeo.tif',
             '201403301509_103001002F462D00.tif',
             '201407311341_1020010031DF9D00.tif'] 
    clims = [[0,255],[100,255],[70,255],[0,255],[0,255]]
  if glacier == 'Kanger':
    xmin = 488000.0
    xmax = 498000.0
    ymin = -2298000.0
    ymax = -2288000.0
    imagefiles = ['20120802134148_LE72300122012215EDC00.tif',
                  'mosaicKang.2014-088.65.37640_1-20mgeo.tif',
                  'mosaicKang.2014-182.163.39074_1-20mgeo.tif',
                  'mosaicKang.2014-292.163.40744_1-20mgeo.tif']
    clims = [[0,3],[50,255],[50,255],[50,255]]
  
  images = []
  images_time = []
  images_type = []
  for file in imagefiles:
    if os.path.isfile(DIRLANDSAT+file):
      images.append(geotifflib.read(DIRLANDSAT+file))
      year,month,day = [int(file[0:4]),int(file[4:6]),int(file[6:8])]
      images_type.append('Landsat-8')
    elif os.path.isfile(DIRWV+file):
      images.append(geotifflib.read(DIRWV+file))
      year,month,day = [int(file[0:4]),int(file[4:6]),int(file[6:8])]
      images_type.append('Worldview')
    else:
      images.append(geotifflib.read(DIRTSX+file))
      if glacier == 'Helheim':
        year,month,day = datelib.doy_to_date(float(file[14:18]),float(file[19:22]))
      elif glacier == 'Kanger':
        year,month,day = datelib.doy_to_date(float(file[11:15]),float(file[16:19]))
      images_type.append('TSX')
    images_time.append([year,month,day,datelib.date_to_fracyear(year,month,day)])
  
  N = len(images) # number of images
  images_labels = ['a','b','c','d','e']
  
  plt.figure(figsize=(5,4))
  gs = matplotlib.gridspec.GridSpec(4,N+1,width_ratios=np.append(np.ones(N),0.1))
  gs.update(left=0.02,wspace=0.02,hspace=0.02,top=0.98,right=0.9)
  for i in range(0,N):
  
    xdem,ydem,zdem,demtime = zslib.grid_near_time(images_time[i][3],glacier)
    xvel,yvel,uvel,vvel,vvel,vtime = vellib.tsx_near_time(images_time[i][3],glacier)
    #if i == 0:
      #xmask,ymask,mask = masklib.load(glacier,np.min(xvel),np.max(xvel),np.min(yvel),np.max(yvel),100)
    xf,yf,zabovefloat = floatlib.extent(xdem,ydem,zdem,demtime,glacier,rho_i=917.0,rho_sw=1020.0,bedsource='cresis',verticaldatum='geoid')
  
    for j in [0,1,3]:
      ax = plt.subplot(gs[j,i])
      plt.imshow(images[i][2],extent=[np.min(images[i][0]),np.max(images[i][0]),np.min(images[i][1]),np.max(images[i][1])],origin='lower',cmap='Greys_r')
      plt.clim(clims[i])
      plt.xlim([xmin,xmax])
      plt.ylim([ymin,ymax])
      plt.xticks([])
      plt.yticks([])
      
    ax = plt.subplot(gs[0,i])
    year,month,day = datelib.fracyear_to_date(images_time[i][3])
    plt.text(xmin+500,ymax-1.25e3,str(year)+'-'+str(month)+'-'+str(int(np.floor(day))),backgroundcolor='w',fontsize=10)
    plt.text(xmin+500,ymin+1e3,images_type[i],backgroundcolor='w',fontsize=10)
    plt.text(xmin+500,ymax-3e3,images_labels[i],fontsize=10,fontweight='bold')
    
    ax = plt.subplot(gs[1,i])
    #vvel = np.ma.masked_array(vvel,mask)
    year,month,day = datelib.fracyear_to_date(vtime)
    xfront,yfront,ftime = icefrontlib.near_time(images_time[i][3],glacier)
    im = plt.imshow(vvel/1e3,extent=[np.min(xvel),np.max(xvel),np.min(yvel),np.max(yvel)],origin='lower')
    plt.plot(xfront,yfront,'k',linewidth=1.5)
    plt.clim([5,10])
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    plt.xticks([])
    plt.yticks([])
    
    if i == N-1:
      ax = plt.subplot(gs[1,N])
      cb = plt.colorbar(im,cax=ax)
      cb.set_ticks(range(4,12))
      cb.ax.tick_params(labelsize=8)
      cb.set_label('Speed (km/yr)',fontsize=10)
    
    ax = plt.subplot(gs[2,i])
    shadeddem = demshadelib.set_shade(zdem,0,220)
    im = plt.imshow(shadeddem,extent=[np.min(xdem),np.max(xdem),np.min(ydem),np.max(ydem)],origin='lower')
    year,month,day = datelib.fracyear_to_date(demtime)
    plt.text(xmin+500,ymax-1.25e3,str(year)+'-'+str(month)+'-'+str(int(np.floor(day))),backgroundcolor='w',fontsize=10)
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    plt.clim([0,220])
    plt.xticks([])
    plt.yticks([])
    
    if i == N-1:
      ax = plt.subplot(gs[2,N])
      cb = plt.colorbar(im,cax=ax)
      cb.set_ticks(np.arange(0,250,50))
      cb.ax.tick_params(labelsize=8)
      cb.set_label('Elevation (m)',fontsize=10)
    
    ax = plt.subplot(gs[3,i])
    im=plt.scatter(xf,yf,c=zabovefloat,lw=0,vmin=-25,vmax=25,cmap='RdBu_r',s=3)
    year,month,day = datelib.fracyear_to_date(demtime)
    plt.text(xmin+500,ymax-1.25e3,str(year)+'-'+str(month)+'-'+str(int(np.floor(day))),backgroundcolor='w',fontsize=10)
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    plt.xticks([])
    plt.yticks([])
    
    if i == N-1:
      ax = plt.subplot(gs[3,N])
      cb = plt.colorbar(im,cax=ax)
      cb.set_ticks(np.arange(-20,40,20))
      cb.ax.tick_params(labelsize=8)
      cb.set_label('Height above \n flotation (m)',fontsize=10)
  
  plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_images.png"),FORMAT='PNG',dpi=800)
  plt.close()
  
  del clims,year,month,day,imagefiles,xmin,xmax,ymin,ymax,xvel,yvel,vvel,vtime,DIRTSX,DIRLANDSAT

###########################################
# Plot velocity vs. terminus through time #
###########################################

if plot_overview == 1:
  plt.figure(figsize=(7.45,7))
  gs = matplotlib.gridspec.GridSpec(7,1)

  # Plot terminus
  plt.subplot(gs[-2, :])
  ax = plt.gca()
  ind = np.where(calvingstyle[:,1] == 'Tabular')[0]
  plt.ylim([-6,6])
  if seasonal:
    xTickPos = np.linspace(np.floor(time1)-0.25,np.ceil(time2)-0.25,(np.ceil(time2)-np.floor(time1))*2+1)
    ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.85','w'],linewidth=0)
  for i in ind:
    if i == ind[1]:
      plt.bar(float(calvingstyle[i,0])-0.025,10.0,width=0.025,color=[0.6,0.6,1],edgecolor='none',label='Tabular',bottom=-5)
    elif i !=0:
      plt.bar(float(calvingstyle[i,0])-0.025,10.0,width=0.025,color=[0.6,0.6,1],edgecolor='none',bottom=-5)
  ind = np.where((calvingstyle[:,1] == 'Mixed'))[0]
  for i in ind:
    if i == ind[1]:
      plt.bar(float(calvingstyle[i,0])-0.025,10.0,width=0.025,color=[0.4,0.8,0.6],edgecolor='none',label='Mixed',bottom=-5)
    elif i != 0:
      plt.bar(float(calvingstyle[i,0])-0.025,10.0,width=0.025,color=[0.4,0.8,0.6],edgecolor='none',bottom=-5)
  ind = np.where((calvingstyle[:,1] == 'Domino'))[0]
  for i in ind:
    if i == ind[1]:
      plt.bar(float(calvingstyle[i,0])-0.025,10.0,width=0.025,color=[1,0.6,0.6],edgecolor='none',label='Nontabular',bottom=-5)
    elif i != 0:
      plt.bar(float(calvingstyle[i,0])-0.025,10.0,width=0.025,color=[1,0.6,0.6],edgecolor='none',bottom=-5)
  if plot_images == 1:
   for i in range(0,len(images)):
      plt.plot([images_time[i][3],images_time[i][3]],ax.get_ylim(),'--',color='0.3')
  plt.legend(loc=2,fontsize=10,numpoints=1,handlelength=0.3,labelspacing=0.05,ncol=3,columnspacing=0.7,handletextpad=0.2,borderpad=0.25)
  nonnan = np.where(~(np.isnan(terminus_val)))[0]
  plt.plot(terminus_time[nonnan],terminus_val[nonnan]/1e3,'ko',linewidth=1,markersize=2)
  x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
  ax.xaxis.set_major_formatter(x_formatter)
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  ax.yaxis.set_minor_locator(AutoMinorLocator(2))
  ax.set_xticklabels([])
  plt.xticks(range(2000,2016),fontsize=10,fontname="Arial")
  plt.xlim([time1,time2])
  plt.yticks(np.arange(-6,8,2),fontsize=10,fontname="Arial")
  plt.ylabel('Terminus \n (km)',fontsize=10,fontname="Arial")
  ax.tick_params('both', length=6, width=1.25, which='major')
  ax.tick_params('both', length=3, width=1, which='minor')
  if glacier == 'Helheim':
    plt.ylim([-3,3])
    plt.text(2008.13,-2.6,'e',fontsize=11,fontname='arial',weight='bold')
  elif glacier == 'Kanger':
    plt.ylim([-4.5,4.5])
    plt.text(2008.13,-3.8,'e',fontsize=11,fontname='arial',weight='bold')
  else:  
    plt.ylim(np.floor((np.min(terminus_val))/1e3),np.ceil((np.max(terminus_val))/1e3))

  # Plot velocities
  ax = plt.subplot(gs[0:-5, :]) 
  #plt.plot([2000,2014],[0,0],'k')
  #coloptions=['k','r','y','g','b']
  coloptions=['r','b','g','limegreen','gold']
  markoptions=['o','o','o','o','o']
  if plot_images == 1:
   for i in range(0,len(images)):
      plt.plot([images_time[i][3],images_time[i][3]],[0,12],'--',color='0.3')
      plt.text(images_time[i][3]+0.05,0.3,images_labels[i],fontsize=10,fontweight='bold')
  if normalized == 1:
    for i in range(0,len(dists_eul)):
      nonnan = np.where(~(np.isnan(vel_val[:,i])))[0]
      plt.plot(vel_time[nonnan],(vel_val[nonnan,i]-vel_val[nonnan[0],i])/1e3,markoptions[i],color=coloptions[i],label=glacier[0]+'{0:02d}'.format(int(abs(dists_eul[i]))),linewidth=1,markersize=3)
    plt.yticks(range(-3,3,1),fontsize=10,fontname='Arial')
    if glacier == 'Kanger':
      plt.ylim([-1.5,3])
  else:
    for i in range(0,len(dists_eul)):
      nonnan = np.where(~(np.isnan(vel_val[:,i])))[0]
      plt.plot(vel_time[nonnan],(vel_val[nonnan,i])/1e3,markoptions[i],color=coloptions[i],label=glacier[0]+'{0:02d}'.format(int(abs(dists_eul[i]))),markersize=3)
    plt.yticks(range(2,12),fontsize=10,fontname="Arial")
    if glacier == 'Helheim':
      plt.ylim([3.5,11])
      plt.text(2008.13,8.0,'a',fontsize=11,fontname='arial',weight='bold')
    elif glacier == 'Kanger':
      plt.ylim([2.5,12])
      plt.text(2008.13,8.2,'a',fontsize=11,fontname='arial',weight='bold')
  plt.ylabel('Glacier speed \n (km yr$^{-1}$)',fontsize=10,fontname="Arial")
  x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
  ax.xaxis.set_major_formatter(x_formatter)
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  ax.yaxis.set_minor_locator(AutoMinorLocator(2))
  plt.legend(loc=2,borderpad=0.3,fontsize=10,numpoints=1,handlelength=0.2,handletextpad=0.5,labelspacing=0.1,ncol=2,columnspacing=0.8)
  matplotlib.rc('font',family="Arial",)
  ax.tick_params('both', length=6, width=1.25, which='major')
  ax.tick_params('both', length=3, width=1, which='minor')
  labels=[]
  plt.xticks(range(2000,2017),fontsize=10,fontname="Arial")
  ax.set_xticklabels([])
  if seasonal:
    xTickPos = np.linspace(np.floor(time1)-0.25,np.ceil(time2)-0.25,(np.ceil(time2)-np.floor(time1))*2+1)
    ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.85','w'],linewidth=0)
  plt.xlim([time1,time2])
 
  # Plot surface elevations
  plt.subplot(gs[-5, :])
  ax = plt.gca()
  plt.ylim([-10,230])
  if plot_images == 1:
   for i in range(0,len(images)):
      plt.plot([images_time[i][3],images_time[i][3]],ax.get_ylim(),'--',color='0.3')
  if seasonal:
    xTickPos = np.linspace(np.floor(time1)-0.25,np.ceil(time2)-0.25,(np.ceil(time2)-np.floor(time1))*2+1)
    plt.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.85','w'],linewidth=0,zorder=1)
  if glacier == 'Helheim':
    plt.fill_between([time1,time2],np.ones(2)*floatlib.height(zb[ind_eul[0]]-50)-floatlib.height(zb[ind_eul[0]]),np.ones(2)*floatlib.height(zb[ind_eul[0]]+50)-floatlib.height(zb[ind_eul[0]]),
    alpha=0.1,facecolor='b',edgecolor='b',antialiased=True,zorder=2)
    plt.plot([time1,time2],[0,0],'b:',linewidth=1.0)
    plt.errorbar(time_wv,zpt_wv[:,1]-floatlib.height(zb[ind_eul[1]]),capsize=1,yerr=zpterror_wv,fmt='o',color=coloptions[1],markersize=3,label='WV')
    plt.errorbar(time_tdm,zpt_tdm[:,1]-floatlib.height(zb[ind_eul[1]]),capsize=1,yerr=zpterror_tdm,fmt='^',color=coloptions[1],markersize=3,label='TDM')
    plt.plot(time_atm,zpt_atm[:,1]-floatlib.height(zb[ind_eul[1]]),'+',color=coloptions[1],markersize=4,label='ATM')
  if glacier == 'Kanger':
    plt.fill_between([time1,time2],np.ones(2)*floatlib.height(zb[ind_eul[2]]-50)-floatlib.height(zb[ind_eul[2]]),np.ones(2)*floatlib.height(zb[ind_eul[2]]+50)-floatlib.height(zb[ind_eul[2]]),
    alpha=0.1,facecolor='g',edgecolor='b',antialiased=True,zorder=2)
    plt.plot([time1,time2],[0,0],'g:',linewidth=1.0)
    plt.errorbar(time_wv,zpt_wv[:,2]-floatlib.height(zb[ind_eul[2]]),capsize=1,yerr=zpterror_wv,fmt='o',color=coloptions[2],markersize=3,label='WV')
    plt.errorbar(time_tdm,zpt_tdm[:,2]-floatlib.height(zb[ind_eul[2]]),capsize=1,yerr=zpterror_tdm,fmt='^',color=coloptions[2],markersize=3,label='TDM')
    plt.errorbar(time_spirit,zpt_spirit[:,2]-floatlib.height(zb[ind_eul[2]]),capsize=1,yerr=zpterror_spirit,fmt='v',color=coloptions[2],markersize=3,label='SPIRIT')
    plt.plot(time_atm,zpt_atm[:,2]-floatlib.height(zb[ind_eul[2]]),'+',color=coloptions[2],markersize=4,label='ATM')
  x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
  ax.xaxis.set_major_formatter(x_formatter)
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  ax.yaxis.set_minor_locator(AutoMinorLocator(2))
  ax.set_xticklabels([])
  plt.xticks(range(2000,2016),fontsize=10,fontname="Arial")
  plt.xlim([time1,time2])
  #plt.ylabel('Height above\n flotation (m)',fontsize=10,fontname="Arial")
  ax.tick_params('both', length=6, width=1.25, which='major')
  ax.tick_params('both', length=3, width=1, which='minor')
  if glacier == 'Helheim':
    plt.yticks(np.arange(80,120,10),fontsize=10,fontname="Arial")
    plt.ylim([70,105])
    plt.text(2008.13,98,'b',fontsize=11,fontname='arial',weight='bold')
  elif glacier == 'Kanger':
    plt.yticks(np.arange(160,250,20),fontsize=10,fontname="Arial")
    plt.ylim([165,230])
    plt.text(2008.13,217,'b',fontsize=11,fontname='arial',weight='bold')
  plt.legend(loc=3,borderpad=0.3,handleheight=0.1,fontsize=10,numpoints=1,handlelength=0.4,labelspacing=0.05,ncol=4,columnspacing=0.7,handletextpad=0.5)

  # Plot surface elevations
  plt.subplot(gs[-4, :])
  ax = plt.gca()
  plt.ylim([-120,120])
  if plot_images == 1:
   for i in range(0,len(images)):
      plt.plot([images_time[i][3],images_time[i][3]],ax.get_ylim(),'--',color='0.3')
  if seasonal:
    xTickPos = np.linspace(np.floor(time1)-0.25,np.ceil(time2)-0.25,(np.ceil(time2)-np.floor(time1))*2+1)
    plt.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.85','w'],linewidth=0,zorder=1)
  if glacier == 'Helheim':
		plt.fill_between([time1,time2],np.ones(2)*floatlib.height(zb[ind_eul[0]]-50)-floatlib.height(zb[ind_eul[0]]),np.ones(2)*floatlib.height(zb[ind_eul[0]]+50)-floatlib.height(zb[ind_eul[0]]),
  			alpha=0.1,facecolor='r',edgecolor='r',antialiased=True,zorder=2)
		plt.plot([time1,time2],[0,0],'r:',linewidth=1.0)
		plt.errorbar(time_wv,zpt_wv[:,0]-floatlib.height(zb[ind_eul[0]]),capsize=1,yerr=zpterror_wv,fmt='o',color=coloptions[0],markersize=3,label='WV')
		plt.errorbar(time_tdm,zpt_tdm[:,0]-floatlib.height(zb[ind_eul[0]]),capsize=1,yerr=zpterror_tdm,fmt='^',color=coloptions[0],markersize=3,label='TDM')
		plt.plot(time_atm,zpt_atm[:,0]-floatlib.height(zb[ind_eul[0]]),'+',color=coloptions[0],markersize=4,label='ATM')
  if glacier == 'Kanger':
		plt.fill_between([time1,time2],np.ones(2)*floatlib.height(zb[ind_eul[1]]-50)-floatlib.height(zb[ind_eul[1]]),np.ones(2)*floatlib.height(zb[ind_eul[1]]+50)-floatlib.height(zb[ind_eul[1]]),
  			alpha=0.1,facecolor='b',edgecolor='b',antialiased=True,zorder=2)
		plt.plot([time1,time2],[0,0],'b:',linewidth=1.0)
		plt.errorbar(time_wv,zpt_wv[:,1]-floatlib.height(zb[ind_eul[1]]),capsize=1,yerr=zpterror_wv,fmt='o',color=coloptions[1],markersize=3,label='WV')
		plt.errorbar(time_tdm,zpt_tdm[:,1]-floatlib.height(zb[ind_eul[1]]),capsize=1,yerr=zpterror_tdm,fmt='^',color=coloptions[1],markersize=3,label='TDM')
		plt.errorbar(time_spirit,zpt_spirit[:,1]-floatlib.height(zb[ind_eul[1]]),capsize=1,yerr=zpterror_spirit,fmt='v',color=coloptions[1],markersize=3,label='SPIRIT')
		plt.plot(time_atm,zpt_atm[:,1]-floatlib.height(zb[ind_eul[1]]),'+',color=coloptions[1],markersize=4,label='ATM')
  x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
  ax.xaxis.set_major_formatter(x_formatter)
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  ax.yaxis.set_minor_locator(AutoMinorLocator(2))
  ax.set_xticklabels([])
  plt.xticks(range(2000,2016),fontsize=10,fontname="Arial")
  plt.xlim([time1,time2])
  plt.ylabel('$\hspace{8}$ Height above flotation \n                          (m)',fontsize=10,fontname="Arial")
  ax.tick_params('both', length=6, width=1.25, which='major')
  ax.tick_params('both', length=3, width=1, which='minor')
  if glacier == 'Helheim':
    plt.yticks(np.arange(-10,40,10),fontsize=10,fontname="Arial")
    plt.ylim([-5,30])
    plt.text(2008.13,23,'c',fontsize=11,fontname='arial',weight='bold')
  elif glacier == 'Kanger':
    plt.yticks(np.arange(-20,60,20),fontsize=10,fontname="Arial")
    plt.ylim([-25,40])
    plt.text(2008.13,30,'c',fontsize=11,fontname='arial',weight='bold')
  plt.legend(loc=3,borderpad=0.3,handleheight=0.1,fontsize=10,numpoints=1,handlelength=0.4,labelspacing=0.05,ncol=4,columnspacing=0.7,handletextpad=0.5)

  
  # Plot thinning rates
  plt.subplot(gs[-3, :])
  ax = plt.gca()
  plt.ylim([-140,140])
  plt.plot([time1,time2],[0,0],'k:',lw=1)
  if plot_images == 1:
   for i in range(0,len(images)):
      plt.plot([images_time[i][3],images_time[i][3]],ax.get_ylim(),'--',color='0.3')
  if seasonal:
    xTickPos = np.linspace(np.floor(time1)-0.25,np.ceil(time2)-0.25,(np.ceil(time2)-np.floor(time1))*2+1)
    ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.85','w'],linewidth=0)
  nonnan = np.where(~(np.isnan(flux_dH[:,0])))[0]
  plt.plot(timeraczs[1:],np.diff(zsrac)/np.diff(timeraczs),c='0.6',label='SMB',lw=1.5,zorder=5)
  plt.errorbar(flux_time[nonnan],flux_dH[nonnan,0],yerr=flux_dH[nonnan,1],fmt='o',color='k',markersize=2.5,label='Flux',zorder=6,capsize=1)
  plt.errorbar(dem_time[:,0],dem_dH[:,0],xerr=dem_time[:,1],yerr=dem_dH[:,1],fmt='^',c='r',markersize=3,label='DEM',capsize=0,zorder=7)
  plt.legend(loc=3,handleheight=0.1,borderpad=0.3,fontsize=10,numpoints=1,handlelength=0.4,labelspacing=0.05,ncol=3,columnspacing=0.7,handletextpad=0.5)
  x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
  ax.xaxis.set_major_formatter(x_formatter)
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  ax.yaxis.set_minor_locator(AutoMinorLocator(2))
  ax.set_xticklabels([])
  plt.xticks(range(2000,2016),fontsize=10,fontname="Arial")
  plt.xlim([time1,time2])
  plt.ylabel('dH/dt (m yr$^{-1}$)',fontsize=10,fontname="Arial")
  ax.tick_params('both', length=6, width=1.25, which='major')
  ax.tick_params('both', length=3, width=1, which='minor')
  if glacier == 'Helheim':
    plt.yticks(np.arange(-100,120,50),fontsize=10,fontname="Arial")
    plt.ylim([-120,80])
    plt.text(2008.13,45,'d',fontsize=11,fontname='arial',weight='bold')
  elif glacier == 'Kanger':
    plt.yticks(np.arange(-100,150,50),fontsize=10,fontname="Arial")
    plt.ylim([-120,80])
    plt.text(2008.13,50,'d',fontsize=11,fontname='arial',weight='bold')

  # Plot presence of rigid melange and/or climate variables
  plt.subplot(gs[-1, :])
  rigid2 = np.zeros(len(velmel_time))
  rigid2[:] = float('nan')
  for i in range(0,len(velmel_time)):
    nonnan = np.where(~(np.isnan(velmel_val[i,:])))[0]
    if len(nonnan) > 0:
      rigid2[i] = 1
  plt.ylabel('Sea Ice Fraction',fontsize=10,fontname="Arial")
  ax = plt.gca()
  plt.ylim([0,1])
  plt.yticks(np.arange(0,1.5,0.5),fontsize=10,fontname='arial')
  x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
  ax.xaxis.set_major_formatter(x_formatter)
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  ax.yaxis.set_minor_locator(AutoMinorLocator(2))
  matplotlib.rc('font',family="Arial",)
  ax.tick_params('both', length=6, width=1.25, which='major')
  ax.tick_params('both', length=3, width=1, which='minor')
  labels=[]
  plt.xticks(range(2000,2017),fontsize=10,fontname="Arial")
  ax.set_xticklabels([])
  if seasonal:
    xTickPos = np.linspace(np.floor(time1)-0.25,np.floor(time2)-0.25,(np.ceil(time2)-np.ceil(time1))*2+1)
    ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.85','w'],linewidth=0)
  for i in range(2000,2017):
    labels.append('Jan \n'+str(i))
  plt.xticks(range(2000,2017))
  ax.set_xticklabels(labels,fontsize=10,fontname='Arial')
  plt.xlim([time1,time2])
  #plt.bar(time1,0.2,bottom=0.4,width=time2-time1,color='0.6',edgecolor='0.6',label=)
  #plt.bar(velmel_time[:,0]-velmel_time[:,1]/2,np.ones(len(velmel_time[:,0]))*0.2,bottom=0.4,width=velmel_time[:,1],color='blue',edgecolor='none',label='Nonrigid')
  #plt.bar(velmel_time[rigid2==1,0]-velmel_time[rigid2==1,1]/2,rigid2[rigid2==1]*0.2,bottom=0.4,width=velmel_time[rigid2==1,1],color='lightblue',edgecolor='none',label='Rigid')
  plt.plot(timesif,sif,'k',lw=1.5,label='SIF')
  ax2 = ax.twinx()
  ax2.plot([0,0],[1,1],'k',lw=1.5,label='SIF')
  ax2.plot(timerac,runrac,'-',c='0.4',dashes=[2,1],lw=1.5,label='Runoff')
  ax2.set_xlim([time1,time2])
  ax2.set_ylim([0,50])
  ax2.set_yticks([0,25,50])
  ax2.set_yticklabels([0,25,50],fontsize=10,fontname='Arial')
  ax2.set_ylabel('Runoff (kg m$^{-2}$ d$^{-1}$)',fontsize=10,fontname='Arial')
  plt.legend(loc=3,borderpad=0.3,fontsize=10,numpoints=1,handlelength=0.7,labelspacing=0.05,ncol=4,columnspacing=0.7,handletextpad=0.5)
  ax.text(2008.13,0.35,'f',fontsize=11,fontname='arial',weight='bold')
  if plot_images == 1:
   for i in range(0,len(images)):
      plt.plot([images_time[i][3],images_time[i][3]],ax.get_ylim(),'--',color='0.3')

  plt.tight_layout()
  plt.subplots_adjust(hspace=0.05,wspace=0) 
  plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_vel_time_"+str(int(time1))+"to"+str(int(time2))+".pdf"),FORMAT='PDF')
  plt.close()

  #########################################
  # Plot overview map for previous figure #
  #########################################
  
  # Image for plotting
  if glacier == "Helheim":
    imagetime = datelib.date_to_fracyear(2014,7,4)
    ximage,yimage,image = geotifflib.readrgb(os.path.join(os.getenv("DATA_HOME"),"Imagery/Landsat/Helheim/TIF/20140704140535_LC82330132014185LGN00.tif"))
  elif glacier == "Kanger":
    imagetime = datelib.date_to_fracyear(2014,7,6)
    ximage,yimage,image = geotifflib.readrgb(os.path.join(os.getenv("DATA_HOME"),"Imagery/Landsat/Kanger/TIF/20140706135251_LC82310122014187LGN00.tif"))

  # Load velocity record
  xvel = np.arange(np.min(ximage),np.max(ximage),100)
  yvel = np.arange(np.min(yimage),np.max(yimage),100)
  vx,vy = vellib.inversion_3D(glacier,xvel,yvel,imagetime,dir_velocity_out='none',blur=False)
  vel = np.sqrt(vx**2+vy**2)
  del vx,vy

  # Load mask
  xmask,ymask,mask = masklib.load_grid(glacier,np.min(xvel),np.max(xvel),np.min(yvel),np.max(yvel),100,icefront_time=datelib.date_to_fracyear(2014,7,4))
  vel_masked = np.ma.masked_array(vel,mask)
  
  fig = plt.figure(figsize=(3,3))

  cx = cubehelix.cmap(start=1.2,rot=-1.1,reverse=True,minLight=0.1,sat=2)
  p=plt.imshow(vel_masked/1e3,extent=[np.min(xvel),np.max(xvel),np.min(yvel),np.max(yvel)],origin='lower',clim=[0,10],cmap=cx)
  ax = plt.gca()
  ax.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
  ax.imshow(vel_masked/1e3,extent=[np.min(xvel),np.max(xvel),np.min(yvel),np.max(yvel)],origin='lower',alpha=0.5,clim=[0,10],cmap=cx)
  ax.axis('equal')
  path = matplotlib.path.Path(np.column_stack([xflux,yflux]))
  patch = matplotlib.patches.PathPatch(path,hatch='xxx',edgecolor='k',facecolor='none',lw=1)
  ax.add_patch(patch)
  #plt.plot(xflux,yflux,'k',linewidth=1)

  if glacier == 'Kanger':
    xmin = 468000.
    xmax = 498000.
    ymin = -2299000.
    ymax = -2264000.
  elif glacier == 'Helheim':
    xmin = 283000.
    xmax = 313000.
    ymin = -2587000.
    ymax = -2552000.
  ax.plot(x[dists>-21e3],y[dists>-21e3],'k',lw=1)
  for i in range(0,len(velocitypoints)):
    ax.plot(velocitypoints[i,0],velocitypoints[i,1],'o',color=coloptions[i],markersize=5)
    if i in [1,2,4]:
      ax.text(velocitypoints[i,0]-3.8e3,velocitypoints[i,1]-1.5e3,glacier[0]+'{0:02d}'.format(int(abs(dists_eul[i]))),fontsize=10,fontname='Arial')
    else:
      ax.text(velocitypoints[i,0],velocitypoints[i,1]+0.75e3,glacier[0]+'{0:02d}'.format(int(abs(dists_eul[i]))),fontsize=10,fontname='Arial')
  ax.set_xticklabels([])
  ax.set_yticklabels([])
  ax.set_yticks([])
  ax.set_xticks([])
  ax.set_ylim([ymin,ymax])
  ax.set_xlim([xmin,xmax])

  plt.tight_layout()
  
  xmin,xmax = plt.xlim()
  ymin,ymax = plt.ylim()
  path = matplotlib.path.Path([[0.57*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
  			[0.99*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
  			[0.99*(xmax-xmin)+xmin,0.72*(ymax-ymin)+ymin],
  			[0.57*(xmax-xmin)+xmin,0.72*(ymax-ymin)+ymin],
  			[0.57*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin]])
  patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1)
  ax.add_patch(patch)
  cbaxes = fig.add_axes([0.59, 0.87, 0.28, 0.02]) 
  cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=[0,5,10]) 
  ax.text(xmin+0.58*(xmax-xmin),ymin+0.81*(ymax-ymin),'Speed (km yr$^{-1}$)',fontsize=10)
  cb.ax.tick_params(labelsize=8)
  ax.plot([xmin+0.61*(xmax-xmin),xmin+0.61*(xmax-xmin)+5e3],[ymin+0.77*(ymax-ymin),ymin+0.77*(ymax-ymin)],'k',linewidth=1.5)
  ax.plot([xmin+0.61*(xmax-xmin),xmin+0.61*(xmax-xmin)],[ymin+0.77*(ymax-ymin),ymin+0.75*(ymax-ymin)],'k',linewidth=1.5)
  ax.plot([xmin+0.61*(xmax-xmin)+5e3,xmin+0.61*(xmax-xmin)+5e3],[ymin+0.77*(ymax-ymin),ymin+0.75*(ymax-ymin)],'k',linewidth=1.5)
  ax.text(xmin+0.64*(xmax-xmin)+5e3,ymin+0.75*(ymax-ymin),'5 km',fontsize=10)
  ax.text(xmin+0.02*(xmax-xmin),ymin+0.93*(ymax-ymin),'a',fontweight='bold',fontsize=11)

  ax2 = fig.add_axes([0.08, 0.08, 0.16, 0.3])
  m = Basemap(width=1600000,height=3000000,
            resolution='l',projection='stere',\
            lat_ts=60,lat_0=72,lon_0=-41.)
  m.drawcoastlines(linewidth=0.75)
  m.drawmapboundary(fill_color='lightblue')
  m.fillcontinents(color='w',lake_color='w')
  if glacier == 'Helheim':
    xg,yg = m([-38.3],[66.1])
  elif glacier == 'Kanger':
    xg,yg = m ([-33.0],[68.63333])
  m.plot(xg,yg,'ro',markersize=4)

  plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_velocity_map.pdf"),FORMAT='PDF',dpi=400)
  plt.close()

  del m,xg,yg,ax,patch,path,xmin,xmax,ymin,ymax,cbaxes,cb,vel_masked,cx,xmask,ymask,mask,xvel,yvel,vel,ximage,yimage,image,imagetime,ax2,fig

##############
# Plot years #
##############

if plot_years == 1:
	years = range(2010,2015)
	plt.figure(figsize=(3.5,5))
	gs = matplotlib.gridspec.GridSpec(4,1)
	colors=['b','g','limegreen','gold','r']
	markers = ['s','o','D','^','v']

	plt.subplot(gs[1])
	ax = plt.gca()
	n=0
	for year in years:
		plt.plot(terminus_time-year,terminus_val/1e3,markers[n]+'-',color=colors[n],markersize=3,lw=0.75)
		n=n+1
	plt.xlim([0,1])
	x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
	ax.xaxis.set_major_formatter(x_formatter)
	ax.xaxis.set_minor_locator(AutoMinorLocator(2))
	ax.yaxis.set_minor_locator(AutoMinorLocator(2))
	ax.set_xticklabels([])
	plt.yticks(np.arange(-3,4,1),fontsize=9,fontname='arial')
	if glacier == 'Helheim':
		plt.ylim([-2.5,2.5])
	elif glacier == 'Kanger':
	  plt.ylim([-4,4])
	plt.ylabel('Terminus (km)',fontsize=10,fontname='arial')
	plt.xticks(np.arange(0,1.1,0.25))

	plt.subplot(gs[0])  
	ax = plt.gca()
	n=0
	for year in years:
		plt.plot(vel_time-year,vel_val[:,1]/1e3,markers[n]+'-',color=colors[n],markersize=3,lw=0.75,label=year)
		n=n+1
	plt.xlim([0,1])
	x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
	ax.xaxis.set_major_formatter(x_formatter)
	ax.xaxis.set_minor_locator(AutoMinorLocator(2))
	ax.yaxis.set_minor_locator(AutoMinorLocator(2))
	ax.set_xticklabels([])
	plt.yticks(range(6,12),fontsize=10,fontname='arial')
	if glacier == 'Helheim':
	  plt.ylim([6,8.5])
	elif glacier == 'Kanger':
	  plt.ylim([7.5,11])
	plt.ylabel('Speed \n (km yr$^{-1}$)',fontsize=10,fontname='arial')
	if glacier == 'Helheim':
	  plt.legend(loc=2,fontsize=9,numpoints=1,ncol=2,handlelength=0.2,handletextpad=0.5,labelspacing=0.2,columnspacing=0.2)
	plt.xticks(np.arange(0,1.1,0.25))
	
	plt.subplot(gs[2])
	ax = plt.gca()
	n=0
	for year in years:
		nonnan = np.where(~(np.isnan(flux_dH[:,0])))[0]
		plt.plot(flux_time[nonnan]-year,flux_dH[nonnan,0],markers[n]+'-',color=colors[n],markersize=3,lw=0.75)
		n=n+1
	plt.xlim([0,1])
	x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
	ax.xaxis.set_major_formatter(x_formatter)
	ax.xaxis.set_minor_locator(AutoMinorLocator(2))
	ax.yaxis.set_minor_locator(AutoMinorLocator(2))
	ax.set_xticklabels([])
	plt.plot([0,1],[0,0],'k',lw=0.75)
	if glacier == 'Helheim':
	  plt.ylim([-110,60])
	  plt.yticks([-100,-50,0,50],fontsize=9,fontname='arial')
	elif glacier == 'Kanger':
	  plt.ylim([-50,50])
	  plt.yticks([-50,-25,0,25,50],fontsize=9,fontname='arial')
	plt.ylabel('dH/dt (m yr$^{-1}$)',fontsize=10,fontname='arial')
	plt.xticks(np.arange(0,1.1,0.25))

	plt.subplot(gs[3])
	ax = plt.gca()
	n=0
	for year in years:

		all = np.r_[zpt_dem[:,1],zpt_atm[:,1]]
		all_time = np.r_[time_dem,time_atm]
		sortind = np.argsort(all_time)
		all = all[sortind]
		all_time = all_time[sortind]
		nonnan = np.where(~(np.isnan(all)))[0]
		plt.plot(all_time[nonnan]-year,all[nonnan],markers[n]+'-',color=colors[n],markersize=3,lw=0.75)
		n=n+1
	plt.xlim([0,1])
	if glacier == 'Helheim':
	  plt.ylim([145,175])
	  plt.yticks([150,160,170],fontsize=9,fontname='arial')
	elif glacier == 'Kanger':
	  plt.ylim([100,160])
	  plt.yticks([100,125,150],fontsize=9,fontname='arial')
	plt.xticks([0,0.25,0.5,0.75,1.],['1-Jan','1-Mar','1-Jun','1-Sept','1-Jan'],fontsize=10) 
	plt.ylabel('Elevation (m)',fontsize=10,fontname='arial')
	#plt.yticks([150,160,170],fontsize=9,fontname='arial')
  
	plt.tight_layout()
	plt.subplots_adjust(hspace=0.05,wspace=0)
	plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_stacked.pdf"),FORMAT='PDF')
	plt.close()
   
###########################
# Plot velocity radargram #
###########################

if plot_radargram == 1:

  flowline_v,flowline_t,termini = vellib.velocity_along_flowline(x,y,glacier,cutoff='terminus',data='TSX')
  flowline_tint=np.linspace(2008.0,2015,1051)
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
  masked_array = np.ma.array (flowline_vint, mask=np.isnan(flowline_vint))
  cmap = matplotlib.cm.jet
  cmap.set_bad('0.65',1.)

  plt.figure(figsize=(6.5,5))
  plt.subplot(121)
  ax=plt.gca()
  plt.imshow(flowline_vint/1e3,extent=[flowline_tint[0],flowline_tint[-1],(dists[-1])/1e3,(dists[0])/1e3],interpolation='nearest',aspect='auto',cmap=cmap)
  plt.fill_between(terminus_time,(terminus_val)/1e3,5,color='w')
  plt.plot(terminus_time,(terminus_val)/1e3,'ko--',markersize=3,linewidth=1)
  ax.get_xaxis().get_major_formatter().set_useOffset(False)
  plt.xticks(range(2000,2017),fontsize=10,fontname="Arial")
  labels=[]
  for i in range(2000,2017):
    labels.append('Jan \n'+str(i))
  ax.set_xticklabels(labels)
  plt.xlim([2008.5,2014.5])
  plt.clim(3,9)
  cb=plt.colorbar(orientation="horizontal",fraction=0.07)
  cb.set_ticks([3,4,5,6,7,8,9])
  cb.set_ticklabels(['3','4','5','6','7','8','9'])
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  cb.set_label('Glacier speed (km/yr)',fontsize=11,fontname="Arial")
  plt.ylabel('Distance along flowline (km)',fontsize=11,fontname="Arial")
  ax.tick_params('both', length=6, width=1.25, which='major')
  ax.tick_params('both', length=3, width=1, which='minor')
  plt.ylim([5,-40])
  plt.yticks(fontsize=10,fontname="Arial")
  #plt.text(2009,50,"a",color='w',fontsize=22,fontweight="bold")

  # Plot percentage
  masked_array = np.ma.array (flowline_vintper, mask=np.isnan(flowline_vintper))
  cmap = matplotlib.cm.bwr
  cmap.set_bad('0.65',1.)

  plt.subplot(122)
  ax=plt.gca()
  plt.imshow(flowline_vintper/1e3,extent=[flowline_tint[0],flowline_tint[-1],(dists[-1])/1e3,(dists[0])/1e3],interpolation='nearest',aspect='auto',cmap=cmap)
  plt.fill_between(terminus_time,(terminus_val)/1e3,5,color='w')
  plt.plot(terminus_time,(terminus_val)/1e3,'ko--',markersize=3,linewidth=1)
  ax.get_xaxis().get_major_formatter().set_useOffset(False)
  plt.xticks(range(2000,2017),fontsize=10,fontname="Arial")
  plt.yticks(fontsize=10,fontname="Arial")
  labels=[]
  for i in range(2000,2017):
    labels.append('Jan \n'+str(i))
  ax.set_xticklabels(labels)
  plt.xlim([2008.5,2014.5])
  plt.clim(-0.5,0.5)
  cb=plt.colorbar(orientation="horizontal",fraction=0.07)
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  cb.set_ticks([-0.4,-0.2,0,0.2,0.4])
  cb.set_label('Change from average (km/yr)',fontsize=11,fontname="Arial")
  ax.set_yticklabels([])
  plt.ylim([5,-40])
  ax.tick_params('both', length=6, width=1.25, which='major')
  ax.tick_params('both', length=3, width=1, which='minor')

  plt.subplots_adjust(hspace=0,wspace=-0.1) 
  plt.tight_layout()

  plt.savefig(os.path.join(os.getenv("HOME"),'Bigtmp/'+glacier+'_velocity_radargram_lines_'+str(time1)+'to'+str(time2)+'.pdf'),FORMAT='PDF')
  plt.close()

######################
# Plot bed, terminus #
######################

if plot_bed == 1:

  # Get morlighem bed
  zb_morlighem = bedlib.morlighem_pts(x,y,verticaldatum='geoid')
  if glacier == 'Helheim':
    ind = np.where(x > 309900)[0]
    zb_morlighem[ind] = 'NaN'
  elif glacier == 'Kanger':
    ind = np.where(x > 490970)[0]
    zb_morlighem[ind] = 'NaN'

  # Plot 1: bed elevation along flowline. Only plot CreSIS data if it is within
  # 500 m of the flowline.

  # Set up plot
  fig = plt.figure(figsize=(5,3))
  ax1 = fig.add_subplot(111)
  ax2 = ax1.twinx()
  ax1.set_zorder(ax2.get_zorder()+1)
  ax1.patch.set_visible(False)

  # Plot morlighem bed
  ax1.plot(dists/1e3,zb_morlighem,label='Morlighem',color='k',linewidth=1.5)

  # Plot Cresis data
  colors=['b','c','g','y','orange','r',glacier[0]]
  if glacier == "Helheim":
    years = ['2001','2006b','2008b','2009','2012','2013','2014']
  elif glacier == "Kanger":
    years = ['2008','2009a','2009b','2012','2013','2014']
  for k in range(0,len(years)):
    bedpts = bedlib.cresis(years[k],glacier,'geoid')
    minind = np.argmin(abs(x-np.min(bedpts[:,0])))
    maxind = np.argmin(abs(x-np.max(bedpts[:,0])))
    ind = range(minind,maxind)
    bed_interp = np.zeros(len(ind))
    bed_interp[:] = 'NaN'
    for i in range(0,len(ind)):
      j = ind[i]
      cind = np.argmin((x[j]-bedpts[:,0])**2+(y[j]-bedpts[:,1])**2)
      if np.sqrt((x[j]-bedpts[cind,0])**2+(y[j]-bedpts[cind,1])**2) < 1500.0:
        bed_interp[i] = bedpts[cind,2]
    ax1.plot((dists[minind:maxind])/1e3,bed_interp,color=colors[k],label=years[k])

  # Set up axes
  ax1.set_xlabel('Distance from terminus (km)',fontsize=10)
  ax1.set_ylabel('Elevation (m asl)',fontsize=10)
  ax1.tick_params(axis='x', labelsize=8)
  ax1.tick_params(axis='y', labelsize=8)
  ax1.set_xlim([-10,5])
  if glacier == 'Helheim':
    ax1.set_ylim([-900,-300])
    ax1.legend(bbox_to_anchor=(0.55,1.0),fontsize=10,ncol=2,labelspacing=0.1,columnspacing=0.2,handletextpad=0.05)
  elif glacier == 'Kanger':
    ax1.set_ylim([-1300,-300])
    ax1.legend(loc=2,fontsize=10,ncol=2,labelspacing=0.1,columnspacing=0.2,handletextpad=0.05)

  # Plot terminus positions 
  ax2.plot(terminus_val/1e3,terminus_time,'.-',color='0.6',markersize=7)
  ind = np.where(calvingstyle[:,1] == 'Tabular')[0]
  for i in ind:
    termpos = terminus_val[np.argmin(abs(float(calvingstyle[i,0])-terminus_time))]
    ax2.plot(termpos/1e3,float(calvingstyle[i,0]),'.',color=[0.6,0.6,1],markersize=7,lw=0)
  ind = np.where(calvingstyle[:,1] == 'Domino')[0]
  for i in ind:
    termpos = terminus_val[np.argmin(abs(float(calvingstyle[i,0])-terminus_time))]
    ax2.plot(termpos/1e3,float(calvingstyle[i,0]),'.',color=[1,0.6,0.6],markersize=7,lw=0)

  ax2.set_ylabel('Time',fontsize=10,color='0.6')
  for tl in ax2.get_yticklabels():
    tl.set_color('0.6')
  ax2.get_yaxis().get_major_formatter().set_useOffset(False)
  ax2.tick_params(axis='y', labelsize=8)
  ax2.set_xlim([-10,5])
  plt.tight_layout()

  # Save figure
  plt.savefig(os.path.join(os.getenv("HOME"),'Bigtmp/'+glacier+'_bed_'+str(time1)+'to'+str(time2)+'.pdf'),FORMAT='PDF')
  plt.close()

  # Plot 2: locations of bed profiles

  # Set up plot
  plt.figure(figsize=(4,4))
  image_extent = [np.min(image[0]),np.max(image[0]),np.max(image[1][np.nonzero(image[1])]),np.min(image[1])]
  plt.imshow(image[2],extent=image_extent,cmap='Greys_r')
  plt.gca().invert_yaxis()
  plt.axis('equal')
  plt.xlim([min(velocitypoints[:,0])-5.0e3,max(velocitypoints[:,0])+5.0e3])
  plt.ylim([min(velocitypoints[:,1])-5.0e3,max(velocitypoints[:,1])+5.0e3])
  ax = plt.gca()

  # Plot bed profiles
  plt.plot(x,y,color='k',linewidth=1.5)
  for i in range(0,len(years)):
    bedpts = bedlib.cresis(years[i],glacier,'geoid')
    plt.plot(bedpts[:,0],bedpts[:,1],label=years[i],color=colors[i])
  plt.legend(fontsize=10,ncol=2)
  ax.set_xticklabels([])
  ax.set_yticklabels([])

  plt.savefig(os.path.join(os.getenv("HOME"),'Bigtmp/'+glacier+'_bed_vel_map.pdf'),FORMAT='PDF')
  plt.close()
# This file compares Helheim's terminus position to its velocity.

import os
import sys
import numpy as np
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"Util/Modules"))
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"BigThreeGlaciers/Tools"))
import velocity, icefronts, bed, glacier_flowline, elevation, fluxgate, flotation
import matplotlib.pyplot as plt
import matplotlib, geotiff, fracyear, dem_shading, icemask
from matplotlib.ticker import AutoMinorLocator
import scipy.signal as signal
import pylab
import matplotlib.cm as cmx
import matplotlib.colors as colors

##########
# Inputs #
##########

# Get arguments
args = sys.argv
glacier = args[1][:] # Options: Kanger, Helheim

# Locations for velocities
dists_eul = -1.0*np.array([2.0,5.0,10.0,20.0,30.0]) # kilometers
dists_mel = np.array([0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0])
#dists_lag = np.array([1.0,5.0,10.0,20.0,30.0])

# Image for plotting
if glacier == "Helheim":
  #image = geotiff.read(os.path.join(os.getenv("HOME"),"Data/Mosaics/Helheim/mosaicHelheim.2014-159.148.38725_1-20mgeo.tif"))
  image = geotiff.readrgb(os.path.join(os.getenv("DATA_HOME"),"Imagery/ASTER/Helheim/20130812141113_AST2125988263_2.tif"))
elif glacier == "Kanger":
  image = geotiff.read(os.path.join(os.getenv("DATA_HOME"),"Mosaics/Kanger/mosaicKang.2014-160.163.38740_1-20mgeo.tif"))
    
################
# Plot options #
################

time1 = 2008.4 #start time for plot
time2 = 2015.75 # end time for plot
seasonal = 1 # plot seasonal bars, to highlight seasonal trends
normalized = 0
lagrangian = 0

# Which plots?
plot_overview = 1
plot_radargram = 0
plot_bed = 1
plot_images = 1

############ 
# Flowline #
############
if glacier == 'Helheim':
  x,y,zb,dists = glacier_flowline.load(glacier,shapefilename='flowline_flightline')
else:
  x,y,zb,dists = glacier_flowline.load(glacier)

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

########################
# Get calving behavior #
########################

calvingstyle = icefronts.calving(glacier)

############################
# Get velocities at points # 
############################

# Find where points are located along flowline
ind_eul=[]
for i in range(0,len(dists_eul)):
  ind_eul.append( (abs(dists - dists_eul[i]*1e3)).argmin() )

# Load velocities for glacier and ice melange
if lagrangian == 1:
  vel_val,vel_time,vel_error,vel_dists,vel_x,vel_y = velocity.velocity_at_lagpoints(x,y,dists_eul*1e3,glacier)
else:
  vel_val,vel_time,vel_error = velocity.velocity_at_eulpoints(x[ind_eul],y[ind_eul],glacier)

velmel_val,velmel_time,velmel_error,velmel_dists,velmel_x,velmel_y = velocity.velocity_at_lagpoints(x,y,dists_mel*1e3,glacier)
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
vel_val[interped < dists_eul[0]*1.0e3,0] = 'NaN'
del interped

###############################################
# Get thinning rates inferred from flux gates #
###############################################

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

xwv,ywv,zwv,timewv = elevation.worldview_grid(glacier,xmin,xmax,ymin,ymax,years='all',verticaldatum='ellipsoid')
dH_time,dH = fluxgate.dem_thinning(glacier,xwv,ywv,zwv,timewv,"fluxgate3",type='absolute')

zpt_atm,time_atm = elevation.atm_at_pts(x[ind_eul],y[ind_eul],glacier,years='all',maxdist=250,verticaldatum='geoid',method='average',cutoff='terminus')
zpt_wv,zpt_std,time_wv = elevation.worldview_at_pts(x[ind_eul],y[ind_eul],glacier,years='all',verticaldatum='geoid',cutoff='terminus',method='average',radius=250)


###############
# Plot images #
###############

if plot_images == 1:
  
  DIRLANDSAT = os.path.join(os.getenv("DATA_HOME"),"Imagery/Landsat/"+glacier+"/TIF/")
  DIRTSX = os.path.join(os.getenv("DATA_HOME"),"Mosaics/"+glacier+"/")
  
  if glacier == 'Helheim':
    xmin = 304000.0
    xmax = 314000.0
    ymin = -2582500.0
    ymax = -2572500.0
    imagefiles = ['mosaicHelheim.2013-216.148.34049_1-20mgeo.tif',
             'mosaicHelheim.2013-304.148.35385_1-20mgeo.tif',
             'mosaicHelheim.2014-027.148.36721_1-20mgeo.tif',
             '20140330140600_LC82330132014089LGN00.tif',
             '20140729135956_LC82320142014210LGN00.tif'] 
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
      images.append(geotiff.read(DIRLANDSAT+file))
      year,month,day = [int(file[0:4]),int(file[4:6]),int(file[6:8])]
      images_type.append('Landsat-8')
    else:
      images.append(geotiff.read(DIRTSX+file))
      if glacier == 'Helheim':
        year,month,day = fracyear.doy_to_date(float(file[14:18]),float(file[19:22]))
      elif glacier == 'Kanger':
        year,month,day = fracyear.doy_to_date(float(file[11:15]),float(file[16:19]))
      images_type.append('TSX')
    images_time.append([year,month,day,fracyear.date_to_fracyear(year,month,day)])
  
  N = len(images) # number of images
  images_labels = ['a','b','c','d','e']
  
  plt.figure(figsize=(5,4))
  gs = matplotlib.gridspec.GridSpec(4,N+1,width_ratios=np.append(np.ones(N),0.1))
  gs.update(left=0.02,wspace=0.02,hspace=0.02,top=0.98,right=0.9)
  for i in range(0,N):
  
    xdem,ydem,zdem,demtime = elevation.grid_near_time(images_time[i][3],glacier)
    xvel,yvel,uvel,vvel,vvel,vtime = velocity.tsx_near_time(images_time[i][3],glacier)
    #if i == 0:
      #xmask,ymask,mask = icemask.load(glacier,np.min(xvel),np.max(xvel),np.min(yvel),np.max(yvel),100)
    xf,yf,zabovefloat = flotation.extent(xdem,ydem,zdem,demtime,glacier,rho_i=917.0,rho_sw=1020.0,bedsource='cresis',verticaldatum='geoid')
  
    for j in [0,1,3]:
      ax = plt.subplot(gs[j,i])
      plt.imshow(images[i][2],extent=[np.min(images[i][0]),np.max(images[i][0]),np.min(images[i][1]),np.max(images[i][1])],origin='lower',cmap='Greys_r')
      plt.clim(clims[i])
      plt.xlim([xmin,xmax])
      plt.ylim([ymin,ymax])
      plt.xticks([])
      plt.yticks([])
      
    ax = plt.subplot(gs[0,i])
    year,month,day = fracyear.fracyear_to_date(images_time[i][3])
    plt.text(xmin+500,ymax-1.25e3,str(year)+'-'+str(month)+'-'+str(int(np.floor(day))),backgroundcolor='w',fontsize=8)
    plt.text(xmin+500,ymin+1e3,images_type[i],backgroundcolor='w',fontsize=8)
    plt.text(xmin+500,ymax-3e3,images_labels[i],fontsize=9,fontweight='bold')
    
    ax = plt.subplot(gs[1,i])
    #vvel = np.ma.masked_array(vvel,mask)
    year,month,day = fracyear.fracyear_to_date(vtime)
    xfront,yfront,ftime = icefronts.near_time(images_time[i][3],glacier)
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
      cb.set_label('Speed (km/yr)',fontsize=8)
    
    ax = plt.subplot(gs[2,i])
    shadeddem = dem_shading.set_shade(zdem,0,220)
    im = plt.imshow(shadeddem,extent=[np.min(xdem),np.max(xdem),np.min(ydem),np.max(ydem)],origin='lower')
    year,month,day = fracyear.fracyear_to_date(demtime)
    plt.text(xmin+500,ymax-1.25e3,str(year)+'-'+str(month)+'-'+str(int(np.floor(day))),backgroundcolor='w',fontsize=8)
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
      cb.set_label('Elevation (m)',fontsize=8)
    
    ax = plt.subplot(gs[3,i])
    im=plt.scatter(xf,yf,c=zabovefloat,lw=0,vmin=-25,vmax=25,cmap='RdBu_r',s=6)
    year,month,day = fracyear.fracyear_to_date(demtime)
    plt.text(xmin+500,ymax-1.25e3,str(year)+'-'+str(month)+'-'+str(int(np.floor(day))),backgroundcolor='w',fontsize=8)
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    plt.xticks([])
    plt.yticks([])
    
    if i == N-1:
      ax = plt.subplot(gs[3,N])
      cb = plt.colorbar(im,cax=ax)
      cb.set_ticks(np.arange(-20,40,20))
      cb.ax.tick_params(labelsize=8)
      cb.set_label('Height above \n flotation (m)',fontsize=8)
  
  plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_images.png"),FORMAT='PNG',dpi=800)
  plt.close()
  
  del clims,year,month,day,imagefiles,xmin,xmax,ymin,ymax,xvel,yvel,vvel,vtime,DIRTSX,DIRLANDSAT

###########################################
# Plot velocity vs. terminus through time #
###########################################

if plot_overview == 1:
  plt.figure(figsize=(6.5,6.5))
  gs = matplotlib.gridspec.GridSpec(6,1)

  # Plot terminus
  plt.subplot(gs[-2, :])
  ax = plt.gca()
  ind = np.where(calvingstyle[:,1] == 'Tabular')[0]
  plt.ylim([-4,4])
  if seasonal:
    xTickPos = np.linspace(np.floor(time1)-0.25,np.ceil(time2)-0.25,(np.ceil(time2)-np.floor(time1))*2+1)
    ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.8','w'],linewidth=0)
  for i in ind:
    if i == ind[1]:
      plt.bar(float(calvingstyle[i,0])-0.03,8.0,width=0.03,color=[0.6,0.6,1],edgecolor='none',label='Tabular',bottom=-4)
    elif i !=0:
      plt.bar(float(calvingstyle[i,0])-0.03,8.0,width=0.03,color=[0.6,0.6,1],edgecolor='none',bottom=-4)
  ind = np.where((calvingstyle[:,1] == 'Mixed'))[0]
  for i in ind:
    if i == ind[1]:
      plt.bar(float(calvingstyle[i,0])-0.03,8.0,width=0.03,color=[0.4,0.8,0.6],edgecolor='none',label='Mixed',bottom=-4)
    elif i != 0:
      plt.bar(float(calvingstyle[i,0])-0.03,8.0,width=0.03,color=[0.4,0.8,0.6],edgecolor='none',bottom=-4)
  ind = np.where((calvingstyle[:,1] == 'Domino'))[0]
  for i in ind:
    if i == ind[1]:
      plt.bar(float(calvingstyle[i,0])-0.03,8.0,width=0.03,color=[1,0.6,0.6],edgecolor='none',label='Nontabular',bottom=-4)
    elif i != 0:
      plt.bar(float(calvingstyle[i,0])-0.03,8.0,width=0.03,color=[1,0.6,0.6],edgecolor='none',bottom=-4)
  if plot_images == 1:
   for i in range(0,len(images)):
      plt.plot([images_time[i][3],images_time[i][3]],ax.get_ylim(),'--',color='0.3')
  plt.legend(loc=2,fontsize=8,numpoints=1,handlelength=0.4,labelspacing=0.05,ncol=3,columnspacing=0.7,handletextpad=0.2)
  nonnan = np.where(~(np.isnan(terminus_val)))[0]
  plt.plot(terminus_time[nonnan],terminus_val[nonnan]/1e3,'ko',linewidth=1,markersize=2)
  x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
  ax.xaxis.set_major_formatter(x_formatter)
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  ax.yaxis.set_minor_locator(AutoMinorLocator(2))
  plt.xticks(range(2000,2016),fontsize=8,fontname="Arial")
  plt.xlim([time1,time2])
  plt.yticks(np.arange(-6,8,2),fontsize=8,fontname="Arial")
  plt.ylabel('Terminus \n (km)',fontsize=8,fontname="Arial")
  ax.tick_params('both', length=8, width=1.5, which='major')
  ax.tick_params('both', length=4, width=1, which='minor')
  if glacier == 'Helheim':
    plt.ylim([-2.2,2.2])
  else:  
    plt.ylim(np.floor((np.min(terminus_val))/1e3),np.ceil((np.max(terminus_val))/1e3))


  # Plot velocities
  ax = plt.subplot(gs[0:-3, :]) 
  #plt.plot([2000,2014],[0,0],'k')
  coloptions=['b','c','g','y','r']
  if plot_images == 1:
   for i in range(0,len(images)):
      plt.plot([images_time[i][3],images_time[i][3]],[0,12],'--',color='0.3')
      plt.text(images_time[i][3]+0.05,0.3,images_labels[i],fontsize=9,fontweight='bold')
  if normalized == 1:
    for i in range(0,len(dists_eul)):
      nonnan = np.where(~(np.isnan(vel_val[:,i])))[0]
      plt.plot(vel_time[nonnan],(vel_val[nonnan,i]-vel_val[nonnan[0],i])/1e3,'o',color=coloptions[i],label=str(dists_eul[i])+' km',linewidth=1,markersize=3)
    plt.yticks(range(-3,3,1),fontsize=8,fontname='Arial')
    if glacier == 'Kanger':
      plt.ylim([-1.5,3])
  else:
    for i in range(0,len(dists_eul)):
      nonnan = np.where(~(np.isnan(vel_val[:,i])))[0]
      plt.plot(vel_time[nonnan],(vel_val[nonnan,i])/1e3,'o',color=coloptions[i],label=str(dists_eul[i])+' km',linewidth=1,markersize=3)
    plt.yticks(range(0,12),fontsize=8,fontname="Arial")
    plt.ylim([0,12])
  plt.ylabel('Glacier speed \n (km/yr)',fontsize=8,fontname="Arial")
  x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
  ax.xaxis.set_major_formatter(x_formatter)
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  ax.yaxis.set_minor_locator(AutoMinorLocator(2))
  plt.legend(loc=2,fontsize=8,numpoints=1,handlelength=0.4,labelspacing=0.1,ncol=2,columnspacing=0.5)
  matplotlib.rc('font',family="Arial",)
  ax.tick_params('both', length=8, width=1.5, which='major')
  ax.tick_params('both', length=4, width=1, which='minor')
  labels=[]
  plt.xticks(range(2000,2017),fontsize=8,fontname="Arial")
  ax.set_xticklabels([])
  if seasonal:
    xTickPos = np.linspace(np.floor(time1)-0.25,np.ceil(time2)-0.25,(np.ceil(time2)-np.floor(time1))*2+1)
    ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.8','w'],linewidth=0)
  plt.xlim([time1,time2])


  # Plot thinning rates from velocities
  plt.subplot(gs[-3, :])
  floatheight=flotation.height(zb[ind_eul])
  ax = plt.gca()
  plt.ylim([-100,100])
  if plot_images == 1:
   for i in range(0,len(images)):
      plt.plot([images_time[i][3],images_time[i][3]],ax.get_ylim(),'--',color='0.3')
  if seasonal:
    xTickPos = np.linspace(np.floor(time1)-0.25,np.ceil(time2)-0.25,(np.ceil(time2)-np.floor(time1))*2+1)
    ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.8','w'],linewidth=0)
    #plt.plot(dH_time,dH[:,0]-dH[0,0],marker='o',linestyle='none',color='r',markersize=3)
  if glacier == 'Helheim':
    plt.plot([time1,time2],[floatheight[0]-zpt_wv[0,0],floatheight[0]-zpt_wv[0,0]],'b:',linewidth=1)
    plt.plot(time_wv,zpt_wv[:,0]-zpt_wv[0,0],'bo',markersize=3)
    plt.plot(time_atm,zpt_atm[:,0]-zpt_wv[0,0],'b+',markersize=4)
    plt.plot(time_wv,zpt_wv[:,1]-zpt_wv[0,1],'co',markersize=3)
    plt.plot(time_atm,zpt_atm[:,1]-zpt_wv[0,1],'c+',markersize=4)
  elif glacier == 'Kanger':
    plt.plot([time1,time2],[floatheight[0]-zpt_atm[4,0],floatheight[0]-zpt_atm[4,0]],'b:',linewidth=1)
    plt.plot(time_wv,zpt_wv[:,0]-zpt_atm[4,0],'bo',markersize=3)
    plt.plot(time_atm,zpt_atm[:,0]-zpt_atm[4,0],'b+',markersize=4)
    plt.plot(time_wv,zpt_wv[:,1]-zpt_atm[4,1],'co',markersize=3)
    plt.plot(time_atm,zpt_atm[:,1]-zpt_atm[4,1],'c+',markersize=4)
  x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
  ax.xaxis.set_major_formatter(x_formatter)
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  ax.yaxis.set_minor_locator(AutoMinorLocator(2))
  ax.set_xticklabels([])
  plt.xticks(range(2000,2016),fontsize=8,fontname="Arial")
  plt.xlim([time1,time2])
  plt.ylabel('Relative elevation (m)',fontsize=8,fontname="Arial")
  ax.tick_params('both', length=8, width=1.5, which='major')
  ax.tick_params('both', length=4, width=1, which='minor')
  if glacier == 'Helheim':
    plt.ylim([-25,5])
    plt.yticks(np.arange(-20,20,10),fontsize=8,fontname="Arial")
  elif glacier == 'Kanger':
    plt.yticks(np.arange(-40,20,20),fontsize=8,fontname="Arial")
    plt.ylim([-50,10])

  # Plot presence of rigid melange
  plt.subplot(gs[-1, :])
  rigid = np.zeros(len(velmel_val[:,0]))
  for i in range(0,len(velmel_val[:,0])):
    nonnan = np.where(~np.isnan(velmel_val[i,:]))[0]
    if len(nonnan) == 1:
      rigid[i] = dists_mel[0]
    elif len(nonnan) > 1:
      ind = np.where(np.diff(nonnan)==1)[0][-1]+1
      rigid[i] = dists_mel[ind]
    else:
      rigid[i] == 0
  plt.ylabel('Rigid melange \n (km)',fontsize=8,fontname="Arial")
  ax = plt.gca()
  plt.ylim([0,5])
  plt.yticks(np.arange(0,5,2),fontsize=8,fontname='arial')
  x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
  ax.xaxis.set_major_formatter(x_formatter)
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  ax.yaxis.set_minor_locator(AutoMinorLocator(2))
  matplotlib.rc('font',family="Arial",)
  ax.tick_params('both', length=8, width=1.5, which='major')
  ax.tick_params('both', length=4, width=1, which='minor')
  labels=[]
  plt.xticks(range(2000,2017),fontsize=8,fontname="Arial")
  ax.set_xticklabels([])
  #if seasonal:
  #  xTickPos = np.linspace(np.floor(time1)-0.25,np.floor(time2)-0.25,(np.ceil(time2)-np.ceil(time1))*2+1)
  #  ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.8','w'],linewidth=0)
  for i in range(2000,2017):
    labels.append('Jan \n'+str(i))
  plt.xticks(range(2000,2017))
  ax.set_xticklabels(labels,fontsize=8,fontname='Arial')
  plt.xlim([time1,time2])
  plt.bar(time1,5,width=time2-time1,color='0.9',edgecolor='k',hatch=4 * "\\")
  plt.bar(velmel_time[:,0]-velmel_time[:,1]/2,np.ones(len(velmel_time[:,0]))*5,width=velmel_time[:,1],color='w',edgecolor='w')
  plt.bar(velmel_time[rigid!=0,0]-velmel_time[rigid!=0,1]/2,rigid[rigid!=0],width=velmel_time[rigid!=0,1],color='k')
  if plot_images == 1:
   for i in range(0,len(images)):
      plt.plot([images_time[i][3],images_time[i][3]],ax.get_ylim(),'--',color='0.3')

  plt.subplots_adjust(hspace=0.05,wspace=0) 
  plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_vel_time_"+str(time1)+"to"+str(time2)+".pdf"),FORMAT='PDF')
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
  for i in range(0,len(velocitypoints)):
    plt.plot(velocitypoints[i,0],velocitypoints[i,1],'o',color=coloptions[i],markersize=7)
  ax.set_xticklabels([])
  ax.set_yticklabels([])

  plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_velocity_map_"+str(time1)+"to"+str(time2)+".pdf"),FORMAT='PDF')
  plt.close()
   
###########################
# Plot velocity radargram #
###########################

if plot_radargram == 1:

  flowline_v,flowline_t,termini = velocity.velocity_along_flowline(x,y,glacier,cutoff='terminus',data='TSX')
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
  plt.xticks(range(2000,2017),fontsize=8,fontname="Arial")
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
  cb.set_label('Glacier speed (km/yr)',fontsize=10,fontname="Arial")
  plt.ylabel('Distance along flowline (km)',fontsize=10,fontname="Arial")
  ax.tick_params('both', length=8, width=1.5, which='major')
  ax.tick_params('both', length=4, width=1, which='minor')
  plt.ylim([5,-40])
  plt.yticks(fontsize=8,fontname="Arial")
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
  ax.tick_params('both', length=8, width=1.5, which='major')
  ax.tick_params('both', length=4, width=1, which='minor')

  plt.subplots_adjust(hspace=0,wspace=-0.1) 
  plt.tight_layout()

  plt.savefig(os.path.join(os.getenv("HOME"),'Bigtmp/'+glacier+'_velocity_radargram_lines_'+str(time1)+'to'+str(time2)+'.pdf'),FORMAT='PDF')
  plt.close()

######################
# Plot bed, terminus #
######################

if plot_bed == 1:

  # Get morlighem bed
  zb_morlighem = bed.morlighem_pts(x,y,verticaldatum='geoid')
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
  colors=['b','c','g','y','orange','r','m']
  if glacier == "Helheim":
    years = ['2001','2006b','2008b','2009','2012','2013','2014']
  elif glacier == "Kanger":
    years = ['2008','2009a','2009b','2012','2013','2014']
  for k in range(0,len(years)):
    bedpts = bed.cresis(years[k],glacier,'geoid')
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
  ax1.set_xlabel('Distance from terminus (km)',fontsize=9)
  ax1.set_ylabel('Elevation (m asl)',fontsize=9)
  ax1.tick_params(axis='x', labelsize=8)
  ax1.tick_params(axis='y', labelsize=8)
  ax1.set_xlim([-10,5])
  if glacier == 'Helheim':
    ax1.set_ylim([-900,-300])
    ax1.legend(bbox_to_anchor=(0.55,1.0),fontsize=9,ncol=2,labelspacing=0.1,columnspacing=0.2,handletextpad=0.05)
  elif glacier == 'Kanger':
    ax1.set_ylim([-1300,-300])
    ax1.legend(loc=2,fontsize=9,ncol=2,labelspacing=0.1,columnspacing=0.2,handletextpad=0.05)

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

  ax2.set_ylabel('Time',fontsize=9,color='0.6')
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
    bedpts = bed.cresis(years[i],glacier,'geoid')
    plt.plot(bedpts[:,0],bedpts[:,1],label=years[i],color=colors[i])
  plt.legend(fontsize=8,ncol=2)
  ax.set_xticklabels([])
  ax.set_yticklabels([])

  plt.savefig(os.path.join(os.getenv("HOME"),'Bigtmp/'+glacier+'_bed_vel_map.pdf'),FORMAT='PDF')
  plt.close()
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
vel_val_Howat_H, vel_time_Howat_H, vel_error_Howat_H = vellib.howat_optical_at_pts(x_H[ind_eul_H],y_H[ind_eul_H],'Helheim')
vel_val_Howat_K, vel_time_Howat_K, vel_error_Howat_K = vellib.howat_optical_at_pts(x_K[ind_eul_K],y_K[ind_eul_K],'Kanger')

# Load elevations
zpt_atm_H,zptstd_atm_H,time_atm_H = zslib.atm_at_pts(x_H[ind_eul_H],y_H[ind_eul_H],'Helheim',years='all',maxdist=100.,verticaldatum='geoid',method='average',cutoff='terminus')
zpt_dem_H,zpterror_dem_H,time_dem_H = zslib.dem_at_pts(x_H[ind_eul_H],y_H[ind_eul_H],'Helheim',years='all',verticaldatum='geoid',cutoff='terminus',method='average',radius=100.)
zpt_atm_K,zptstd_atm_K,time_atm_K = zslib.atm_at_pts(x_K[ind_eul_K],y_K[ind_eul_K],'Kanger',years='all',maxdist=100.,verticaldatum='geoid',method='average',cutoff='terminus')
zpt_dem_K,zpterror_dem_K,time_dem_K = zslib.dem_at_pts(x_K[ind_eul_K],y_K[ind_eul_K],'Kanger',years='all',verticaldatum='geoid',cutoff='terminus',method='average',radius=100.)

# Color options for plotting
coloptions=['r','b','g','limegreen','gold']


###############
# Main timeseries figure #
###############

for glacier in ['Kanger','Helheim']:
  if glacier == 'Helheim':
    x = x_H; y = y_H; zb = zb_H; dists = dists_H
    vel_val = vel_val_H; vel_time = vel_time_H
    vel_val_Howat = vel_val_Howat_H; vel_time_Howat = vel_time_Howat_H
    terminus_val = terminus_val_H; terminus_time = terminus_time_H
    time_dem = time_dem_H; zpt_dem = zpt_dem_H; zpterror_dem = zpterror_dem_H
    time_atm = time_atm_H; zpt_atm = zpt_atm_H
    ind_eul = ind_eul_H
  elif glacier == 'Kanger':  
    x = x_K; y = y_K; zb = zb_K; dists = dists_K
    vel_val = vel_val_K; vel_time = vel_time_K
    vel_val_Howat = vel_val_Howat_K; vel_time_Howat = vel_time_Howat_K
    terminus_val = terminus_val_K; terminus_time = terminus_time_K
    time_dem = time_dem_K; zpt_dem = zpt_dem_K; zpterror_dem = zpterror_dem_K
    time_atm = time_atm_K; zpt_atm = zpt_atm_K
    ind_eul = ind_eul_K
    
  plt.figure(figsize=(8,8))
  time1 = 2000.0; time2 = 2016.75; 
  years = np.arange(np.floor(time1),np.ceil(time2)+1) 
  gs = matplotlib.gridspec.GridSpec(6,1)
  matplotlib.rc('font',family='Arial')

  # Plot terminus
  plt.subplot(gs[0, :])
  ax = plt.gca()
  plt.plot(terminus_time,terminus_val/1e3,'ko',markersize=5)
  x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
  ax.xaxis.set_major_formatter(x_formatter)
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  plt.xticks(range(2000,2017),fontsize=14,fontname="Arial")
  ax.set_xticklabels([])
  plt.xlim([time1,time2])
  plt.yticks(np.arange(-6,8,3),fontsize=14,fontname="Arial")
  plt.ylabel('Terminus \n (km)',fontsize=14,fontname="Arial")
  ax.tick_params('both', length=6, width=1.25, which='major')
  ax.tick_params('x', length=6, width=1.25, which='minor')
  plt.ylim([-5,6])

  # Plot velocities
  ax = plt.subplot(gs[1:4, :]) 
  coloptions=['limegreen','r','b','g','gold']
  x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
  ax.xaxis.set_major_formatter(x_formatter)
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  ax.yaxis.set_minor_locator(AutoMinorLocator(2))
  matplotlib.rc('font',family="Arial",)
  ax.tick_params('both', length=6, width=1.25, which='major')
  ax.tick_params('both', length=6, width=1.25, which='minor')
  for i in range(1,len(dists_eul)):
    nonnan = np.where(~(np.isnan(vel_val[:,i])))[0]
    plt.plot(vel_time[nonnan],(vel_val[nonnan,i])/1e3,'ko',markerfacecolor=coloptions[i],label=glacier[0]+'{0:02d}'.format(int(abs(dists_eul[i]))),markersize=6)
    plt.plot(vel_time_Howat[vel_time_Howat[:,0] < 2009, 0], vel_val_Howat[vel_time_Howat[:,0] < 2009,i]/1e3,'ko',markerfacecolor=coloptions[i], markersize=6)
  plt.legend(loc=2,borderpad=0.3,fontsize=14,numpoints=1,handlelength=0.2,handletextpad=0.5,labelspacing=0.1,ncol=2,columnspacing=0.8)
  plt.yticks(np.arange(2,16,2),fontsize=14,fontname="Arial")
  plt.ylabel('Glacier speed \n (km/yr)',fontsize=14,fontname="Arial")
  plt.xticks(np.arange(2000,2017,2),fontsize=14,fontname="Arial")
  ax.set_xticklabels([])
  plt.xlim([time1,time2])
  if glacier == 'Helheim':
    plt.ylim([3,11])
  elif glacier == 'Kanger':
    plt.ylim([1.5,14.5])

  # Plot surface elevations
  plt.subplot(gs[4:, :])
  ax = plt.gca()
  if glacier == 'Helheim':
    plt.plot(time_atm,zpt_atm[:,1],'ko',markerfacecolor=coloptions[1],markersize=6)    
    plt.plot(time_dem,zpt_dem[:,1],'ko',markerfacecolor=coloptions[1],markersize=6)
  elif glacier == 'Kanger':
    plt.plot(time_dem,zpt_dem[:,1],'ko',markerfacecolor=coloptions[1],markersize=6)
    plt.plot(time_atm,zpt_atm[:,1],'ko',markerfacecolor=coloptions[1],markersize=6)
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
  ax.set_xticklabels(labels,fontsize=14,fontname='Arial')
  plt.xlim([time1,time2])
  if glacier == 'Helheim':
    plt.ylabel('Surface elevation \n (m)',fontsize=14,fontname="Arial")
    plt.yticks(np.arange(0,350,50),fontsize=14,fontname="Arial")
    plt.ylim([70,290])
  elif glacier == 'Kanger':
    plt.ylabel('Surface elevation \n (m)',fontsize=14,fontname="Arial")
    plt.yticks(np.arange(50,400,50),fontsize=14,fontname="Arial")
    plt.ylim([90,340])

  plt.tight_layout()
  plt.subplots_adjust(hspace=0.0,wspace=0,top=0.98,right=0.98,left=0.12,bottom=0.08) 
  plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Defense_"+glacier+"_all.pdf"),FORMAT='PDF',dpi=400)
  plt.close()

del i, ax, labels, time1, time2

#####################
# Make overview map #
#####################

coast=np.loadtxt(os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/Coastlines/greenland_coast_wgs84.txt"))
xvel,yvel,vel=geotifflib.read(os.path.join(os.getenv("DATA_HOME"),"Velocity/Random/Greenland/LandsatAll/landsatall_v_wgs84.tif"))

# Set up basemap
fig=plt.figure(figsize=(3.3,6))
ax = plt.gca()
ax.tick_params('font',family='Arial',color='k',fontsize=12)
m = Basemap(width=1600000,height=2800000,
            resolution='l',projection='stere',\
            lat_ts=60,lat_0=72,lon_0=-41.)

nx = int((m.xmax-m.xmin)/1000.)+1; ny = int((m.ymax-m.ymin)/1000.)+1
velgrid = m.transform_scalar(vel,xvel,yvel,nx,ny)
xcoast,ycoast = m( coast[:,0],coast[:,1] )

clip = Path(np.column_stack([xcoast,ycoast]))
clip = PathPatch(clip,transform=ax.transData,facecolor='k')

cx = cubehelix.cmap(start=1.2,rot=-1.1,reverse=True,minLight=0.22,sat=3)
ax.add_patch(clip)
cs = m.imshow(velgrid,cmap=cx,vmin=1,vmax=1000,zorder=2,norm=matplotlib.colors.LogNorm())
cb = m.colorbar(cs,location='bottom',ticks=[1,10,100,1000],format='%d')
plt.setp(plt.getp(cb.ax.axes, 'xticklabels'), color='k')
cb.ax.xaxis.set_tick_params(labelsize=12)
cb.set_label('Glacier speed (m/yr)',fontsize=16,fontname='Arial',color='k')
cb.update_ticks()

cs.set_clip_path(clip)

m.plot(xcoast,ycoast,'k')

xg,yg = m (-33.0,68.63333)
#m.plot(xg,yg,'ko',markerfacecolor='r',markersize=10,zorder=3)

xg,yg = m([-38.3],[66.4])
#m.plot(xg,yg,'ko',markerfacecolor='r',markersize=10,zorder=3)

ax.axis('off')
plt.tight_layout()
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Defense_overviewmap.pdf"),format='PDF',dpi=400)
plt.close()

############################
# Make basemap for Helheim #
############################

imagetime = datelib.date_to_fracyear(2014,7,4)
ximage,yimage,image = geotifflib.readrgb(os.path.join(os.getenv("DATA_HOME"),"Imagery/Landsat/Helheim/TIF/20140704140535_LC82330132014185LGN00.tif"))

# Load velocity record
xvel = np.arange(np.min(ximage),np.max(ximage),100)
yvel = np.arange(np.min(yimage),np.max(yimage),100)
vx,vy = vellib.inversion_3D('Helheim',xvel,yvel,imagetime,dir_velocity_out='none',blur=False)
vel = np.sqrt(vx**2+vy**2)
del vx,vy

# Load mask
xmask,ymask,mask = masklib.load_grid('Helheim',np.min(xvel),np.max(xvel),np.min(yvel),np.max(yvel),100,icefront_time=datelib.date_to_fracyear(2014,7,4))
vel_masked = np.ma.masked_array(vel,mask)

fig = plt.figure(figsize=(2.5,2.5))
matplotlib.rc('font',family='Arial')

cx = cubehelix.cmap(start=1.2,rot=-1.1,reverse=True,minLight=0.1,sat=2)
p=plt.imshow(vel_masked/1e3,extent=[np.min(xvel),np.max(xvel),np.min(yvel),np.max(yvel)],origin='lower',clim=[0,10],cmap=cx)
ax = plt.gca()
ax.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
ax.imshow(vel_masked/1e3,extent=[np.min(xvel),np.max(xvel),np.min(yvel),np.max(yvel)],origin='lower',alpha=0.7,clim=[0,10],cmap=cx)
ax.axis('equal')

xmin = 283000.
xmax = 313000.
ymin = -2587000.
ymax = -2552000.

ax.set_xticklabels([]); ax.set_yticklabels([])
ax.set_yticks([]); ax.set_xticks([])
ax.set_ylim([ymin,ymax]); ax.set_xlim([xmin,xmax])
plt.tight_layout()

xmin,xmax = plt.xlim()
ymin,ymax = plt.ylim()
path = matplotlib.path.Path([[0.5*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.99*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.99*(xmax-xmin)+xmin,0.67*(ymax-ymin)+ymin],
                        [0.5*(xmax-xmin)+xmin,0.67*(ymax-ymin)+ymin],
                        [0.5*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin]])
patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1)
ax.add_patch(patch)
cbaxes = fig.add_axes([0.54, 0.91, 0.35, 0.02])
cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=[0,5,10])
ax.text(xmin+0.53*(xmax-xmin),ymin+0.78*(ymax-ymin),'Speed (km/yr)',fontsize=11,fontname='Arial')
cb.ax.tick_params(labelsize=11)
ax.plot([xmin+0.54*(xmax-xmin),xmin+0.54*(xmax-xmin)+5e3],[ymin+0.72*(ymax-ymin),ymin+0.72*(ymax-ymin)],'k',linewidth=1.5)
ax.plot([xmin+0.54*(xmax-xmin),xmin+0.54*(xmax-xmin)],[ymin+0.72*(ymax-ymin),ymin+0.70*(ymax-ymin)],'k',linewidth=1.5)
ax.plot([xmin+0.54*(xmax-xmin)+5e3,xmin+0.54*(xmax-xmin)+5e3],[ymin+0.72*(ymax-ymin),ymin+0.70*(ymax-ymin)],'k',linewidth=1.5)
ax.text(xmin+0.58*(xmax-xmin)+5e3,ymin+0.69*(ymax-ymin),'5 km',fontsize=11,fontname='Arial')

plt.subplots_adjust(wspace=0,hspace=0,top=0.98,right=0.98,left=0.02,bottom=0.02)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Defense_Helheim_map.pdf"),FORMAT='PDF',dpi=600,transparent=True)

for i in range(1,len(ind_eul_H)):
  ax.plot(x_H[ind_eul_H[i]],y_H[ind_eul_H[i]],'ko',markerfacecolor=coloptions[i],markersize=5)
  if i in [1,2,4]:
    ax.text(x_H[ind_eul_H[i]]-4e3,y_H[ind_eul_H[i]]-1.5e3,'H'+'{0:02d}'.format(int(abs(dists_eul[i]))),fontsize=9,fontname='Arial')
  else:
    ax.text(x_H[ind_eul_H[i]],y_H[ind_eul_H[i]]+0.75e3,'H'+'{0:02d}'.format(int(abs(dists_eul[i]))),fontsize=9,fontname='Arial')

plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Defense_Helheim_map_points.pdf"),FORMAT='PDF',dpi=600,transparent=True)
plt.close()

###########################
# Make basemap for Kanger #
###########################

imagetime = datelib.date_to_fracyear(2014,7,6)
ximage,yimage,image = geotifflib.readrgb(os.path.join(os.getenv("DATA_HOME"),"Imagery/Landsat/Kanger/TIF/20140706135251_LC82310122014187LGN00.tif"))

# Load velocity record
xvel = np.arange(np.min(ximage),np.max(ximage),100)
yvel = np.arange(np.min(yimage),np.max(yimage),100)
vx,vy = vellib.inversion_3D('Kanger',xvel,yvel,imagetime,dir_velocity_out='none',blur=False)
vel = np.sqrt(vx**2+vy**2)
del vx,vy

# Load mask
xmask,ymask,mask = masklib.load_grid('Kanger',np.min(xvel),np.max(xvel),np.min(yvel),np.max(yvel),100,icefront_time=datelib.date_to_fracyear(2014,7,4))
vel_masked = np.ma.masked_array(vel,mask)

fig = plt.figure(figsize=(2.5,2.5))
matplotlib.rc('font',family='Arial')

cx = cubehelix.cmap(start=1.2,rot=-1.1,reverse=True,minLight=0.1,sat=2)
p=plt.imshow(vel_masked/1e3,extent=[np.min(xvel),np.max(xvel),np.min(yvel),np.max(yvel)],origin='lower',clim=[0,10],cmap=cx)
ax = plt.gca()
ax.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
ax.imshow(vel_masked/1e3,extent=[np.min(xvel),np.max(xvel),np.min(yvel),np.max(yvel)],origin='lower',alpha=0.7,clim=[0,10],cmap=cx)
ax.axis('equal')

xmin = 468000.
xmax = 498000.
ymin = -2299000.
ymax = -2264000.

ax.set_xticklabels([]); ax.set_yticklabels([])
ax.set_yticks([]); ax.set_xticks([])
ax.set_ylim([ymin,ymax]); ax.set_xlim([xmin,xmax])
plt.tight_layout()

xmin,xmax = plt.xlim()
ymin,ymax = plt.ylim()
path = matplotlib.path.Path([[0.5*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.99*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.99*(xmax-xmin)+xmin,0.67*(ymax-ymin)+ymin],
                        [0.5*(xmax-xmin)+xmin,0.67*(ymax-ymin)+ymin],
                        [0.5*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin]])
patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1)
ax.add_patch(patch)
cbaxes = fig.add_axes([0.54, 0.91, 0.35, 0.02])
cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=[0,5,10])
ax.text(xmin+0.53*(xmax-xmin),ymin+0.78*(ymax-ymin),'Speed (km/yr)',fontsize=11,fontname='Arial')
cb.ax.tick_params(labelsize=11)
ax.plot([xmin+0.54*(xmax-xmin),xmin+0.54*(xmax-xmin)+5e3],[ymin+0.72*(ymax-ymin),ymin+0.72*(ymax-ymin)],'k',linewidth=1.5)
ax.plot([xmin+0.54*(xmax-xmin),xmin+0.54*(xmax-xmin)],[ymin+0.72*(ymax-ymin),ymin+0.70*(ymax-ymin)],'k',linewidth=1.5)
ax.plot([xmin+0.54*(xmax-xmin)+5e3,xmin+0.54*(xmax-xmin)+5e3],[ymin+0.72*(ymax-ymin),ymin+0.70*(ymax-ymin)],'k',linewidth=1.5)
ax.text(xmin+0.58*(xmax-xmin)+5e3,ymin+0.69*(ymax-ymin),'5 km',fontsize=11,fontname='Arial')

plt.subplots_adjust(wspace=0,hspace=0,top=0.98,right=0.98,left=0.02,bottom=0.02)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Defense_Kanger_map.pdf"),FORMAT='PDF',dpi=600,transparent=True)

for i in range(1,len(ind_eul_K)):
  ax.plot(x_K[ind_eul_K[i]],y_K[ind_eul_K[i]],'ko',markerfacecolor=coloptions[i],markersize=5)
  if i in [1,2,4]:
    ax.text(x_K[ind_eul_K[i]]-4e3,y_K[ind_eul_K[i]]-1.5e3,'K'+'{0:02d}'.format(int(abs(dists_eul[i]))),fontsize=9,fontname='Arial')
  else:
    ax.text(x_K[ind_eul_K[i]],y_K[ind_eul_K[i]]+0.75e3,'K'+'{0:02d}'.format(int(abs(dists_eul[i]))),fontsize=9,fontname='Arial')

plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Defense_Kanger_map_points.pdf"),FORMAT='PDF',dpi=600,transparent=True)
ax.add_patch(patch)
plt.close()

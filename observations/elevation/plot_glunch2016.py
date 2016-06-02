# Make some plots for my Glunch 2016 presentation

import datelib, geotifflib, glaclib, zslib, fluxlib, vellib, masklib, icefrontlib, climlib, floatlib, zslib
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

# Flux gates
xgate_H,ygate_H = fluxlib.fluxbox_geometry('Helheim',"fluxgate1")
xgate_K,ygate_K = fluxlib.fluxbox_geometry('Kanger',"fluxgate1")

# Load DEMs	
xdem_H,ydem_H,zdem_H,timedem_H,errordem_H = zslib.dem_grid('Helheim',285000.0,320000.0,-2588000.0,-2566000.0,years='all',verticaldatum='ellipsoid',method='nearest',return_error=True)
xdem_K,ydem_K,zdem_K,timedem_K,errordem_K = zslib.dem_grid('Kanger',449800.0,503000.0,-2302000.0,-2266000.0,years='all',verticaldatum='ellipsoid',method='nearest',return_error=True)


# Thinning rates
time_H,dH_H,Q_H,hbar_H,ubar_H,smb_H,width_H,area_H = fluxlib.fluxgate_thinning('Helheim','fluxgate1','cresis',10.,timing='velocity')
time_K,dH_K,Q_K,hbar_K,ubar_K,smb_K,width_K,area_K = fluxlib.fluxgate_thinning('Kanger','fluxgate1','cresis',10.,timing='velocity')
dem_time_H,dem_dH_H = fluxlib.dem_thinning('Helheim',xdem_H,ydem_H,zdem_H,timedem_H,errordem_H,"fluxgate1")
dem_time_K,dem_dH_K = fluxlib.dem_thinning('Kanger',xdem_K,ydem_K,zdem_K,timedem_K,errordem_K,"fluxgate1")

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

# Filter length for 
filt_len = 11.

# Load smb from RACMO
xrac_H,yrac_H,smbrac_H,timerac_H = climlib.racmo_at_pts(np.mean(xgate_H),np.mean(ygate_H),'smb',filt_len=filt_len)
xrac_K,yrac_K,smbrac_K,timerac_K = climlib.racmo_at_pts(np.mean(xgate_K),np.mean(ygate_K),'smb',filt_len=filt_len)

xrac_H,yrac_H,runrac_H,timerac_H = climlib.racmo_at_pts(np.mean(xgate_H),np.mean(ygate_H),'runoff',filt_len=filt_len)
xrac_K,yrac_K,runrac_K,timerac_K = climlib.racmo_at_pts(np.mean(xgate_K),np.mean(ygate_K),'runoff',filt_len=filt_len)

ind = np.argmin(abs(dists_H-np.nanmax(terminus_val_H)))
xsst_H,ysst_H,sst_H,timesst_H = climlib.SIF_at_pts(np.mean(x_H[ind]),np.mean(y_H[ind]),filt_len=filt_len,variable='sst')
ind = np.argmin(abs(dists_K-np.nanmax(terminus_val_K)))
xsst_K,ysst_K,sst_K,timesst_K = climlib.SIF_at_pts(np.mean(x_K[ind]),np.mean(y_K[ind]),filt_len=filt_len,variable='sst')

ind = np.argmin(abs(dists_H-np.nanmax(terminus_val_H)))
xsif_H,ysif_H,sif_H,timesif_H = climlib.SIF_at_pts(np.mean(x_H[ind]),np.mean(y_H[ind]),filt_len=filt_len)
ind = np.argmin(abs(dists_K-np.nanmax(terminus_val_K)))
xsif_K,ysif_K,sif_K,timesif_K = climlib.SIF_at_pts(np.mean(x_K[ind]),np.mean(y_K[ind]),filt_len=filt_len)

yH_runoff,day1H_runoff,day2H_runoff,meltlengthH_runoff,totalH_runoff = climlib.seasonlength(timerac_H,runrac_H,'runoff')
yH_smb,day1H_smb,day2_smb,meltlengthH_smb,totalH_smb = climlib.seasonlength(timerac_H,smbrac_H,'smb')

yK_runoff,day1K_runoff,day2K_runoff,meltlengthK_runoff,totalK_runoff = climlib.seasonlength(timerac_K,runrac_K,'runoff')
yK_smb,day1K_smb,day2_smb,meltlengthK_smb,totalK_smb = climlib.seasonlength(timerac_K,smbrac_K,'smb')

# Color options for plotting
coloptions=['r','b','g','limegreen','gold']

#####################
# Make overview map #
##################### 

coast=np.loadtxt(os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/Coastlines/greenland_coast_wgs84.txt"))
xvel,yvel,vel=geotifflib.read(os.path.join(os.getenv("DATA_HOME"),"Velocity/Random/Greenland/LandsatAll/landsatall_v_wgs84.tif"))
#xvel,yvel,vel=geotifflib.read(os.path.join(os.getenv("DATA_HOME"),"Velocity/Random/Greenland/track-07to10/vel_07to10_v_wgs84.tif"))

# Set up basemap 
fig=plt.figure(figsize=(3.3,6))
ax = plt.gca()
ax.tick_params('font',family='Arial',color='w',fontsize=12)
m = Basemap(width=1600000,height=2800000,
            resolution='l',projection='stere',\
            lat_ts=60,lat_0=72,lon_0=-41.)

nx = int((m.xmax-m.xmin)/1000.)+1; ny = int((m.ymax-m.ymin)/1000.)+1
velgrid = m.transform_scalar(vel,xvel,yvel,nx,ny)
xcoast,ycoast = m( coast[:,0],coast[:,1] )

clip = Path(np.column_stack([xcoast,ycoast]))
clip = PathPatch(clip,transform=ax.transData,facecolor='w')

cx = cubehelix.cmap(start=1.2,rot=-1.1,reverse=True,minLight=0.22,sat=3)
ax.add_patch(clip)
cs = m.imshow(velgrid,cmap=cx,vmin=0,vmax=200,zorder=2)
#cs.cmap.set_over('mediumvioletred')
#cbaxes = fig.add_axes([0.58, 0.2, 0.25, 0.02]) 
cb = m.colorbar(cs,location='bottom',ticks=[0,100,200]) 
plt.setp(plt.getp(cb.ax.axes, 'xticklabels'), color='w')
cb.ax.xaxis.set_tick_params(labelsize=12)
cb.set_label('Glacier speed (m/yr)',fontsize=16,fontname='Arial',color='w',fontweight='bold') 
cb.update_ticks()
 
cs.set_clip_path(clip)

m.plot(xcoast,ycoast,'k')

xg,yg = m (-33.0,68.63333)
m.plot(xg,yg,'ro',markersize=5,zorder=3)

xg,yg = m([-38.3],[66.1])
m.plot(xg,yg,'ro',markersize=5,zorder=3)

ax.axis('off')
plt.tight_layout()
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_overviewmap.pdf"),format='PDF',transparent=True,dpi=600)
plt.close()

del m,ax,xg,yg,cs,cb,clip,xcoast,ycoast,nx,ny,fig,xvel,yvel,vel,coast,velgrid

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
ax.imshow(vel_masked/1e3,extent=[np.min(xvel),np.max(xvel),np.min(yvel),np.max(yvel)],origin='lower',alpha=0.5,clim=[0,10],cmap=cx)
ax.axis('equal')

xmin = 283000.
xmax = 313000.
ymin = -2587000.
ymax = -2552000.

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
ax.plot(x_H[ind_eul_H[1]],y_H[ind_eul_H[1]],'ro',markersize=5)

plt.subplots_adjust(wspace=0,hspace=0,top=0.98,right=0.98,left=0.02,bottom=0.02)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_Helheim_map.pdf"),FORMAT='PDF',dpi=600,transparent=True)

ax.plot(x_H[dists_H>-21e3],y_H[dists_H>-21e3],'k',lw=1.5)
path = matplotlib.path.Path(np.column_stack([xgate_H,ygate_H]))
patch = matplotlib.patches.PathPatch(path,hatch='xxx',edgecolor='k',facecolor='none',lw=1)
ax.add_patch(patch)
for i in range(0,len(ind_eul_H)):
  ax.plot(x_H[ind_eul_H[i]],y_H[ind_eul_H[i]],'o',color=coloptions[i],markersize=5)
  if i in [1,2,4]:
    ax.text(x_H[ind_eul_H[i]]-4e3,y_H[ind_eul_H[i]]-1.5e3,'H'+'{0:02d}'.format(int(abs(dists_eul[i]))),fontsize=9,fontname='Arial')
  else:
    ax.text(x_H[ind_eul_H[i]],y_H[ind_eul_H[i]]+0.75e3,'H'+'{0:02d}'.format(int(abs(dists_eul[i]))),fontsize=9,fontname='Arial')


plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_Helheim_map1.pdf"),FORMAT='PDF',dpi=600,transparent=True)
plt.close()

del ax,patch,path,xmin,xmax,ymin,ymax,cbaxes,cb,vel_masked,cx,xmask,ymask,mask,xvel,yvel,vel,ximage,yimage,image,imagetime,ax2,fig

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
ax.imshow(vel_masked/1e3,extent=[np.min(xvel),np.max(xvel),np.min(yvel),np.max(yvel)],origin='lower',alpha=0.5,clim=[0,10],cmap=cx)
ax.axis('equal')

xmin = 468000.
xmax = 498000.
ymin = -2299000.
ymax = -2264000.


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
ax.plot(x_K[ind_eul_K[1]],y_K[ind_eul_K[1]],'bo',markersize=5)

plt.subplots_adjust(wspace=0,hspace=0,top=0.98,right=0.98,left=0.02,bottom=0.02)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_Kanger_map.pdf"),FORMAT='PDF',dpi=600,transparent=True)

ax.plot(x_K[dists_K>-21e3],y_K[dists_K>-21e3],'k',lw=1.5)
path = matplotlib.path.Path(np.column_stack([xgate_K,ygate_K]))
patch = matplotlib.patches.PathPatch(path,hatch='xxx',edgecolor='k',facecolor='none',lw=1)
ax.add_patch(patch)
for i in range(0,len(ind_eul_K)):
  ax.plot(x_K[ind_eul_K[i]],y_K[ind_eul_K[i]],'o',color=coloptions[i],markersize=5)
  if i in [1,2,4]:
    ax.text(x_K[ind_eul_K[i]]-4e3,y_K[ind_eul_K[i]]-1.5e3,'K'+'{0:02d}'.format(int(abs(dists_eul[i]))),fontsize=9,fontname='Arial')
  else:
    ax.text(x_K[ind_eul_K[i]],y_K[ind_eul_K[i]]+0.75e3,'K'+'{0:02d}'.format(int(abs(dists_eul[i]))),fontsize=9,fontname='Arial')

plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_Kanger_map1.pdf"),FORMAT='PDF',dpi=600,transparent=True)
plt.close()

del ax,patch,path,xmin,xmax,ymin,ymax,cbaxes,cb,vel_masked,cx,xmask,ymask,mask,xvel,yvel,vel,ximage,yimage,image,imagetime,ax2,fig

############################
# Plot behavior since 2000 #
############################

# Helheim ice-front position
plt.figure(figsize=(9,3))
ax=plt.gca()
for child in ax.get_children():
  if isinstance(child, matplotlib.spines.Spine):
    child.set_color('k')
ax.tick_params(axis='both',colors='k',which='both')
ax.arrow(2007,3.7,0,0.6,head_width=0.1,fc='k',ec='k',lw=1.5,head_length=0.1)
ax.text(2007.2,3.7,'Advance',fontsize=16)
ax.arrow(2007,-2.7,0,-0.6,head_width=0.1,fc='k',ec='k',lw=1.5,head_length=0.1)
ax.text(2007.2,-3.3,'Retreat',fontsize=16)
plt.yticks(np.arange(-4,6,2),fontsize=12)
plt.ylim([-4,5])
plt.xlim([2000,2016])
plt.plot(terminus_time_H,terminus_val_H/1e3,'ro',markersize=4)
plt.ylabel('Ice-front position (km)',color='k',fontsize=14)
x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
ax.xaxis.set_major_formatter(x_formatter)
ax.text(2000.3,-3.4,'Helheim',fontsize=20)
plt.tight_layout()
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_Helheim_retreat2000.pdf"),format="PDF")
plt.close()

# Kanger ice-front position
plt.figure(figsize=(9,3))
ax=plt.gca()
for child in ax.get_children():
  if isinstance(child, matplotlib.spines.Spine):
    child.set_color('k')
ax.tick_params(axis='both',colors='k',which='both')
ax.arrow(2007,3.7,0,0.6,head_width=0.1,fc='k',ec='k',lw=1.5,head_length=0.1)
ax.text(2007.2,3.7,'Advance',fontsize=16)
ax.arrow(2007,-2.7,0,-0.6,head_width=0.1,fc='k',ec='k',lw=1.5,head_length=0.1)
ax.text(2007.2,-3.3,'Retreat',fontsize=16)
plt.yticks(np.arange(-4,6,2),fontsize=12)
plt.yticks(np.arange(-4,6,2),fontsize=12)
plt.ylim([-4,5])
plt.xlim([2000,2016])
plt.plot(terminus_time_K,terminus_val_K/1e3,'bo',markersize=4)
plt.ylabel('Ice-front position (km)',color='k',fontsize=14)
x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
ax.xaxis.set_major_formatter(x_formatter)
ax.text(2000.3,-3.4,'Kangerdlugssuaq',fontsize=20)
plt.tight_layout()
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_Kanger_retreat2000.pdf"),format="PDF")
plt.close()

# Helheim speedup
plt.figure(figsize=(9,3))
ax=plt.gca()
for child in ax.get_children():
  if isinstance(child, matplotlib.spines.Spine):
    child.set_color('k')
ax.tick_params(axis='both',colors='k',which='both')
ax.arrow(2007,11.2,0,0.6,head_width=0.1,fc='k',ec='k',lw=1.5,head_length=0.1)
ax.text(2007.2,11.2,'Speedup',fontsize=16)
ax.arrow(2007,5.7,0,-0.6,head_width=0.1,fc='k',ec='k',lw=1.5,head_length=0.1)
ax.text(2007.2,5.1,'Slowdown',fontsize=16)
plt.xlim([2000,2016])
plt.plot(vel_time_H,vel_val_H[:,1]/1e3,'ro',markersize=4)
plt.ylabel('Glacier speed (km/yr)',color='k',fontsize=14)
plt.ylim([4.5,12.5])
x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
ax.xaxis.set_major_formatter(x_formatter)
ax.text(2000.3,11,'Helheim',fontsize=20)
plt.tight_layout()
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_Helheim_speedup2000.pdf"),format="PDF")
plt.close()

# Kanger speedup
plt.figure(figsize=(9,3))
ax=plt.gca()
for child in ax.get_children():
  if isinstance(child, matplotlib.spines.Spine):
    child.set_color('k')
ax.tick_params(axis='both',colors='k',which='both')
ax.arrow(2007,11.2,0,0.6,head_width=0.1,fc='k',ec='k',lw=1.5,head_length=0.1)
ax.text(2007.2,11.2,'Speedup',fontsize=16)
ax.arrow(2007,5.7,0,-0.6,head_width=0.1,fc='k',ec='k',lw=1.5,head_length=0.1)
ax.text(2007.2,5.1,'Slowdown',fontsize=16)
#plt.yticks(np.arange(-4,6,2),fontsize=12)
#plt.ylim([-4,5])
plt.xlim([2000,2016])
plt.plot(vel_time_K,vel_val_K[:,1]/1e3,'bo',markersize=4)
plt.ylabel('Glacier speed (km/yr)',color='k',fontsize=14)
plt.ylim([4.5,12.5])
x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
ax.xaxis.set_major_formatter(x_formatter)
ax.text(2000.3,11,'Kangerdlugssuaq',fontsize=20)
plt.tight_layout()
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_Kanger_speedup2000.pdf"),format="PDF")
plt.close()

# Helheim thinning
plt.figure(figsize=(9,3))
ax=plt.gca()
for child in ax.get_children():
  if isinstance(child, matplotlib.spines.Spine):
    child.set_color('k')
ax.tick_params(axis='both',colors='k',which='both')
ax.arrow(2007,295,0,15,head_width=0.1,fc='k',ec='k',lw=1.5,head_length=3)
ax.text(2007.2,295,'Thickening',fontsize=16)
ax.arrow(2007,125,0,-15,head_width=0.1,fc='k',ec='k',lw=1.5,head_length=3)
ax.text(2007.2,110,'Thinning',fontsize=16)
#plt.yticks(np.arange(-4,6,2),fontsize=12)
#plt.ylim([-4,5])
plt.xlim([2000,2016])
plt.plot(time_atm_H,zpt_atm_H[:,1],'ro',markersize=4)
plt.plot(time_dem_H,zpt_dem_H[:,1],'ro',markersize=4)
plt.ylabel('Surface elevation (m asl)',color='k',fontsize=14)
#plt.ylim([4.5,12.5])
plt.ylim([90,330])
x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
ax.xaxis.set_major_formatter(x_formatter)
ax.text(2000.3,107,'Helheim',fontsize=20)
plt.tight_layout()
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_Helheim_thinning2000.pdf"),format="PDF")
plt.close()

# Kanger thinning
plt.figure(figsize=(9,3))
ax=plt.gca()
for child in ax.get_children():
  if isinstance(child, matplotlib.spines.Spine):
    child.set_color('k')
ax.tick_params(axis='both',colors='k',which='both')
ax.arrow(2007,295,0,15,head_width=0.1,fc='k',ec='k',lw=1.5,head_length=3)
ax.text(2007.2,295,'Thickening',fontsize=16)
ax.arrow(2007,125,0,-15,head_width=0.1,fc='k',ec='k',lw=1.5,head_length=3)
ax.text(2007.2,110,'Thinning',fontsize=16)
#plt.yticks(np.arange(-4,6,2),fontsize=12)
#plt.ylim([-4,5])
plt.xlim([2000,2016])
plt.plot(time_atm_K,zpt_atm_K[:,1],'bo',markersize=4)
plt.plot(time_dem_K,zpt_dem_K[:,1],'bo',markersize=4)
plt.ylim([90,330])
plt.ylabel('Surface elevation (m asl)',color='k',fontsize=14)
#plt.ylim([4.5,12.5])
x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
ax.xaxis.set_major_formatter(x_formatter)
ax.text(2000.3,107,'Kangerdlugssuaq',fontsize=20)
plt.tight_layout()
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_Kanger_thinning2000.pdf"),format="PDF")
plt.close()

#######################
# Plot bed elevations #
#######################

for glacier in ['Helheim','Kanger']:
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
  elif glacier == 'Kanger':
    plt.plot(dists/1e3,atm_data['20010520'][:,2],'b',label='20 May 2001',lw=1.2)
    ind = np.where((terminus_time_K > 2008.) & (terminus_time_K < 2016.))
    plt.plot([np.min(terminus_val_K[ind])/1e3,np.min(terminus_val_K[ind])/1e3],[-1500,1500],'0.5',lw=1.5)
    plt.plot([np.max(terminus_val_K[ind])/1e3,np.max(terminus_val_K[ind])/1e3],[-1500,1500],'0.5')
    plt.plot(dists/1e3,zsdem[0,:].T,color='r',linewidth=1.2,label='14 July 2013')
    ind = np.argmin(abs(dists_K--5e3))
    plt.plot(dists_K[0:ind]/1e3,floatlib.height(zb_K[0:ind]),'k',linewidth=1.2,label='Flotation',dashes=[2,2,2,2])
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
  plt.savefig(os.path.join(os.getenv("HOME"),'Bigtmp/Glunch_'+glacier+'_zs_flowline.pdf'),FORMAT='PDF',dpi=600)
  plt.close()

del zsdem,atm_data,x,y,dists,ax,cresis,zcresis,dcresis,fig,gs,junk,ind,tcresis,cutoff


###############
# Main figure #
###############

for glacier in ['Helheim','Kanger']:

  vbars = 'runoff'
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
    
    demtimes = ['20130209','20130508','20130804','20131031','20140127','20140330']
  elif glacier == 'Kanger':
    x = x_K; y = y_K; zb = zb_K; dists = dists_K
    vel_val = vel_val_K; vel_time = vel_time_K
    terminus_val = terminus_val_K; terminus_time = terminus_time_K
    time_dem = time_dem_K; zpt_dem = zpt_dem_K; zpterror_dem = zpterror_dem_K
    time_atm = time_atm_K; zpt_atm = zpt_atm_K
    day1_runoff = day1H_runoff; day2_runoff = day2H_runoff; year_runoff = yH_runoff
    ind_eul = ind_eul_K
    smb = smb_K; dem_dH = dem_dH_K; dem_time = dem_time_K;
    time = time_K; dH = dH_K
    
    demtimes = ['20100807','20110307','20110501','20110826','20111106','20120213']
      
  calvingstyle = icefrontlib.calving(glacier)

  for k in range(0,4+len(demtimes)):
    
    if k > 3: 
      if glacier == 'Helheim':
        plt.figure(figsize=(3.5,6.5))  
        time1 = 2013.0; time2 = 2014.5
      elif glacier == 'Kanger':  
        plt.figure(figsize=(3.5,6.5))  
        time1 = 2010.0; time2 = 2012.5
    else:
      plt.figure(figsize=(7.45,6.5))
      time1 = 2008.0; time2 = 2016.0;  
 
    gs = matplotlib.gridspec.GridSpec(7,1)
    matplotlib.rc('font',family='Arial')

    # Plot terminus
    plt.subplot(gs[0, :])
    ax = plt.gca()
    ind = np.where(calvingstyle[:,1] == 'Tabular')[0]
    plt.ylim([-6,6])

    if vbars == 'runoff':
      for i in range(0,len(year_runoff)):
        path = matplotlib.path.Path([[day1_runoff[i],-5],[day1_runoff[i],5],[day2_runoff[i],5],[day2_runoff[i],-5]])
        patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
        ax.add_patch(patch)
    if plot_calving:
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
          plt.bar(float(calvingstyle[i,0])-0.025,10.0,width=0.025,color=[1,0.6,0.6],edgecolor='none',label='Non-tabular',bottom=-5)
        elif i != 0:
          plt.bar(float(calvingstyle[i,0])-0.025,10.0,width=0.025,color=[1,0.6,0.6],edgecolor='none',bottom=-5)
      plt.legend(loc=2,fontsize=11,numpoints=1,handlelength=0.3,labelspacing=0.05,ncol=3,columnspacing=0.7,handletextpad=0.2,borderpad=0.25)
    nonnan = np.where(~(np.isnan(terminus_val)))[0]
    plt.plot(terminus_time[nonnan],terminus_val[nonnan]/1e3,'ko',linewidth=1,markersize=2.5)
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
    ax.set_xticklabels(labels,fontsize=12,fontname='Arial')
    plt.xlim([time1,time2])
    plt.yticks(np.arange(-6,8,2),fontsize=12,fontname="Arial")
    plt.ylabel('Terminus (km)',fontsize=12,fontname="Arial")
    ax.tick_params('both', length=6, width=1.25, which='major')
    ax.tick_params('both', length=3, width=1, which='minor')
    if k > 3:
      demtime = datelib.date_to_fracyear(int(demtimes[k-4][0:4]),int(demtimes[k-4][4:6]),int(demtimes[k-4][6:8]))
      plt.plot([demtime,demtime],[-5,5],'k--',lw=1.5)
    if glacier == 'Helheim':
      plt.ylim([-3,3])
    elif glacier == 'Kanger':
      plt.ylim([-4.5,4.5])
    else:  
      plt.ylim(np.floor((np.min(terminus_val))/1e3),np.ceil((np.max(terminus_val))/1e3))


    # Plot velocities
    ax = plt.subplot(gs[1:-4, :]) 
    coloptions=['r','b','g','limegreen','gold']
    markoptions=['o','o','o','o','o','o']

    if k > 0:
      for i in range(0,len(dists_eul)):
        nonnan = np.where(~(np.isnan(vel_val[:,i])))[0]
        plt.plot(vel_time[nonnan],(vel_val[nonnan,i])/1e3,markoptions[i],color=coloptions[i],label=glacier[0]+'{0:02d}'.format(int(abs(dists_eul[i]))),markersize=3.5)
      if k < 4:
        plt.legend(loc=2,borderpad=0.3,fontsize=11,numpoints=1,handlelength=0.2,handletextpad=0.5,labelspacing=0.1,ncol=2,columnspacing=0.8)
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
    plt.xlim([time1,time2])
    if k > 3:
      plt.plot([demtime,demtime],[0,12],'k--',lw=1.5)
    if glacier == 'Helheim':
      plt.ylim([3.5,11])
    elif glacier == 'Kanger':
      plt.ylim([2.5,12])


    # Plot surface elevations
    plt.subplot(gs[-4:-2, :])
    ax = plt.gca()
    # Set up points for legend
    plt.errorbar(0,0,capsize=1,yerr=0.5,fmt='o',color='k',markersize=3.5,label='DEM')
    plt.plot(0,0,'+',color='k',markersize=5,label='ATM')
    if glacier == 'Helheim':
      plt.plot(0,0,'rs',label='H02',markersize=3.5)
      plt.plot(0,0,'bs',label='H05',markersize=3.5)
    elif glacier == 'Kanger':
      plt.plot(0,0,'bs',label='K05',markersize=3.5)
      plt.plot(0,0,'gs',label='K10',markersize=3.5)
    if vbars == 'runoff':
      for i in range(0,len(year_runoff)):
        path = matplotlib.path.Path([[day1_runoff[i],-50],[day1_runoff[i],240],[day2_runoff[i],240],[day2_runoff[i],-50]])
        patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
        ax.add_patch(patch)
    if glacier == 'Helheim' and (k > 1):
		  plt.fill_between([time1,time2],np.ones(2)*floatlib.height(zb[ind_eul[0]]-50)-floatlib.height(zb[ind_eul[0]]),np.ones(2)*floatlib.height(zb[ind_eul[0]]+50)-floatlib.height(zb[ind_eul[0]]),
  			alpha=0.1,facecolor='r',edgecolor='r',antialiased=True,zorder=2)
		  plt.plot([time1,time2],[0,0],'r:',linewidth=1.0)
		  plt.errorbar(time_dem,zpt_dem[:,0]-floatlib.height(zb[ind_eul[0]]),capsize=1,yerr=zpterror_dem,fmt='o',color=coloptions[0],markersize=3.5)
		  plt.plot(time_atm,zpt_atm[:,0]-floatlib.height(zb[ind_eul[0]]),'+',color=coloptions[0],markersize=5)
    elif glacier == 'Kanger' and (k > 1):
		  plt.fill_between([time1,time2],np.ones(2)*floatlib.height(zb[ind_eul[1]]-50)-floatlib.height(zb[ind_eul[1]]),np.ones(2)*floatlib.height(zb[ind_eul[1]]+50)-floatlib.height(zb[ind_eul[1]]),
  			alpha=0.1,facecolor='b',edgecolor='b',antialiased=True,zorder=2)
		  plt.plot([time1,time2],[0,0],'b:',linewidth=1.0)
		  plt.errorbar(time_dem,zpt_dem[:,1]-floatlib.height(zb[ind_eul[1]]),capsize=1,yerr=zpterror_dem,fmt='o',color=coloptions[1],markersize=3.5)
		  plt.plot(time_atm,zpt_atm[:,1]-floatlib.height(zb[ind_eul[1]]),'+',color=coloptions[1],markersize=5)
    x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax.xaxis.set_major_formatter(x_formatter)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.set_xticklabels([])
    plt.xticks(range(2000,2016),fontsize=12,fontname="Arial")
    plt.xlim([time1,time2])
    ax.tick_params('both', length=6, width=1.25, which='major')
    ax.tick_params('both', length=3, width=1, which='minor')
    if k > 3:
      plt.plot([demtime,demtime],[-30,100],'k--',lw=1.5)  
    if glacier == 'Helheim':
      plt.ylabel('  Height above flotation \n at H02 (m)',fontsize=12,fontname="Arial")
      plt.yticks(np.arange(-10,50,10),fontsize=12,fontname="Arial")
      plt.ylim([-5,45])
    elif glacier == 'Kanger':
      plt.ylabel('  Height above flotation \n at K05 (m)',fontsize=12,fontname="Arial")
      plt.yticks(np.arange(-20,80,20),fontsize=12,fontname="Arial")
      plt.ylim([-25,62])
    if (k > 1) and (k < 4):
      ind = [0,2,3,1]
      handles, labels = ax.get_legend_handles_labels()
      labels = [labels[i] for i in ind]
      handles = [handles[i] for i in ind]
      plt.legend(handles,labels,loc=3,borderpad=0.3,handleheight=0.1,fontsize=11,numpoints=1,handlelength=0.4,labelspacing=0.05,ncol=2,columnspacing=0.7,handletextpad=0.5)

    ax2 = ax.twinx()
    if glacier == 'Helheim' and (k>1):
      ax2.errorbar(time_dem,zpt_dem[:,1]-floatlib.height(zb[ind_eul[1]]),capsize=1,yerr=zpterror_dem,fmt='o',color=coloptions[1],markersize=3.5,label='WV')
      ax2.plot(time_atm,zpt_atm[:,1]-floatlib.height(zb[ind_eul[1]]),'+',color=coloptions[1],markersize=5,label='ATM')
    elif glacier == 'Kanger' and (k>1):
      ax2.errorbar(time_dem,zpt_dem[:,2]-floatlib.height(zb[ind_eul[2]]),capsize=1,yerr=zpterror_dem,fmt='o',color=coloptions[2],markersize=3.5,label='WV')
      ax2.plot(time_atm,zpt_atm[:,2]-floatlib.height(zb[ind_eul[2]]),'+',color=coloptions[2],markersize=5,label='ATM')
    ax2.set_xlim([time1,time2])
    ax2.tick_params('both', length=6, width=1.25, which='major')
    ax2.tick_params('both', length=3, width=1, which='minor')
    ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
    if glacier == 'Helheim':
      ax2.set_yticks(np.arange(60,120,10))
      ax2.set_ylim([55,105])
      ax2.set_yticklabels([60,70,80,90,100],fontsize=12,fontname='arial')
      ax2.set_ylabel('    Height above flotation at H05 (m)',fontsize=12,fontname='Arial')
    elif glacier == 'Kanger':
      ax2.set_yticks(np.arange(150,250,20))
      ax2.set_ylim([145,232])
      ax2.set_yticklabels(np.arange(150,250,10),fontsize=12,fontname='arial')
      ax2.set_ylabel('   Height above flotation at K10 (m)',fontsize=12,fontname='Arial')


    plt.subplot(gs[-2:,:])
    ax = plt.gca()
    plt.plot([time1,time2],[0,0],'k')
    if k > 2:
      plt.plot(timerac_H,smbrac_H*365.25/900.,'0.4',lw=1.2,label='SMB')
      plt.errorbar(dem_time[:,0],dem_dH[:,0],xerr=dem_time[:,1],yerr=dem_dH[:,1],fmt='bo',markersize=3,capsize=1,lw=0.5,label='DEM')
      plt.errorbar(time,dH[:,0]+smb,fmt='rs',markersize=3,capsize=1,lw=1.5,label='Flux-gate')
      if k < 4:
        plt.legend(loc=2,numpoints=1,ncol=3,handletextpad=0.2,fontsize=11,columnspacing=0.05,markerscale=1,handlelength=1.2,borderpad=0.2)
    plt.ylabel(r'dh/dt (m/yr)',fontname='Arial',fontsize=12)
    plt.yticks(np.arange(-160,160,40),fontsize=12,fontname='Arial')
    if k > 3:
      plt.plot([demtime,demtime],[-200,150],'k--',lw=1.5)  
    plt.ylim([-165,115])
    for i in range(0,len(yH_runoff)):
      path = matplotlib.path.Path([[day1_runoff[i],-165],[day1_runoff[i],160],[day2_runoff[i],160],[day2_runoff[i],-165]])
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

    plt.tight_layout()
    if k > 3:
      if glacier == 'Helheim':
        plt.subplots_adjust(hspace=0.03,wspace=0,top=0.93,right=0.83,left=0.2,bottom=0.06) 
      else:
        plt.subplots_adjust(hspace=0.03,wspace=0,top=0.93,right=0.83,left=0.24,bottom=0.06)
    else:
      plt.subplots_adjust(hspace=0.03,wspace=0,top=0.93,right=0.92,left=0.11,bottom=0.06) 
    plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_"+glacier+"_all"+str(k)+".pdf"),FORMAT='PDF',dpi=600)
    plt.close()

  del k,i, demtime, path, patch, ax, labels, handles, time1, time2

  ###################################
  # Plot individual DEMs and images #
  ###################################

  if glacier == 'Helheim':
    imagefiles = ['mosaicHelheim.2013-040.148.31377_1-20mgeo.tif',\
                  'mosaicHelheim.2013-128.148.32713_1-20mgeo.tif',\
                  'mosaicHelheim.2013-216.148.34049_1-20mgeo.tif',\
                  'mosaicHelheim.2013-304.148.35385_1-20mgeo.tif',\
                  'mosaicHelheim.2014-027.148.36721_1-20mgeo.tif',\
                  'mosaicHelheim.2014-093.148.37723_1-20mgeo.tif',]
    xmin = 304000.0
    xmax = 313600.0
    ymin = -2582000.0
    ymax = -2573000.0
    
    ind = np.where(np.isnan(zb))[0]
    ind2 = np.where(~(np.isnan(zb)))[0]
    zb[ind] = np.interp(dists[ind],dists[ind2],zb[ind2])
    
  elif glacier == 'Kanger':
    imagefiles = ['mosaicKang.2010-235.163.17698_1-20mgeo.tif',\
                  'mosaicKang.2011-068.163.3974_1-20mgeo.tif',\
                  'mosaicKang.2011-123.163.4809_1-20mgeo.tif',\
                  'mosaicKang.2011-233.163.23209_1-20mgeo.tif',\
                  'mosaicKang.2011-310.163.24378_1-20mgeo.tif',\
                  'mosaicKang.2012-044.163.25881_1-20mgeo.tif']
    
    xmin = 483700.0
    xmax = 497200.0
    ymin = -2297500.0
    ymax = -2283500.0
    
    ind = np.argmin(abs(dists--4.9e3))
    zb[ind:]=float('nan')

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


  for k in range(0,len(demtimes)):
    
    # Plot glacier cross section
    zs_demtime,demtime = zslib.dem_along_flowline(x,y,glacier,years=demtimes[k],verticaldatum='geoid',filt_len=100.)

    fig = plt.figure(figsize=(4,2.3))
    plt.plot(dcresis/1e3,zcresis,'.',c='0.7',markersize=4)
    plt.plot(dists/1e3,floatlib.height(zb),'k:',lw=2)
    plt.plot(dists/1e3,zs_demtime[0,:],'b',lw=1.5)
    plt.plot(dists/1e3,floatlib.icebottom(zb,zs_demtime[0,:]),'b',lw=1.5)
    ind = np.where(~(np.isnan(zs_demtime[0,:])))[0][-1]
    plt.plot([dists[ind]/1e3,dists[ind]/1e3],[floatlib.icebottom(zb,zs_demtime[0,:])[ind],zs_demtime[0,ind]],'b',lw=1.5)
    plt.plot(dists/1e3,zb,'k',lw=1.5)
    plt.xlabel('Distance along flowline (km)',fontsize=10,fontname='Arial')
    plt.ylabel('Elevation (m asl)',fontsize=10,fontname='Arial')
    plt.xticks(np.arange(-20,10,5),fontsize=10,fontname='Arial')
    plt.xlim([-20,3])
    plt.yticks(np.arange(-1000,1000,500),fontsize=10,fontname='Arial')
    if glacier == 'Helheim':
      plt.ylim([-1200,700])
    elif glacier == 'Kanger':
      plt.ylim([-1300,700])
      #plt.text(-2.5,280,demtimes[k][0:4]+' '+datelib.month_to_string(int(demtimes[k][4:6]))+' '+demtimes[k][6:8],fontsize=11,fontname='Arial')
    plt.text(-4,500,demtimes[k][0:4]+' '+datelib.month_to_string(int(demtimes[k][4:6]))+' '+demtimes[k][6:8],fontsize=11,fontname='Arial')
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.03,wspace=0,top=0.97,right=0.97,left=0.18,bottom=0.18) 
    plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_"+glacier+"_dem"+str(k)+".pdf"),FORMAT='PDF',dpi=600)
    plt.close()
    
    
    # Plot satellite image
    
    DIR = os.path.join(os.getenv("DATA_HOME"),"Mosaics/"+glacier+"/")
    file = geotifflib.read(DIR+imagefiles[k])
    
    fig = plt.figure(figsize=(3.4,3))
    plt.imshow(file[2],extent=[file[0][0],file[0][-1],file[1][0],file[1][-1]],origin='lower',cmap='Greys_r')
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    plt.xticks([])
    plt.yticks([])
    plt.plot(x,y,'k',lw=1.5)
    if glacier == 'Helheim':
      year,month,day = datelib.doy_to_date(int(imagefiles[k][14:18]),int(imagefiles[k][19:22]))
    elif glacier == 'Kanger':
      year,month,day = datelib.doy_to_date(int(imagefiles[k][11:15]),int(imagefiles[k][16:19]))      
    plt.text(xmin+0.04*(xmax-xmin),0.91*(ymax-ymin)+ymin,str(year)+' '+datelib.month_to_string(month)+' '+str(int(np.round(day))),fontsize=15,fontname='Arial')
    for i in range(0,3):
      plt.plot(x[ind_eul[i]],y[ind_eul[i]],'o',color=coloptions[i],markersize=10)
      plt.text(x[ind_eul[i]]+300,y[ind_eul[i]]+300,glacier[0]+'{0:02d}'.format(int(abs(dists_eul[i]))),fontsize=14,fontname='Arial')

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.03,wspace=0,top=0.99,right=0.99,left=0.01,bottom=0.01) 
    plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_"+glacier+"_image"+str(k)+".pdf"),FORMAT='PDF',dpi=600,transparent=True)
    plt.close()
  
  del cresis,dcresis,zcresis

######################
# Compare ice fronts #
######################

for glacier in ['Helheim','Kanger']:

  vbars = 'runoff'

  if glacier == 'Helheim':
    x = x_H; y = y_H; zb = zb_H; dists = dists_H
    terminus_val = terminus_val_H; terminus_time = terminus_time_H
    day1_runoff = day1H_runoff; day2_runoff = day2H_runoff; year_runoff = yH_runoff
    ind_eul = ind_eul_H
    time = time_H; dH = dH_H
    timesif = timesif_H; sif = sif_H
    timerac = timerac_H; runrac = runrac_H
    colorpt = 'r'
    
  elif glacier == 'Kanger':
    x = x_K; y = y_K; zb = zb_K; dists = dists_K
    terminus_val = terminus_val_K; terminus_time = terminus_time_K
    day1_runoff = day1H_runoff; day2_runoff = day2H_runoff; year_runoff = yH_runoff
    ind_eul = ind_eul_K
    timesif = timesif_K; sif = sif_K
    timerac = timerac_K; runrac = runrac_K
    colorpt = 'b'
      
  calvingstyle = icefrontlib.calving(glacier)

  for k in range(0,2):
    
    if k == 0:
      plot_calving=0
    elif k == 1:
      plot_calving=1

    plt.figure(figsize=(4.5,4.5))
    time1 = 2008.0; time2 = 2016.0;  
 
    gs = matplotlib.gridspec.GridSpec(4,1)
    matplotlib.rc('font',family='Arial')

    # Plot terminus
    plt.subplot(gs[0:2, :])
    ax = plt.gca()
    ind = np.where(calvingstyle[:,1] == 'Tabular')[0]
    plt.ylim([-6,6])

    if vbars == 'runoff':
      for i in range(0,len(year_runoff)):
        path = matplotlib.path.Path([[day1_runoff[i],-5],[day1_runoff[i],5],[day2_runoff[i],5],[day2_runoff[i],-5]])
        patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
        ax.add_patch(patch)
    if plot_calving:
      for i in ind:
        if i == ind[1]:
          plt.bar(float(calvingstyle[i,0])-0.05,10.0,width=0.05,color=[0.6,0.6,1],edgecolor='none',label='Tabular',bottom=-5)
        elif i !=0:
          plt.bar(float(calvingstyle[i,0])-0.05,10.0,width=0.05,color=[0.6,0.6,1],edgecolor='none',bottom=-5)
      ind = np.where((calvingstyle[:,1] == 'Mixed'))[0]
      for i in ind:
        if i == ind[1]:
          plt.bar(float(calvingstyle[i,0])-0.05,10.0,width=0.05,color=[0.4,0.8,0.6],edgecolor='none',label='Mixed',bottom=-5)
        elif i != 0:
          plt.bar(float(calvingstyle[i,0])-0.05,10.0,width=0.05,color=[0.4,0.8,0.6],edgecolor='none',bottom=-5)
      ind = np.where((calvingstyle[:,1] == 'Domino'))[0]
      for i in ind:
        if i == ind[1]:
          plt.bar(float(calvingstyle[i,0])-0.05,10.0,width=0.05,color=[1,0.6,0.6],edgecolor='none',label='Non-tabular',bottom=-5)
        elif i != 0:
          plt.bar(float(calvingstyle[i,0])-0.05,10.0,width=0.05,color=[1,0.6,0.6],edgecolor='none',bottom=-5)
      plt.legend(loc=2,fontsize=11,numpoints=1,handlelength=0.3,labelspacing=0.05,ncol=3,columnspacing=0.7,handletextpad=0.2,borderpad=0.25)
    nonnan = np.where(~(np.isnan(terminus_val)))[0]
    plt.plot(terminus_time[nonnan],terminus_val[nonnan]/1e3,'o',linewidth=1,markersize=3,color=colorpt)
    x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax.xaxis.set_major_formatter(x_formatter)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    plt.xticks(range(2000,2017))
    plt.xlim([time1,time2])
    ax.set_xticklabels([])
    plt.yticks(np.arange(-6,8,2),fontsize=11,fontname="Arial")
    plt.ylabel('Terminus (km)',fontsize=11,fontname="Arial")
    ax.tick_params('both', length=6, width=1.25, which='major')
    ax.tick_params('both', length=3, width=1, which='minor')
    plt.ylim([-4.5,4.5])


    # Plot SIF
    ax = plt.subplot(gs[2, :]) 
    plt.plot(timesif,sif,lw=1,color=colorpt)
    x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax.xaxis.set_major_formatter(x_formatter)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.set_xticklabels([])
    plt.xticks(range(2000,2016),fontsize=11,fontname="Arial")
    plt.xlim([time1,time2])
    ax.tick_params('both', length=6, width=1.25, which='major')
    ax.tick_params('both', length=3, width=1, which='minor')
    plt.ylim([0,1])
    plt.yticks([0,1.1,0.5])
    plt.ylabel('Sea ice fraction',fontsize=11,fontname='Arial')
    for i in range(0,len(yH_runoff)):
      path = matplotlib.path.Path([[day1_runoff[i],-165],[day1_runoff[i],160],[day2_runoff[i],160],[day2_runoff[i],-165]])
      patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
      ax.add_patch(patch) 

    # Plot surface runoff
    plt.subplot(gs[-1:,:])
    plt.plot(timerac,runrac,color=colorpt,lw=1)
    ax = plt.gca()
    plt.ylabel('Runoff \n'+r'(kg/m$^2$yr',fontname='Arial',fontsize=11)
    plt.plot
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
    ax.set_xticklabels(labels,fontsize=11,fontname='Arial')
    plt.xlim([time1,time2])
    plt.ylim([0,59])
    plt.yticks(np.arange(0,60,20))
    for i in range(0,len(yH_runoff)):
      path = matplotlib.path.Path([[day1_runoff[i],-165],[day1_runoff[i],160],[day2_runoff[i],160],[day2_runoff[i],-165]])
      patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
      ax.add_patch(patch) 

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.03,wspace=0,top=0.98,right=0.95,left=0.16,bottom=0.09) 
    plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_"+glacier+"_icefront"+str(k)+".pdf"),FORMAT='PDF',dpi=600)
    plt.close()

  del k,i, path, patch, ax, labels, time1, time2

#################
# Plot ice flux #
#################


for glacier in ['Helheim','Kanger']:
  if glacier == 'Helheim':
    balance = 31.9 # from Howat 2010
    time = time_H
    Q = Q_H
    colorpt = 'r'
  elif glacier == 'Kanger':
    balance = 22.5 # from Howat 2010
    time = time_K
    Q = Q_K
    colorpt = 'b'
  fig = plt.figure(figsize=(5,4))
  ax = plt.gca()
  plt.plot([2008,2016.5],[balance,balance],'k--',lw=2)
  inttime = np.arange(2008,2016,0.01)
  ind = np.where(~(np.isnan(Q[:,1])))[0]
  meanQ = np.mean(np.interp(inttime,time[ind],Q[ind,1]))
  #plt.plot([2008,2016],[meanQ/1e9,meanQ/1e9],'--',lw=2,color=colorpt)
  plt.plot(time,Q[:,1]/1e9,'o',markersize=4,color=colorpt)
  plt.ylabel('Flux (km$^3$/yr)',fontsize=13,fontname='Arial')
  plt.xlim([2008,2016])
  x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
  ax.xaxis.set_major_formatter(x_formatter)
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  ax.yaxis.set_minor_locator(AutoMinorLocator(2))
  ax.tick_params('both', length=6, width=1.25, which='major')
  ax.tick_params('both', length=3, width=1, which='minor')
  plt.xticks(np.arange(2008,2017),fontsize=11)
  plt.xlim([2008,2016.3])
  plt.ylim([21,41])
  labels=[]
  for i in range(2008,2017):
    labels.append('Jan \n'+str(i))
  ax.set_xticklabels(labels,fontsize=12,fontname='Arial')
  if glacier == 'Helheim':
    plt.text(2011,32.5,'Balance flux',fontsize=14)
  elif glacier == 'Kanger':
    plt.text(2008.2,23.2,'Balance flux',fontsize=14)

  plt.tight_layout()
  plt.subplots_adjust(hspace=0.03,wspace=0,top=0.97,right=0.95,left=0.12,bottom=0.11) 
  plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_"+glacier+"_flux.pdf"),FORMAT='PDF',dpi=600)
  plt.close()

del inttime, time, Q, balance fig, ax, labels, x_formatter, colorpt, meanQ

################################
# Ice-front position vs. speed #
################################

fig = plt.figure(figsize=(4,4))
ind = np.where((vel_time_H > 2008.) & (vel_time_H < 2016.) & ~(np.isnan(vel_val_H[:,1])))[0]
tint = np.interp(vel_time_H,terminus_time_H,terminus_val_H)
slope,intercept,r,p,std_err = stats.linregress(tint[ind],vel_val_H[ind,1])
plt.plot(tint[ind]/1e3,(tint[ind]*slope+intercept)/1e3,'k',lw=1.5)
plt.text(0.8,10,'slope = {0:.3f}'.format(slope),fontsize=14,fontname='Arial')
plt.plot(tint[ind]/1e3,vel_val_H[ind,1]/1e3,'ro')
plt.xticks(np.arange(-4,6,2),fontsize=14)
plt.yticks(np.arange(6,12,1),fontsize=14)
plt.xlim([-4,4])
plt.ylim([5.5,10.5])
plt.xlabel('Terminus (km)',fontsize=14,fontname='Arial')
plt.ylabel('Glacier speed at H05 (km/yr)',fontsize=14,fontname='Arial')
plt.tight_layout()
plt.subplots_adjust(hspace=0.03,wspace=0,top=0.97,right=0.95,left=0.14,bottom=0.125) 
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_Helheim_H05.pdf"),FORMAT='PDF',dpi=600)
plt.close()

fig = plt.figure(figsize=(4,4))
ind = np.where((vel_time_K > 2008.) & (vel_time_K < 2016.)  & ~(np.isnan(vel_val_K[:,1])))[0]
tint = np.interp(vel_time_K,terminus_time_K,terminus_val_K)
slope,intercept,r,p,std_err = stats.linregress(tint[ind],vel_val_K[ind,1])
plt.plot(tint[ind]/1e3,(tint[ind]*slope+intercept)/1e3,'k',lw=1.5)
plt.text(0.8,10,'slope = {0:.3f}'.format(slope),fontsize=14,fontname='Arial')
plt.plot(tint[ind]/1e3,vel_val_K[ind,1]/1e3,'bo')
plt.xticks(np.arange(-4,6,2),fontsize=14)
plt.yticks(np.arange(6,12,1),fontsize=14)
plt.xlim([-4,4])
plt.ylim([5.5,10.5])
plt.xlabel('Terminus (km)',fontsize=14,fontname='Arial')
plt.ylabel('Glacier speed at K05 (km/yr)',fontsize=14,fontname='Arial')
plt.tight_layout()
plt.subplots_adjust(hspace=0.03,wspace=0,top=0.97,right=0.95,left=0.14,bottom=0.125) 
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_Kanger_K05.pdf"),FORMAT='PDF',dpi=600)
plt.close()

fig = plt.figure(figsize=(4,4))
ind = np.where((vel_time_H > 2008.) & (vel_time_H < 2016.) & ~(np.isnan(vel_val_H[:,3])))[0]
tint = np.interp(vel_time_H,terminus_time_H,terminus_val_H)
slope,intercept,r,p,std_err = stats.linregress(tint[ind],vel_val_H[ind,3])
plt.plot(tint[ind]/1e3,(tint[ind]*slope+intercept)/1e3,'k',lw=1.5)
plt.text(0.8,7.5,'slope = {0:.3f}'.format(slope),fontsize=14,fontname='Arial')
plt.plot(tint[ind]/1e3,vel_val_H[ind,3]/1e3,'ro')
plt.xticks(np.arange(-4,6,2),fontsize=14)
plt.yticks(np.arange(3,12,1),fontsize=14)
plt.xlim([-4,4])
plt.ylim([3,8])
plt.xlabel('Terminus (km)',fontsize=14,fontname='Arial')
plt.ylabel('Glacier speed at H15 (km/yr)',fontsize=14,fontname='Arial')
plt.tight_layout()
plt.subplots_adjust(hspace=0.03,wspace=0,top=0.97,right=0.95,left=0.14,bottom=0.125) 
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_Helheim_H15.pdf"),FORMAT='PDF',dpi=600)
plt.close()

fig = plt.figure(figsize=(4,4))
ind = np.where((vel_time_K > 2008.) & (vel_time_K < 2016.) & ~(np.isnan(vel_val_K[:,3])))[0]
tint = np.interp(vel_time_K,terminus_time_K,terminus_val_K)
slope,intercept,r,p,std_err = stats.linregress(tint[ind],vel_val_K[ind,3])
plt.plot(tint[ind]/1e3,(tint[ind]*slope+intercept)/1e3,'k',lw=1.5)
plt.text(0.8,7.5,'slope = {0:.3f}'.format(slope),fontsize=14,fontname='Arial')
plt.plot(tint[ind]/1e3,vel_val_K[ind,3]/1e3,'bo')
plt.xticks(np.arange(-4,6,2),fontsize=14)
plt.yticks(np.arange(3,12,1),fontsize=14)
plt.xlim([-4,4])
plt.ylim([3,8])
plt.xlabel('Terminus (km)',fontsize=14,fontname='Arial')
plt.ylabel('Glacier speed at K15 (km/yr)',fontsize=14,fontname='Arial')
plt.tight_layout()
plt.subplots_adjust(hspace=0.03,wspace=0,top=0.97,right=0.95,left=0.14,bottom=0.125) 
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_Kanger_K15.pdf"),FORMAT='PDF',dpi=600)
plt.close()

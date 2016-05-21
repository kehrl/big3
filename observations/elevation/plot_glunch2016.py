# Make some plots for my Glunch 2016 presentation

import datelib, geotifflib, glaclib, zslib, fluxlib, vellib, masklib, icefrontlib, climlib
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
terminus_val_H, terminus_time_H = icefrontlib.distance_along_flowline(x_H,y_H,dists_H,'Helheim',type='icefront',time1=2000.,time2=2016.)
terminus_val_K, terminus_time_K = icefrontlib.distance_along_flowline(x_K,y_K,dists_K,'Kanger',type='icefront',time1=2000.,time2=2016.)

# Locations where we want velocities and surface elevations
dists_eul = -1*np.array([5.,10.,15.])
ind_eul_H=[]
for i in range(0,len(dists_eul)):
  ind_eul_H.append( (abs(dists_H - dists_eul[i]*1e3)).argmin() )
ind_eul_K=[]
for i in range(0,len(dists_eul)):
  ind_eul_K.append( (abs(dists_K - dists_eul[i]*1e3)).argmin() )

# Load velocities
vel_val_H,vel_time_H,vel_error_H = vellib.velocity_at_eulpoints(x_H[ind_eul_H],y_H[ind_eul_H],'Helheim',data='TSX')
vel_val_K,vel_time_K,vel_error_K = vellib.velocity_at_eulpoints(x_K[ind_eul_K],y_K[ind_eul_K],'Kanger',data='TSX')

# Load elevations
zpt_atm_H,zptstd_atm_H,time_atm_H = zslib.atm_at_pts(x_H[ind_eul_H],y_H[ind_eul_H],'Helheim',years='all',maxdist=500.,verticaldatum='geoid',method='average',cutoff='terminus')
zpt_dem_H,zpterror_dem_H,time_dem_H = zslib.dem_at_pts(x_H[ind_eul_H],y_H[ind_eul_H],'Helheim',years='all',verticaldatum='geoid',cutoff='terminus',method='average',radius=500.)
zpt_atm_K,zptstd_atm_K,time_atm_K = zslib.atm_at_pts(x_K[ind_eul_K],y_K[ind_eul_K],'Kanger',years='all',maxdist=500.,verticaldatum='geoid',method='average',cutoff='terminus')
zpt_dem_K,zpterror_dem_K,time_dem_K = zslib.dem_at_pts(x_K[ind_eul_K],y_K[ind_eul_K],'Kanger',years='all',verticaldatum='geoid',cutoff='terminus',method='average',radius=500.)

# Filter length for 
filt_len = 11.

# Load smb from RACMO
xrac_H,yrac_H,smbrac_H,timerac_H = climlib.racmo_at_pts(np.mean(xgate_H),np.mean(ygate_H),'smb',filt_len=filt_len)
xrac_K,yrac_K,smbrac_K,timerac_K = climlib.racmo_at_pts(np.mean(xgate_K),np.mean(ygate_K),'smb',filt_len=filt_len)

xrac_H,yrac_H,runrac_H,timerac_H = climlib.racmo_at_pts(np.mean(xgate_H),np.mean(ygate_H),'runoff',filt_len=filt_len)
xrac_K,yrac_K,runrac_K,timerac_K = climlib.racmo_at_pts(np.mean(xgate_K),np.mean(ygate_K),'runoff',filt_len=filt_len)

xrac_H,yrac_H,zsrac_H,timezsrac_H = climlib.racmo_at_pts(np.mean(xgate_H)-10e3,np.mean(ygate_H)+10e3,'zs',filt_len='none')
xrac_K,yrac_K,zsrac_K,timezsrac_K = climlib.racmo_at_pts(np.mean(xgate_K),np.mean(ygate_K),'zs',filt_len='none')

yH_runoff,day1H_runoff,day2H_runoff,meltlengthH_runoff,totalH_runoff = climlib.seasonlength(timerac_H,runrac_H,'runoff')
yH_smb,day1H_smb,day2_smb,meltlengthH_smb,totalH_smb = climlib.seasonlength(timerac_H,smbrac_H,'smb')

yK_runoff,day1K_runoff,day2K_runoff,meltlengthK_runoff,totalK_runoff = climlib.seasonlength(timerac_K,runrac_K,'runoff')
yK_smb,day1K_smb,day2_smb,meltlengthK_smb,totalK_smb = climlib.seasonlength(timerac_K,smbrac_K,'smb')

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

cx = cubehelix.cmap(start=1.2,rot=-1.1,reverse=True,minLight=0.1,sat=3)
ax.add_patch(clip)
cs = m.imshow(velgrid,cmap=cx,vmin=0,vmax=200,zorder=2)
#cbaxes = fig.add_axes([0.58, 0.2, 0.25, 0.02]) 
cb =m.colorbar(cs,location='bottom',ticks=[0,100,200]) 
plt.setp(plt.getp(cb.ax.axes, 'xticklabels'), color='w')
cb.ax.xaxis.set_tick_params(labelsize=12)
cb.set_label('Velocity (m/yr)',fontsize=16,fontname='Arial',color='w',fontweight='bold') 
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

plt.subplots_adjust(wspace=0,hspace=0,top=0.98,right=0.98,left=0.02,bottom=0.02)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_Helheim_map.pdf"),FORMAT='PDF',dpi=600,transparent=True)

path = matplotlib.path.Path(np.column_stack([xgate_H,ygate_H]))
patch = matplotlib.patches.PathPatch(path,hatch='xxx',edgecolor='r',facecolor='none',lw=1)
ax.add_patch(patch)
ax.plot(x_H,y_H,'r',lw=1)
zeroind = np.argmin(abs(dists_H--5e3))
ax.plot(x_H[zeroind],y_H[zeroind],'ro',markersize=6)
#ax.text(x_H[zeroind]+0.2e3,y_H[zeroind]-1.7e3,'H05',fontsize=8,fontname='Arial')

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

plt.subplots_adjust(wspace=0,hspace=0,top=0.98,right=0.98,left=0.02,bottom=0.02)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_Kanger_map.pdf"),FORMAT='PDF',dpi=600,transparent=True)

path = matplotlib.path.Path(np.column_stack([xgate_K,ygate_K]))
patch = matplotlib.patches.PathPatch(path,hatch='xxx',edgecolor='r',facecolor='none',lw=1)
ax.add_patch(patch)

ax.plot(x_K,y_K,'r',lw=1)
zeroind = np.argmin(abs(dists_K--5e3))
ax.plot(x_K[zeroind],y_K[zeroind],'ro',markersize=6)
#ax.text(x_K[zeroind]+0.2e3,y_K[zeroind]-1.7e3,'K05',fontsize=8,fontname='Arial')

plt.close()

del ax,patch,path,xmin,xmax,ymin,ymax,cbaxes,cb,vel_masked,cx,xmask,ymask,mask,xvel,yvel,vel,ximage,yimage,image,imagetime,ax2,fig

###########################
# Plot terminus positions #
###########################

plt.figure(figsize=(8,3))
ax=plt.gca()
for child in ax.get_children():
  if isinstance(child, matplotlib.spines.Spine):
    child.set_color('w')
ax.tick_params(axis='both',colors='w',which='both')
plt.plot(terminus_time_H,terminus_val_H,'wo',markersize=5)
x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
ax.xaxis.set_major_formatter(x_formatter)

plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_Helheim_retreat2000.pdf"),format="PDF",transparent=True)









##############
# Big figure #
##############

plt.figure(figsize=(10,10))
gs = matplotlib.gridspec.GridSpec(5,1)
matplotlib.rc('font',family='Arial')
time1=2008.
time2=2016.

plt.subplot(gs[0, :])
ax = plt.gca()
plt.plot(0,0,'k^',markersize=7,label='Helheim')
plt.plot(0,0,'ro',markersize=7,label='Kanger')
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1H_runoff[i],-120],[day1H_runoff[i],100],[day2H_runoff[i],100],[day2H_runoff[i],-120]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
  ax.add_patch(patch)
plt.xticks(range(2000,2017),fontsize=20,fontname="Arial")
ax.set_xticklabels([])
plt.xlim([time1,time2])
plt.yticks(np.arange(-4,6,2),fontsize=18)
plt.ylim([-4,4])
plt.ylabel('Terminus\n(km)',fontsize=20)
ax.tick_params('both', length=8, width=1.25, which='major')
plt.legend(loc=4,fontsize=18,numpoints=1,ncol=2,borderpad=0.2,handletextpad=0.2,columnspacing=0.5,handlelength=0.5)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.tick_params('x', length=4, width=1.25, which='minor')

plt.subplot(gs[1, :])
ax = plt.gca()
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1H_runoff[i],-120],[day1H_runoff[i],100],[day2H_runoff[i],100],[day2H_runoff[i],-120]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
  ax.add_patch(patch)
plt.xticks(range(2000,2017),fontsize=20,fontname="Arial")
ax.set_xticklabels([])
plt.xlim([time1,time2])
plt.ylabel('Glacier speed\n'+r'(km yr$^{-1}$)',fontsize=20)
plt.yticks(np.arange(6,12,2),fontsize=18)
plt.ylim([5.5,10.7])
ax.tick_params('both', length=8, width=1.25, which='major')
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.tick_params('x', length=4, width=1.25, which='minor')

#plt.subplot(gs[2, :])
#ax = plt.gca()
#for i in range(0,len(yH_runoff)):
#  path = matplotlib.path.Path([[day1H_runoff[i],-120],[day1H_runoff[i],100],[day2H_runoff[i],100],[day2H_runoff[i],-120]])
#  patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
#  ax.add_patch(patch)
#plt.xticks(range(2000,2017),fontsize=20,fontname="Arial")
#ax.set_xticklabels([])
#plt.xlim([time1,time2])
#plt.ylabel('Strain rate\n'+r'(yr$^{-1}$)',fontsize=20)
#plt.yticks(np.arange(0.1,0.5,0.1),fontsize=18)
#plt.ylim([0.1,0.4])
#ax.tick_params('both', length=8, width=1.25, which='major')
#ax.xaxis.set_minor_locator(AutoMinorLocator(2))
#ax.tick_params('x', length=4, width=1.25, which='minor')

plt.subplot(gs[2, :])
ax = plt.gca()
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1H_runoff[i],-120],[day1H_runoff[i],100],[day2H_runoff[i],100],[day2H_runoff[i],-120]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
  ax.add_patch(patch)
plt.xticks(range(2000,2017),fontsize=20,fontname="Arial")
ax.set_xticklabels([])
plt.xlim([time1,time2])
plt.ylabel('dH/dt\n'+r'(m yr$^{-1}$)',fontsize=20)
plt.yticks(np.arange(-80,120,40),fontsize=18)
plt.ylim([-90,70])
plt.plot([time1,time2],[0,0],'k',lw=1.25)
ax.tick_params('both', length=8, width=1.25, which='major')
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.tick_params('x', length=4, width=1.25, which='minor')

plt.subplot(gs[3, :])
ax = plt.gca()
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1H_runoff[i],-120],[day1H_runoff[i],100],[day2H_runoff[i],100],[day2H_runoff[i],-120]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
  ax.add_patch(patch)
plt.xticks(range(2000,2017),fontsize=20,fontname="Arial")
ax.set_xticklabels([])
plt.xlim([time1,time2])
plt.ylabel('Elevation\n(m)',fontsize=20)
plt.yticks(np.arange(-60,40,20),fontsize=18)
plt.ylim([-50,15])
ax.tick_params('both', length=8, width=1.25, which='major')
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.tick_params('x', length=4, width=1.25, which='minor')

plt.subplot(gs[4, :])
ax = plt.gca()
ax.plot(timerac_H,runrac_H,color='k',lw=2,label='Runoff')
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1H_runoff[i],-120],[day1H_runoff[i],100],[day2H_runoff[i],100],[day2H_runoff[i],-120]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
  ax.add_patch(patch)
x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
ax.xaxis.set_major_formatter(x_formatter)
ax.tick_params('both', length=8, width=1.25, which='major')
labels=[]
for i in range(2000,2017):
  labels.append('Jan \n'+str(i))
plt.xticks(range(2000,2017),fontsize=20,fontname='Arial')
ax.set_xticklabels(labels,fontsize=20,fontname='Arial')
plt.xlim([time1,time2])
plt.ylabel('Runoff\n'+r'(kg m$^{-2}$ d$^{-1}$)',fontsize=20)
plt.yticks([0,20,40],fontsize=18)
plt.ylim([0,50])
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.tick_params('x', length=4, width=1.25, which='minor')

plt.tight_layout()
plt.subplots_adjust(hspace=0.05,wspace=0,top=0.98,right=0.94,left=0.12,bottom=0.06) 

plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_main_1.pdf"),DPI=600,format='PDF')

plt.subplot(gs[0,:])
plt.plot(terminus_time_H,terminus_val_H/1e3,'^',color='k',label='Helheim',markersize=6)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_main_2.pdf"),DPI=600,format='PDF')

plt.subplot(gs[1,:])
plt.plot(vel_time_H,vel_val_H[:,0]/1e3,'^',color='k',markersize=6)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_main_3.pdf"),DPI=600,format='PDF')

#plt.subplot(gs[2,:])
#plt.plot(vel_time_H,(vel_val_H[:,1]-vel_val_H[:,2])/5e3,'^',color='k',markersize=6)
#plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_main_4.pdf"),DPI=600,format='PDF')

plt.subplot(gs[2,:])
plt.plot(timezsrac_H[1:],np.convolve(np.diff(zsrac_H)/np.diff(timezsrac_H),np.ones((11,))/11.,mode='same'),c='0.4',lw=2,label='SMB')
plt.plot(time_H,dH_H[:,0]+smb_H,'^',color='k',markersize=6)
plt.plot(dem_time_H,dem_dH_H[:,0],'^',color='k',markersize=6)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_main_4.pdf"),DPI=600,format='PDF')

plt.subplot(gs[3,:])
ind = np.argmin(abs(time_atm_H-2008.5))
plt.plot(time_atm_H,zpt_atm_H[:,0]-zpt_atm_H[ind,0],'^',color='k',markersize=6)
plt.plot(time_dem_H,zpt_dem_H[:,0]-zpt_atm_H[ind,0],'^',color='k',markersize=6)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_main_5.pdf"),DPI=600,format='PDF')

plt.subplot(gs[-1,:])
ax.plot(timerac_K,runrac_K,color='r',lw=2,label='Runoff')
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_main_6.pdf"),DPI=600,format='PDF')

plt.subplot(gs[0,:])
plt.plot(terminus_time_K,terminus_val_K/1e3,'o',color='r',label='Kanger',markersize=6)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_main_7.pdf"),DPI=600,format='PDF')

plt.subplot(gs[1,:])
plt.plot(vel_time_K,vel_val_K[:,0]/1e3,'o',color='r',markersize=6)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_main_8.pdf"),DPI=600,format='PDF')

#plt.subplot(gs[2,:])
#plt.plot(vel_time_K,(vel_val_K[:,1]-vel_val_K[:,2])/5e3,'o',color='r',markersize=6)
#plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_main_10.pdf"),DPI=600,format='PDF')

plt.subplot(gs[2,:])
#plt.plot(timezsrac_H[1:],np.convolve(np.diff(zsrac_H)/np.diff(timezsrac_H),np.ones((11,))/11.,mode='same'),c='0.4',lw=2,label='SMB')
plt.plot(time_K,dH_K[:,0]+smb_K,'o',color='r',markersize=6)
plt.plot(dem_time_K,dem_dH_K[:,0],'o',color='r',markersize=6)

plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_main_9.pdf"),DPI=600,format='PDF')

plt.subplot(gs[3,:])
ind = np.argmin(abs(time_atm_K-2008.5))
#plt.plot(time_atm_K,zpt_atm_K[:,0]-zpt_atm_K[ind,0],'o',color='r',markersize=6)
plt.plot(time_dem_K,zpt_dem_K[:,0]-zpt_atm_K[ind,0],'o',color='r',markersize=6)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Glunch_main_10.pdf"),DPI=600,format='PDF')

plt.close()
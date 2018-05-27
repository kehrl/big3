# Make overview figure for Helheim/Kanger inversions.
#
# Laura Kehrl, 10 April 2018

import os
import numpy as np
import cubehelix, matplotlib
import masklib, inverselib, geotifflib, datelib, glaclib, vellib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# Cutoff for velocities
velocity_cutoff = 1000 #m/yr

# Image for plotting
imagetime_HG = datelib.date_to_fracyear(2014,7,4)
ximage_HG,yimage_HG,image_HG = geotifflib.readrgb(os.path.join(os.getenv("DATA_HOME"),\
        "Imagery/Landsat/Helheim/TIF/20140704140535_LC82330132014185LGN00.tif"))

imagetime_KG = datelib.date_to_fracyear(2014,7,6)
ximage_KG,yimage_KG,image_KG = geotifflib.readrgb(os.path.join(os.getenv("DATA_HOME"),\
        "Imagery/Landsat/Kanger/TIF/20140706135251_LC82310122014187LGN00.tif"))

# Glacier extents for inversions
extent_HG = np.loadtxt('/Users/kehrl/Models/Helheim/3D/INV_SSA_ModelT/'+
        'DEM20120316_modelT_Lcurve/inputs/mesh_extent.dat')
hole1_HG = np.loadtxt('/Users/kehrl/Models/Helheim/3D/INV_SSA_ModelT/'+
        'DEM20120316_modelT_Lcurve/inputs/mesh_hole1.dat')
hole2_HG = np.loadtxt('/Users/kehrl/Models/Helheim/3D/INV_SSA_ModelT/'+
        'DEM20120316_modelT_Lcurve/inputs/mesh_hole2.dat')
extent_KG = np.loadtxt('/Users/kehrl/Models/Kanger/3D/INV_SSA_ModelT/'+
        'DEM20120522_modelT_Lcurve/inputs/mesh_extent.dat')

# Load velocity records
xvel_HG = np.arange(np.min(ximage_HG),np.max(ximage_HG),100)
yvel_HG = np.arange(np.min(yimage_HG),np.max(yimage_HG),100)
vx,vy = vellib.inversion_3D('Helheim',xvel_HG,yvel_HG,imagetime_HG,\
        dir_velocity_out='none',blur=False)
vel_HG = np.sqrt(vx**2+vy**2)

xvel_KG = np.arange(np.min(ximage_KG),np.max(ximage_KG),100)
yvel_KG = np.arange(np.min(yimage_KG),np.max(yimage_KG),100)
vx,vy = vellib.inversion_3D('Kanger',xvel_KG,yvel_KG,imagetime_KG,\
        dir_velocity_out='none',blur=False)
vel_KG = np.sqrt(vx**2+vy**2)
del vx,vy

# Load masks
xmask,ymask,mask = masklib.load_grid('Helheim',np.min(xvel_HG),np.max(xvel_HG),\
        np.min(yvel_HG),np.max(yvel_HG),100,icefront_time=imagetime_HG)
mask[:,xmask > 277000+38e3] = 1
vel_masked_HG = np.ma.masked_array(vel_HG,mask)

xmask,ymask,mask = masklib.load_grid('Kanger',np.min(xvel_KG),np.max(xvel_KG),\
        np.min(yvel_KG),np.max(yvel_KG),100,icefront_time=imagetime_KG)
mask[:,xmask > 462000+38e3] = 1
vel_masked_KG = np.ma.masked_array(vel_KG,mask)
del vel_KG, vel_HG, xmask, ymask, mask

# Load min velocities for showing region where velocity is above cutoff
x_cutoff_KG, y_cutoff_KG, vsurfini_cutoff_KG, ind_cutoff_grid, ind_cutoff  = inverselib.get_velocity_cutoff('Kanger',\
        velocity_cutoff=velocity_cutoff,model_dir='INV_SSA_ModelT')
x_cutoff_HG, y_cutoff_HG, vsurfini_cutoff_HG, ind_cutoff_grid, ind_cutoff  = inverselib.get_velocity_cutoff('Helheim',\
        velocity_cutoff=velocity_cutoff,model_dir='INV_SSA_ModelT')

# Figure
fig = plt.figure(figsize=(3.2,1.8))
gs = matplotlib.gridspec.GridSpec(1,2)

ax = plt.subplot(gs[0])
cx = cubehelix.cmap(start=1.2,rot=-1.1,reverse=True,minLight=0.1,sat=2)
p=plt.imshow(vel_masked_KG/1e3,extent=[np.min(xvel_KG),np.max(xvel_KG),np.min(yvel_KG),\
        np.max(yvel_KG)],origin='lower',clim=[0,8],cmap=cx)
ax = plt.gca()
ax.imshow(image_KG[:,:,0],extent=[np.min(ximage_KG),np.max(ximage_KG),np.min(yimage_KG),\
        np.max(yimage_KG)],cmap='Greys_r',origin='lower',clim=[0,0.6])
ax.imshow(vel_masked_KG/1e3,extent=[np.min(xvel_KG),np.max(xvel_KG),np.min(yvel_KG),\
        np.max(yvel_KG)],origin='lower',alpha=0.7,clim=[0,8],cmap=cx)
plt.contour(x_cutoff_KG,y_cutoff_KG,vsurfini_cutoff_KG,[velocity_cutoff],colors='k',linestyles='dashed',linewidths=1)
plt.contour(x_cutoff_KG,y_cutoff_KG,vsurfini_cutoff_KG,[4000],colors='w',linestyles='dashed',linewidths=1)
ax.axis('equal')
ax.plot(np.r_[extent_KG[:,0],extent_KG[0,0]],np.r_[extent_KG[:,1],extent_KG[0,1]],'k',lw=1)
plt.xlim([462000,462000+41e3])
plt.ylim([-2306000,-2306000+44e3])#[-2299000,-2264000])
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_yticks([])
ax.set_xticks([])
xmin,xmax = plt.xlim()
ymin,ymax = plt.ylim()
ax.text(xmin+0.04*(xmax-xmin),ymin+0.93*(ymax-ymin),'(a)',fontsize=8,fontweight='bold')
#ax.text(xmin+0.08*(xmax-xmin),ymin+0.95*(ymax-ymin),'Kangerlussuaq Glacier',fontsize=8)
ax.text(xmin+0.15*(xmax-xmin),ymin+0.93*(ymax-ymin),'KG',fontsize=8)

path = matplotlib.path.Path([[0.02*(xmax-xmin)+xmin,0.28*(ymax-ymin)+ymin],
                        [0.45*(xmax-xmin)+xmin,0.28*(ymax-ymin)+ymin],
                        [0.45*(xmax-xmin)+xmin,-0.01*(ymax-ymin)+ymin],
                        [0.02*(xmax-xmin)+xmin,-0.01*(ymax-ymin)+ymin],
                        [0.02*(xmax-xmin)+xmin,0.28*(ymax-ymin)+ymin]])
patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1)
ax.add_patch(patch)
cbaxes = fig.add_axes([0.045, 0.225, 0.17, 0.03])
cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=[0,4,8],extend='max')
ax.text(xmin+0.11*(xmax-xmin),ymin+0.02*(ymax-ymin),'km yr$^{-1}$',fontsize=8)
cb.ax.tick_params(labelsize=8)
ax.plot([xmin+0.5*(xmax-xmin),xmin+0.5*(xmax-xmin)+5e3],\
        [ymin+0.045*(ymax-ymin),ymin+0.045*(ymax-ymin)],'k',lw=1)
ax.plot([xmin+0.5*(xmax-xmin),xmin+0.5*(xmax-xmin)],\
        [ymin+0.045*(ymax-ymin),ymin+0.025*(ymax-ymin)],'k',lw=1)
ax.plot([xmin+0.5*(xmax-xmin)+5e3,xmin+0.5*(xmax-xmin)+5e3],\
        [ymin+0.045*(ymax-ymin),ymin+0.025*(ymax-ymin)],'k',lw=1)
ax.text(xmin+0.52*(xmax-xmin)+5e3,ymin+0.015*(ymax-ymin),'5 km',fontsize=8)

ax2 = fig.add_axes([0.295, 0.515, 0.25, 0.45])
m = Basemap(width=1600000,height=3000000,
            resolution='l',projection='stere',\
            lat_ts=60,lat_0=72,lon_0=-41.)
m.drawcoastlines(linewidth=0.75)
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='w',lake_color='w')
xg,yg = m ([-33.0],[68.63333])
m.plot(xg,yg,'ro',mec='k',mew=0.5,markersize=5)

ax = plt.subplot(gs[1])
p=plt.imshow(vel_masked_HG/1e3,extent=[np.min(xvel_HG),np.max(xvel_HG),np.min(yvel_HG),\
        np.max(yvel_HG)],origin='lower',clim=[0,10],cmap=cx)
ax = plt.gca()
ax.imshow(image_HG[:,:,0],extent=[np.min(ximage_HG),np.max(ximage_HG),np.min(yimage_HG),\
        np.max(yimage_HG)],cmap='Greys_r',origin='lower',clim=[0,0.6])
ax.imshow(vel_masked_HG/1e3,extent=[np.min(xvel_HG),np.max(xvel_HG),np.min(yvel_HG),\
        np.max(yvel_HG)],origin='lower',alpha=0.7,clim=[0,8],cmap=cx)
plt.contour(x_cutoff_HG,y_cutoff_HG,vsurfini_cutoff_HG,[velocity_cutoff],colors='k',linestyles='dashed',linewidths=1)
plt.contour(x_cutoff_HG,y_cutoff_HG,vsurfini_cutoff_HG,[4000],colors='w',linestyles='dashed',linewidths=1)
ax.plot(np.r_[extent_HG[:,0],extent_HG[0,0]],np.r_[extent_HG[:,1],extent_HG[0,1]],'k',lw=1)
ax.plot(np.r_[hole1_HG[:,0],hole1_HG[0,0]],np.r_[hole1_HG[:,1],hole1_HG[0,1]],'k',lw=1)
ax.plot(np.r_[hole2_HG[:,0],hole2_HG[0,0]],np.r_[hole2_HG[:,1],hole2_HG[0,1]],'k',lw=1)
ax.axis('equal')
plt.xlim([277000,277000+41e3])
plt.ylim([-2585500,-2585500+44e3])
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_yticks([])
ax.set_xticks([])
xmin,xmax = plt.xlim()
ymin,ymax = plt.ylim()
ax.text(xmin+0.04*(xmax-xmin),ymin+0.93*(ymax-ymin),'(b)',fontsize=8,fontweight='bold')
ax.text(xmin+0.15*(xmax-xmin),ymin+0.93*(ymax-ymin),'HG',fontsize=8)

ax2 = fig.add_axes([0.78, 0.515, 0.25, 0.45])
m = Basemap(width=1600000,height=3000000,
            resolution='l',projection='stere',\
            lat_ts=60,lat_0=72,lon_0=-41.)
m.drawcoastlines(linewidth=0.75)
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents(color='w',lake_color='w')
xg,yg = m([-38.3],[66.1])
m.plot(xg,yg,'ro',mec='k',mew=0.5,markersize=5)

plt.tight_layout
plt.subplots_adjust(wspace=0.0,hspace=0.0,top=0.99,right=0.99,left=0.01,bottom=0.01)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Overview_maps.pdf"),\
          FORMAT='PDF',DPI=300)
plt.close()


import sys
import os
import matplotlib.pyplot as plt
from pylab import *
import argparse
import scipy
import elmerreadlib, glaclib, colorlib, geotifflib, meshlib
import cubehelix
from matplotlib.path import Path
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

##########
# Inputs #
##########

args = sys.argv

# Get inputs to file
parser = argparse.ArgumentParser()
parser.add_argument("-glacier", dest="glacier",required = True,
                help = "Glacier name.")
parser.add_argument("-date",dest="date",required = False, default = 0,
                help = "Date for inversion.")

args, _ = parser.parse_known_args(sys.argv)
glacier = args.glacier
date = args.date

# Set default date if none is given
if date == 0:
 if glacier == 'Helheim':
   date = '20120316'
 elif glacier == 'Kanger':
   date = '20120213'

# Directory
DIR_FS_ModelT = os.path.join(os.getenv("MODEL_HOME"),glacier+"/Results/INV_FS_ModelT/")
DIR_FS_ConstantT = os.path.join(os.getenv("MODEL_HOME"),glacier+"/Results/INV_FS_ConstantT/")
DIR_SSA_ModelT = os.path.join(os.getenv("MODEL_HOME"),glacier+"/Results/INV_SSA_ModelT/")
DIR_SSA_ConstantT = os.path.join(os.getenv("MODEL_HOME"),glacier+"/Results/INV_SSA_ConstantT/")

###############################
# Load shapefiles for regions #
###############################

regions = True 
if regions == True:
  extent_region_U = meshlib.shp_to_xy(os.path.join(os.getenv("DATA_HOME"),\
          "ShapeFiles/Glaciers/3D/"+glacier+"/glacier_region1.shp"))
  extent_region_U = np.row_stack([extent_region_U,extent_region_U[0,:]])

  extent_region_M = meshlib.shp_to_xy(os.path.join(os.getenv("DATA_HOME"),\
                    "ShapeFiles/Glaciers/3D/"+glacier+"/glacier_region5.shp"))
  extent_region_M = np.row_stack([extent_region_M,extent_region_M[0,:]])

  extent_region_L = meshlib.shp_to_xy(os.path.join(os.getenv("DATA_HOME"),\
          "ShapeFiles/Glaciers/3D/"+glacier+"/glacier_region2.shp"))
  extent_region_L = np.row_stack([extent_region_L,extent_region_L[0,:]])
  
  extent_region_S1 = meshlib.shp_to_xy(os.path.join(os.getenv("DATA_HOME"),\
          "ShapeFiles/Glaciers/3D/"+glacier+"/glacier_region3.shp"))
  extent_region_S1 = np.row_stack([extent_region_S1,extent_region_S1[0,:]])
  
  extent_region_S2 = meshlib.shp_to_xy(os.path.join(os.getenv("DATA_HOME"),\
          "ShapeFiles/Glaciers/3D/"+glacier+"/glacier_region4.shp"))
  extent_region_S2 = np.row_stack([extent_region_S2,extent_region_S2[0,:]])
  
  region_U = Path(np.column_stack((extent_region_U[:,0],extent_region_U[:,1])))
  region_M = Path(np.column_stack((extent_region_M[:,0],extent_region_M[:,1])))  
  region_L = Path(np.column_stack((extent_region_L[:,0],extent_region_L[:,1])))
  region_S1 = Path(np.column_stack((extent_region_S1[:,0],extent_region_S1[:,1])))
  region_S2 = Path(np.column_stack((extent_region_S2[:,0],extent_region_S2[:,1])))

################
# Load Results #
################

# Load TIF files
x,y,taub_FS_ModelT = geotifflib.read(DIR_FS_ModelT+'DEM'+date+'_bed_mod_taub.tif') 
x,y,taub_FS_ConstantT = geotifflib.read(DIR_FS_ConstantT+'DEM'+date+'_bed_mod_taub.tif')
x,y,taub_SSA_ModelT = geotifflib.read(DIR_SSA_ModelT+'DEM'+date+'_bed_mod_taub.tif')
x,y,taub_SSA_ConstantT = geotifflib.read(DIR_SSA_ConstantT+'DEM'+date+'_bed_mod_taub.tif')

x,y,ub_FS_ModelT = geotifflib.read(DIR_FS_ModelT+'DEM'+date+'_bed_mod_ub.tif')
x,y,vb_FS_ModelT = geotifflib.read(DIR_FS_ModelT+'DEM'+date+'_bed_mod_vb.tif')
x,y,us_FS_ModelT = geotifflib.read(DIR_FS_ModelT+'DEM'+date+'_surf_mod_us.tif')
x,y,vs_FS_ModelT = geotifflib.read(DIR_FS_ModelT+'DEM'+date+'_surf_mod_vs.tif')

x,y,ub_FS_ConstantT = geotifflib.read(DIR_FS_ConstantT+'DEM'+date+'_bed_mod_ub.tif')
x,y,vb_FS_ConstantT = geotifflib.read(DIR_FS_ConstantT+'DEM'+date+'_bed_mod_vb.tif')
x,y,us_FS_ConstantT = geotifflib.read(DIR_FS_ConstantT+'DEM'+date+'_surf_mod_us.tif')
x,y,vs_FS_ConstantT = geotifflib.read(DIR_FS_ConstantT+'DEM'+date+'_surf_mod_vs.tif')

x,y,zs = geotifflib.read(DIR_FS_ConstantT+'DEM'+date+'_surf_mea_zs.tif')
x,y,zb = geotifflib.read(DIR_FS_ConstantT+'DEM'+date+'_bed_mod_zb.tif')
x,y,umea = geotifflib.read(DIR_FS_ConstantT+'DEM'+date+'_surf_mea_us.tif')
x,y,vmea = geotifflib.read(DIR_FS_ConstantT+'DEM'+date+'_surf_mea_vs.tif')
x,y,depthT = geotifflib.read(DIR_FS_ModelT+'DEM'+date+'_depthT.tif')

################
# Plot options #
################

if glacier == 'Kanger':
  xmin = 468000.
  xmax = 500800.
  ymin = -2300000.
  ymax = -2260200.
elif glacier == 'Helheim':
  xmin = 280200.
  xmax = 313000.
  ymin = -2585100.
  ymax = -2545300.

# Image for plotting
try:
  if glacier == "Helheim":
    ximage,yimage,image = geotifflib.readrgb(os.path.join(os.getenv("DATA_HOME"),"Imagery/Landsat/Helheim/TIF/20140704140535_LC82330132014185LGN00.tif"))
  elif glacier == "Kanger":
    ximage,yimage,image = geotifflib.readrgb(os.path.join(os.getenv("DATA_HOME"),"Imagery/Landsat/Kanger/TIF/20140706135251_LC82310122014187LGN00.tif"))
  sat_image = bool(True)
except:
  print "No satellite image for background"
  sat_image = bool(False)

############################
# Calculate driving stress #
############################

# Parameters
rho_i = 917.
g = 9.8

# Calculate grad(h) & H for driving stress
dhdx = np.zeros_like(zs)
dhdy = np.zeros_like(zs)
dhdx[:,:] = float('nan')
dhdy[:,:] = float('nan')

H = zs-zb

dhdx[1:-1,1:-1] = (zs[1:-1,2:]-zs[1:-1,0:-2])/(x[2:]-x[0:-2])
dhdy[1:-1,1:-1] = ((zs[2:,1:-1]-zs[0:-2,1:-1]).T/(y[2:]-y[0:-2])).T
#tauds = rho_i*g*H*np.sqrt((dhdx)**2+(dhdy)**2)
tauds_flow = -rho_i*g*H*(dhdx*umea+dhdy*vmea)/np.sqrt(umea**2+vmea**2)

del dhdx,dhdy,rho_i,g,H

#######################
# Plot driving stress #
#######################

fig = plt.figure(figsize=(2.4,3))
ax = plt.gca()
plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
p = plt.imshow(tauds_flow/1e3,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',vmin=-500,vmax=500,cmap='RdBu_r')
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
plt.xticks([])
plt.yticks([])
xmin,xmax = plt.xlim()
path = matplotlib.path.Path([[0.47*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.98*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.98*(xmax-xmin)+xmin,0.74*(ymax-ymin)+ymin],
                        [0.47*(xmax-xmin)+xmin,0.74*(ymax-ymin)+ymin],
                        [0.47*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin]])
patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1)
ax.add_patch(patch)
cbaxes = fig.add_axes([0.60, 0.905, 0.25, 0.035])
cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=np.arange(-500,600,500))
cb.ax.tick_params(labelsize=9)
ax.text(xmin+0.51*(xmax-xmin),ymin+0.805*(ymax-ymin),r'$\tau_d \cdot u / ||u||$ (kPa)',fontsize=9,fontname='Arial')
ax.plot([xmin+0.63*(xmax-xmin),xmin+0.63*(xmax-xmin)+5e3],[ymin+0.78*(ymax-ymin),ymin+0.78*(ymax-ymin)],'k',linewidth=1.5)
ax.plot([xmin+0.63*(xmax-xmin),xmin+0.63*(xmax-xmin)],[ymin+0.78*(ymax-ymin),ymin+0.76*(ymax-ymin)],'k',linewidth=1.5)
ax.plot([xmin+0.63*(xmax-xmin)+5e3,xmin+0.63*(xmax-xmin)+5e3],[ymin+0.78*(ymax-ymin),ymin+0.76*(ymax-ymin)],'k',linewidth=1.5)
ax.text(xmin+0.67*(xmax-xmin)+5e3,ymin+0.76*(ymax-ymin),'5 km',fontsize=9,fontname='Arial')
plt.tight_layout()
plt.subplots_adjust(hspace=0.03,wspace=0.03,top=0.99,right=0.99,left=0.01,bottom=0.01)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+'_inversion_'+date+'_taud.pdf'),FORMAT='PDF',dpi=600)
plt.close()

#################################################################################
# Plot basal shear stress for different forward models and temperature profiles #
#################################################################################

fig = plt.figure(figsize=(4.8,5.8))
matplotlib.rc('font',family='Arial')
gs = matplotlib.gridspec.GridSpec(2,2)

plt.subplot(gs[0,0])
ax = plt.gca()
plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
p=plt.imshow(taub_FS_ConstantT*1e3,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',vmin=0,vmax=500)
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
plt.xticks([])
plt.yticks([])
path = matplotlib.path.Path([[0.57*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.98*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.98*(xmax-xmin)+xmin,0.71*(ymax-ymin)+ymin],
                        [0.57*(xmax-xmin)+xmin,0.71*(ymax-ymin)+ymin],
                        [0.57*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin]])
patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1)
ax.add_patch(patch)
cbaxes = fig.add_axes([0.3, 0.945, 0.17, 0.025])
cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=np.arange(0,600,200))
ax.text(xmin+0.68*(xmax-xmin),ymin+0.78*(ymax-ymin),r'$\tau_b$ (kPa)',fontsize=10)
cb.ax.tick_params(labelsize=10)
ax.plot([xmin+0.63*(xmax-xmin),xmin+0.63*(xmax-xmin)+5e3],[ymin+0.75*(ymax-ymin),ymin+0.75*(ymax-ymin)],'k',linewidth=1.5)
ax.plot([xmin+0.63*(xmax-xmin),xmin+0.63*(xmax-xmin)],[ymin+0.75*(ymax-ymin),ymin+0.73*(ymax-ymin)],'k',linewidth=1.5)
ax.plot([xmin+0.63*(xmax-xmin)+5e3,xmin+0.63*(xmax-xmin)+5e3],[ymin+0.75*(ymax-ymin),ymin+0.73*(ymax-ymin)],'k',linewidth=1.5)
ax.text(xmin+0.66*(xmax-xmin)+5e3,ymin+0.725*(ymax-ymin),'5 km',fontsize=10)
ax.text(xmin+0.03*(xmax-xmin),ymin+0.92*(ymax-ymin),'a',fontweight='bold',fontsize=10)

plt.subplot(gs[0,1])
ax = plt.gca()
plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
p=plt.imshow(taub_FS_ModelT*1e3,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',vmin=0,vmax=500)
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
plt.xticks([])
plt.yticks([])
path = matplotlib.path.Path([[0.57*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.98*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.98*(xmax-xmin)+xmin,0.77*(ymax-ymin)+ymin],
                        [0.57*(xmax-xmin)+xmin,0.77*(ymax-ymin)+ymin],
                        [0.57*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin]])
patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1)
ax.add_patch(patch)
cbaxes = fig.add_axes([0.795, 0.945, 0.17, 0.025])
cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=np.arange(0,600,200))
ax.text(xmin+0.68*(xmax-xmin),ymin+0.79*(ymax-ymin),r'$\tau_b$ (kPa)',fontsize=10)
cb.ax.tick_params(labelsize=10)
ax.text(xmin+0.03*(xmax-xmin),ymin+0.92*(ymax-ymin),'b',fontweight='bold',fontsize=10)

plt.subplot(gs[1,0])
ax = plt.gca()
plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
p=plt.imshow(taub_SSA_ConstantT*1e3,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',vmin=0,vmax=500)
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
plt.xticks([])
plt.yticks([])
path = matplotlib.path.Path([[0.57*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.98*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.98*(xmax-xmin)+xmin,0.77*(ymax-ymin)+ymin],
                        [0.57*(xmax-xmin)+xmin,0.77*(ymax-ymin)+ymin],
                        [0.57*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin]])
patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1)
ax.add_patch(patch)
cbaxes = fig.add_axes([0.3, 0.45, 0.17, 0.025])
cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=np.arange(0,600,200))
ax.text(xmin+0.68*(xmax-xmin),ymin+0.79*(ymax-ymin),r'$\tau_b$ (kPa)',fontsize=10)
cb.ax.tick_params(labelsize=10)
ax.text(xmin+0.03*(xmax-xmin),ymin+0.92*(ymax-ymin),'c',fontweight='bold',fontsize=10)

plt.subplot(gs[1,1])
ax = plt.gca()
plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
p=plt.imshow(taub_SSA_ModelT*1e3,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',vmin=0,vmax=500)
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
plt.xticks([])
plt.yticks([])
path = matplotlib.path.Path([[0.57*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.98*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.98*(xmax-xmin)+xmin,0.77*(ymax-ymin)+ymin],
                        [0.57*(xmax-xmin)+xmin,0.77*(ymax-ymin)+ymin],
                        [0.57*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin]])
patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1)
ax.add_patch(patch)
cbaxes = fig.add_axes([0.795, 0.45, 0.17, 0.025])
cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=np.arange(0,600,200))
ax.text(xmin+0.68*(xmax-xmin),ymin+0.79*(ymax-ymin),r'$\tau_b$ (kPa)',fontsize=10)
cb.ax.tick_params(labelsize=10)
ax.text(xmin+0.03*(xmax-xmin),ymin+0.92*(ymax-ymin),'d',fontweight='bold',fontsize=10)

plt.tight_layout()
plt.subplots_adjust(hspace=0.03,wspace=0.03,top=0.99,right=0.99,left=0.01,bottom=0.01)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_inversion_"+date"_temperature_taub.pdf"),FORMAT='PDF',dpi=600)
plt.close()

############################################################
# Plot depth-averaged temperatures and temperature profile #
############################################################

if glacier == 'Helheim':
  pts = np.array([[301875,-2576310],[290102,-2562290]])
  labels=['A','B']
  colorlabels = ['b','r','g']
  styles = ['-','--',':']
elif glacier == 'Kanger':
  pts = np.array([[486354,-2289805],[478226,-2276526]])
  labels=['C','D']
  colorlabels = ['b','r','g']
  styles = ['-','--',':']

fig = plt.figure(figsize=(3.75,3))
matplotlib.rc('font',family='Arial')
gs = matplotlib.gridspec.GridSpec(1,3)

plt.subplot(gs[0:2])
ax = plt.gca()
plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
p=plt.imshow(depthT,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',vmin=-15,vmax=0,cmap='jet')
plt.xticks([])
plt.yticks([])
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
for i in range(0,len(pts[:,0])):
  plt.plot(pts[i,0],pts[i,1],'ko',markerfacecolor=colorlabels[i],lw=0.5)
  plt.text(pts[i,0]+1000,pts[i,1]-800,labels[i],fontsize=9,fontname='Arial')

xmin,xmax = plt.xlim()
ymin,ymax = plt.ylim()
path = matplotlib.path.Path([[0.42*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
  			[0.99*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
  			[0.99*(xmax-xmin)+xmin,0.73*(ymax-ymin)+ymin],
  			[0.42*(xmax-xmin)+xmin,0.73*(ymax-ymin)+ymin],
  			[0.42*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin]])
patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1)
ax.add_patch(patch)
cbaxes = fig.add_axes([0.28, 0.91, 0.26, 0.025]) 
cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=np.arange(-15,5,5)) 
ax.text(xmin+0.48*(xmax-xmin),ymin+0.81*(ymax-ymin),'Temperature ($^o$C)',fontsize=8)
cb.ax.tick_params(labelsize=8)
ax.plot([xmin+0.6*(xmax-xmin),xmin+0.58*(xmax-xmin)+5e3],[ymin+0.78*(ymax-ymin),ymin+0.78*(ymax-ymin)],'k',linewidth=1.5)
ax.plot([xmin+0.58*(xmax-xmin),xmin+0.58*(xmax-xmin)],[ymin+0.78*(ymax-ymin),ymin+0.76*(ymax-ymin)],'k',linewidth=1.5)
ax.plot([xmin+0.58*(xmax-xmin)+5e3,xmin+0.58*(xmax-xmin)+5e3],[ymin+0.78*(ymax-ymin),ymin+0.76*(ymax-ymin)],'k',linewidth=1.5)
ax.text(xmin+0.62*(xmax-xmin)+5e3,ymin+0.76*(ymax-ymin),'5 km',fontsize=8)
ax.text(xmin+0.03*(xmax-xmin),ymin+0.93*(ymax-ymin),'a',fontweight='bold',fontsize=10)

plt.subplot(gs[2])
ax = plt.gca()
for i in range(0,len(pts[:,0])):
  tempvtu = elmerreadlib.pvtu_file(DIR_FS_ModelT+'adjoint_beta0001.pvtu',['constant temperature'])
  column = elmerreadlib.values_in_column(tempvtu,pts[i,0],pts[i,1])
  plt.plot(column['constant temperature'],np.max(column['z'])-column['z'],styles[i],lw=1.5,label=labels[i],color=colorlabels[i])
ax.yaxis.tick_right()
ax.yaxis.set_label_position('right')
plt.legend(loc=3,borderpad=0.3,fontsize=8,numpoints=1,handlelength=2.0,handletextpad=0.5,labelspacing=0.1,ncol=1,columnspacing=0.8)
ax.invert_yaxis()
plt.ylabel('Depth (m)',fontsize=9,fontname='Arial')
plt.xlabel(r'Temperature ($^o$C)',fontsize=9,fontname='Arial')
plt.xticks(np.arange(-15,5,5),fontname='Arial',fontsize=8)
plt.yticks(np.arange(0,1800,200),fontname='Arial',fontsize=8)
plt.text(-15,120,'b',fontsize=10,fontweight='bold')

plt.xlim([-16,0])
plt.ylim([1600,0])

plt.tight_layout()
plt.subplots_adjust(hspace=0.02,wspace=0.02,top=0.97,right=0.86,left=0.01,bottom=0.13)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_temperature.pdf"),FORMAT='PDF',dpi=600)
plt.close()

######################
# Plot sliding ratio #
######################

fig = plt.figure(figsize=(4.9,3))
matplotlib.rc('font',family='Arial')
gs = matplotlib.gridspec.GridSpec(1,2)

cx = cubehelix.cmap(start=1.2,rot=-1.1,reverse=True,minLight=0.1,sat=2)

plt.subplot(gs[0])
ax = plt.gca()
plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
p=plt.imshow(np.sqrt(ub_FS_ConstantT**2+vb_FS_ConstantT**2)/np.sqrt(us_FS_ConstantT**2+vs_FS_ConstantT**2),extent=[x[0],x[-1],y[0],y[-1]],origin='lower',vmin=0,vmax=1)
plt.xticks([])
plt.yticks([])
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
#for i in range(0,len(pts[:,0])):
#  plt.plot(pts[i,0],pts[i,1],'o',color=colorlabels[i])
#  plt.text(pts[i,0]+1000,pts[i,1]-800,labels[i],fontsize=9,fontname='Arial')

xmin,xmax = plt.xlim()
ymin,ymax = plt.ylim()
path = matplotlib.path.Path([[0.42*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.99*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.99*(xmax-xmin)+xmin,0.73*(ymax-ymin)+ymin],
                        [0.42*(xmax-xmin)+xmin,0.73*(ymax-ymin)+ymin],
                        [0.42*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin]])
patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1)
ax.add_patch(patch)
cbaxes = fig.add_axes([0.24, 0.91, 0.22, 0.025])
cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=np.arange(0,1.5,0.5))
ax.text(xmin+0.48*(xmax-xmin),ymin+0.81*(ymax-ymin),'Sliding ratio',fontsize=8)
cb.ax.tick_params(labelsize=8)
ax.plot([xmin+0.63*(xmax-xmin),xmin+0.63*(xmax-xmin)+5e3],[ymin+0.78*(ymax-ymin),ymin+0.78*(ymax-ymin)],'k',linewidth=1.5)
ax.plot([xmin+0.63*(xmax-xmin),xmin+0.63*(xmax-xmin)],[ymin+0.78*(ymax-ymin),ymin+0.76*(ymax-ymin)],'k',linewidth=1.5)
ax.plot([xmin+0.63*(xmax-xmin)+5e3,xmin+0.63*(xmax-xmin)+5e3],[ymin+0.78*(ymax-ymin),ymin+0.76*(ymax-ymin)],'k',linewidth=1.5)
ax.text(xmin+0.67*(xmax-xmin)+5e3,ymin+0.76*(ymax-ymin),'5 km',fontsize=8)
ax.text(xmin+0.03*(xmax-xmin),ymin+0.93*(ymax-ymin),'a',fontweight='bold',fontsize=10)

plt.subplot(gs[1])
ax = plt.gca()
plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
p=plt.imshow(np.sqrt(ub_FS_ConstantT**2+vb_FS_ConstantT**2)/np.sqrt(us_FS_ConstantT**2+vs_FS_ConstantT**2),extent=[x[0],x[-1],y[0],y[-1]],origin='lower',vmin=0,vmax=1)
plt.xticks([])
plt.yticks([])
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
#for i in range(0,len(pts[:,0])):
#  plt.plot(pts[i,0],pts[i,1],'o',color=colorlabels[i])
#  plt.text(pts[i,0]+1000,pts[i,1]-800,labels[i],fontsize=9,fontname='Arial')

xmin,xmax = plt.xlim()
ymin,ymax = plt.ylim()
path = matplotlib.path.Path([[0.42*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.99*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.99*(xmax-xmin)+xmin,0.78*(ymax-ymin)+ymin],
                        [0.42*(xmax-xmin)+xmin,0.78*(ymax-ymin)+ymin],
                        [0.42*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin]])
patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1)
ax.add_patch(patch)
cbaxes = fig.add_axes([0.74, 0.91, 0.22, 0.025])
cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=np.arange(0,1.5,0.5))
ax.text(xmin+0.48*(xmax-xmin),ymin+0.81*(ymax-ymin),'Sliding ratio',fontsize=8)
cb.ax.tick_params(labelsize=8)
ax.text(xmin+0.03*(xmax-xmin),ymin+0.93*(ymax-ymin),'b',fontweight='bold',fontsize=10)

plt.tight_layout()
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_inversion_slidingratio.pdf"),FORMAT='PDF',dpi=600)
plt.close()

###############################################################
# Compare results from different models for different regions #
###############################################################

if regions == True:
  
  # Percentile range for calculating statistics
  q=50.
  
  # Get averages for Full Stokes
  ind_U_x = []
  ind_U_y = []
  ind_M_x = []
  ind_M_y = []
  ind_L_x = []
  ind_L_y = []
  ind_S_x = []
  ind_S_y = []
  for i in range(0,len(x)):
    for j in range(0,len(y)):
      if region_U.contains_point([x[i],y[j]]):
        ind_U_x.append(i)
        ind_U_y.append(j)
      if region_M.contains_point([x[i],y[j]]):
        ind_M_x.append(i)
        ind_M_y.append(j)
      if region_L.contains_point([x[i],y[j]]):
        ind_L_x.append(i)
        ind_L_y.append(j) 
      if (region_S1.contains_point([x[i],y[j]])) or (region_S2.contains_point([x[i],y[j]])):
        ind_S_x.append(i)
        ind_S_y.append(j)
  
  # Get values for FS-ConstantT
  region_U_taub_FS_ConstantT = np.array([np.mean(taub_FS_ConstantT[ind_U_y,ind_U_x]),\
          np.percentile(taub_FS_ConstantT[ind_U_y,ind_U_x],q/2),\
          np.percentile(taub_FS_ConstantT[ind_U_y,ind_U_x],100-q/2)])*1e3
  region_M_taub_FS_ConstantT = np.array([np.mean(taub_FS_ConstantT[ind_M_y,ind_M_x]),\
          np.percentile(taub_FS_ConstantT[ind_M_y,ind_M_x],q/2),\
          np.percentile(taub_FS_ConstantT[ind_M_y,ind_M_x],100-q/2)])*1e3
  region_L_taub_FS_ConstantT = np.array([np.mean(taub_FS_ConstantT[ind_L_y,ind_L_x]),\
          np.percentile(taub_FS_ConstantT[ind_L_y,ind_L_x],q/2),\
          np.percentile(taub_FS_ConstantT[ind_L_y,ind_L_x],100-q/2)])*1e3
  region_S_taub_FS_ConstantT = np.array([np.mean(taub_FS_ConstantT[ind_S_y,ind_S_x]),\
          np.percentile(taub_FS_ConstantT[ind_S_y,ind_S_x],q/2),\
          np.percentile(taub_FS_ConstantT[ind_S_y,ind_S_x],100-q/2)])*1e3

  # Get values for FS-ModelT
  region_U_taub_FS_ModelT = np.array([np.mean(taub_FS_ModelT[ind_U_y,ind_U_x]),\
          np.percentile(taub_FS_ModelT[ind_U_y,ind_U_x],q/2),\
          np.percentile(taub_FS_ModelT[ind_U_y,ind_U_x],100-q/2)])*1e3
  region_M_taub_FS_ModelT = np.array([np.mean(taub_FS_ModelT[ind_M_y,ind_M_x]),\
          np.percentile(taub_FS_ModelT[ind_M_y,ind_M_x],q/2),\
          np.percentile(taub_FS_ModelT[ind_M_y,ind_M_x],100-q/2)])*1e3
  region_L_taub_FS_ModelT = np.array([np.mean(taub_FS_ModelT[ind_L_y,ind_L_x]),\
          np.percentile(taub_FS_ModelT[ind_L_y,ind_L_x],q/2),\
          np.percentile(taub_FS_ModelT[ind_L_y,ind_L_x],100-q/2)])*1e3
  region_S_taub_FS_ModelT = np.array([np.mean(taub_FS_ModelT[ind_S_y,ind_S_x]),\
          np.percentile(taub_FS_ModelT[ind_S_y,ind_S_x],q/2),\
          np.percentile(taub_FS_ModelT[ind_S_y,ind_S_x],100-q/2)])*1e3
  
  # Get values for SSA-ConstantT
  region_U_taub_SSA_ConstantT = np.array([np.mean(taub_SSA_ConstantT[ind_U_y,ind_U_x]),\
          np.percentile(taub_SSA_ConstantT[ind_U_y,ind_U_x],q/2),\
          np.percentile(taub_SSA_ConstantT[ind_U_y,ind_U_x],100-q/2)])*1e3
  region_M_taub_SSA_ConstantT = np.array([np.mean(taub_SSA_ConstantT[ind_M_y,ind_M_x]),\
          np.percentile(taub_SSA_ConstantT[ind_M_y,ind_M_x],q/2),\
          np.percentile(taub_SSA_ConstantT[ind_M_y,ind_M_x],100-q/2)])*1e3
  region_L_taub_SSA_ConstantT = np.array([np.mean(taub_SSA_ConstantT[ind_L_y,ind_L_x]),\
          np.percentile(taub_SSA_ConstantT[ind_L_y,ind_L_x],q/2),\
          np.percentile(taub_SSA_ConstantT[ind_L_y,ind_L_x],100-q/2)])*1e3
  region_S_taub_SSA_ConstantT = np.array([np.mean(taub_SSA_ConstantT[ind_S_y,ind_S_x]),\
          np.percentile(taub_SSA_ConstantT[ind_S_y,ind_S_x],q/2),\
          np.percentile(taub_SSA_ConstantT[ind_S_y,ind_S_x],100-q/2)])*1e3
  
  # Get values for SSA-ModelT
  region_U_taub_SSA_ModelT = np.array([np.mean(taub_SSA_ModelT[ind_U_y,ind_U_x]),\
          np.percentile(taub_SSA_ModelT[ind_U_y,ind_U_x],q/2),\
          np.percentile(taub_SSA_ModelT[ind_U_y,ind_U_x],100-q/2)])*1e3
  region_M_taub_SSA_ModelT = np.array([np.mean(taub_SSA_ModelT[ind_M_y,ind_M_x]),\
          np.percentile(taub_SSA_ModelT[ind_M_y,ind_M_x],q/2),\
          np.percentile(taub_SSA_ModelT[ind_M_y,ind_M_x],100-q/2)])*1e3
  region_L_taub_SSA_ModelT = np.array([np.mean(taub_SSA_ModelT[ind_L_y,ind_L_x]),\
          np.percentile(taub_SSA_ModelT[ind_L_y,ind_L_x],q/2),\
          np.percentile(taub_SSA_ModelT[ind_L_y,ind_L_x],100-q/2)])*1e3
  region_S_taub_SSA_ModelT = np.array([np.mean(taub_SSA_ModelT[ind_S_y,ind_S_x]),\
          np.percentile(taub_SSA_ModelT[ind_S_y,ind_S_x],q/2),\
          np.percentile(taub_SSA_ModelT[ind_S_y,ind_S_x],100-q/2)])*1e3
  
  region_U_taud = np.mean(tauds_flow[ind_U_y,ind_U_x])/1e3
  region_M_taud = np.mean(tauds_flow[ind_M_y,ind_M_x])/1e3  
  region_L_taud = np.mean(tauds_flow[ind_L_y,ind_L_x])/1e3
  region_S_taud = np.mean(tauds_flow[ind_S_y,ind_S_x])/1e3
  
  ###################
  # Plot comparison #
  ###################

  # Plot overview map
  fig = plt.figure(figsize=(2.4,3))
  ax = plt.gca()
  plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
  #p = plt.imshow(tauds_flow/1e3,extent=[xb[0],xb[-1],yb[0],yb[-1]],origin='lower',vmin=-500,vmax=500,cmap='RdBu_r')
  plt.xlim([xmin,xmax])
  plt.ylim([ymin,ymax])
  plt.xticks([])
  plt.yticks([])
  xmin,xmax = plt.xlim()
  ymin,ymax = plt.ylim()
  #plt.plot(extent[:,0],extent[:,1],'k',lw=2)
  patch = matplotlib.patches.PathPatch(region_L,edgecolor='k',facecolor='b',lw=2,alpha=0.5)
  ax.add_patch(patch)
  patch = matplotlib.patches.PathPatch(region_S1,edgecolor='k',facecolor='k',lw=2,alpha=0.5)
  ax.add_patch(patch)
  patch = matplotlib.patches.PathPatch(region_S2,edgecolor='k',facecolor='k',lw=2,alpha=0.5)
  ax.add_patch(patch)
  patch = matplotlib.patches.PathPatch(region_U,edgecolor='k',facecolor='r',lw=2,alpha=0.5)
  ax.add_patch(patch)
  patch = matplotlib.patches.PathPatch(region_M,edgecolor='k',facecolor='w',lw=2,alpha=0.5)
  ax.add_patch(patch)
  path = matplotlib.path.Path([[0.6*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.98*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.98*(xmax-xmin)+xmin,0.90*(ymax-ymin)+ymin],
                        [0.6*(xmax-xmin)+xmin,0.90*(ymax-ymin)+ymin],
                        [0.6*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin]])
  patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1)
  ax.add_patch(patch)
  #cbaxes = fig.add_axes([0.60, 0.905, 0.25, 0.035])
  #cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=np.arange(-500,600,500))
  #ax.text(xmin+0.51*(xmax-xmin),ymin+0.805*(ymax-ymin),r'$\tau_d \cdot u / ||u||$ (kPa)',fontsize=10)
  #cb.ax.tick_params(labelsize=10)
  ax.plot([xmin+0.63*(xmax-xmin),xmin+0.63*(xmax-xmin)+5e3],[ymin+0.94*(ymax-ymin),ymin+0.94*(ymax-ymin)],'k',linewidth=1.5)
  ax.plot([xmin+0.63*(xmax-xmin),xmin+0.63*(xmax-xmin)],[ymin+0.94*(ymax-ymin),ymin+0.92*(ymax-ymin)],'k',linewidth=1.5)
  ax.plot([xmin+0.63*(xmax-xmin)+5e3,xmin+0.63*(xmax-xmin)+5e3],[ymin+0.94*(ymax-ymin),ymin+0.92*(ymax-ymin)],'k',linewidth=1.5)
  ax.text(xmin+0.65*(xmax-xmin)+5e3,ymin+0.92*(ymax-ymin),'5 km',fontsize=10)
  
  ax.text(xmin+0.63*(xmax-xmin),ymin+0.2*(ymax-ymin),'L',fontweight='bold',fontsize=10,backgroundcolor='w')
  ax.text(xmin+0.18*(xmax-xmin),ymin+0.4*(ymax-ymin),'S',fontweight='bold',fontsize=10,backgroundcolor='w')
  ax.text(xmin+0.35*(xmax-xmin),ymin+0.8*(ymax-ymin),'S',fontweight='bold',fontsize=10,backgroundcolor='w')
  ax.text(xmin+0.23*(xmax-xmin),ymin+0.6*(ymax-ymin),'U',fontweight='bold',fontsize=10,backgroundcolor='w')
  ax.text(xmin+0.5*(xmax-xmin),ymin+0.33*(ymax-ymin),'M',fontweight='bold',fontsize=10,backgroundcolor='w') 

  plt.tight_layout()
  plt.subplots_adjust(hspace=0.03,wspace=0.03,top=0.99,right=0.99,left=0.01,bottom=0.01)
  plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_inversion_regions_map.pdf"),FORMAT='PDF',dpi=600)
  plt.close()
  
  
  # Plot bar graph
  fig = plt.figure(figsize=(3,3))
  ax = plt.gca()
  width=0.15
  
  ind = 0.7+np.arange(4)
  plt.plot([0.5,5],[1,1],'k--',lw=2)

  # FS-CT
  ax.bar(ind[0],region_L_taub_FS_ConstantT[0]/region_L_taud,width,color='r',edgecolor='k',\
       yerr=[[region_L_taub_FS_ConstantT[0]/region_L_taud-region_L_taub_FS_ConstantT[1]/region_L_taud],\
       [region_L_taub_FS_ConstantT[2]/region_L_taud-region_L_taub_FS_ConstantT[0]/region_L_taud]],\
       error_kw=dict(ecolor='k',capsize=2),label='FS-CT')
  ax.bar(ind[1],region_M_taub_FS_ConstantT[0]/region_M_taud,width,color='r',edgecolor='k',\
       yerr=[[region_M_taub_FS_ConstantT[0]/region_M_taud-region_M_taub_FS_ConstantT[1]/region_M_taud],\
       [region_M_taub_FS_ConstantT[2]/region_M_taud-region_M_taub_FS_ConstantT[0]/region_M_taud]],\
       error_kw=dict(ecolor='k',capsize=2))
  ax.bar(ind[2],region_U_taub_FS_ConstantT[0]/region_U_taud,width,color='r',edgecolor='k',\
       yerr=[[region_U_taub_FS_ConstantT[0]/region_U_taud-region_U_taub_FS_ConstantT[1]/region_U_taud],\
       [region_U_taub_FS_ConstantT[2]/region_U_taud-region_U_taub_FS_ConstantT[0]/region_U_taud]],\
       error_kw=dict(ecolor='k',capsize=2))
  ax.bar(ind[3],region_S_taub_FS_ConstantT[0]/region_S_taud,width,color='r',edgecolor='k',\
       yerr=[[region_S_taub_FS_ConstantT[0]/region_S_taud-region_S_taub_FS_ConstantT[1]/region_S_taud],\
       [region_S_taub_FS_ConstantT[2]/region_S_taud-region_S_taub_FS_ConstantT[0]/region_S_taud]],\
       error_kw=dict(ecolor='k',capsize=2))

  # FS-MT
  ax.bar(ind[0]+width,region_L_taub_FS_ModelT[0]/region_L_taud,width,color=(0.9,0.7,0.7),edgecolor='k',\
       yerr=[[region_L_taub_FS_ModelT[0]/region_L_taud-region_L_taub_FS_ModelT[1]/region_L_taud],\
       [region_L_taub_FS_ModelT[2]/region_L_taud-region_L_taub_FS_ModelT[0]/region_L_taud]],\
       error_kw=dict(ecolor='k',capsize=2),label='FS-MT')
  ax.bar(ind[1]+width,region_M_taub_FS_ModelT[0]/region_M_taud,width,color=(0.9,0.7,0.7),edgecolor='k',\
       yerr=[[region_M_taub_FS_ModelT[0]/region_M_taud-region_M_taub_FS_ModelT[1]/region_M_taud],\
       [region_M_taub_FS_ModelT[2]/region_M_taud-region_M_taub_FS_ModelT[0]/region_M_taud]],\
       error_kw=dict(ecolor='k',capsize=2))
  ax.bar(ind[2]+width,region_U_taub_FS_ModelT[0]/region_U_taud,width,color=(0.9,0.7,0.7),edgecolor='k',\
       yerr=[[region_U_taub_FS_ModelT[0]/region_U_taud-region_U_taub_FS_ModelT[1]/region_U_taud],\
       [region_U_taub_FS_ModelT[2]/region_U_taud-region_U_taub_FS_ModelT[0]/region_U_taud]],\
       error_kw=dict(ecolor='k',capsize=2))
  ax.bar(ind[3]+width,region_S_taub_FS_ModelT[0]/region_S_taud,width,color=(0.9,0.7,0.7),edgecolor='k',\
       yerr=[[region_S_taub_FS_ModelT[0]/region_S_taud-region_S_taub_FS_ModelT[1]/region_S_taud],\
       [region_S_taub_FS_ModelT[2]/region_S_taud-region_S_taub_FS_ModelT[0]/region_S_taud]],\
       error_kw=dict(ecolor='k',capsize=2))

  # SSA-CT
  ax.bar(ind[0]+2*width,region_L_taub_SSA_ConstantT[0]/region_L_taud,width,color='b',edgecolor='k',\
       yerr=[[region_L_taub_SSA_ConstantT[0]/region_L_taud-region_L_taub_SSA_ConstantT[1]/region_L_taud],\
       [region_L_taub_SSA_ConstantT[2]/region_L_taud-region_L_taub_SSA_ConstantT[0]/region_L_taud]],\
       error_kw=dict(ecolor='k',capsize=2),label='SSA-CT')
  ax.bar(ind[1]+2*width,region_M_taub_SSA_ConstantT[0]/region_M_taud,width,color='b',edgecolor='k',\
       yerr=[[region_M_taub_SSA_ConstantT[0]/region_M_taud-region_M_taub_SSA_ConstantT[1]/region_M_taud],\
       [region_M_taub_SSA_ConstantT[2]/region_M_taud-region_M_taub_SSA_ConstantT[0]/region_L_taud]],\
       error_kw=dict(ecolor='k',capsize=2))
  ax.bar(ind[2]+2*width,region_U_taub_SSA_ConstantT[0]/region_U_taud,width,color='b',edgecolor='k',\
       yerr=[[region_U_taub_SSA_ConstantT[0]/region_U_taud-region_U_taub_SSA_ConstantT[1]/region_U_taud],\
       [region_U_taub_SSA_ConstantT[2]/region_U_taud-region_U_taub_SSA_ConstantT[0]/region_U_taud]],\
       error_kw=dict(ecolor='k',capsize=2))
  ax.bar(ind[3]+2*width,region_S_taub_SSA_ConstantT[0]/region_S_taud,width,color='b',edgecolor='k',\
       yerr=[[region_S_taub_SSA_ConstantT[0]/region_S_taud-region_S_taub_SSA_ConstantT[1]/region_S_taud],
       [region_S_taub_SSA_ConstantT[2]/region_S_taud-region_S_taub_SSA_ConstantT[0]/region_S_taud]],\
       error_kw=dict(ecolor='k',capsize=2))

  # SSA-MT
  ax.bar(ind[0]+3*width,region_L_taub_SSA_ModelT[0]/region_L_taud,width,color=(0.7,0.7,0.9),edgecolor='k',\
       yerr=[[region_L_taub_SSA_ModelT[0]/region_L_taud-region_L_taub_SSA_ModelT[1]/region_L_taud],\
       [region_L_taub_SSA_ModelT[2]/region_L_taud-region_L_taub_SSA_ModelT[0]/region_L_taud]],\
       error_kw=dict(ecolor='k',capsize=2),label='SSA-MT')
  ax.bar(ind[1]+3*width,region_M_taub_SSA_ModelT[0]/region_M_taud,width,color=(0.7,0.7,0.9),edgecolor='k',\
       yerr=[[region_M_taub_SSA_ModelT[0]/region_M_taud-region_M_taub_SSA_ModelT[1]/region_M_taud],\
       [region_M_taub_SSA_ModelT[2]/region_M_taud-region_M_taub_SSA_ModelT[0]/region_M_taud]],\
       error_kw=dict(ecolor='k',capsize=2))
  ax.bar(ind[2]+3*width,region_U_taub_SSA_ModelT[0]/region_U_taud,width,color=(0.7,0.7,0.9),edgecolor='k',\
       yerr=[[region_U_taub_SSA_ModelT[0]/region_U_taud-region_U_taub_SSA_ModelT[1]/region_U_taud],\
       [region_U_taub_SSA_ModelT[2]/region_U_taud-region_U_taub_SSA_ModelT[0]/region_U_taud]],\
       error_kw=dict(ecolor='k',capsize=2))
  ax.bar(ind[3]+3*width,region_S_taub_SSA_ModelT[0]/region_S_taud,width,color=(0.7,0.7,0.9),edgecolor='k',\
       yerr=[[region_S_taub_SSA_ModelT[0]/region_S_taud-region_S_taub_SSA_ModelT[1]/region_S_taud],\
       [region_S_taub_SSA_ModelT[2]/region_S_taud-region_S_taub_SSA_ModelT[0]/region_S_taud]],\
       error_kw=dict(ecolor='k',capsize=2))
  
  plt.xlim([0.5,4.5])
  plt.yticks(np.arange(0.0,1.3,0.2),fontsize=12,fontname='Arial')
  plt.ylabel(r'$\tau_b/\tau_d$',fontname='Arial',fontsize=10)
  plt.xticks([1,2,3],fontsize=12,fontname='Arial')
  ax.set_xticklabels(('L','U','S'))
  plt.legend(loc=2,fontsize=10,numpoints=1,ncol=1,handlelength=0.4,handletextpad=0.2,labelspacing=0.2,\
       columnspacing=0.2,framealpha=1)
  
  plt.subplots_adjust(hspace=0.02,wspace=0.02,top=0.96,right=0.98,left=0.2,bottom=0.09)
  plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_inversion_regions_"+date+"_stress.pdf"),FORMAT='PDF',dpi=600)
  plt.close()



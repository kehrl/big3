import sys
import os
import matplotlib.pyplot as plt
from pylab import *
import argparse
import scipy
import datelib, inverselib, elmerreadlib, glaclib, geotifflib, meshlib
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

args, _ = parser.parse_known_args(sys.argv)
glacier = args.glacier

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
  if glacier == 'Helheim':
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

  if glacier == 'Kanger':
    extent_region_L = meshlib.shp_to_xy(os.path.join(os.getenv("DATA_HOME"),\
          "ShapeFiles/Glaciers/3D/"+glacier+"/glacier_region2.shp"))
    extent_region_L = np.row_stack([extent_region_L,extent_region_L[0,:]])
    extent_region_M = meshlib.shp_to_xy(os.path.join(os.getenv("DATA_HOME"),\
          "ShapeFiles/Glaciers/3D/"+glacier+"/glacier_region3.shp"))
    extent_region_M = np.row_stack([extent_region_M,extent_region_M[0,:]])
    extent_region_U = meshlib.shp_to_xy(os.path.join(os.getenv("DATA_HOME"),\
          "ShapeFiles/Glaciers/3D/"+glacier+"/glacier_region4.shp"))
    extent_region_U = np.row_stack([extent_region_U,extent_region_U[0,:]])
    extent_region_S = meshlib.shp_to_xy(os.path.join(os.getenv("DATA_HOME"),\
          "ShapeFiles/Glaciers/3D/"+glacier+"/glacier_region1.shp"))
    extent_region_S = np.row_stack([extent_region_S,extent_region_S[0,:]])
    
    region_S = Path(np.column_stack((extent_region_S[:,0],extent_region_S[:,1])))
    region_L = Path(np.column_stack((extent_region_L[:,0],extent_region_L[:,1])))
    region_M = Path(np.column_stack((extent_region_M[:,0],extent_region_M[:,1])))
    region_U = Path(np.column_stack((extent_region_U[:,0],extent_region_U[:,1])))

# Get mesh extent
if glacier == 'Kanger':
    extent = np.loadtxt(os.path.join(os.getenv("MODEL_HOME"),"Kanger/3D/INV_SSA_ModelT/DEM20120522_modelT/inputs/mesh_extent.dat"))
elif glacier == 'Helheim':
    extent = np.loadtxt(os.path.join(os.getenv("MODEL_HOME"),"Helheim/3D/INV_SSA_ModelT/DEM20040803_modelT/inputs/mesh_extent.dat"))

# Get indices where velocity remains above a certain cutoff 
cutoff = 1000
x_cutoff,y_cutoff,vsurfini_cutoff,ind_cutoff_grid,ind_cutoff  = inverselib.get_velocity_cutoff(glacier,velocity_cutoff=cutoff,\
        SSA=False,model_dir='INV_FS_ModelT')

################
# Load Results #
################

# Get inversion dates
dates = []
files = os.listdir(DIR_FS_ModelT)
for file in files:
  if file.endswith('_bed_mod_taub.tif'):
      dates.append(file[3:11])

# Get dimensions of TIF files
x,y,test = geotifflib.read(DIR_SSA_ModelT+'DEM'+dates[0]+'_bed_mod_taub.tif')

nt = len(dates)
nx = len(x)
ny = len(y)
del test, x, y

# Convert dates to time
times = np.zeros([nt,1])
for i in range(0,nt):
    times[i] = datelib.date_to_fracyear(float(dates[i][0:4]),float(dates[i][4:6]),float(dates[i][6:]))

# Set up output variables
taub_FS_ModelT = np.zeros([ny,nx,nt])
taub_FS_ConstantT = np.zeros([ny,nx,nt])
taub_SSA_ModelT = np.zeros([ny,nx,nt])
taub_SSA_ConstantT = np.zeros([ny,nx,nt])

beta_FS_ModelT = np.zeros([ny,nx,nt])
beta_FS_ConstantT = np.zeros([ny,nx,nt])
beta_SSA_ModelT = np.zeros([ny,nx,nt])
beta_SSA_ConstantT = np.zeros([ny,nx,nt])

ub_FS_ModelT = np.zeros([ny,nx,nt])
vb_FS_ModelT = np.zeros([ny,nx,nt])
us_FS_ModelT = np.zeros([ny,nx,nt])
vs_FS_ModelT = np.zeros([ny,nx,nt])

ub_FS_ConstantT = np.zeros([ny,nx,nt])
vb_FS_ConstantT = np.zeros([ny,nx,nt])
us_FS_ConstantT = np.zeros([ny,nx,nt])
vs_FS_ConstantT = np.zeros([ny,nx,nt])

zs = np.zeros([ny,nx,nt])
zb = np.zeros([ny,nx,nt])
umea = np.zeros([ny,nx,nt])
vmea = np.zeros([ny,nx,nt])

# Load TIF files
i = 0
for date in dates:
    x,y,taub_FS_ModelT[:,:,i] = geotifflib.read(DIR_FS_ModelT+'DEM'+date+'_bed_mod_taub.tif') 
    x,y,taub_FS_ConstantT[:,:,i] = geotifflib.read(DIR_FS_ConstantT+'DEM'+date+'_bed_mod_taub.tif')
    x,y,taub_SSA_ModelT[:,:,i] = geotifflib.read(DIR_SSA_ModelT+'DEM'+date+'_bed_mod_taub.tif')
    x,y,taub_SSA_ConstantT[:,:,i] = geotifflib.read(DIR_SSA_ConstantT+'DEM'+date+'_bed_mod_taub.tif')

    x,y,beta_FS_ModelT[:,:,i] = geotifflib.read(DIR_FS_ModelT+'DEM'+date+'_bed_mod_beta.tif')
    x,y,beta_FS_ConstantT[:,:,i] = geotifflib.read(DIR_FS_ConstantT+'DEM'+date+'_bed_mod_beta.tif')
    x,y,beta_SSA_ModelT[:,:,i] = geotifflib.read(DIR_SSA_ModelT+'DEM'+date+'_bed_mod_beta.tif')
    x,y,beta_SSA_ConstantT[:,:,i] = geotifflib.read(DIR_SSA_ConstantT+'DEM'+date+'_bed_mod_beta.tif')

    x,y,ub_FS_ModelT[:,:,i] = geotifflib.read(DIR_FS_ModelT+'DEM'+date+'_bed_mod_ub.tif')
    x,y,vb_FS_ModelT[:,:,i] = geotifflib.read(DIR_FS_ModelT+'DEM'+date+'_bed_mod_vb.tif')
    x,y,us_FS_ModelT[:,:,i] = geotifflib.read(DIR_FS_ModelT+'DEM'+date+'_surf_mod_us.tif')
    x,y,vs_FS_ModelT[:,:,i] = geotifflib.read(DIR_FS_ModelT+'DEM'+date+'_surf_mod_vs.tif')

    x,y,ub_FS_ConstantT[:,:,i] = geotifflib.read(DIR_FS_ConstantT+'DEM'+date+'_bed_mod_ub.tif')
    x,y,vb_FS_ConstantT[:,:,i] = geotifflib.read(DIR_FS_ConstantT+'DEM'+date+'_bed_mod_vb.tif')
    x,y,us_FS_ConstantT[:,:,i] = geotifflib.read(DIR_FS_ConstantT+'DEM'+date+'_surf_mod_us.tif')
    x,y,vs_FS_ConstantT[:,:,i] = geotifflib.read(DIR_FS_ConstantT+'DEM'+date+'_surf_mod_vs.tif')

    x,y,zs[:,:,i] = geotifflib.read(DIR_FS_ConstantT+'DEM'+date+'_surf_mea_zs.tif')
    x,y,zb[:,:,i] = geotifflib.read(DIR_FS_ConstantT+'DEM'+date+'_bed_mod_zb.tif')
    x,y,umea[:,:,i] = geotifflib.read(DIR_FS_ConstantT+'DEM'+date+'_surf_mea_us.tif')
    x,y,vmea[:,:,i] = geotifflib.read(DIR_FS_ConstantT+'DEM'+date+'_surf_mea_vs.tif')

    i = i+1

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
        ximage,yimage,image = geotifflib.readrgb(os.path.join(os.getenv("DATA_HOME"),\
                "Imagery/Landsat/Helheim/TIF/20140704140535_LC82330132014185LGN00.tif"))
    elif glacier == "Kanger":
        ximage,yimage,image = geotifflib.readrgb(os.path.join(os.getenv("DATA_HOME"),\
                "Imagery/Landsat/Kanger/TIF/20140706135251_LC82310122014187LGN00.tif"))
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
dhdx = np.zeros([ny,nx,nt])
dhdy = np.zeros([ny,nx,nt])
dhdx[:,:,:] = float('nan')
dhdy[:,:,:] = float('nan')

H = zs-zb

for i in range(0,nt):
    dhdx[1:-1,1:-1,i] = (zs[1:-1,2:,i]-zs[1:-1,0:-2,i])/(x[2:]-x[0:-2])
    dhdy[1:-1,1:-1,i] = ((zs[2:,1:-1,i]-zs[0:-2,1:-1,i]).T/(y[2:]-y[0:-2])).T
#tauds = rho_i*g*H*np.sqrt((dhdx)**2+(dhdy)**2)
tauds_flow = -rho_i*g*H*(dhdx*umea+dhdy*vmea)/np.sqrt(umea**2+vmea**2)

del dhdx,dhdy,rho_i,g,H

#######################
# Plot driving stress #
#######################

for i in range(0,nt):
    date = dates[i]

    fig = plt.figure(figsize=(2.4,3))
    ax = plt.gca()
  
    plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],\
            cmap='Greys_r',origin='lower',clim=[0,0.6])
    p = plt.imshow(tauds_flow[:,:,i]/1e3,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',vmin=-500,vmax=500,cmap='RdBu_r')
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
    plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+'_inversion_'+date+'_taud.pdf'),FORMAT='PDF',dpi=300)
    plt.close()

#################################################################################
# Plot basal shear stress for different forward models and temperature profiles #
#################################################################################

for i in range(0,nt):
    date = dates[i]

    fig = plt.figure(figsize=(4.8,5.8))
    matplotlib.rc('font',family='Arial')
    gs = matplotlib.gridspec.GridSpec(2,2)

    plt.subplot(gs[0,0])
    ax = plt.gca()
    plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],\
            cmap='Greys_r',origin='lower',clim=[0,0.6])
    p=plt.imshow(taub_FS_ConstantT[:,:,i]*1e3,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',vmin=0,vmax=500)
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
    ax.text(xmin+0.68*(xmax-xmin),ymin+0.78*(ymax-ymin),r'$\tau_b$ (kPa)',fontsize=8)
    cb.ax.tick_params(labelsize=8)
    ax.plot([xmin+0.63*(xmax-xmin),xmin+0.63*(xmax-xmin)+5e3],[ymin+0.75*(ymax-ymin),ymin+0.75*(ymax-ymin)],'k',linewidth=1.5)
    ax.plot([xmin+0.63*(xmax-xmin),xmin+0.63*(xmax-xmin)],[ymin+0.75*(ymax-ymin),ymin+0.73*(ymax-ymin)],'k',linewidth=1.5)
    ax.plot([xmin+0.63*(xmax-xmin)+5e3,xmin+0.63*(xmax-xmin)+5e3],[ymin+0.75*(ymax-ymin),ymin+0.73*(ymax-ymin)],'k',linewidth=1.5)
    ax.text(xmin+0.66*(xmax-xmin)+5e3,ymin+0.725*(ymax-ymin),'5 km',fontsize=8)
    ax.text(xmin+0.03*(xmax-xmin),ymin+0.92*(ymax-ymin),'a FS-CT',fontweight='bold',fontsize=8)

    plt.subplot(gs[0,1])
    ax = plt.gca()
    plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],\
            cmap='Greys_r',origin='lower',clim=[0,0.6])
    p=plt.imshow(taub_FS_ModelT[:,:,i]*1e3,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',vmin=0,vmax=500)
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
    ax.text(xmin+0.68*(xmax-xmin),ymin+0.79*(ymax-ymin),r'$\tau_b$ (kPa)',fontsize=8)
    cb.ax.tick_params(labelsize=8)
    ax.text(xmin+0.03*(xmax-xmin),ymin+0.92*(ymax-ymin),'b FS-MT',fontweight='bold',fontsize=8)

    plt.subplot(gs[1,0])
    ax = plt.gca()
    plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],\
          cmap='Greys_r',origin='lower',clim=[0,0.6])
    p=plt.imshow(taub_SSA_ConstantT[:,:,i]*1e3,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',vmin=0,vmax=500)
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
    ax.text(xmin+0.68*(xmax-xmin),ymin+0.79*(ymax-ymin),r'$\tau_b$ (kPa)',fontsize=8)
    cb.ax.tick_params(labelsize=8)
    ax.text(xmin+0.03*(xmax-xmin),ymin+0.92*(ymax-ymin),'c SSA-CT',fontweight='bold',fontsize=8)

    plt.subplot(gs[1,1])
    ax = plt.gca()
    plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],\
            cmap='Greys_r',origin='lower',clim=[0,0.6])
    p=plt.imshow(taub_SSA_ModelT[:,:,i]*1e3,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',vmin=0,vmax=500)
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
    ax.text(xmin+0.68*(xmax-xmin),ymin+0.79*(ymax-ymin),r'$\tau_b$ (kPa)',fontsize=8)
    cb.ax.tick_params(labelsize=8)
    ax.text(xmin+0.03*(xmax-xmin),ymin+0.92*(ymax-ymin),'d SSA-MT',fontweight='bold',fontsize=8)

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.03,wspace=0.03,top=0.99,right=0.99,left=0.01,bottom=0.01)
    plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_inversion_"+date+"_taub.pdf"),FORMAT='PDF',dpi=300)
    plt.close()

############################################################
# Plot depth-averaged temperatures and temperature profile #
############################################################

if glacier == 'Helheim':
    pts = np.array([[301875,-2576310],[290102,-2562290]])
    labels=['C','D']
    colorlabels = ['b','r','g']
    styles = ['-','--',':']
elif glacier == 'Kanger':
    pts = np.array([[486354,-2289805],[478226,-2276526]])
    labels=['A','B']
    colorlabels = ['b','r','g']
    styles = ['-','--',':']

fig = plt.figure(figsize=(3.75,3))
matplotlib.rc('font',family='Arial')
gs = matplotlib.gridspec.GridSpec(1,3)

plt.subplot(gs[0:2])
ax = plt.gca()
plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],\
        cmap='Greys_r',origin='lower',clim=[0,0.6])
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
ax.plot([xmin+0.58*(xmax-xmin),xmin+0.58*(xmax-xmin)+5e3],[ymin+0.78*(ymax-ymin),ymin+0.78*(ymax-ymin)],'k',linewidth=1.5)
ax.plot([xmin+0.58*(xmax-xmin),xmin+0.58*(xmax-xmin)],[ymin+0.78*(ymax-ymin),ymin+0.76*(ymax-ymin)],'k',linewidth=1.5)
ax.plot([xmin+0.58*(xmax-xmin)+5e3,xmin+0.58*(xmax-xmin)+5e3],[ymin+0.78*(ymax-ymin),ymin+0.76*(ymax-ymin)],'k',linewidth=1.5)
ax.text(xmin+0.6*(xmax-xmin)+5e3,ymin+0.75*(ymax-ymin),'5 km',fontsize=8)
if glacier == 'Helheim':
  ax.text(xmin+0.03*(xmax-xmin),ymin+0.93*(ymax-ymin),'(c)',fontweight='bold',fontsize=8)
elif glacier == 'Kanger':
  ax.text(xmin+0.03*(xmax-xmin),ymin+0.93*(ymax-ymin),'(a)',fontweight='bold',fontsize=8)

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
plt.ylabel('Depth (m)',fontsize=8,fontname='Arial')
plt.xlabel(r'Temperature ($^o$C)',fontsize=8,fontname='Arial')
plt.xticks(np.arange(-15,5,5),fontname='Arial',fontsize=8)
plt.yticks(np.arange(0,1600,250),fontname='Arial',fontsize=8)
if glacier == 'Helheim':
  plt.text(-15,120,'(d)',fontsize=8,fontweight='bold')
elif glacier == 'Kanger':
  plt.text(-15,120,'(b)',fontsize=8,fontweight='bold')

plt.xlim([-16,0])
plt.ylim([1600,0])

plt.tight_layout()
plt.subplots_adjust(hspace=0.02,wspace=0.02,top=0.97,right=0.86,left=0.01,bottom=0.13)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_temperature.pdf"),FORMAT='PDF',dpi=300)
plt.close()

######################
# Plot sliding ratio #
######################

if glacier == 'Kanger':
    xmin = 465000; xmax = 464000+33e3
    ymin = -2296000; ymax = -2296000+32e3
    top = 1.65
    topd = 0
    labels = ['(a)','(b)']
elif glacier == 'Helheim':
    xmin = 277000; xmax = 277000+33e3
    ymin = -2585500; ymax = -2585500+42e3
    top = 2.0
    topd = 0.01
    labels = ['(c)','(d)']

for i in range(0,nt):
    date = dates[i]

    # Make plot
    fig = plt.figure(figsize=(3.2,top))
    matplotlib.rc('font',family='Arial')
    gs = matplotlib.gridspec.GridSpec(1,2)
    plt.subplot(gs[0])
    ax = plt.gca()
    plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],\
            cmap='Greys_r',origin='lower',clim=[0,0.6])
    p=plt.imshow(np.sqrt(ub_FS_ConstantT[:,:,i]**2+vb_FS_ConstantT[:,:,i]**2)/np.sqrt(us_FS_ConstantT[:,:,i]**2+\
            vs_FS_ConstantT[:,:,i]**2),extent=[x[0],x[-1],y[0],y[-1]],origin='lower',vmin=0,vmax=1,cmap='viridis_r')
    plt.plot(np.r_[extent[:,0],extent[0,0]],np.r_[extent[:,1],extent[0,1]],'k')
    plt.xticks([])
    plt.yticks([])
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    if glacier == 'Kanger':
        path = matplotlib.path.Path([[0.52*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.98*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.98*(xmax-xmin)+xmin,(0.62+topd*4)*(ymax-ymin)+ymin],
                        [0.52*(xmax-xmin)+xmin,(0.62+topd*4)*(ymax-ymin)+ymin],
                        [0.52*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin]])
        patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1)
        ax.add_patch(patch)
        cbaxes = fig.add_axes([0.30, 0.90+topd, 0.15, 0.03])
        cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=np.arange(0,1.5,0.5))
        ax.text(xmin+0.65*(xmax-xmin),ymin+(0.725+topd*2)*(ymax-ymin),r'$u_b/u_s$',fontsize=8)
        cb.ax.tick_params(labelsize=8)
        ax.plot([xmin+0.58*(xmax-xmin),xmin+0.58*(xmax-xmin)+5e3],\
                [ymin+(0.68+topd*4)*(ymax-ymin),ymin+(0.68+topd*4)*(ymax-ymin)],'k',linewidth=1.5)
        ax.plot([xmin+0.58*(xmax-xmin),xmin+0.58*(xmax-xmin)],\
                [ymin+(0.68+topd*4)*(ymax-ymin),ymin+(0.66+topd*4)*(ymax-ymin)],'k',linewidth=1.5)
        ax.plot([xmin+0.58*(xmax-xmin)+5e3,xmin+0.58*(xmax-xmin)+5e3],\
                [ymin+(topd*4+0.68)*(ymax-ymin),ymin+(topd*4+0.66)*(ymax-ymin)],'k',linewidth=1.5)
        ax.text(xmin+0.62*(xmax-xmin)+5e3,ymin+(topd*4+0.65)*(ymax-ymin),'5 km',fontsize=8)
    ax.text(xmin+0.03*(xmax-xmin),ymin+0.92*(ymax-ymin),labels[0],fontweight='bold',fontsize=8)
    ax.text(xmin+0.15*(xmax-xmin),ymin+0.92*(ymax-ymin),'FS-CT',fontsize=8)
    plt.subplot(gs[1])
    ax = plt.gca()
    plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],\
            cmap='Greys_r',origin='lower',clim=[0,0.6])
    p=plt.imshow(np.sqrt(ub_FS_ModelT[:,:,i]**2+vb_FS_ModelT[:,:,i]**2)/np.sqrt(us_FS_ModelT[:,:,i]**2+\
            vs_FS_ModelT[:,:,i]**2),extent=[x[0],x[-1],y[0],y[-1]],origin='lower',vmin=0,vmax=1,cmap='viridis_r')
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    plt.xticks([])
    plt.yticks([])
    plt.plot(np.r_[extent[:,0],extent[0,0]],np.r_[extent[:,1],extent[0,1]],'k')
    ax.text(xmin+0.03*(xmax-xmin),ymin+0.92*(ymax-ymin),labels[1],fontweight='bold',fontsize=8)
    ax.text(xmin+0.15*(xmax-xmin),ymin+0.92*(ymax-ymin),'FS-MT',fontsize=8)
    plt.subplots_adjust(wspace=0.02,hspace=0.02,top=0.99,right=0.99,left=0.01,bottom=0.01)
    plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_inversion_"+date+"_slidingratio.pdf"),FORMAT='PDF',dpi=300)
    plt.close()

####################################################################################
# Bin taub to see if there is a distance that minimizes differences between models #
####################################################################################

bin_taub = False
if bin_taub:
    # First let's just try the full stokes results
    bins = np.r_[2,np.arange(5,51,5)]

    # Set up output variables
    taub_FS_ConstantT_bins = np.zeros([ny,nx,nt,len(bins)])
    taub_FS_ConstantT_bins[:,:,:,:] = np.float('NaN')
    taub_FS_ModelT_bins = np.zeros([ny,nx,nt,len(bins)])
    taub_FS_ModelT_bins[:,:,:,:] = np.float('NaN')
    taub_SSA_ConstantT_bins = np.zeros([ny,nx,nt,len(bins)])
    taub_SSA_ConstantT_bins[:,:,:,:] = np.float('NaN')
    taub_SSA_ModelT_bins = np.zeros([ny,nx,nt,len(bins)])
    taub_SSA_ModelT_bins[:,:,:,:] = np.float('NaN')

    taud_bins = np.zeros([ny,nx,nt,len(bins)])
    taud_bins[:,:,:,:] = np.float('NaN')

    taub_FSCT_FSMT = np.zeros([nt,len(bins)])
    taub_FSCT_SSACT = np.zeros([nt,len(bins)])
    taub_FSCT_SSAMT = np.zeros([nt,len(bins)])
    taub_FSMT_SSACT = np.zeros([nt,len(bins)])
    taub_FSMT_SSAMT = np.zeros([nt,len(bins)])
    taub_SSACT_SSAMT = np.zeros([nt,len(bins)])

    # Compute spatial averages
    for k in range(0,len(bins)):
        for l in range(0,nt):
            kernel = np.ones([1+2*bins[k],1+2*bins[k]])/((1+2*bins[k])*(1+2*bins[k]))
            taub_FS_ConstantT_bins[:,:,l,k] = scipy.ndimage.convolve(taub_FS_ConstantT[:,:,l],kernel,\
                    mode='constant',cval=np.float('nan'))  
            taub_FS_ModelT_bins[:,:,l,k] = scipy.ndimage.convolve(taub_FS_ModelT[:,:,l],kernel,\
                    mode='constant',cval=np.float('nan'))
            taub_SSA_ConstantT_bins[:,:,l,k] = scipy.ndimage.convolve(taub_SSA_ConstantT[:,:,l],\
                    kernel,mode='constant',cval=np.float('nan'))
            taub_SSA_ModelT_bins[:,:,l,k] = scipy.ndimage.convolve(taub_SSA_ModelT[:,:,l],kernel,\
                    mode='constant',cval=np.float('nan'))
            taud_bins[:,:,l,k] = scipy.ndimage.convolve(tauds_flow[:,:,l],kernel,\
                    mode='constant',cval=np.float('nan'))

    for j in range(0,nt):
        nonnan = np.where(~(np.isnan(taub_FS_ConstantT_bins[:,:,0,-1])))
        for i in range(0,len(bins)):
            taub_FSCT_FSMT[j,i] = np.sqrt(np.mean(abs(taub_FS_ConstantT_bins[nonnan[0],nonnan[1],j,i]-\
                    taub_FS_ModelT_bins[nonnan[0],nonnan[1],j,i])**2)) 
            taub_FSCT_SSACT[j,i] = np.sqrt(np.mean(abs(taub_FS_ConstantT_bins[nonnan[0],nonnan[1],j,i]-\
                    taub_SSA_ConstantT_bins[nonnan[0],nonnan[1],j,i])**2))
            taub_FSCT_SSAMT[j,i] = np.sqrt(np.mean(abs(taub_FS_ConstantT_bins[nonnan[0],nonnan[1],j,i]-\
                    taub_SSA_ModelT_bins[nonnan[0],nonnan[1],j,i])**2))
            taub_FSMT_SSACT[j,i] = np.sqrt(np.mean(abs(taub_FS_ModelT_bins[nonnan[0],nonnan[1],j,i]-\
                    taub_SSA_ConstantT_bins[nonnan[0],nonnan[1],j,i])**2))
            taub_FSMT_SSAMT[j,i] = np.sqrt(np.mean(abs(taub_FS_ModelT_bins[nonnan[0],nonnan[1],j,i]-\
                    taub_SSA_ModelT_bins[nonnan[0],nonnan[1],j,i])**2))
            taub_SSACT_SSAMT[j,i] = np.sqrt(np.mean(abs(taub_SSA_ConstantT_bins[nonnan[0],nonnan[1],j,i]-\
                    taub_SSA_ModelT_bins[nonnan[0],nonnan[1],j,i])**2))

    # Plot RMS vs bin size
    plt.figure(figsize=(3.5,3.5))
    dx = x[1]-x[0]
    plt.plot(bins*dx*2,np.mean(taub_FSCT_FSMT,axis=0)*1e3,'k',label='FS-CT,   FS-MT')
    plt.plot(bins*dx*2,np.mean(taub_FSCT_SSACT,axis=0)*1e3,'r',label='FS-CT,   SSA-CT')
    plt.plot(bins*dx*2,np.mean(taub_FSCT_SSAMT,axis=0)*1e3,'b',label='FS-CT,   SSA-MT')
    plt.plot(bins*dx*2,np.mean(taub_FSMT_SSACT,axis=0)*1e3,'g',label='FS-MT,   SSA-CT')
    plt.plot(bins*dx*2,np.mean(taub_FSMT_SSAMT,axis=0)*1e3,'m',label='FS-MT,   SSA-MT')
    plt.plot(bins*dx*2,np.mean(taub_SSACT_SSAMT,axis=0)*1e3,'c',label='SSA-CT, SSA-MT')
    p=plt.legend(loc=1,ncol=1,fontsize=8)
    plt.ylabel('RMS Difference (kPa)')
    plt.xlabel('Bin size (m)')
    plt.xlim([0,5000])
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.03,wspace=0.03,top=0.99,right=0.95,left=0.16,bottom=0.13)
    plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_inversion_bindists.pdf"),FORMAT='PDF',dpi=300)
    plt.close()

    # Plot areas
    #i = 150
    #j = 460
    #bin = 8
    #plt.figure(figsize=(3.5,3.5))
    #plt.plot(times,,'ko-',label=r'$\tau_d$')
    #plt.plot(times,taub_FS_ConstantT_bins[j,i,:,bin]*1e3,'ro-',label='FS-CT')
    #plt.plot(times,taub_FS_ModelT_bins[j,i,:,bin]*1e3,'r^--',label='FS-MT')
    #plt.plot(times,taub_SSA_ConstantT_bins[j,i,:,bin]*1e3,'bo-',label='SSA-CT')
    #plt.plot(times,taub_SSA_ModelT_bins[j,i,:,bin]*1e3,'b^--',label='SSA-MT')



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
            if glacier == 'Helheim':
                if (region_S1.contains_point([x[i],y[j]])) or (region_S2.contains_point([x[i],y[j]])):
                    ind_S_x.append(i)
                    ind_S_y.append(j)
            elif glacier == 'Kanger':
                if (region_S.contains_point([x[i],y[j]])):
                    ind_S_x.append(i)
                    ind_S_y.append(j)

    # Set up output variables
    region_U_taub_FS_ConstantT = np.zeros([nt,3])
    region_M_taub_FS_ConstantT = np.zeros([nt,3])
    region_L_taub_FS_ConstantT = np.zeros([nt,3])
    region_S_taub_FS_ConstantT = np.zeros([nt,3])
    region_U_taub_FS_ModelT = np.zeros([nt,3])
    region_M_taub_FS_ModelT = np.zeros([nt,3])
    region_L_taub_FS_ModelT = np.zeros([nt,3])
    region_S_taub_FS_ModelT = np.zeros([nt,3])
    region_U_taub_SSA_ConstantT = np.zeros([nt,3])
    region_M_taub_SSA_ConstantT = np.zeros([nt,3])
    region_L_taub_SSA_ConstantT = np.zeros([nt,3])
    region_S_taub_SSA_ConstantT = np.zeros([nt,3])
    region_U_taub_SSA_ModelT = np.zeros([nt,3])
    region_M_taub_SSA_ModelT = np.zeros([nt,3])
    region_L_taub_SSA_ModelT = np.zeros([nt,3])
    region_S_taub_SSA_ModelT = np.zeros([nt,3])

    region_U_beta_FS_ConstantT = np.zeros([nt,3])
    region_M_beta_FS_ConstantT = np.zeros([nt,3])
    region_L_beta_FS_ConstantT = np.zeros([nt,3])
    region_S_beta_FS_ConstantT = np.zeros([nt,3])
    region_U_beta_FS_ModelT = np.zeros([nt,3])
    region_M_beta_FS_ModelT = np.zeros([nt,3])
    region_L_beta_FS_ModelT = np.zeros([nt,3])
    region_S_beta_FS_ModelT = np.zeros([nt,3])
    region_U_beta_SSA_ConstantT = np.zeros([nt,3])
    region_M_beta_SSA_ConstantT = np.zeros([nt,3])
    region_L_beta_SSA_ConstantT = np.zeros([nt,3])
    region_S_beta_SSA_ConstantT = np.zeros([nt,3])
    region_U_beta_SSA_ModelT = np.zeros([nt,3])
    region_M_beta_SSA_ModelT = np.zeros([nt,3])
    region_L_beta_SSA_ModelT = np.zeros([nt,3])
    region_S_beta_SSA_ModelT = np.zeros([nt,3])

    region_U_taud = np.zeros([nt,1])
    region_M_taud = np.zeros([nt,1])
    region_L_taud = np.zeros([nt,1])
    region_S_taud = np.zeros([nt,1])

    for i in range(0,nt):      
        # Get values for FS-ConstantT
        region_U_taub_FS_ConstantT[i,:] = np.array([np.mean(taub_FS_ConstantT[ind_U_y,ind_U_x,i]),\
            np.percentile(taub_FS_ConstantT[ind_U_y,ind_U_x,i],q/2),\
            np.percentile(taub_FS_ConstantT[ind_U_y,ind_U_x,i],100-q/2)])*1e3
        region_U_beta_FS_ConstantT[i,:] = np.array([np.mean(beta_FS_ConstantT[ind_U_y,ind_U_x,i]),\
            np.percentile(beta_FS_ConstantT[ind_U_y,ind_U_x,i],q/2),\
            np.percentile(beta_FS_ConstantT[ind_U_y,ind_U_x,i],100-q/2)])
        region_M_taub_FS_ConstantT[i,:] = np.array([np.mean(taub_FS_ConstantT[ind_M_y,ind_M_x,i]),\
            np.percentile(taub_FS_ConstantT[ind_M_y,ind_M_x,i],q/2),\
            np.percentile(taub_FS_ConstantT[ind_M_y,ind_M_x,i],100-q/2)])*1e3
        region_M_beta_FS_ConstantT[i,:] = np.array([np.mean(beta_FS_ConstantT[ind_M_y,ind_M_x,i]),\
            np.percentile(beta_FS_ConstantT[ind_M_y,ind_M_x,i],q/2),\
            np.percentile(beta_FS_ConstantT[ind_M_y,ind_M_x,i],100-q/2)])
        region_L_taub_FS_ConstantT[i,:] = np.array([np.mean(taub_FS_ConstantT[ind_L_y,ind_L_x,i]),\
            np.percentile(taub_FS_ConstantT[ind_L_y,ind_L_x,i],q/2),\
            np.percentile(taub_FS_ConstantT[ind_L_y,ind_L_x,i],100-q/2)])*1e3
        region_L_beta_FS_ConstantT[i,:] = np.array([np.mean(beta_FS_ConstantT[ind_L_y,ind_L_x,i]),\
            np.percentile(beta_FS_ConstantT[ind_L_y,ind_L_x,i],q/2),\
            np.percentile(beta_FS_ConstantT[ind_L_y,ind_L_x,i],100-q/2)])
        region_S_taub_FS_ConstantT[i,:] = np.array([np.mean(taub_FS_ConstantT[ind_S_y,ind_S_x,i]),\
            np.percentile(taub_FS_ConstantT[ind_S_y,ind_S_x,i],q/2),\
            np.percentile(taub_FS_ConstantT[ind_S_y,ind_S_x,i],100-q/2)])*1e3
        region_S_beta_FS_ConstantT[i,:] = np.array([np.mean(beta_FS_ConstantT[ind_S_y,ind_S_x,i]),\
            np.percentile(beta_FS_ConstantT[ind_S_y,ind_S_x,i],q/2),\
            np.percentile(beta_FS_ConstantT[ind_S_y,ind_S_x,i],100-q/2)])

        # Get values for FS-ModelT
        region_U_taub_FS_ModelT[i,:] = np.array([np.mean(taub_FS_ModelT[ind_U_y,ind_U_x,i]),\
            np.percentile(taub_FS_ModelT[ind_U_y,ind_U_x,i],q/2),\
            np.percentile(taub_FS_ModelT[ind_U_y,ind_U_x,i],100-q/2)])*1e3
        region_U_beta_FS_ModelT[i,:] = np.array([np.mean(beta_FS_ModelT[ind_U_y,ind_U_x,i]),\
            np.percentile(beta_FS_ModelT[ind_U_y,ind_U_x,i],q/2),\
            np.percentile(beta_FS_ModelT[ind_U_y,ind_U_x,i],100-q/2)])
        region_M_taub_FS_ModelT[i,:] = np.array([np.mean(taub_FS_ModelT[ind_M_y,ind_M_x,i]),\
            np.percentile(taub_FS_ModelT[ind_M_y,ind_M_x,i],q/2),\
            np.percentile(taub_FS_ModelT[ind_M_y,ind_M_x,i],100-q/2)])*1e3
        region_M_beta_FS_ModelT[i,:] = np.array([np.mean(beta_FS_ModelT[ind_M_y,ind_M_x,i]),\
            np.percentile(beta_FS_ModelT[ind_M_y,ind_M_x,i],q/2),\
            np.percentile(beta_FS_ModelT[ind_M_y,ind_M_x,i],100-q/2)])
        region_L_taub_FS_ModelT[i,:] = np.array([np.mean(taub_FS_ModelT[ind_L_y,ind_L_x,i]),\
            np.percentile(taub_FS_ModelT[ind_L_y,ind_L_x,i],q/2),\
            np.percentile(taub_FS_ModelT[ind_L_y,ind_L_x,i],100-q/2)])*1e3
        region_L_beta_FS_ModelT[i,:] = np.array([np.mean(beta_FS_ModelT[ind_L_y,ind_L_x,i]),\
            np.percentile(beta_FS_ModelT[ind_L_y,ind_L_x,i],q/2),\
            np.percentile(beta_FS_ModelT[ind_L_y,ind_L_x,i],100-q/2)])
        region_S_taub_FS_ModelT[i,:] = np.array([np.mean(taub_FS_ModelT[ind_S_y,ind_S_x,i]),\
            np.percentile(taub_FS_ModelT[ind_S_y,ind_S_x,i],q/2),\
            np.percentile(taub_FS_ModelT[ind_S_y,ind_S_x,i],100-q/2)])*1e3
        region_S_beta_FS_ModelT[i,:] = np.array([np.mean(beta_FS_ModelT[ind_S_y,ind_S_x,i]),\
            np.percentile(beta_FS_ModelT[ind_S_y,ind_S_x,i],q/2),\
            np.percentile(beta_FS_ModelT[ind_S_y,ind_S_x,i],100-q/2)])

        # Get values for SSA-ConstantT
        region_U_taub_SSA_ConstantT[i,:] = np.array([np.mean(taub_SSA_ConstantT[ind_U_y,ind_U_x,i]),\
            np.percentile(taub_SSA_ConstantT[ind_U_y,ind_U_x,i],q/2),\
            np.percentile(taub_SSA_ConstantT[ind_U_y,ind_U_x,i],100-q/2)])*1e3
        region_U_beta_SSA_ConstantT[i,:] = np.array([np.mean(beta_SSA_ConstantT[ind_U_y,ind_U_x,i]),\
            np.percentile(beta_SSA_ConstantT[ind_U_y,ind_U_x,i],q/2),\
            np.percentile(beta_SSA_ConstantT[ind_U_y,ind_U_x,i],100-q/2)])
        region_M_taub_SSA_ConstantT[i,:] = np.array([np.mean(taub_SSA_ConstantT[ind_M_y,ind_M_x,i]),\
            np.percentile(taub_SSA_ConstantT[ind_M_y,ind_M_x,i],q/2),\
            np.percentile(taub_SSA_ConstantT[ind_M_y,ind_M_x,i],100-q/2)])*1e3
        region_M_beta_SSA_ConstantT[i,:] = np.array([np.mean(beta_SSA_ConstantT[ind_M_y,ind_M_x,i]),\
            np.percentile(beta_SSA_ConstantT[ind_M_y,ind_M_x,i],q/2),\
            np.percentile(beta_SSA_ConstantT[ind_M_y,ind_M_x,i],100-q/2)])
        region_L_taub_SSA_ConstantT[i,:] = np.array([np.mean(taub_SSA_ConstantT[ind_L_y,ind_L_x,i]),\
            np.percentile(taub_SSA_ConstantT[ind_L_y,ind_L_x,i],q/2),\
            np.percentile(taub_SSA_ConstantT[ind_L_y,ind_L_x,i],100-q/2)])*1e3
        region_L_beta_SSA_ConstantT[i,:] = np.array([np.mean(beta_SSA_ConstantT[ind_L_y,ind_L_x,i]),\
            np.percentile(beta_SSA_ConstantT[ind_L_y,ind_L_x,i],q/2),\
            np.percentile(beta_SSA_ConstantT[ind_L_y,ind_L_x,i],100-q/2)])
        region_S_taub_SSA_ConstantT[i,:] = np.array([np.mean(taub_SSA_ConstantT[ind_S_y,ind_S_x,i]),\
            np.percentile(taub_SSA_ConstantT[ind_S_y,ind_S_x,i],q/2),\
            np.percentile(taub_SSA_ConstantT[ind_S_y,ind_S_x,i],100-q/2)])*1e3
        region_S_beta_SSA_ConstantT[i,:] = np.array([np.mean(beta_SSA_ConstantT[ind_S_y,ind_S_x,i]),\
            np.percentile(beta_SSA_ConstantT[ind_S_y,ind_S_x,i],q/2),\
            np.percentile(beta_SSA_ConstantT[ind_S_y,ind_S_x,i],100-q/2)])

        # Get values for SSA-ModelT
        region_U_taub_SSA_ModelT[i,:] = np.array([np.mean(taub_SSA_ModelT[ind_U_y,ind_U_x,i]),\
            np.percentile(taub_SSA_ModelT[ind_U_y,ind_U_x,i],q/2),\
            np.percentile(taub_SSA_ModelT[ind_U_y,ind_U_x,i],100-q/2)])*1e3
        region_U_beta_SSA_ModelT[i,:] = np.array([np.mean(beta_SSA_ModelT[ind_U_y,ind_U_x,i]),\
            np.percentile(beta_SSA_ModelT[ind_U_y,ind_U_x,i],q/2),\
            np.percentile(beta_SSA_ModelT[ind_U_y,ind_U_x,i],100-q/2)])
        region_M_taub_SSA_ModelT[i,:] = np.array([np.mean(taub_SSA_ModelT[ind_M_y,ind_M_x,i]),\
            np.percentile(taub_SSA_ModelT[ind_M_y,ind_M_x,i],q/2),\
            np.percentile(taub_SSA_ModelT[ind_M_y,ind_M_x,i],100-q/2)])*1e3
        region_M_beta_SSA_ModelT[i,:] = np.array([np.mean(beta_SSA_ModelT[ind_M_y,ind_M_x,i]),\
            np.percentile(beta_SSA_ModelT[ind_M_y,ind_M_x,i],q/2),\
            np.percentile(beta_SSA_ModelT[ind_M_y,ind_M_x,i],100-q/2)])
        region_L_taub_SSA_ModelT[i,:] = np.array([np.mean(taub_SSA_ModelT[ind_L_y,ind_L_x,i]),\
            np.percentile(taub_SSA_ModelT[ind_L_y,ind_L_x,i],q/2),\
            np.percentile(taub_SSA_ModelT[ind_L_y,ind_L_x,i],100-q/2)])*1e3
        region_L_beta_SSA_ModelT[i,:] = np.array([np.mean(beta_SSA_ModelT[ind_L_y,ind_L_x,i]),\
            np.percentile(beta_SSA_ModelT[ind_L_y,ind_L_x,i],q/2),\
            np.percentile(beta_SSA_ModelT[ind_L_y,ind_L_x,i],100-q/2)])
        region_S_taub_SSA_ModelT[i,:] = np.array([np.mean(taub_SSA_ModelT[ind_S_y,ind_S_x,i]),\
            np.percentile(taub_SSA_ModelT[ind_S_y,ind_S_x,i],q/2),\
            np.percentile(taub_SSA_ModelT[ind_S_y,ind_S_x,i],100-q/2)])*1e3
        region_S_beta_SSA_ModelT[i,:] = np.array([np.mean(beta_SSA_ModelT[ind_S_y,ind_S_x,i]),\
            np.percentile(beta_SSA_ModelT[ind_S_y,ind_S_x,i],q/2),\
            np.percentile(beta_SSA_ModelT[ind_S_y,ind_S_x,i],100-q/2)])

        region_U_taud[i] = np.mean(tauds_flow[ind_U_y,ind_U_x,i])/1e3
        region_M_taud[i] = np.mean(tauds_flow[ind_M_y,ind_M_x,i])/1e3  
        region_L_taud[i] = np.mean(tauds_flow[ind_L_y,ind_L_x,i])/1e3
        region_S_taud[i] = np.mean(tauds_flow[ind_S_y,ind_S_x,i])/1e3
  
    ###################
    # Plot comparison #
    ###################

    # Plot overview map
    fig = plt.figure(figsize=(2.4,3))
    ax = plt.gca()
    plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],\
            cmap='Greys_r',origin='lower',clim=[0,0.6])
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    plt.xticks([])
    plt.yticks([])
    xmin,xmax = plt.xlim()
    ymin,ymax = plt.ylim()
    patch = matplotlib.patches.PathPatch(region_L,edgecolor='k',facecolor='c',lw=2,alpha=0.5)
    ax.add_patch(patch)
    if glacier == 'Helheim':
        patch = matplotlib.patches.PathPatch(region_S1,edgecolor='k',facecolor='k',lw=2,alpha=0.5)
        ax.add_patch(patch)
        patch = matplotlib.patches.PathPatch(region_S2,edgecolor='k',facecolor='k',lw=2,alpha=0.5)
        ax.add_patch(patch)
    elif glacier == 'Kanger':
        patch = matplotlib.patches.PathPatch(region_S,edgecolor='k',facecolor='k',lw=2,alpha=0.5)
        ax.add_patch(patch)
    patch = matplotlib.patches.PathPatch(region_U,edgecolor='k',facecolor='g',lw=2,alpha=0.5)
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
    #ax.text(xmin+0.51*(xmax-xmin),ymin+0.805*(ymax-ymin),r'$\tau_d \cdot u / ||u||$ (kPa)',fontsize=8)
    #cb.ax.tick_params(labelsize=8)
    ax.plot([xmin+0.63*(xmax-xmin),xmin+0.63*(xmax-xmin)+5e3],[ymin+0.94*(ymax-ymin),ymin+0.94*(ymax-ymin)],'k',linewidth=1.5)
    ax.plot([xmin+0.63*(xmax-xmin),xmin+0.63*(xmax-xmin)],[ymin+0.94*(ymax-ymin),ymin+0.92*(ymax-ymin)],'k',linewidth=1.5)
    ax.plot([xmin+0.63*(xmax-xmin)+5e3,xmin+0.63*(xmax-xmin)+5e3],[ymin+0.94*(ymax-ymin),ymin+0.92*(ymax-ymin)],'k',linewidth=1.5)
    ax.text(xmin+0.65*(xmax-xmin)+5e3,ymin+0.92*(ymax-ymin),'5 km',fontsize=8)
  
    if glacier == 'Helheim': 
        ax.text(xmin+0.63*(xmax-xmin),ymin+0.2*(ymax-ymin),'L',fontweight='bold',fontsize=8,bbox=dict(boxstyle='square,pad=0.1', fc='white', ec='none'))
        ax.text(xmin+0.18*(xmax-xmin),ymin+0.4*(ymax-ymin),'S',fontweight='bold',fontsize=8,bbox=dict(boxstyle='square,pad=0.1', fc='white', ec='none'))
        ax.text(xmin+0.35*(xmax-xmin),ymin+0.8*(ymax-ymin),'S',fontweight='bold',fontsize=8,bbox=dict(boxstyle='square,pad=0.1', fc='white', ec='none'))
        ax.text(xmin+0.23*(xmax-xmin),ymin+0.6*(ymax-ymin),'U',fontweight='bold',fontsize=8,bbox=dict(boxstyle='square,pad=0.1', fc='white', ec='none'))
        ax.text(xmin+0.5*(xmax-xmin),ymin+0.33*(ymax-ymin),'M',fontweight='bold',fontsize=8,bbox=dict(boxstyle='square,pad=0.1', fc='white', ec='none')) 
    elif glacier == 'Kanger':
        ax.text(xmin+0.51*(xmax-xmin),ymin+0.27*(ymax-ymin),'L',fontweight='bold',fontsize=8,bbox=dict(boxstyle='square,pad=0.1', fc='white', ec='none'))
        ax.text(xmin+0.17*(xmax-xmin),ymin+0.4*(ymax-ymin),'S',fontweight='bold',fontsize=8,bbox=dict(boxstyle='square,pad=0.1', fc='white', ec='none'))
        ax.text(xmin+0.3*(xmax-xmin),ymin+0.55*(ymax-ymin),'U',fontweight='bold',fontsize=8,bbox=dict(boxstyle='square,pad=0.1', fc='white', ec='none'))
        ax.text(xmin+0.39*(xmax-xmin),ymin+0.38*(ymax-ymin),'M',fontweight='bold',fontsize=8,bbox=dict(boxstyle='square,pad=0.1', fc='white', ec='none'))

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.03,wspace=0.03,top=0.99,right=0.99,left=0.01,bottom=0.01)
    plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_inversion_regions_map.pdf"),FORMAT='PDF',dpi=300)
    plt.close()
  
    # Plot bar graphs for taub
    for i in range(0,nt):
        date = dates[i]

        fig = plt.figure(figsize=(3,3))
        ax = plt.gca()
        width=0.15
  
        ind = 0.7+np.arange(4)
        plt.plot([0.5,5],[1,1],'k--',lw=2)

        # FS-CT
        ax.bar(ind[0],region_L_taub_FS_ConstantT[i,0]/region_L_taud[i],width,color='r',edgecolor='k',\
            yerr=[[region_L_taub_FS_ConstantT[i,0]/region_L_taud[i]-region_L_taub_FS_ConstantT[i,1]/region_L_taud[i]],\
            [region_L_taub_FS_ConstantT[i,2]/region_L_taud[i]-region_L_taub_FS_ConstantT[i,0]/region_L_taud[i]]],\
            error_kw=dict(ecolor='k',capsize=2),label='FS-CT')
        ax.bar(ind[1],region_M_taub_FS_ConstantT[i,0]/region_M_taud[i],width,color='r',edgecolor='k',\
            yerr=[[region_M_taub_FS_ConstantT[i,0]/region_M_taud[i]-region_M_taub_FS_ConstantT[i,1]/region_M_taud[i]],\
            [region_M_taub_FS_ConstantT[i,2]/region_M_taud[i]-region_M_taub_FS_ConstantT[i,0]/region_M_taud[i]]],\
            error_kw=dict(ecolor='k',capsize=2))
        ax.bar(ind[2],region_U_taub_FS_ConstantT[i,0]/region_U_taud[i],width,color='r',edgecolor='k',\
            yerr=[[region_U_taub_FS_ConstantT[i,0]/region_U_taud[i]-region_U_taub_FS_ConstantT[i,1]/region_U_taud[i]],\
            [region_U_taub_FS_ConstantT[i,2]/region_U_taud[i]-region_U_taub_FS_ConstantT[i,0]/region_U_taud[i]]],\
            error_kw=dict(ecolor='k',capsize=2))
        ax.bar(ind[3],region_S_taub_FS_ConstantT[i,0]/region_S_taud[i],width,color='r',edgecolor='k',\
            yerr=[[region_S_taub_FS_ConstantT[i,0]/region_S_taud[i]-region_S_taub_FS_ConstantT[i,1]/region_S_taud[i]],\
            [region_S_taub_FS_ConstantT[i,2]/region_S_taud[i]-region_S_taub_FS_ConstantT[i,0]/region_S_taud[i]]],\
            error_kw=dict(ecolor='k',capsize=2))

        # FS-MT
        ax.bar(ind[0]+width,region_L_taub_FS_ModelT[i,0]/region_L_taud[i],width,color=(0.9,0.7,0.7),edgecolor='k',\
            yerr=[[region_L_taub_FS_ModelT[i,0]/region_L_taud[i]-region_L_taub_FS_ModelT[i,1]/region_L_taud[i]],\
            [region_L_taub_FS_ModelT[i,2]/region_L_taud[i]-region_L_taub_FS_ModelT[i,0]/region_L_taud[i]]],\
            error_kw=dict(ecolor='k',capsize=2),label='FS-MT')
        ax.bar(ind[1]+width,region_M_taub_FS_ModelT[i,0]/region_M_taud[i],width,color=(0.9,0.7,0.7),edgecolor='k',\
            yerr=[[region_M_taub_FS_ModelT[i,0]/region_M_taud[i]-region_M_taub_FS_ModelT[i,1]/region_M_taud[i]],\
            [region_M_taub_FS_ModelT[i,2]/region_M_taud[i]-region_M_taub_FS_ModelT[i,0]/region_M_taud[i]]],\
            error_kw=dict(ecolor='k',capsize=2))
        ax.bar(ind[2]+width,region_U_taub_FS_ModelT[i,0]/region_U_taud[i],width,color=(0.9,0.7,0.7),edgecolor='k',\
            yerr=[[region_U_taub_FS_ModelT[i,0]/region_U_taud[i]-region_U_taub_FS_ModelT[i,1]/region_U_taud[i]],\
            [region_U_taub_FS_ModelT[i,2]/region_U_taud[i]-region_U_taub_FS_ModelT[i,0]/region_U_taud[i]]],\
            error_kw=dict(ecolor='k',capsize=2))
        ax.bar(ind[3]+width,region_S_taub_FS_ModelT[i,0]/region_S_taud[i],width,color=(0.9,0.7,0.7),edgecolor='k',\
            yerr=[[region_S_taub_FS_ModelT[i,0]/region_S_taud[i]-region_S_taub_FS_ModelT[i,1]/region_S_taud[i]],\
            [region_S_taub_FS_ModelT[i,2]/region_S_taud[i]-region_S_taub_FS_ModelT[i,0]/region_S_taud[i]]],\
            error_kw=dict(ecolor='k',capsize=2))

        # SSA-CT
        ax.bar(ind[0]+2*width,region_L_taub_SSA_ConstantT[i,0]/region_L_taud[i],width,color='b',edgecolor='k',\
            yerr=[[region_L_taub_SSA_ConstantT[i,0]/region_L_taud[i]-region_L_taub_SSA_ConstantT[i,1]/region_L_taud[i]],\
            [region_L_taub_SSA_ConstantT[i,2]/region_L_taud[i]-region_L_taub_SSA_ConstantT[i,0]/region_L_taud[i]]],\
            error_kw=dict(ecolor='k',capsize=2),label='SSA-CT')
        ax.bar(ind[1]+2*width,region_M_taub_SSA_ConstantT[i,0]/region_M_taud[i],width,color='b',edgecolor='k',\
            yerr=[[region_M_taub_SSA_ConstantT[i,0]/region_M_taud[i]-region_M_taub_SSA_ConstantT[i,1]/region_M_taud[i]],\
            [region_M_taub_SSA_ConstantT[i,2]/region_M_taud[i]-region_M_taub_SSA_ConstantT[i,0]/region_L_taud[i]]],\
            error_kw=dict(ecolor='k',capsize=2))
        ax.bar(ind[2]+2*width,region_U_taub_SSA_ConstantT[i,0]/region_U_taud[i],width,color='b',edgecolor='k',\
            yerr=[[region_U_taub_SSA_ConstantT[i,0]/region_U_taud[i]-region_U_taub_SSA_ConstantT[i,1]/region_U_taud[i]],\
            [region_U_taub_SSA_ConstantT[i,2]/region_U_taud[i]-region_U_taub_SSA_ConstantT[i,0]/region_U_taud[i]]],\
            error_kw=dict(ecolor='k',capsize=2))
        ax.bar(ind[3]+2*width,region_S_taub_SSA_ConstantT[i,0]/region_S_taud[i],width,color='b',edgecolor='k',\
            yerr=[[region_S_taub_SSA_ConstantT[i,0]/region_S_taud[i]-region_S_taub_SSA_ConstantT[i,1]/region_S_taud[i]],
            [region_S_taub_SSA_ConstantT[i,2]/region_S_taud[i]-region_S_taub_SSA_ConstantT[i,0]/region_S_taud[i]]],\
            error_kw=dict(ecolor='k',capsize=2))

        # SSA-MT
        ax.bar(ind[0]+3*width,region_L_taub_SSA_ModelT[i,0]/region_L_taud[i],width,color=(0.7,0.7,0.9),edgecolor='k',\
            yerr=[[region_L_taub_SSA_ModelT[i,0]/region_L_taud[i]-region_L_taub_SSA_ModelT[i,1]/region_L_taud[i]],\
            [region_L_taub_SSA_ModelT[i,2]/region_L_taud[i]-region_L_taub_SSA_ModelT[i,0]/region_L_taud[i]]],\
            error_kw=dict(ecolor='k',capsize=2),label='SSA-MT')
        ax.bar(ind[1]+3*width,region_M_taub_SSA_ModelT[i,0]/region_M_taud[i],width,color=(0.7,0.7,0.9),edgecolor='k',\
            yerr=[[region_M_taub_SSA_ModelT[i,0]/region_M_taud[i]-region_M_taub_SSA_ModelT[i,1]/region_M_taud[i]],\
            [region_M_taub_SSA_ModelT[i,2]/region_M_taud[i]-region_M_taub_SSA_ModelT[i,0]/region_M_taud[i]]],\
            error_kw=dict(ecolor='k',capsize=2))
        ax.bar(ind[2]+3*width,region_U_taub_SSA_ModelT[i,0]/region_U_taud[i],width,color=(0.7,0.7,0.9),edgecolor='k',\
            yerr=[[region_U_taub_SSA_ModelT[i,0]/region_U_taud[i]-region_U_taub_SSA_ModelT[i,1]/region_U_taud[i]],\
            [region_U_taub_SSA_ModelT[i,2]/region_U_taud[i]-region_U_taub_SSA_ModelT[i,0]/region_U_taud[i]]],\
            error_kw=dict(ecolor='k',capsize=2))
        ax.bar(ind[3]+3*width,region_S_taub_SSA_ModelT[i,0]/region_S_taud[i],width,color=(0.7,0.7,0.9),edgecolor='k',\
            yerr=[[region_S_taub_SSA_ModelT[i,0]/region_S_taud[i]-region_S_taub_SSA_ModelT[i,1]/region_S_taud[i]],\
            [region_S_taub_SSA_ModelT[i,2]/region_S_taud[i]-region_S_taub_SSA_ModelT[i,0]/region_S_taud[i]]],\
            error_kw=dict(ecolor='k',capsize=2))
  
        plt.xlim([0.5,4.5])
        plt.yticks(np.arange(0.0,1.3,0.2),fontsize=8,fontname='Arial')
        plt.ylim([0,1.2])
        plt.ylabel(r'$\tau_b/\tau_d$',fontname='Arial',fontsize=8)
        plt.xticks([1,2,3,4],fontsize=8,fontname='Arial')
        ax.set_xticklabels(('L','M','U','S'))
        plt.legend(loc=2,fontsize=8,numpoints=1,ncol=1,handlelength=0.4,handletextpad=0.2,labelspacing=0.2,\
            columnspacing=0.2,framealpha=1)
  
        plt.subplots_adjust(hspace=0.02,wspace=0.02,top=0.96,right=0.98,left=0.2,bottom=0.09)
        plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_inversion_regions_"+date+"_stress.pdf"),FORMAT='PDF',dpi=300)
        plt.close()

    # Plot bar graphs for beta
    for i in range(0,nt):
        date = dates[i]

        fig = plt.figure(figsize=(3,3))
        ax = plt.gca()
        width=0.15

        nd = 0.7+np.arange(4)

        # FS-CT
        ax.bar(ind[0],region_L_beta_FS_ConstantT[i,0],width,color='r',edgecolor='k',\
            yerr=[[region_L_beta_FS_ConstantT[i,0]-region_L_beta_FS_ConstantT[i,1]],\
            [region_L_beta_FS_ConstantT[i,2]-region_L_beta_FS_ConstantT[i,0]]],\
            error_kw=dict(ecolor='k',capsize=2),label='FS-CT')
        ax.bar(ind[1],region_M_beta_FS_ConstantT[i,0],width,color='r',edgecolor='k',\
            yerr=[[region_M_beta_FS_ConstantT[i,0]-region_M_beta_FS_ConstantT[i,1]],\
            [region_M_beta_FS_ConstantT[i,2]-region_M_beta_FS_ConstantT[i,0]]],\
            error_kw=dict(ecolor='k',capsize=2))
        ax.bar(ind[2],region_U_beta_FS_ConstantT[i,0],width,color='r',edgecolor='k',\
            yerr=[[region_U_beta_FS_ConstantT[i,0]-region_U_beta_FS_ConstantT[i,1]],\
            [region_U_beta_FS_ConstantT[i,2]-region_U_beta_FS_ConstantT[i,0]]],\
            error_kw=dict(ecolor='k',capsize=2))
        ax.bar(ind[3],region_S_beta_FS_ConstantT[i,0],width,color='r',edgecolor='k',\
            yerr=[[region_S_beta_FS_ConstantT[i,0]-region_S_beta_FS_ConstantT[i,0]],\
            [region_S_beta_FS_ConstantT[i,2]-region_S_beta_FS_ConstantT[i,0]]],\
            error_kw=dict(ecolor='k',capsize=2))

        # FS-MT 
        ax.bar(ind[0]+width,region_L_beta_FS_ModelT[i,0],width,color=(0.9,0.7,0.7),edgecolor='k',\
            yerr=[[region_L_beta_FS_ModelT[i,0]-region_L_beta_FS_ModelT[i,1]],\
            [region_L_beta_FS_ModelT[i,2]-region_L_beta_FS_ModelT[i,0]]],\
            error_kw=dict(ecolor='k',capsize=2),label='FS-MT')
        ax.bar(ind[1]+width,region_M_beta_FS_ModelT[i,0],width,color=(0.9,0.7,0.7),edgecolor='k',\
            yerr=[[region_M_beta_FS_ModelT[i,0]-region_M_beta_FS_ModelT[i,1]],\
            [region_M_beta_FS_ModelT[i,2]-region_M_beta_FS_ModelT[i,0]]],\
            error_kw=dict(ecolor='k',capsize=2))
        ax.bar(ind[2]+width,region_U_beta_FS_ModelT[i,0],width,color=(0.9,0.7,0.7),edgecolor='k',\
            yerr=[[region_U_beta_FS_ModelT[i,0]-region_U_beta_FS_ModelT[i,1]],\
            [region_U_beta_FS_ModelT[i,2]-region_U_beta_FS_ModelT[i,0]]],\
            error_kw=dict(ecolor='k',capsize=2))
        ax.bar(ind[3]+width,region_S_beta_FS_ModelT[i,0],width,color=(0.9,0.7,0.7),edgecolor='k',\
            yerr=[[region_S_beta_FS_ModelT[i,0]-region_S_beta_FS_ModelT[i,1]],\
            [region_S_beta_FS_ModelT[i,2]-region_S_beta_FS_ModelT[i,0]]],\
            error_kw=dict(ecolor='k',capsize=2))

        # SSA-CT 
        ax.bar(ind[0]+2*width,region_L_beta_SSA_ConstantT[i,0],width,color='b',edgecolor='k',\
            yerr=[[region_L_beta_SSA_ConstantT[i,0]-region_L_beta_SSA_ConstantT[i,1]],\
            [region_L_beta_SSA_ConstantT[i,2]-region_L_beta_SSA_ConstantT[i,0]]],\
            error_kw=dict(ecolor='k',capsize=2),label='SSA-CT')
        ax.bar(ind[1]+2*width,region_M_beta_SSA_ConstantT[i,0],width,color='b',edgecolor='k',\
            yerr=[[region_M_beta_SSA_ConstantT[i,0]-region_M_beta_SSA_ConstantT[i,1]],\
            [region_M_beta_SSA_ConstantT[i,2]-region_M_beta_SSA_ConstantT[i,0]]],\
            error_kw=dict(ecolor='k',capsize=2))
        ax.bar(ind[2]+2*width,region_U_beta_SSA_ConstantT[i,0],width,color='b',edgecolor='k',\
            yerr=[[region_U_beta_SSA_ConstantT[i,0]-region_U_beta_SSA_ConstantT[i,1]],\
            [region_U_beta_SSA_ConstantT[i,2]-region_U_beta_SSA_ConstantT[i,0]]],\
            error_kw=dict(ecolor='k',capsize=2))
        ax.bar(ind[3]+2*width,region_S_beta_SSA_ConstantT[i,0],width,color='b',edgecolor='k',\
            yerr=[[region_S_beta_SSA_ConstantT[i,0]-region_S_beta_SSA_ConstantT[i,1]],\
            [region_S_beta_SSA_ConstantT[i,2]-region_S_beta_SSA_ConstantT[i,0]]],\
            error_kw=dict(ecolor='k',capsize=2))

        # SSA-MT
        ax.bar(ind[0]+3*width,region_L_beta_SSA_ModelT[i,0],width,color=(0.7,0.7,0.9),edgecolor='k',\
            yerr=[[region_L_beta_SSA_ModelT[i,0]-region_L_beta_SSA_ModelT[i,1]],\
            [region_L_beta_SSA_ModelT[i,2]-region_L_beta_SSA_ModelT[i,0]]],\
            error_kw=dict(ecolor='k',capsize=2),label='SSA-MT')
        ax.bar(ind[1]+3*width,region_M_beta_SSA_ModelT[i,0],width,color=(0.7,0.7,0.9),edgecolor='k',\
            yerr=[[region_M_beta_SSA_ModelT[i,0]-region_M_beta_SSA_ModelT[i,1]],\
            [region_M_beta_SSA_ModelT[i,2]-region_M_beta_SSA_ModelT[i,0]]],\
            error_kw=dict(ecolor='k',capsize=2))
        ax.bar(ind[2]+3*width,region_U_beta_SSA_ModelT[i,0],width,color=(0.7,0.7,0.9),edgecolor='k',\
            yerr=[[region_U_beta_SSA_ModelT[i,0]-region_U_beta_SSA_ModelT[i,1]],\
            [region_U_beta_SSA_ModelT[i,2]-region_U_beta_SSA_ModelT[i,0]]],\
            error_kw=dict(ecolor='k',capsize=2))
        ax.bar(ind[3]+3*width,region_S_beta_SSA_ModelT[i,0],width,color=(0.7,0.7,0.9),edgecolor='k',\
            yerr=[[region_S_beta_SSA_ModelT[i,0]-region_S_beta_SSA_ModelT[i,1]],\
            [region_S_beta_SSA_ModelT[i,2]-region_S_beta_SSA_ModelT[i,0]]],\
            error_kw=dict(ecolor='k',capsize=2))

        plt.xlim([0.5,4.5])
        #plt.yticks(np.arange(0.0,1.3,0.2),fontsize=8,fontname='Arial')
        #plt.ylim([0,1.2])
        plt.ylabel(r'$\beta^2$ (MPa yr$^{-1}$ m)',fontname='Arial',fontsize=8)
        plt.xticks([1,2,3,4],fontsize=8,fontname='Arial')
        ax.set_xticklabels(('L','M','U','S'))
        plt.legend(loc=0,fontsize=8,numpoints=1,ncol=2,handlelength=0.4,handletextpad=0.2,labelspacing=0.2,\
            columnspacing=0.5,framealpha=1)

        plt.subplots_adjust(hspace=0.02,wspace=0.02,top=0.96,right=0.98,left=0.26,bottom=0.09)
        plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_inversion_regions_"+date+"_beta.pdf"),FORMAT='PDF',dpi=300)
        plt.close()
  
    #########################################
    # Plot timeseries for different regions #
    #########################################

    fig = plt.figure(figsize=(5,6))
    matplotlib.rc('font',family='Arial')
    gs = matplotlib.gridspec.GridSpec(4,1)

    n = 1
    if glacier == 'Helheim':
        tlim = 20
        m = range(1,len(times))
    elif glacier == 'Kanger':
        tlim = 30
        m = range(1,len(times)-1)
    plt.subplot(gs[0])
    ax = plt.gca()
    plt.plot(times[m],region_L_taud[m]-region_L_taud[n],'ko-',label=r'$\tau_d$')  
    plt.plot(times[m],region_L_taub_FS_ConstantT[m,0]-region_L_taub_FS_ConstantT[n,0],'ro-')
    plt.plot(times[m],region_L_taub_FS_ModelT[m,0]-region_L_taub_FS_ModelT[n,0],'r^--')  
    plt.plot(times[m],region_L_taub_SSA_ConstantT[m,0]-region_L_taub_SSA_ConstantT[n,0],'bo-')
    plt.plot(times[m],region_L_taub_SSA_ModelT[m,0]-region_L_taub_SSA_ModelT[n,0],'b^--')
    ax.set_xticks(range(2000,2017))
    ax.set_xticklabels([])
    plt.xlim([2011,2015])
    plt.ylim([-1*tlim,tlim])
    plt.yticks(np.arange(-1*tlim+5,30,tlim-5),fontname='Arial',fontsize=8)
    ax.tick_params(axis='x',which='both',direction='in')
    plt.text(2011.05,13.5,'L',fontweight='bold',fontsize=8,fontname='Arial')

    plt.subplot(gs[1])
    ax = plt.gca()
    plt.plot(times[m],region_M_taud[m]-region_M_taud[n],'ko-',label=r'$\tau_d$')
    plt.plot(times[m],region_M_taub_FS_ConstantT[m,0]-region_M_taub_FS_ConstantT[n,0],'ro-')
    plt.plot(times[m],region_M_taub_FS_ModelT[m,0]-region_M_taub_FS_ModelT[n,0],'r^--')
    plt.plot(times[m],region_M_taub_SSA_ConstantT[m,0]-region_M_taub_SSA_ConstantT[n,0],'bo-')
    plt.plot(times[m],region_M_taub_SSA_ModelT[m,0]-region_M_taub_SSA_ModelT[n,0],'b^--')
    ax.set_xticks(range(2000,2017))
    ax.set_xticklabels([])
    plt.xlim([2011,2015])
    plt.ylim([-1*tlim,tlim])
    ax.tick_params(axis='x',which='both',direction='in')
    plt.yticks(np.arange(-1*tlim+5,30,tlim-5),fontname='Arial',fontsize=8)
    plt.text(2011.05,13.5,'M',fontweight='bold',fontsize=8,fontname='Arial')

    plt.subplot(gs[2])
    ax = plt.gca()
    plt.plot(times[m],region_U_taud[m]-region_U_taud[n],'ko-',label=r'$\tau_d$')  
    plt.plot(times[m],region_U_taub_FS_ConstantT[m,0]-region_U_taub_FS_ConstantT[n,0],'ro-',label='FS-CT')
    plt.plot(times[m],region_U_taub_FS_ModelT[m,0]-region_U_taub_FS_ModelT[n,0],'r^--',label='FS-MT')
    plt.plot(times[m],region_U_taub_SSA_ConstantT[m,0]-region_U_taub_SSA_ConstantT[n,0],'bo-',label='SSA-CT')
    plt.plot(times[m],region_U_taub_SSA_ModelT[m,0]-region_U_taub_SSA_ModelT[n,0],'b^--',label='SSA-MT')
    ax.set_xticks(range(2000,2017))
    ax.set_xticklabels([])
    plt.xlim([2011,2015])
    if glacier == 'Helheim':
        plt.ylim([-1*tlim,tlim])
    elif glacier == 'Kanger':
        plt.ylim([-1*tlim+10,tlim+10])
    ax.tick_params(axis='x',which='both',direction='in')
    plt.yticks(np.arange(-1*tlim+5,30,tlim-5),fontname='Arial',fontsize=8)
    plt.ylabel(r'$\tau_b$ Deviation (kPa)',fontsize=8,fontname='Arial')
    plt.text(2011.05,13.5,'U',fontweight='bold',fontsize=8,fontname='Arial')
    if glacier == 'Kanger':
        plt.legend(loc=4,fontsize=8,numpoints=1,handlelength=2.5,labelspacing=0.05,ncol=5,columnspacing=0.7,\
                handletextpad=0.3,borderpad=0.25)
  
    plt.subplot(gs[3])
    ax = plt.gca()
    plt.plot(times[m],region_S_taud[m]-region_S_taud[n],'ko-',label=r'$\tau_d$')
    plt.plot(times[m],region_S_taub_FS_ConstantT[m,0]-region_S_taub_FS_ConstantT[n,0],'ro-',label='FS-CT')
    plt.plot(times[m],region_S_taub_FS_ModelT[m,0]-region_S_taub_FS_ModelT[n,0],'r^--',label='FS-MT')
    plt.plot(times[m],region_S_taub_SSA_ConstantT[m,0]-region_S_taub_SSA_ConstantT[n,0],'bo-',label='SSA-CT')
    plt.plot(times[m],region_S_taub_SSA_ModelT[m,0]-region_S_taub_SSA_ModelT[n,0],'b^--',label='SSA-MT')
    plt.ylim([-8,8])
    plt.yticks(np.arange(-5,8,5),fontsize=8,fontname='Arial')
    labels=[]
    for i in range(2000,2017):
        labels.append('Jan \n'+str(i))
    ax.set_xticks(range(2000,2017))
    ax.set_xticklabels(labels,fontsize=8,fontname='Arial')
    ax.tick_params(axis='x',which='both',direction='in')
    plt.xlim([2011,2015])
    plt.text(2011.05,5.5,'S',fontweight='bold',fontsize=8,fontname='Arial')
    if glacier == 'Helheim':
        plt.legend(loc=1,fontsize=8,numpoints=1,handlelength=2.5,labelspacing=0.05,ncol=5,columnspacing=0.7,\
                handletextpad=0.3,borderpad=0.25)

    plt.subplots_adjust(hspace=0.05,wspace=0.05,top=0.96,right=0.96,left=0.15,bottom=0.09)
    plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_inversion_regions_timeseries_stress.pdf"),FORMAT='PDF',dpi=300)
    plt.close()

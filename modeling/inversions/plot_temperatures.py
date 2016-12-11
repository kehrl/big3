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

args, _ = parser.parse_known_args(sys.argv)
glacier = args.glacier

# Model Parameters
method = 'adjoint'
regpar_fs = '5e11'
regpar_ssa = '1e10'

# Boundaries
bbed = 4
bsur = 5

# Directory
DIR_FS = os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/Temperatures_FS/")
DIR_SSA = os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/Temperatures_SSA/")

###############################
# Load shapefiles for regions #
###############################

if glacier == 'Helheim':
  extent_region_U = meshlib.shp_to_xy(os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/Glaciers/3D/"+glacier+"/glacier_region1.shp"))
  extent_region_U = np.row_stack([extent_region_U,extent_region_U[0,:]])
  extent_region_L = meshlib.shp_to_xy(os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/Glaciers/3D/"+glacier+"/glacier_region2.shp"))
  extent_region_L = np.row_stack([extent_region_L,extent_region_L[0,:]])
  extent_region_M1 = meshlib.shp_to_xy(os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/Glaciers/3D/"+glacier+"/glacier_region3.shp"))
  extent_region_M1 = np.row_stack([extent_region_M1,extent_region_M1[0,:]])
  extent_region_M2 = meshlib.shp_to_xy(os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/Glaciers/3D/"+glacier+"/glacier_region4.shp"))
  extent_region_M2 = np.row_stack([extent_region_M2,extent_region_M2[0,:]])
  
  region_U = Path(np.column_stack((extent_region_U[:,0],extent_region_U[:,1])))
  region_L = Path(np.column_stack((extent_region_L[:,0],extent_region_L[:,1])))
  region_M1 = Path(np.column_stack((extent_region_M1[:,0],extent_region_M1[:,1])))
  region_M2 = Path(np.column_stack((extent_region_M2[:,0],extent_region_M2[:,1])))

####################
# Load SSA Results #
####################

# Get directories
dirs = os.listdir(DIR_SSA)
for dir in dirs:
  if dir.endswith('ModelT'):
    dirs2 = os.listdir(DIR_SSA+dir+'/mesh2d/inversion_'+method)
    for dir2 in dirs2:
      if dir2.startswith('lambda_'+regpar_ssa) and not(dir2.endswith('.pdf')):
        modelT_ssa = elmerreadlib.result_file(DIR_SSA+dir+'/mesh2d/','inversion_'+method+'/'+dir2+'/adjoint_beta_ssa.result',['ssavelocity 1','ssavelocity 2','beta','vsurfini 1','vsurfini 2'])
        # Get bed,surface from inputs to calculate driving stress
        #xb,yb,zbgrid = elmerreadlib.input_file(DIR_SSA+dir+"/inputs/zbdem.xy")
        #xb,yb,zsgrid = elmerreadlib.input_file(DIR_SSA+dir+"/inputs/zsdem.xy")
        #xgrid,ygrid,ugrid = elmerreadlib.input_file(DIR_SSA+dir+"/inputs/udem.xy")
        #xgrid,ygrid,vgrid = elmerreadlib.input_file(DIR_SSA+dir+"/inputs/vdem.xy")

    # Mesh boundaries
    #extent = np.loadtxt(DIR_SSA+dir+"/inputs/mesh_extent.dat")
    #try:
    #  hole1 = np.loadtxt(DIR_SSA+dir+"/inputs/mesh_hole1.dat")
    #  hole2 = np.loadtxt(DIR_SSA+dir+"/inputs/mesh_hole2.dat")
    #  holes=[hole1,hole2]
    #except:
    #  holes = []

  #elif dir.endswith('HighT'):
  #  dirs2 = os.listdir(DIR_SSA+dir+'/mesh2d/inversion_'+method)
  #  for dir2 in dirs2:
  #    if dir2.startswith('lambda_'+regpar_ssa) and not(dir2.endswith('.pdf')):
  #      highT_bed = elmerreadlib.saveline_boundary(DIR_SSA+dir+'/mesh2d/inversion_'+method+'/'+dir2+'/',method+"_beta.dat",bbed,['velocity','beta'])
  #      highT_sur = elmerreadlib.saveline_boundary(DIR_SSA+dir+'/mesh2d/inversion_'+method+'/'+dir2+'/',method+"_beta.dat",bsur,['velocity','vsurfini'])
  elif dir.endswith('LowT'):
    dirs2 = os.listdir(DIR_SSA+dir+'/mesh2d/inversion_'+method)
    for dir2 in dirs2:
      if dir2.startswith('lambda_'+regpar_ssa) and not(dir2.endswith('.pdf')):
       lowT_ssa = elmerreadlib.result_file(DIR_SSA+dir+'/mesh2d/','inversion_'+method+'/'+dir2+'/adjoint_beta_ssa.result',['ssavelocity 1','ssavelocity 2','beta','vsurfini 1','vsurfini 2'])


############################
# Load Full-Stokes Results #
############################

# Get directories
dirs = os.listdir(DIR_FS)
for dir in dirs:
  if dir.endswith('ModelT'):
    dirs2 = os.listdir(DIR_FS+dir+'/mesh2d/inversion_'+method)
    for dir2 in dirs2:
      if dir2.startswith('lambda_'+regpar_fs) and not(dir2.endswith('.pdf')):
        modelT_bed_fs = elmerreadlib.saveline_boundary(DIR_FS+dir+'/mesh2d/inversion_'+method+'/'+dir2+'/',method+"_beta.dat",bbed,['velocity','beta'])
        modelT_sur_fs = elmerreadlib.saveline_boundary(DIR_FS+dir+'/mesh2d/inversion_'+method+'/'+dir2+'/',method+"_beta.dat",bsur,['velocity','vsurfini'])
        modelT_fs = elmerreadlib.pvtu_file(DIR_FS+dir+'/mesh2d/inversion_'+method+'/'+dir2+'/'+method+'_beta0010.pvtu',['velocity','constant temperature'])
        # Get bed,surface from inputs to calculate driving stress
        xb,yb,zbgrid = elmerreadlib.input_file(DIR_FS+dir+"/inputs/zbdem.xy")
        xb,yb,zsgrid = elmerreadlib.input_file(DIR_FS+dir+"/inputs/zsdem.xy")
        xgrid,ygrid,ugrid = elmerreadlib.input_file(DIR_FS+dir+"/inputs/udem.xy")
        xgrid,ygrid,vgrid = elmerreadlib.input_file(DIR_FS+dir+"/inputs/vdem.xy")
    modelT_depth_fs = elmerreadlib.depth_averaged_variables(modelT_fs)

    # Mesh boundaries
    extent = np.loadtxt(DIR_FS+dir+"/inputs/mesh_extent.dat")
    try:
      hole1 = np.loadtxt(DIR_FS+dir+"/inputs/mesh_hole1.dat")
      hole2 = np.loadtxt(DIR_FS+dir+"/inputs/mesh_hole2.dat")
      holes=[hole1,hole2]
    except:
      holes = []

  #elif dir.endswith('HighT'):
  #  dirs2 = os.listdir(DIR_FS+dir+'/mesh2d/inversion_'+method)
  #  for dir2 in dirs2:
  #    if dir2.startswith('lambda_'+regpar_fs) and not(dir2.endswith('.pdf')):
  #      highT_bed_fs = elmerreadlib.saveline_boundary(DIR_FS+dir+'/mesh2d/inversion_'+method+'/'+dir2+'/',method+"_beta.dat",bbed,['velocity','beta'])
  #      highT_sur_fs = elmerreadlib.saveline_boundary(DIR_FS+dir+'/mesh2d/inversion_'+method+'/'+dir2+'/',method+"_beta.dat",bsur,['velocity','vsurfini'])
  elif dir.endswith('LowT'):
    dirs2 = os.listdir(DIR_FS+dir+'/mesh2d/inversion_'+method)
    for dir2 in dirs2:
      if dir2.startswith('lambda_'+regpar_fs) and not(dir2.endswith('.pdf')):
        lowT_bed_fs = elmerreadlib.saveline_boundary(DIR_FS+dir+'/mesh2d/inversion_'+method+'/'+dir2+'/',method+"_beta.dat",bbed,['velocity','beta'])
        lowT_sur_fs = elmerreadlib.saveline_boundary(DIR_FS+dir+'/mesh2d/inversion_'+method+'/'+dir2+'/',method+"_beta.dat",bsur,['velocity','vsurfini'])
        lowT_fs = elmerreadlib.pvtu_file(DIR_FS+dir+'/mesh2d/inversion_'+method+'/'+dir2+'/'+method+'_beta0010.pvtu',['velocity','constant temperature'])
    lowT_depth_fs = elmerreadlib.depth_averaged_variables(lowT_fs)

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
ny,nx = np.shape(zsgrid)

xbgrid,ybgrid = np.meshgrid(xb,yb)
xb_flat = xbgrid.flatten()
yb_flat = ybgrid.flatten()

# Make mask
mask = np.ones([ny,nx],dtype='bool')
if len(extent) > 0:
  exterior = Path(np.column_stack((extent[:,0],extent[:,1])))
  pts = exterior.contains_points(np.column_stack((xb_flat,yb_flat)))
  bool = ~(np.reshape(pts,(ny,nx)))
  mask[bool==0] = 0

if len(holes) > 0:
  # Mask out points inside holes
  for i in range(0,len(holes)):
    hole = Path(np.column_stack((holes[i][:,0],holes[i][:,1])))
    pts = hole.contains_points(np.column_stack((xb_flat,yb_flat)))
    bool = (np.reshape(pts,(ny,nx)))
    mask[bool==1] = 1

# Interpolate velocity to zb,zs grid so that we can calculate the along-flow component of the driving stress
fu = scipy.interpolate.RegularGridInterpolator((ygrid,xgrid),ugrid)
fv = scipy.interpolate.RegularGridInterpolator((ygrid,xgrid),vgrid)
u = np.reshape(fu((yb_flat,xb_flat)),(ny,nx))
v = np.reshape(fv((yb_flat,xb_flat)),(ny,nx))
del xb_flat,yb_flat,fu,fv,xbgrid,ybgrid

# Calculate grad(h) & H for driving stress
dhdx = np.zeros([ny,nx])
dhdy = np.zeros([ny,nx])
dhdx[:,:] = float('nan')
dhdy[:,:] = float('nan')

H = zsgrid-zbgrid

dhdx[1:-1,1:-1] = (zsgrid[1:-1,2:]-zsgrid[1:-1,0:-2])/(xb[2:]-xb[0:-2])
dhdy[1:-1,1:-1] = ((zsgrid[2:,1:-1]-zsgrid[0:-2,1:-1]).T/(yb[2:]-yb[0:-2])).T
#tauds = rho_i*g*H*np.sqrt((dhdx)**2+(dhdy)**2)
tauds_flow = -rho_i*g*H*(dhdx*u+dhdy*v)/np.sqrt(u**2+v**2)

# Interpolate tauds_flow onto grid
grid_modelT = elmerreadlib.grid3d(modelT_ssa,'taub',holes,extent)
xm = grid_modelT[0]
ym = grid_modelT[1]
xmesh,ymesh = np.meshgrid(xm,ym)

f = scipy.interpolate.RegularGridInterpolator([yb,xb],tauds_flow,method='linear')
tauds_flow_grid = np.reshape(f((ymesh.flatten(),xmesh.flatten())),(len(ym),len(xm)))

del xmesh, ymesh, grid_modelT

tauds_flow = np.ma.masked_array(tauds_flow,mask)

del dhdx,dhdy,rho_i,g,H,ny,nx

#######################
# Plot driving stress #
#######################

fig = plt.figure(figsize=(2.4,3))
ax = plt.gca()
plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
p = plt.imshow(tauds_flow/1e3,extent=[xb[0],xb[-1],yb[0],yb[-1]],origin='lower',vmin=-500,vmax=500,cmap='RdBu_r')
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
ax.text(xmin+0.51*(xmax-xmin),ymin+0.805*(ymax-ymin),r'$\tau_d \cdot u / ||u||$ (kPa)',fontsize=11)
cb.ax.tick_params(labelsize=11)
ax.plot([xmin+0.63*(xmax-xmin),xmin+0.63*(xmax-xmin)+5e3],[ymin+0.78*(ymax-ymin),ymin+0.78*(ymax-ymin)],'k',linewidth=1.5)
ax.plot([xmin+0.63*(xmax-xmin),xmin+0.63*(xmax-xmin)],[ymin+0.78*(ymax-ymin),ymin+0.76*(ymax-ymin)],'k',linewidth=1.5)
ax.plot([xmin+0.63*(xmax-xmin)+5e3,xmin+0.63*(xmax-xmin)+5e3],[ymin+0.78*(ymax-ymin),ymin+0.76*(ymax-ymin)],'k',linewidth=1.5)
ax.text(xmin+0.67*(xmax-xmin)+5e3,ymin+0.76*(ymax-ymin),'5 km',fontsize=10)
plt.tight_layout()
plt.subplots_adjust(hspace=0.03,wspace=0.03,top=0.99,right=0.99,left=0.01,bottom=0.01)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_inversions_temperature_taud.pdf"),FORMAT='PDF',dpi=600)
plt.close()

#################################################################################
# Plot basal shear stress for different forward models and temperature profiles #
#################################################################################

fig = plt.figure(figsize=(4.8,5.8))
matplotlib.rc('font',family='Arial')
gs = matplotlib.gridspec.GridSpec(2,2)

plt.subplot(gs[0,0])
ax = plt.gca()
grid = elmerreadlib.grid3d(lowT_bed_fs,'taub',holes,extent)
plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
p=plt.imshow(grid[2]*1e3,extent=[grid[0][0],grid[0][-1],grid[1][0],grid[1][-1]],origin='lower',vmin=0,vmax=500)
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
ax.text(xmin+0.68*(xmax-xmin),ymin+0.78*(ymax-ymin),r'$\tau_b$ (kPa)',fontsize=11)
cb.ax.tick_params(labelsize=11)
ax.plot([xmin+0.63*(xmax-xmin),xmin+0.63*(xmax-xmin)+5e3],[ymin+0.75*(ymax-ymin),ymin+0.75*(ymax-ymin)],'k',linewidth=1.5)
ax.plot([xmin+0.63*(xmax-xmin),xmin+0.63*(xmax-xmin)],[ymin+0.75*(ymax-ymin),ymin+0.73*(ymax-ymin)],'k',linewidth=1.5)
ax.plot([xmin+0.63*(xmax-xmin)+5e3,xmin+0.63*(xmax-xmin)+5e3],[ymin+0.75*(ymax-ymin),ymin+0.73*(ymax-ymin)],'k',linewidth=1.5)
ax.text(xmin+0.66*(xmax-xmin)+5e3,ymin+0.725*(ymax-ymin),'5 km',fontsize=10)
ax.text(xmin+0.03*(xmax-xmin),ymin+0.92*(ymax-ymin),'a',fontweight='bold',fontsize=16)

plt.subplot(gs[0,1])
ax = plt.gca()
grid = elmerreadlib.grid3d(modelT_bed_fs,'taub',holes,extent)
plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
p=plt.imshow(grid[2]*1e3,extent=[grid[0][0],grid[0][-1],grid[1][0],grid[1][-1]],origin='lower',vmin=0,vmax=500)
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
ax.text(xmin+0.68*(xmax-xmin),ymin+0.79*(ymax-ymin),r'$\tau_b$ (kPa)',fontsize=11)
cb.ax.tick_params(labelsize=11)
ax.text(xmin+0.03*(xmax-xmin),ymin+0.92*(ymax-ymin),'b',fontweight='bold',fontsize=16)

plt.subplot(gs[1,0])
ax = plt.gca()
grid = elmerreadlib.grid3d(lowT_ssa,'taub',holes,extent)
plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
p=plt.imshow(grid[2]*1e3,extent=[grid[0][0],grid[0][-1],grid[1][0],grid[1][-1]],origin='lower',vmin=0,vmax=500)
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
ax.text(xmin+0.68*(xmax-xmin),ymin+0.79*(ymax-ymin),r'$\tau_b$ (kPa)',fontsize=11)
cb.ax.tick_params(labelsize=11)
ax.text(xmin+0.03*(xmax-xmin),ymin+0.92*(ymax-ymin),'c',fontweight='bold',fontsize=16)

plt.subplot(gs[1,1])
ax = plt.gca()
grid = elmerreadlib.grid3d(modelT_ssa,'taub',holes,extent)
plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
p=plt.imshow(grid[2]*1e3,extent=[grid[0][0],grid[0][-1],grid[1][0],grid[1][-1]],origin='lower',vmin=0,vmax=500)
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
ax.text(xmin+0.68*(xmax-xmin),ymin+0.79*(ymax-ymin),r'$\tau_b$ (kPa)',fontsize=11)
cb.ax.tick_params(labelsize=11)
ax.text(xmin+0.03*(xmax-xmin),ymin+0.92*(ymax-ymin),'d',fontweight='bold',fontsize=16)

plt.tight_layout()
plt.subplots_adjust(hspace=0.03,wspace=0.03,top=0.99,right=0.99,left=0.01,bottom=0.01)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_inversions_temperature_taub.pdf"),FORMAT='PDF',dpi=600)
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
grid_temp = elmerreadlib.grid3d(modelT_depth_fs,'constant temperature',holes,extent)
plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
p=plt.imshow(grid_temp[2],extent=[grid_temp[0][0],grid_temp[0][-1],grid_temp[1][0],grid_temp[1][-1]],origin='lower',vmin=-15,vmax=0)
plt.xticks([])
plt.yticks([])
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
for i in range(0,len(pts[:,0])):
  plt.plot(pts[i,0],pts[i,1],'o',color=colorlabels[i])
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
cbaxes = fig.add_axes([0.29, 0.91, 0.26, 0.025]) 
cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=np.arange(-15,5,5)) 
ax.text(xmin+0.48*(xmax-xmin),ymin+0.81*(ymax-ymin),'Temperature ($^o$C)',fontsize=8)
cb.ax.tick_params(labelsize=8)
ax.plot([xmin+0.63*(xmax-xmin),xmin+0.63*(xmax-xmin)+5e3],[ymin+0.78*(ymax-ymin),ymin+0.78*(ymax-ymin)],'k',linewidth=1.5)
ax.plot([xmin+0.63*(xmax-xmin),xmin+0.63*(xmax-xmin)],[ymin+0.78*(ymax-ymin),ymin+0.76*(ymax-ymin)],'k',linewidth=1.5)
ax.plot([xmin+0.63*(xmax-xmin)+5e3,xmin+0.63*(xmax-xmin)+5e3],[ymin+0.78*(ymax-ymin),ymin+0.76*(ymax-ymin)],'k',linewidth=1.5)
ax.text(xmin+0.67*(xmax-xmin)+5e3,ymin+0.76*(ymax-ymin),'5 km',fontsize=8)
ax.text(xmin+0.03*(xmax-xmin),ymin+0.93*(ymax-ymin),'a',fontweight='bold',fontsize=10)


plt.subplot(gs[2])
ax = plt.gca()
for i in range(0,len(pts[:,0])):
  column = elmerreadlib.values_in_column(modelT_fs,pts[i,0],pts[i,1])
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

plt.tight_layout()
plt.subplots_adjust(hspace=0.02,wspace=0.02,top=0.97,right=0.86,left=0.01,bottom=0.12)
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
grid_vsur = elmerreadlib.grid3d(lowT_sur_fs,'velocity',holes,extent)
grid_vbed = elmerreadlib.grid3d(lowT_bed_fs,'velocity',holes,extent)
plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
p=plt.imshow(grid_vbed[2]/grid_vsur[2],extent=[grid_vsur[0][0],grid_vsur[0][-1],grid_vsur[1][0],grid_vsur[1][-1]],origin='lower',vmin=0,vmax=1)
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
grid_vsur = elmerreadlib.grid3d(modelT_sur_fs,'velocity',holes,extent)
grid_vbed = elmerreadlib.grid3d(modelT_bed_fs,'velocity',holes,extent)
plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
p=plt.imshow(grid_vbed[2]/grid_vsur[2],extent=[grid_vsur[0][0],grid_vsur[0][-1],grid_vsur[1][0],grid_vsur[1][-1]],origin='lower',vmin=0,vmax=1)
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
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_temperature_slidingratio.pdf"),FORMAT='PDF',dpi=600)
plt.close()

###############################################################
# Compare results from different models for different regions #
###############################################################

if glacier == 'Helheim':
  
  # Percentile range for calculating statistics
  q=50.
  
  # Get averages for Full Stokes
  ind_U_x = []
  ind_U_y = []
  ind_L_x = []
  ind_L_y = []
  ind_M_x = []
  ind_M_y = []
  grid_lowT = elmerreadlib.grid3d(lowT_bed_fs,'taub',holes,extent)
  grid_modelT = elmerreadlib.grid3d(modelT_bed_fs,'taub',holes,extent)
  for i in range(0,len(grid_lowT[0])):
    for j in range(0,len(grid_lowT[1])):
      if region_U.contains_point([grid_lowT[0][i],grid_lowT[1][j]]):
        ind_U_x.append(i)
        ind_U_y.append(j)
      if region_L.contains_point([grid_lowT[0][i],grid_lowT[1][j]]):
        ind_L_x.append(i)
        ind_L_y.append(j) 
      if (region_M1.contains_point([grid_lowT[0][i],grid_lowT[1][j]])) or (region_M2.contains_point([grid_lowT[0][i],grid_lowT[1][j]])):
        ind_M_x.append(i)
        ind_M_y.append(j)
  
  region_U_taub_lowT_fs = np.array([np.mean(grid_lowT[2][ind_U_y,ind_U_x]),np.percentile(grid_lowT[2][ind_U_y,ind_U_x],q/2),\
                        np.percentile(grid_lowT[2][ind_U_y,ind_U_x],100-q/2)])*1e3
  region_L_taub_lowT_fs = np.array([np.mean(grid_lowT[2][ind_L_y,ind_L_x]),np.percentile(grid_lowT[2][ind_L_y,ind_L_x],q/2),\
                        np.percentile(grid_lowT[2][ind_L_y,ind_L_x],100-q/2)])*1e3
  region_M_taub_lowT_fs = np.array([np.mean(grid_lowT[2][ind_M_y,ind_M_x]),np.percentile(grid_lowT[2][ind_M_y,ind_M_x],q/2),\
                        np.percentile(grid_lowT[2][ind_M_y,ind_M_x],100-q/2)])*1e3
  
  region_U_taub_modelT_fs = np.array([np.mean(grid_modelT[2][ind_U_y,ind_U_x]),np.percentile(grid_modelT[2][ind_U_y,ind_U_x],q/2),\
                          np.percentile(grid_modelT[2][ind_U_y,ind_U_x],100-q/2)])*1e3
  region_L_taub_modelT_fs = np.array([np.mean(grid_modelT[2][ind_L_y,ind_L_x]),np.percentile(grid_modelT[2][ind_L_y,ind_L_x],q/2),\
                          np.percentile(grid_modelT[2][ind_L_y,ind_L_x],100-q/2)])*1e3
  region_M_taub_modelT_fs = np.array([np.mean(grid_modelT[2][ind_M_y,ind_M_x]),np.percentile(grid_modelT[2][ind_M_y,ind_M_x],q/2),\
                          np.percentile(grid_modelT[2][ind_M_y,ind_M_x],100-q/2)])*1e3
  
  # Get averages for SSA
  ind_U_x = []
  ind_U_y = []
  ind_L_x = []
  ind_L_y = []
  ind_M_x = []
  ind_M_y = []
  grid_lowT = elmerreadlib.grid3d(lowT_ssa,'taub',holes,extent)
  grid_modelT = elmerreadlib.grid3d(modelT_ssa,'taub',holes,extent)
  for i in range(0,len(grid_lowT[0])):
    for j in range(0,len(grid_lowT[1])):
      if region_U.contains_point([grid_lowT[0][i],grid_lowT[1][j]]):
        ind_U_x.append(i)
        ind_U_y.append(j)
      if region_L.contains_point([grid_lowT[0][i],grid_lowT[1][j]]):
        ind_L_x.append(i)
        ind_L_y.append(j)
      if (region_M1.contains_point([grid_lowT[0][i],grid_lowT[1][j]])) or (region_M2.contains_point([grid_lowT[0][i],grid_lowT[1][j]])):
        ind_M_x.append(i)
        ind_M_y.append(j)
  
  region_U_taub_lowT_ssa = np.array([np.mean(grid_lowT[2][ind_U_y,ind_U_x]),np.percentile(grid_lowT[2][ind_U_y,ind_U_x],q/2),\
                         np.percentile(grid_lowT[2][ind_U_y,ind_U_x],100-q/2)])*1e3
  region_L_taub_lowT_ssa = np.array([np.mean(grid_lowT[2][ind_L_y,ind_L_x]),np.percentile(grid_lowT[2][ind_L_y,ind_L_x],q/2),\
                         np.percentile(grid_lowT[2][ind_L_y,ind_L_x],100-q/2)])*1e3
  region_M_taub_lowT_ssa = np.array([np.mean(grid_lowT[2][ind_M_y,ind_M_x]),np.percentile(grid_lowT[2][ind_M_y,ind_M_x],q/2),\
                         np.percentile(grid_lowT[2][ind_M_y,ind_M_x],100-q/2)])*1e3
  
  region_U_taub_modelT_ssa = np.array([np.mean(grid_modelT[2][ind_U_y,ind_U_x]),np.percentile(grid_modelT[2][ind_U_y,ind_U_x],q/2),\
                           np.percentile(grid_modelT[2][ind_U_y,ind_U_x],100-q/2)])*1e3
  region_L_taub_modelT_ssa = np.array([np.mean(grid_modelT[2][ind_L_y,ind_L_x]),np.percentile(grid_modelT[2][ind_L_y,ind_L_x],q/2),\
                           np.percentile(grid_modelT[2][ind_L_y,ind_L_x],100-q/2)])*1e3
  region_M_taub_modelT_ssa = np.array([np.mean(grid_modelT[2][ind_M_y,ind_M_x]),np.percentile(grid_modelT[2][ind_M_y,ind_M_x],q/2),\
                           np.percentile(grid_modelT[2][ind_M_y,ind_M_x],100-q/2)])*1e3
  
  # Get average driving stresses
  ind_U_x = []
  ind_U_y = []
  ind_L_x = []
  ind_L_y = []
  ind_M_x = []
  ind_M_y = []
  for i in range(0,len(xb)):
    for j in range(0,len(yb)):
      if region_U.contains_point([xb[i],yb[j]]):
        ind_U_x.append(i)
        ind_U_y.append(j)
      if region_L.contains_point([xb[i],yb[j]]):
        ind_L_x.append(i)
        ind_L_y.append(j)
      if (region_M1.contains_point([xb[i],yb[j]])) or (region_M2.contains_point([xb[i],yb[j]])):
        ind_M_x.append(i)
        ind_M_y.append(j)
  
  region_U_taud = np.mean(tauds_flow[ind_U_y,ind_U_x])/1e3
  region_L_taud = np.mean(tauds_flow[ind_L_y,ind_L_x])/1e3
  region_M_taud = np.mean(tauds_flow[ind_M_y,ind_M_x])/1e3
  
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
  patch = matplotlib.patches.PathPatch(region_M1,edgecolor='k',facecolor='k',lw=2,alpha=0.5)
  ax.add_patch(patch)
  patch = matplotlib.patches.PathPatch(region_M2,edgecolor='k',facecolor='k',lw=2,alpha=0.5)
  ax.add_patch(patch)
  patch = matplotlib.patches.PathPatch(region_U,edgecolor='k',facecolor='r',lw=2,alpha=0.5)
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
  #ax.text(xmin+0.51*(xmax-xmin),ymin+0.805*(ymax-ymin),r'$\tau_d \cdot u / ||u||$ (kPa)',fontsize=11)
  #cb.ax.tick_params(labelsize=11)
  ax.plot([xmin+0.63*(xmax-xmin),xmin+0.63*(xmax-xmin)+5e3],[ymin+0.94*(ymax-ymin),ymin+0.94*(ymax-ymin)],'k',linewidth=1.5)
  ax.plot([xmin+0.63*(xmax-xmin),xmin+0.63*(xmax-xmin)],[ymin+0.94*(ymax-ymin),ymin+0.92*(ymax-ymin)],'k',linewidth=1.5)
  ax.plot([xmin+0.63*(xmax-xmin)+5e3,xmin+0.63*(xmax-xmin)+5e3],[ymin+0.94*(ymax-ymin),ymin+0.92*(ymax-ymin)],'k',linewidth=1.5)
  ax.text(xmin+0.65*(xmax-xmin)+5e3,ymin+0.92*(ymax-ymin),'5 km',fontsize=10)
  
  ax.text(xmin+0.63*(xmax-xmin),ymin+0.2*(ymax-ymin),'L',fontweight='bold',fontsize=11,backgroundcolor='w')
  ax.text(xmin+0.18*(xmax-xmin),ymin+0.4*(ymax-ymin),'M',fontweight='bold',fontsize=11,backgroundcolor='w')
  ax.text(xmin+0.35*(xmax-xmin),ymin+0.8*(ymax-ymin),'M',fontweight='bold',fontsize=11,backgroundcolor='w')
  ax.text(xmin+0.23*(xmax-xmin),ymin+0.6*(ymax-ymin),'U',fontweight='bold',fontsize=11,backgroundcolor='w')
  
  plt.tight_layout()
  plt.subplots_adjust(hspace=0.03,wspace=0.03,top=0.99,right=0.99,left=0.01,bottom=0.01)
  plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_inversions_temperature_regions.pdf"),FORMAT='PDF',dpi=600)
  plt.close()
  
  
  # Plot bar graph
  fig = plt.figure(figsize=(3,3))
  ax = plt.gca()
  width=0.15
  
  ind = 0.7+np.arange(3)
  plt.plot([0.5,4],[1,1],'k--',lw=2)
  ax.bar(ind[0],region_L_taub_lowT_fs[0]/region_L_taud,width,color='r',edgecolor='k',\
       yerr=[[region_L_taub_lowT_fs[0]/region_L_taud-region_L_taub_lowT_fs[1]/region_L_taud],\
       [region_L_taub_lowT_fs[2]/region_L_taud-region_L_taub_lowT_fs[0]/region_L_taud]],error_kw=dict(ecolor='k',capsize=2),label='a: FS,CT')
  ax.bar(ind[1],region_U_taub_lowT_fs[0]/region_U_taud,width,color='r',\
       yerr=[[region_U_taub_lowT_fs[0]/region_U_taud-region_U_taub_lowT_fs[1]/region_U_taud],\
       [region_U_taub_lowT_fs[2]/region_U_taud-region_U_taub_lowT_fs[0]/region_U_taud]],error_kw=dict(ecolor='k',capsize=2))
  ax.bar(ind[2],region_M_taub_lowT_fs[0]/region_M_taud,width,color='r',\
       yerr=[[region_M_taub_lowT_fs[0]/region_M_taud-region_M_taub_lowT_fs[1]/region_M_taud],\
       [region_M_taub_lowT_fs[2]/region_M_taud-region_M_taub_lowT_fs[0]/region_M_taud]],error_kw=dict(ecolor='k',capsize=2))
  ax.bar(ind[0]+width,region_L_taub_modelT_fs[0]/region_L_taud,width,color=(0.9,0.65,0.65),\
       yerr=[[region_L_taub_modelT_fs[0]/region_L_taud-region_L_taub_modelT_fs[1]/region_L_taud],\
       [region_L_taub_modelT_fs[2]/region_L_taud-region_L_taub_modelT_fs[0]/region_L_taud]],error_kw=dict(ecolor='k',capsize=2),label='b: FS,VT')
  ax.bar(ind[1]+width,region_U_taub_modelT_fs[0]/region_U_taud,width,color=(0.9,0.7,0.7),\
       yerr=[[region_U_taub_modelT_fs[0]/region_U_taud-region_U_taub_modelT_fs[1]/region_U_taud],\
       [region_U_taub_modelT_fs[2]/region_U_taud-region_U_taub_modelT_fs[0]/region_U_taud]],error_kw=dict(ecolor='k',capsize=2))
  ax.bar(ind[2]+width,region_M_taub_modelT_fs[0]/region_M_taud,width,color=(0.9,0.7,0.7),\
       yerr=[[region_M_taub_modelT_fs[0]/region_M_taud-region_M_taub_modelT_fs[1]/region_M_taud],\
       [region_M_taub_modelT_fs[2]/region_M_taud-region_M_taub_modelT_fs[0]/region_M_taud]],error_kw=dict(ecolor='k',capsize=2))
  ax.bar(ind[0]+2*width,region_L_taub_lowT_ssa[0]/region_L_taud,width,color='b',\
       yerr=[[region_L_taub_lowT_ssa[0]/region_L_taud-region_L_taub_lowT_ssa[1]/region_L_taud],\
       [region_L_taub_lowT_ssa[2]/region_L_taud-region_L_taub_lowT_ssa[0]/region_L_taud]],error_kw=dict(ecolor='k',capsize=2),label='c: SSA,CT')
  ax.bar(ind[1]+2*width,region_U_taub_lowT_ssa[0]/region_U_taud,width,color='b',\
       yerr=[[region_U_taub_lowT_ssa[0]/region_U_taud-region_U_taub_lowT_ssa[1]/region_U_taud],\
       [region_U_taub_lowT_ssa[2]/region_U_taud-region_U_taub_lowT_ssa[0]/region_U_taud]],error_kw=dict(ecolor='k',capsize=2))
  ax.bar(ind[2]+2*width,region_M_taub_lowT_ssa[0]/region_M_taud,width,color='b',\
       yerr=[[region_M_taub_lowT_ssa[0]/region_M_taud-region_M_taub_lowT_ssa[1]/region_M_taud],\
       [region_M_taub_lowT_ssa[2]/region_M_taud-region_M_taub_lowT_ssa[0]/region_M_taud]],error_kw=dict(ecolor='k',capsize=2))
  ax.bar(ind[0]+3*width,region_L_taub_modelT_ssa[0]/region_L_taud,width,color=(0.7,0.7,0.9),\
       yerr=[[region_L_taub_modelT_ssa[0]/region_L_taud-region_L_taub_modelT_ssa[1]/region_L_taud],\
       [region_L_taub_modelT_ssa[2]/region_L_taud-region_L_taub_modelT_ssa[0]/region_L_taud]],error_kw=dict(ecolor='k',capsize=2),label='d: SSA,VT')
  ax.bar(ind[1]+3*width,region_U_taub_modelT_ssa[0]/region_U_taud,width,color=(0.7,0.7,0.9),\
       yerr=[[region_U_taub_modelT_ssa[0]/region_U_taud-region_U_taub_modelT_ssa[1]/region_U_taud],\
       [region_U_taub_modelT_ssa[2]/region_U_taud-region_U_taub_modelT_ssa[0]/region_U_taud]],error_kw=dict(ecolor='k',capsize=2))
  ax.bar(ind[2]+3*width,region_M_taub_modelT_ssa[0]/region_M_taud,width,color=(0.7,0.7,0.9),\
       yerr=[[region_M_taub_modelT_ssa[0]/region_M_taud-region_M_taub_modelT_ssa[1]/region_M_taud],\
       [region_M_taub_modelT_ssa[2]/region_M_taud-region_M_taub_modelT_ssa[0]/region_M_taud]],error_kw=dict(ecolor='k',capsize=2))
  
  plt.xlim([0.5,3.5])
  plt.yticks(np.arange(0.0,1.3,0.2),fontsize=12,fontname='Arial')
  plt.ylabel(r'$\tau_b/\tau_d$',fontname='Arial',fontsize=16)
  plt.xticks([1,2,3],fontsize=12,fontname='Arial')
  ax.set_xticklabels(('L','U','M'))
  plt.legend(loc=2,fontsize=10,numpoints=1,ncol=1,handlelength=0.4,handletextpad=0.2,labelspacing=0.2,columnspacing=0.2)
  
  plt.subplots_adjust(hspace=0.02,wspace=0.02,top=0.96,right=0.98,left=0.2,bottom=0.09)
  plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_inversions_temperature_stress.pdf"),FORMAT='PDF',dpi=600)
  plt.close()



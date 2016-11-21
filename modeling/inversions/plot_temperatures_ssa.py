import sys
import os
import matplotlib.pyplot as plt
from pylab import *
import argparse
import scipy
import elmerreadlib,glaclib,colorlib,geotifflib
import cubehelix
from matplotlib.path import Path
import numpy as np
import matplotlib
impot matplotlib.pyplot as plt

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
regpar = '1e10'

# Directory
DIR = os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/Temperatures_SSA/")

#############
# Load data #
#############

# Get directories
dirs = os.listdir(DIR)
for dir in dirs:
  if dir.endswith('ModelT'):
    dirs2 = os.listdir(DIR+dir+'/mesh2d/inversion_'+method)
    for dir2 in dirs2:
      if dir2.startswith('lambda_'+regpar) and not(dir2.endswith('.pdf')):
        modelT = elmerreadlib.result_file(DIR+dir+'/mesh2d/','inversion_'+method+'/'+dir2+'/adjoint_beta_ssa.result',['ssavelocity 1','ssavelocity 2','beta','vsurfini 1','vsurfini 2'])
        # Get bed,surface from inputs to calculate driving stress
        xb,yb,zbgrid = elmerreadlib.input_file(DIR+dir+"/inputs/zbdem.xy")
        xb,yb,zsgrid = elmerreadlib.input_file(DIR+dir+"/inputs/zsdem.xy")
        xgrid,ygrid,ugrid = elmerreadlib.input_file(DIR+dir+"/inputs/udem.xy")
        xgrid,ygrid,vgrid = elmerreadlib.input_file(DIR+dir+"/inputs/vdem.xy")

    # Mesh boundaries
    extent = np.loadtxt(DIR+dir+"/inputs/mesh_extent.dat")
    try:
      hole1 = np.loadtxt(DIR+dir+"/inputs/mesh_hole1.dat")
      hole2 = np.loadtxt(DIR+dir+"/inputs/mesh_hole2.dat")
      holes=[hole1,hole2]
    except:
      holes = []

  #elif dir.endswith('HighT'):
  #  dirs2 = os.listdir(DIR+dir+'/mesh2d/inversion_'+method)
  #  for dir2 in dirs2:
  #    if dir2.startswith('lambda_'+regpar) and not(dir2.endswith('.pdf')):
  #      highT_bed = elmerreadlib.saveline_boundary(DIR+dir+'/mesh2d/inversion_'+method+'/'+dir2+'/',method+"_beta.dat",bbed,['velocity','beta'])
  #      highT_sur = elmerreadlib.saveline_boundary(DIR+dir+'/mesh2d/inversion_'+method+'/'+dir2+'/',method+"_beta.dat",bsur,['velocity','vsurfini'])
  elif dir.endswith('LowT'):
    dirs2 = os.listdir(DIR+dir+'/mesh2d/inversion_'+method)
    for dir2 in dirs2:
      if dir2.startswith('lambda_'+regpar) and not(dir2.endswith('.pdf')):
       lowT = elmerreadlib.result_file(DIR+dir+'/mesh2d/','inversion_'+method+'/'+dir2+'/adjoint_beta_ssa.result',['ssavelocity 1','ssavelocity 2','beta','vsurfini 1','vsurfini 2'])

# Plot extent
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

tauds_flow = np.ma.masked_array(tauds_flow,mask)

del dhdx,dhdy,rho_i,g,H,ny,nx

######################################################
# Make plot of inversions for different temperatures #
######################################################

fig = plt.figure(figsize=(7.2,3))
matplotlib.rc('font',family='Arial')
gs = matplotlib.gridspec.GridSpec(1,3)

plt.subplot(gs[0])
ax = plt.gca()
plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
p = plt.imshow(tauds_flow/1e3,extent=[xb[0],xb[-1],yb[0],yb[-1]],origin='lower',vmin=-500,vmax=500,cmap='RdBu_r')
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
plt.xticks([])
plt.yticks([])
xmin,xmax = plt.xlim()
ymin,ymax = plt.ylim()
path = matplotlib.path.Path([[0.47*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.98*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.98*(xmax-xmin)+xmin,0.77*(ymax-ymin)+ymin],
                        [0.47*(xmax-xmin)+xmin,0.77*(ymax-ymin)+ymin],
                        [0.47*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin]])
patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1)
ax.add_patch(patch)
cbaxes = fig.add_axes([0.19, 0.91, 0.11, 0.025])
cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=np.arange(-500,600,500))
ax.text(xmin+0.56*(xmax-xmin),ymin+0.83*(ymax-ymin),r'$\tau_d \cdot u / ||u||$ (kPa)',fontsize=8)
cb.ax.tick_params(labelsize=8)
ax.plot([xmin+0.63*(xmax-xmin),xmin+0.63*(xmax-xmin)+5e3],[ymin+0.81*(ymax-ymin),ymin+0.81*(ymax-ymin)],'k',linewidth=1.5)
ax.plot([xmin+0.63*(xmax-xmin),xmin+0.63*(xmax-xmin)],[ymin+0.81*(ymax-ymin),ymin+0.79*(ymax-ymin)],'k',linewidth=1.5)
ax.plot([xmin+0.63*(xmax-xmin)+5e3,xmin+0.63*(xmax-xmin)+5e3],[ymin+0.81*(ymax-ymin),ymin+0.79*(ymax-ymin)],'k',linewidth=1.5)
ax.text(xmin+0.67*(xmax-xmin)+5e3,ymin+0.79*(ymax-ymin),'5 km',fontsize=8)
ax.text(xmin+0.03*(xmax-xmin),ymin+0.93*(ymax-ymin),'a',fontweight='bold',fontsize=10)

plt.subplot(gs[1])
ax = plt.gca()
grid = elmerreadlib.grid3d(lowT,'taub',holes,extent)
plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
p=plt.imshow(grid[2]*1e3,extent=[grid[0][0],grid[0][-1],grid[1][0],grid[1][-1]],origin='lower',vmin=0,vmax=500)
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
plt.xticks([])
plt.yticks([])
path = matplotlib.path.Path([[0.57*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.98*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.98*(xmax-xmin)+xmin,0.81*(ymax-ymin)+ymin],
                        [0.57*(xmax-xmin)+xmin,0.81*(ymax-ymin)+ymin],
                        [0.57*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin]])
patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1)
ax.add_patch(patch)
cbaxes = fig.add_axes([0.535, 0.91, 0.11, 0.025])
cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=np.arange(0,600,200))
ax.text(xmin+0.68*(xmax-xmin),ymin+0.83*(ymax-ymin),r'$\tau_b$ (kPa)',fontsize=8)
cb.ax.tick_params(labelsize=8)
ax.text(xmin+0.03*(xmax-xmin),ymin+0.93*(ymax-ymin),'b',fontweight='bold',fontsize=10)

plt.subplot(gs[2])
ax = plt.gca()
grid = elmerreadlib.grid3d(modelT,'taub',holes,extent)
plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
p=plt.imshow(grid[2]*1e3,extent=[grid[0][0],grid[0][-1],grid[1][0],grid[1][-1]],origin='lower',vmin=0,vmax=500)
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
plt.xticks([])
plt.yticks([])
path = matplotlib.path.Path([[0.57*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.98*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.98*(xmax-xmin)+xmin,0.81*(ymax-ymin)+ymin],
                        [0.57*(xmax-xmin)+xmin,0.81*(ymax-ymin)+ymin],
                        [0.57*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin]])
patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1)
ax.add_patch(patch)
cbaxes = fig.add_axes([0.865, 0.91, 0.11, 0.025])
cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=np.arange(0,600,200))
ax.text(xmin+0.68*(xmax-xmin),ymin+0.83*(ymax-ymin),r'$\tau_b$ (kPa)',fontsize=8)
cb.ax.tick_params(labelsize=8)
ax.text(xmin+0.03*(xmax-xmin),ymin+0.93*(ymax-ymin),'c',fontweight='bold',fontsize=10)

plt.tight_layout()
plt.subplots_adjust(hspace=0.02,wspace=0.02,top=0.99,right=0.99,left=0.01,bottom=0.01)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_inversions_ssa_temperature.pdf"),FORMAT='PDF',dpi=600)
plt.close()


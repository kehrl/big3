import sys
import os
import matplotlib.pyplot as plt
from pylab import *
import argparse
import scipy
import datelib, elmerreadlib, glaclib, geotifflib, meshlib
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
parser.add_argument("-mesh",dest="mesh",required = True,
                help = "Mesh directory name.")

args, _ = parser.parse_known_args(sys.argv)
glacier = args.glacier
mesh = args.mesh

# Directory
DIR = os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/"+mesh)

###############################
# Load shapefiles for regions #
###############################

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

elif glacier == 'Kanger':
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

# Plot options
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

######################
# Get input datasets #
######################

x,y,umea = elmerreadlib.input_file(DIR+'/inputs/udem.xy')
x,y,vmea = elmerreadlib.input_file(DIR+'/inputs/vdem.xy')
x,y,zs = elmerreadlib.input_file(DIR+'/inputs/zsdem.xy')
x,y,zb = elmerreadlib.input_file(DIR+'/inputs/zbdem.xy')

nx = len(x)
ny = len(y)

############################
# Calculate driving stress #
############################

# Parameters
rho_i = 917.
g = 9.8

# Calculate grad(h) & H for driving stress
dhdx = np.zeros([ny,nx])
dhdy = np.zeros([ny,nx])
dhdx[:,:] = float('nan')
dhdy[:,:] = float('nan')

H = zs-zb

dhdx[1:-1,1:-1] = (zs[1:-1,2:]-zs[1:-1,0:-2])/(x[2:]-x[0:-2])
dhdy[1:-1,1:-1] = ((zs[2:,1:-1]-zs[0:-2,1:-1]).T/(y[2:]-y[0:-2])).T
tauds_flow = -rho_i*g*H*(dhdx*umea+dhdy*vmea)/np.sqrt(umea**2+vmea**2)

del dhdx,dhdy,rho_i,g,H

#################################
# Get regularization parameters #
#################################

regpars = []
fid = open(DIR+'/mesh2d/inversion_adjoint/summary.dat','r')
lines = fid.readlines()
for line in lines[2:]:
    p = line.split()
    regpars.append(p[0])
fid.close()

regpars = ['1e12','2e12','3e12','4e12','5e12','6e12','7e12','8e12','9e12','1e13']
nt = len(regpars)

###############################################################
# Compare results from different models for different regions #
###############################################################

# Get one of the models to get node coordinates
dirs = os.listdir(DIR+'/mesh2d/inversion_adjoint/')
for dir in dirs:
    if dir.startswith('lambda_'+regpars[0]) and not(dir.endswith('.pdf')):
        try:
            vtu = elmerreadlib.pvtu_file(DIR+'/mesh2d/inversion_adjoint/'+dir+'/adjoint_beta_ssa0001.pvtu',\
                    ['ssavelocity','beta','vsurfini'])
        except:
            vtu = elmerreadlib.pvtu_file(DIR+'/mesh2d/inversion_adjoint/'+dir+'/adjoint_beta0001.pvtu',\
                    ['velocity','beta','vsurfini'])
bed = elmerreadlib.values_in_layer(vtu,'bed')

# Percentile range for calculating statistics
q = 50.
  
# Get indices for averages
ind_U = []
ind_M = []
ind_L = []
ind_S = []
for i in range(0,len(bed['x'])):
    if region_U.contains_point([bed['x'][i],bed['y'][i]]):
        ind_U.append(i)
    if region_M.contains_point([bed['x'][i],bed['y'][i]]):
        ind_M.append(i)
    if region_L.contains_point([bed['x'][i],bed['y'][i]]):
        ind_L.append(i)
    if glacier == 'Helheim':
        if (region_S1.contains_point([bed['x'][i],bed['y'][i]])) or (region_S2.contains_point([bed['x'][i],bed['y'][i]])):
            ind_S.append(i)
    elif glacier == 'Kanger':
        if (region_S.contains_point([bed['x'][i],bed['y'][i]])):
            ind_S.append(i)

# Set up output variables
region_U_taub = np.zeros([nt,3])
region_M_taub = np.zeros([nt,3])
region_L_taub = np.zeros([nt,3])
region_S_taub = np.zeros([nt,3])

region_U_beta = np.zeros([nt,3])
region_M_beta = np.zeros([nt,3])    
region_L_beta = np.zeros([nt,3])
region_S_beta = np.zeros([nt,3])
    
region_U_taud = np.zeros([nt,1])
region_M_taud = np.zeros([nt,1])
region_L_taud = np.zeros([nt,1])
region_S_taud = np.zeros([nt,1])

for i in range(0,nt):     
    regpar = regpars[i]
    del bed
    for dir in dirs:
        if dir.startswith('lambda_'+regpar) and not(dir.endswith('.pdf')):
            try:
                vtu = elmerreadlib.pvtu_file(DIR+'/mesh2d/inversion_adjoint/'+dir+'/adjoint_beta_ssa0001.pvtu',\
                    ['ssavelocity','beta','vsurfini'])
            except:
                vtu = elmerreadlib.pvtu_file(DIR+'/mesh2d/inversion_adjoint/'+dir+'/adjoint_beta0001.pvtu',\
                    ['velocity','beta','vsurfini'])
    bed = elmerreadlib.values_in_layer(vtu,'bed')

    region_U_taub[i,:] = np.array([np.mean(bed['taub'][ind_U]),\
            np.percentile(bed['taub'][ind_U],q/2),\
            np.percentile(bed['taub'][ind_U],100-q/2)])*1e3
    region_U_beta[i,:] = np.array([np.mean(bed['beta'][ind_U]),\
            np.percentile(bed['beta'][ind_U],q/2),\
            np.percentile(bed['beta'][ind_U],100-q/2)])
    region_M_taub[i,:] = np.array([np.mean(bed['taub'][ind_M]),\
            np.percentile(bed['taub'][ind_M],q/2),\
            np.percentile(bed['taub'][ind_M],100-q/2)])*1e3
    region_M_beta[i,:] = np.array([np.mean(bed['beta'][ind_M]),\
            np.percentile(bed['beta'][ind_M],q/2),\
            np.percentile(bed['beta'][ind_M],100-q/2)])
    region_L_taub[i,:] = np.array([np.mean(bed['taub'][ind_L]),\
            np.percentile(bed['taub'][ind_L],q/2),\
            np.percentile(bed['taub'][ind_L],100-q/2)])*1e3
    region_L_beta[i,:] = np.array([np.mean(bed['beta'][ind_L]),\
            np.percentile(bed['beta'][ind_L],q/2),\
            np.percentile(bed['beta'][ind_L],100-q/2)])
    region_S_taub[i,:] = np.array([np.mean(bed['taub'][ind_S]),\
            np.percentile(bed['taub'][ind_S],q/2),\
            np.percentile(bed['taub'][ind_S],100-q/2)])*1e3
    region_S_beta[i,:] = np.array([np.mean(bed['beta'][ind_S]),\
            np.percentile(bed['beta'][ind_S],q/2),\
            np.percentile(bed['beta'][ind_S],100-q/2)])

    #region_U_taud[i] = np.mean(tauds_flow[ind_U_y,ind_U_x,i])/1e3
    #region_M_taud[i] = np.mean(tauds_flow[ind_M_y,ind_M_x,i])/1e3  
    #region_L_taud[i] = np.mean(tauds_flow[ind_L_y,ind_L_x,i])/1e3
    #region_S_taud[i] = np.mean(tauds_flow[ind_S_y,ind_S_x,i])/1e3
  
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
#plt.plot(extent[:,0],extent[:,1],'k',lw=2)
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
#ax.text(xmin+0.51*(xmax-xmin),ymin+0.805*(ymax-ymin),r'$\tau_d \cdot u / ||u||$ (kPa)',fontsize=10)
#cb.ax.tick_params(labelsize=10)
ax.plot([xmin+0.63*(xmax-xmin),xmin+0.63*(xmax-xmin)+5e3],[ymin+0.94*(ymax-ymin),ymin+0.94*(ymax-ymin)],'k',linewidth=1.5)
ax.plot([xmin+0.63*(xmax-xmin),xmin+0.63*(xmax-xmin)],[ymin+0.94*(ymax-ymin),ymin+0.92*(ymax-ymin)],'k',linewidth=1.5)
ax.plot([xmin+0.63*(xmax-xmin)+5e3,xmin+0.63*(xmax-xmin)+5e3],[ymin+0.94*(ymax-ymin),ymin+0.92*(ymax-ymin)],'k',linewidth=1.5)
ax.text(xmin+0.65*(xmax-xmin)+5e3,ymin+0.92*(ymax-ymin),'5 km',fontsize=10)
  
if glacier == 'Helheim': 
    ax.text(xmin+0.63*(xmax-xmin),ymin+0.2*(ymax-ymin),'L',fontweight='bold',fontsize=10,\
            bbox=dict(boxstyle='square,pad=0.1', fc='white', ec='none'))
    ax.text(xmin+0.18*(xmax-xmin),ymin+0.4*(ymax-ymin),'S',fontweight='bold',fontsize=10,\
            bbox=dict(boxstyle='square,pad=0.1', fc='white', ec='none'))
    ax.text(xmin+0.35*(xmax-xmin),ymin+0.8*(ymax-ymin),'S',fontweight='bold',fontsize=10,\
            bbox=dict(boxstyle='square,pad=0.1', fc='white', ec='none'))
    ax.text(xmin+0.23*(xmax-xmin),ymin+0.6*(ymax-ymin),'U',fontweight='bold',fontsize=10,\
            bbox=dict(boxstyle='square,pad=0.1', fc='white', ec='none'))
    ax.text(xmin+0.5*(xmax-xmin),ymin+0.33*(ymax-ymin),'M',fontweight='bold',fontsize=10,\
            bbox=dict(boxstyle='square,pad=0.1', fc='white', ec='none')) 
elif glacier == 'Kanger':
    ax.text(xmin+0.51*(xmax-xmin),ymin+0.27*(ymax-ymin),'L',fontweight='bold',fontsize=10,\
            bbox=dict(boxstyle='square,pad=0.1', fc='white', ec='none'))
    ax.text(xmin+0.17*(xmax-xmin),ymin+0.4*(ymax-ymin),'S',fontweight='bold',fontsize=10,\
            bbox=dict(boxstyle='square,pad=0.1', fc='white', ec='none'))
    ax.text(xmin+0.3*(xmax-xmin),ymin+0.55*(ymax-ymin),'U',fontweight='bold',fontsize=10,\
            bbox=dict(boxstyle='square,pad=0.1', fc='white', ec='none'))
    ax.text(xmin+0.39*(xmax-xmin),ymin+0.38*(ymax-ymin),'M',fontweight='bold',fontsize=10,\
            bbox=dict(boxstyle='square,pad=0.1', fc='white', ec='none'))

plt.tight_layout()
plt.subplots_adjust(hspace=0.03,wspace=0.03,top=0.99,right=0.99,left=0.01,bottom=0.01)
plt.savefig(os.path.join(os.getenv("HOME"),DIR+"/figures/Regions_map.pdf"),FORMAT='PDF',dpi=600)
plt.close()
  
# Plot bar graphs for taub
fig = plt.figure(figsize=(4,6))
gs = matplotlib.gridspec.GridSpec(4,1)
width=0.15

ax = plt.subplot(gs[0])
r = (np.max(region_L_taub[:,0])+np.min(region_L_taub[:,0]))/2
plt.bar(np.arange(0,len(regpars)),r-region_L_taub[:,0],bottom=r,color='r',edgecolor='k')#,\
plt.ylim([r-4,r+4])
xmin,xmax = ax.get_xlim()
plt.plot([xmin,xmax],[r,r],'k',zorder=0)
plt.xlim([xmin,xmax])
ax.set_xticklabels([])
plt.text(-0.5,r+2.9,'L',fontweight='bold')

ax = plt.subplot(gs[1])
r = (np.max(region_M_taub[:,0])+np.min(region_M_taub[:,0]))/2
plt.bar(np.arange(0,len(regpars)),r-region_M_taub[:,0],bottom=r,color='r',edgecolor='k')
plt.ylim([r-4,r+4])
plt.plot([xmin,xmax],[r,r],'k',zorder=0)
plt.xlim([xmin,xmax])
ax.set_xticklabels([])
plt.text(-0.5,r+2.9,'M',fontweight='bold')

ax = plt.subplot(gs[2])
r = (np.max(region_U_taub[:,0])+np.min(region_U_taub[:,0]))/2
plt.bar(np.arange(0,len(regpars)),r-region_U_taub[:,0],bottom=r,color='r',edgecolor='k')
plt.ylim([r-4,r+4])
plt.plot([xmin,xmax],[r,r],'k',zorder=0)
plt.xlim([xmin,xmax])
ax.set_xticklabels([])
plt.ylabel(r'            $\tau_b$ (kPa)',fontname='arial')
plt.text(-0.5,r+2.9,'U',fontweight='bold')

ax = plt.subplot(gs[3])
r = (np.max(region_S_taub[:,0])+np.min(region_S_taub[:,0]))/2
plt.bar(np.arange(0,len(regpars)),r-region_S_taub[:,0],bottom=r,color='r',edgecolor='k')
plt.ylim([r-4,r+4])
plt.plot([xmin,xmax],[r,r],'k',zorder=0)
plt.xlim([xmin,xmax])
ax.set_xticks(np.arange(0,len(regpars)))
ax.set_xticklabels(regpars,rotation=45)
plt.text(-0.5,r+2.9,'S',fontweight='bold')
plt.xlabel(r'$\lambda$',fontname='arial')
  
plt.subplots_adjust(hspace=0.03,wspace=0.03,top=0.96,right=0.98,left=0.2,bottom=0.1)
plt.savefig(os.path.join(os.getenv("HOME"),DIR+"/figures/Regions_taub_regpar.pdf"),FORMAT='PDF',dpi=600)
plt.close()

import sys
import os
import matplotlib.pyplot as plt
from pylab import *
import argparse
import masklib, meshlib, elmerreadlib, glaclib, colorlib, datelib, icefrontlib
import scipy
from matplotlib.path import Path
from matplotlib.ticker import AutoMinorLocator

##########
# Inputs #
##########

args = sys.argv

# Get inputs to file
parser = argparse.ArgumentParser()
parser.add_argument("-dir", dest="dir", required = True,
        help = "Directory where all the model results are located.")
parser.add_argument("-glacier", dest="glacier",required = True,
        help = "Glacier name.")
parser.add_argument("-method",dest="method",required=False,default='adjoint', 
        help = "Inverse method (adjoint or robin).")
parser.add_argument("-regpar",dest="regpar",required=False,default='5e11',
        help = "Regularization parameter.")

args, _ = parser.parse_known_args(sys.argv)
DIR = args.dir
glacier = args.glacier
method = args.method
regpar = args.regpar

# Surface boundary numbers for saveline
bbed = 4 
bsur = 5

# Parameters
g = 9.8
rho_i = 917.01

# Times
time1 = 2011.
time2 = 2015.

# Input directory
DIR = os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/"+DIR+"/")

# Load some data for plots
x,y,zb,dists = glaclib.load_flowline(glacier)
terminus_val,terminus_time = icefrontlib.distance_along_flowline(x,y,dists,glacier)

################
# Load results #
################

# Get inversion directories and sort them by date
DEMDIRs = os.listdir(DIR)
argsort = np.argsort(DEMDIRs)

# Load taub, beta for each inversion
nt = len(DEMDIRs)
n = 0
for k in argsort:
  DEMDIR = DEMDIRs[k]

  # Mesh boundaries
  extent = np.loadtxt(DIR+DEMDIR+"/inputs/mesh_extent.dat")
  # extent = meshlib.shp_to_xy(os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/Glaciers/3D/"+glacier+"/glacier_extent_inversion_small_front.shp"))
  try:
    hole1 = np.loadtxt(DIR+DEMDIR+"/inputs/mesh_hole1.dat")
    hole2 = np.loadtxt(DIR+DEMDIR+"/inputs/mesh_hole2.dat")
    holes=[hole1,hole2]  
  except:
    holes = []

  # Get results from inversion
  invdir = DIR+DEMDIR+"/mesh2d/inversion_"+method+"/"
  regpardirs = os.listdir(invdir)
  for regpardir in regpardirs:
    if (regpardir.startswith('lambda_'+regpar)) and not(regpardir.endswith('pdf')):
      finaldir = invdir+regpardir+'/'
  bed = elmerreadlib.saveline_boundary(finaldir,method+"_beta.dat",bbed,['velocity','beta'])
  surf = elmerreadlib.saveline_boundary(finaldir,method+"_beta.dat",bsur,['velocity','vsurfini'])

  # Get bed,surface from inputs to calculate driving stress
  xgrid,ygrid,zbgrid = elmerreadlib.input_file(DIR+DEMDIR+"/inputs/zbdem.xy")  
  xgrid,ygrid,zsgrid = elmerreadlib.input_file(DIR+DEMDIR+"/inputs/zsdem.xy")

  if n == 0:
    nx = len(xgrid)
    ny = len(ygrid)
    zb = np.zeros([ny,nx,nt])
    zs = np.zeros([ny,nx,nt])
    times = np.zeros([nt])
    taubs = np.zeros([ny,nx,nt])
    tauds = np.zeros([ny,nx,nt])
    tauds_flow = np.zeros([ny,nx,nt])
    betas = np.zeros([ny,nx,nt])
    vsurf = np.zeros([ny,nx,nt])
    vbed = np.zeros([ny,nx,nt])
    xb = xgrid
    yb = ygrid
    xm = bed['x']
    ym = bed['y']
    xb_grid,yb_grid = np.meshgrid(xb,yb)
    xb_flat = xb_grid.flatten()
    yb_flat = yb_grid.flatten()
    del xb_grid,yb_grid
    
    # Make masks
    masks = np.ones([ny,nx,nt],dtype='bool')
    if len(extent) > 0:
      exterior = Path(np.column_stack((extent[:,0],extent[:,1])))
      pts = exterior.contains_points(np.column_stack((xb_flat,yb_flat)))
      bool = ~(np.reshape(pts,(ny,nx)))
      masks[bool==0,:] = 0
  
    if len(holes) > 0:
      # Mask out points inside holes
      for i in range(0,len(holes)):
        hole = Path(np.column_stack((holes[i][:,0],holes[i][:,1])))
        pts = hole.contains_points(np.column_stack((xb_flat,yb_flat)))
        bool = (np.reshape(pts,(ny,nx)))
        masks[bool==1,:] = 1

  taubs[:,:,n] = np.reshape(scipy.interpolate.griddata((bed['x'],bed['y']),bed['taub'],(xb_flat,yb_flat)),(ny,nx))
  betas[:,:,n] = np.reshape(scipy.interpolate.griddata((bed['x'],bed['y']),bed['beta']**2,(xb_flat,yb_flat)),(ny,nx))
  vbed[:,:,n] = np.reshape(scipy.interpolate.griddata((bed['x'],bed['y']),bed['velocity'],(xb_flat,yb_flat)),(ny,nx))
  vsurf[:,:,n] = np.reshape(scipy.interpolate.griddata((surf['x'],surf['y']),surf['vsurfini'],(xb_flat,yb_flat)),(ny,nx))
  times[n] = datelib.date_to_fracyear(int(DEMDIR[3:7]),int(DEMDIR[7:9]),int(DEMDIR[9:11]))

  # Get velocities and interpolate to grid
  xgrid,ygrid,ugrid = elmerreadlib.input_file(DIR+DEMDIR+"/inputs/udem.xy")
  xgrid,ygrid,vgrid = elmerreadlib.input_file(DIR+DEMDIR+"/inputs/vdem.xy")
  fu = scipy.interpolate.RegularGridInterpolator((ygrid,xgrid),ugrid)
  fv = scipy.interpolate.RegularGridInterpolator((ygrid,xgrid),vgrid)
  u = np.reshape(fu((yb_flat,xb_flat)),(ny,nx))
  v = np.reshape(fv((yb_flat,xb_flat)),(ny,nx))

  zs[:,:,n] = zsgrid
  zb[:,:,n] = zbgrid

  # Calculate grad(h) & H for driving stress
  dhdx = np.zeros([ny,nx])
  dhdy = np.zeros([ny,nx])
  H = zsgrid-zbgrid

  dhdx[:,:] = float('nan')
  dhdy[:,:] = float('nan')

  dhdx[1:-1,1:-1] = (zsgrid[1:-1,2:]-zsgrid[1:-1,0:-2])/(xb[2:]-xb[0:-2])
  dhdy[1:-1,1:-1] = ((zsgrid[2:,1:-1]-zsgrid[0:-2,1:-1]).T/(yb[2:]-yb[0:-2])).T
  tauds[:,:,n] = rho_i*g*H*np.sqrt((dhdx)**2+(dhdy)**2)
  tauds_flow[:,:,n] = -rho_i*g*H*(dhdx*u+dhdy*v)/np.sqrt(u**2+v**2)

  n = n+1

# Mask variables outside of mesh extent
vsurf = np.ma.masked_array(vsurf,masks)
vbed = np.ma.masked_array(vbed,masks)
taubs = np.ma.masked_array(taubs,masks)
tauds = np.ma.masked_array(tauds,masks)
betas = np.ma.masked_array(betas,masks)
tauds_flow = np.ma.masked_array(tauds_flow,masks)

################
# Plot options #
################

# Extent
xmin = np.min(extent[:,0])
xmax = np.max(extent[:,0])
ymin = np.min(extent[:,1])
ymax = np.max(extent[:,1])

############################################
# Plot means and variability for variables #
############################################

plt.figure(figsize=(7,5))
gs = matplotlib.gridspec.GridSpec(2,3)
matplotlib.rc('font',family='Arial')

plt.subplot(gs[0,0])
plt.imshow(np.mean(tauds_flow,axis=2)/1e3,extent=[xb[0],xb[-1],yb[0],yb[-1]],origin='lower',vmin=-500,vmax=500,cmap='RdBu_r')
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
plt.xticks([])
plt.yticks([])
cb = plt.colorbar(ticks=np.arange(-400,600,200))
cb.set_label(r'$\tau_d*u/||u||$ (kPa)',size=8,fontname='arial')
cb.ax.tick_params(labelsize=8)

plt.subplot(gs[0,1])
plt.imshow(np.mean(taubs,axis=2)*1e3,extent=[xb[0],xb[-1],yb[0],yb[-1]],origin='lower',vmin=0,vmax=500)
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
plt.xticks([])
plt.yticks([])
cb = plt.colorbar(ticks=np.arange(0,600,100))
cb.set_label(r'$\tau_b$ (kPa)',size=8,fontname='arial')
cb.ax.tick_params(labelsize=8)

plt.subplot(gs[0,2])
plt.imshow(np.mean(vsurf,axis=2)/1e3,extent=[xb[0],xb[-1],yb[0],yb[-1]],origin='lower',vmin=0,vmax=8)
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
plt.xticks([])
plt.yticks([])
cb = plt.colorbar(ticks=np.arange(0,10,2))
cb.set_label('Velocity (km/yr)',size=8,fontname='arial')
cb.ax.tick_params(labelsize=8)

plt.subplot(gs[1,0])
plt.imshow(np.std(tauds_flow,axis=2)/1e3,extent=[xb[0],xb[-1],yb[0],yb[-1]],origin='lower',vmin=0,vmax=100)
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
plt.xticks([])
plt.yticks([])
cb = plt.colorbar(ticks=np.arange(0,120,20))
cb.set_label(r'$\tau_d$ variability (kPa)',size=8,fontname='arial')
cb.ax.tick_params(labelsize=8)

plt.subplot(gs[1,1])
plt.imshow(np.std(taubs,axis=2)*1e3,extent=[xb[0],xb[-1],yb[0],yb[-1]],origin='lower',vmin=0,vmax=100)
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
plt.xticks([])
plt.yticks([])
cb = plt.colorbar(ticks=np.arange(0,120,20))
cb.set_label(r'$\tau_b$ variability (kPa)',size=8,fontname='arial')
cb.ax.tick_params(labelsize=8)

plt.subplot(gs[1,2])
plt.imshow(np.std(vsurf,axis=2)/1e3,extent=[xb[0],xb[-1],yb[0],yb[-1]],origin='lower',vmin=0,vmax=0.500)
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
plt.xticks([])
plt.yticks([])
cb = plt.colorbar(ticks=np.arange(0,1,0.2))
cb.set_label(r'Velocity variability (km/yr)',size=8,fontname='arial')
cb.ax.tick_params(labelsize=8)

plt.tight_layout()
plt.subplots_adjust(hspace=0.15,wspace=0.1,right=0.95,left=0.01,bottom=0.02,top=0.98)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_inversion_variability.pdf"),FORMAT='PDF',dpi=600)
plt.close()

#########################################################
# Calculate averages for different regions through time #
#########################################################
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

  taubs_U = np.zeros([len(times),3])
  taubs_L = np.zeros([len(times),3])
  taubs_M = np.zeros([len(times),3])
  
  tauds_U = np.zeros([len(times),3])
  tauds_L = np.zeros([len(times),3])
  tauds_M = np.zeros([len(times),3])

  betas_U = np.zeros([len(times),3])
  betas_L = np.zeros([len(times),3])
  betas_M = np.zeros([len(times),3])

  vbeds_U = np.zeros([len(times),3])
  vbeds_L = np.zeros([len(times),3])
  vbeds_M = np.zeros([len(times),3])

  vsurf_U = np.zeros([len(times),3])
  vsurf_L = np.zeros([len(times),3])
  vsurf_M = np.zeros([len(times),3])

  # Percentile range for calculating statistics
  q=50.

  # Get averages for Full Stokes
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

  for i in range(0,len(times)):
    tauds_U[i,0] = np.nanmean(tauds_flow[ind_U_y,ind_U_x,i])/1e3
    tauds_U[i,1] = np.percentile(tauds_flow[ind_U_y,ind_U_x,i],q/2)/1e3
    tauds_U[i,2] = np.percentile(tauds_flow[ind_U_y,ind_U_x,i],100-q/2)/1e3
    taubs_U[i,0] = np.nanmean(taubs[ind_U_y,ind_U_x,i])*1e3
    taubs_U[i,1] = np.percentile(taubs[ind_U_y,ind_U_x,i],q/2)*1e3
    taubs_U[i,2] = np.percentile(taubs[ind_U_y,ind_U_x,i],100-q/2)*1e3 
    betas_U[i,0] = np.nanmean(betas[ind_U_y,ind_U_x,i])
    betas_U[i,1] = np.percentile(betas[ind_U_y,ind_U_x,i],q/2)
    betas_U[i,2] = np.percentile(betas[ind_U_y,ind_U_x,i],100-q/2)    
    vbeds_U[i,0] = np.nanmean(vbed[ind_U_y,ind_U_x,i])
    vbeds_U[i,1] = np.percentile(vbed[ind_U_y,ind_U_x,i],q/2)
    vbeds_U[i,2] = np.percentile(vbed[ind_U_y,ind_U_x,i],100-q/2)   
    vsurf_U[i,0] = np.nanmean(vsurf[ind_U_y,ind_U_x,i])
    vsurf_U[i,1] = np.percentile(vsurf[ind_U_y,ind_U_x,i],q/2)
    vsurf_U[i,2] = np.percentile(vsurf[ind_U_y,ind_U_x,i],100-q/2)

    tauds_L[i,0] = np.nanmean(tauds_flow[ind_L_y,ind_L_x,i])/1e3
    tauds_L[i,1] = np.percentile(tauds_flow[ind_L_y,ind_L_x,i],q/2)/1e3
    tauds_L[i,2] = np.percentile(tauds_flow[ind_L_y,ind_L_x,i],100-q/2)/1e3
    taubs_L[i,0] = np.nanmean(taubs[ind_L_y,ind_L_x,i])*1e3
    taubs_L[i,1] = np.percentile(taubs[ind_L_y,ind_L_x,i],q/2)*1e3
    taubs_L[i,2] = np.percentile(taubs[ind_L_y,ind_L_x,i],100-q/2)*1e3     
    betas_L[i,0] = np.nanmean(betas[ind_L_y,ind_L_x,i])
    betas_L[i,1] = np.percentile(betas[ind_L_y,ind_L_x,i],q/2)
    betas_L[i,2] = np.percentile(betas[ind_L_y,ind_L_x,i],100-q/2)  
    vbeds_L[i,0] = np.nanmean(vbed[ind_L_y,ind_L_x,i])
    vbeds_L[i,1] = np.percentile(vbed[ind_L_y,ind_L_x,i],q/2)
    vbeds_L[i,2] = np.percentile(vbed[ind_L_y,ind_L_x,i],100-q/2)
    vsurf_L[i,0] = np.nanmean(vsurf[ind_L_y,ind_L_x,i])
    vsurf_L[i,1] = np.percentile(vsurf[ind_L_y,ind_L_x,i],q/2)
    vsurf_L[i,2] = np.percentile(vsurf[ind_L_y,ind_L_x,i],100-q/2)   

    tauds_M[i,0] = np.nanmean(tauds_flow[ind_M_y,ind_M_x,i])/1e3
    tauds_M[i,1] = np.percentile(tauds_flow[ind_M_y,ind_M_x,i],q/2)/1e3
    tauds_M[i,2] = np.percentile(tauds_flow[ind_M_y,ind_M_x,i],100-q/2)/1e3
    taubs_M[i,0] = np.nanmean(taubs[ind_M_y,ind_M_x,i])*1e3
    taubs_M[i,1] = np.percentile(taubs[ind_M_y,ind_M_x,i],q/2)*1e3  
    taubs_M[i,2] = np.percentile(taubs[ind_M_y,ind_M_x,i],100-q/2)*1e3
    betas_M[i,0] = np.nanmean(betas[ind_M_y,ind_M_x,i])
    betas_M[i,1] = np.percentile(betas[ind_M_y,ind_M_x,i],q/2)
    betas_M[i,2] = np.percentile(betas[ind_M_y,ind_M_x,i],100-q/2)
    vbeds_M[i,0] = np.nanmean(vbed[ind_M_y,ind_M_x,i])
    vbeds_M[i,1] = np.percentile(vbed[ind_M_y,ind_M_x,i],q/2)
    vbeds_M[i,2] = np.percentile(vbed[ind_M_y,ind_M_x,i],100-q/2)
    vsurf_M[i,0] = np.nanmean(vsurf[ind_M_y,ind_M_x,i])
    vsurf_M[i,1] = np.percentile(vsurf[ind_M_y,ind_M_x,i],q/2)
    vsurf_M[i,2] = np.percentile(vsurf[ind_M_y,ind_M_x,i],100-q/2)


  ############################# 
  # Plot averages for regions #
  #############################

  # Make plot for lower glacier
  
  plt.figure(figsize=(3.75,3))
  gs = matplotlib.gridspec.GridSpec(3,1)
  matplotlib.rc('font',family='Arial') 

  plt.subplot(gs[0])
  ax = plt.gca()
  plt.plot(times,vbeds_L[:,0],'ko--',markersize=5)
  plt.ylabel(r'$u_b$ (km/yr)',fontsize=12,fontname='Arial')
  plt.xlim([time1,time2])
  plt.xticks(np.arange(time1,time2+1,1))
  ax.set_xticklabels([])
  #plt.yticks(np.arange(4500,6500,500))
  #plt.ylim([4500,6000])
  plt.yticks(np.arange(3000,5500,500))
  plt.ylim([3500,5000])
  xTickPos = np.linspace(np.floor(time1)+0.25,np.ceil(time2)+0.25,(np.ceil(time2)-np.floor(time1))*2+1)
  ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.85','w'],linewidth=0)
  ax.text(time1+0.1,3650,'a',fontsize=12,fontname='Arial',fontweight='bold')
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  plt.text(time1+0.1,4650,'Region L',fontsize=12,fontname='Arial',fontweight='bold')
  #plt.text(time1+0.1,5650,'Region L',fontsize=12,fontname='Arial',fontweight='bold')
  ax.text(time1+0.1,3650,'a',fontsize=12,fontname='Arial',fontweight='bold')

  plt.subplot(gs[1])
  ax = plt.gca()
  plt.plot(times,tauds_L[:,0],'bo--',markersize=5,label=r'$\tau_d$')
  plt.plot(0,0,'rs--',markersize=5,label=r'$\tau_b$')
  plt.ylabel(r'$\tau_d$ (kPa)',fontsize=12,fontname='Arial')
  ax.set_ylim([179,191])
  ax.set_yticks(np.arange(180,192,4))
  xTickPos = np.linspace(np.floor(time1)+0.25,np.ceil(time2)+0.25,(np.ceil(time2)-np.floor(time1))*2+1)
  ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.85','w'],linewidth=0)
  ax.set_xlim([time1,time2])
  ax.set_xticklabels([])
  ax.text(time1+0.1,180,'b',fontsize=12,fontname='Arial',fontweight='bold')
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  
  plt.subplot(gs[2])
  ax = plt.gca()
  plt.plot(times,taubs_L[:,0],'ro--',markersize=5)
  plt.ylabel(r'$\tau_b$ (kPa)',fontsize=12,fontname='Arial')
  #ax.set_yticks(np.arange(22,35,4))
  #ax.set_ylim([21,35])
  #ax.text(time1+0.1,22,'c',fontsize=12,fontname='Arial',fontweight='bold')
  ax.set_ylim([59,71])
  ax.set_yticks(np.arange(60,72,4))
  plt.xticks(np.arange(time1,time2+1,1))
  ax.text(time1+0.1,60,'c',fontsize=12,fontname='Arial',fontweight='bold')
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  labels=[]
  plt.xticks(np.arange(2000,2017))
  for i in range(2000,2017):
    labels.append('Jan \n'+str(i))
  plt.xlim([time1,time2])
  xTickPos = np.linspace(np.floor(time1)+0.25,np.ceil(time2)+0.25,(np.ceil(time2)-np.floor(time1))*2+1)
  ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.85','w'],linewidth=0)  
  ax.set_xticklabels(labels,fontsize=12,fontname='Arial')

  plt.tight_layout()
  plt.subplots_adjust(hspace=0.03,wspace=0.03,right=0.94,left=0.2,bottom=0.15,top=0.95)
  plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_inversion_timeseries_lower.pdf"),FORMAT='PDF',dpi=600)
  plt.close()

  # Make plot for upper glacier
  
  plt.figure(figsize=(3.75,3))
  gs = matplotlib.gridspec.GridSpec(3,1)
  matplotlib.rc('font',family='Arial')

  plt.subplot(gs[0])
  ax = plt.gca()
  plt.plot(times,vbeds_U[:,0],'ko--',markersize=5)
  plt.ylabel(r'$u_b$ (m/yr)',fontsize=12,fontname='Arial')
  plt.xlim([time1,time2])
  plt.xticks(np.arange(time1,time2+1,1))
  ax.set_xticklabels([])
  #plt.yticks(np.arange(1950,2155,100))
  #plt.ylim([1900,2100])
  plt.yticks(np.arange(800,1005,100))
  plt.ylim([800,1000])
  xTickPos = np.linspace(np.floor(time1)+0.25,np.ceil(time2)+0.25,(np.ceil(time2)-np.floor(time1))*2+1)
  ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.85','w'],linewidth=0)
  ax.text(time1+0.1,820,'a',fontsize=12,fontname='Arial',fontweight='bold')
  ax.text(time1+0.5,960,'Region U',fontsize=12,fontname='Arial',fontweight='bold')
  ax.text(time1+0.5,2060,'Region U',fontsize=12,fontname='Arial',fontweight='bold')  
  ax.text(time1+0.1,1920,'a',fontsize=12,fontname='Arial',fontweight='bold')  
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))

  plt.subplot(gs[1])
  ax = plt.gca()
  plt.plot(times,tauds_U[:,0],'bo--',markersize=5,label=r'$\tau_d$')
  plt.plot(0,0,'rs--',markersize=5,label=r'$\tau_b$')
  ax.set_ylim([339,351])
  plt.xlim([time1,time2])
  ax.set_yticks(np.arange(340,352,4))
  ax.set_ylabel(r'$\tau_d$ (kPa)',fontsize=12,fontname='Arial')
  xTickPos = np.linspace(np.floor(time1)+0.25,np.ceil(time2)+0.25,(np.ceil(time2)-np.floor(time1))*2+1)
  ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.85','w'],linewidth=0)  
  ax.set_xticklabels([])
  ax.text(time1+0.1,340,'b',fontsize=12,fontname='Arial',fontweight='bold')

  plt.subplot(gs[2])
  ax = plt.gca()
  plt.plot(times,taubs_U[:,0],'ro--',markersize=5)
  #ax.set_ylim([192,207])
  #ax.set_yticks(np.arange(192,208,4))
  ax.set_ylim([231,243])
  ax.set_yticks(np.arange(232,244,4))
  ax.set_ylabel(r'$\tau_b$ (kPa)',fontsize=12,fontname='Arial')
  plt.xticks(np.arange(time1,time2+1,1))
  ax.legend(loc=1,borderpad=0.3,fontsize=10,numpoints=1,handlelength=0.7,labelspacing=0.05,ncol=4,columnspacing=0.7,handletextpad=0.5)
  ax.text(time1+0.1,232,'c',fontsize=12,fontname='Arial',fontweight='bold')
  #ax.text(time1+0.1,200,'c',fontsize=12,fontname='Arial',fontweight='bold')  
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  labels=[]
  plt.xticks(np.arange(2000,2017))
  for i in range(2000,2017):
    labels.append('Jan \n'+str(i))
  plt.xlim([time1,time2])
  xTickPos = np.linspace(np.floor(time1)+0.25,np.ceil(time2)+0.25,(np.ceil(time2)-np.floor(time1))*2+1)
  ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.85','w'],linewidth=0)
  ax.set_xticklabels(labels,fontsize=12,fontname='Arial')

  plt.tight_layout()
  plt.subplots_adjust(hspace=0.03,wspace=0.03,right=0.94,left=0.20,bottom=0.15,top=0.95)
  plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_inversion_timeseries_upper.pdf"),FORMAT='PDF',dpi=600)
  plt.close()

  # Make plot for margins 

  plt.figure(figsize=(3.75,3))
  gs = matplotlib.gridspec.GridSpec(3,1)
  matplotlib.rc('font',family='Arial')

  plt.subplot(gs[0])
  ax = plt.gca()
  plt.plot(times,vbeds_M[:,0],'ko--',markersize=5)
  plt.ylabel(r'$u_b$ (m/yr)',fontsize=12,fontname='Arial')
  plt.xlim([time1,time2])
  plt.xticks(np.arange(time1,time2+1,1))
  ax.set_xticklabels([])
  #plt.yticks(np.arange(300,400,20))
  #plt.ylim([300,380])
  plt.yticks(np.arange(100,200,20))
  plt.ylim([100,180])
  xTickPos = np.linspace(np.floor(time1)+0.25,np.ceil(time2)+0.25,(np.ceil(time2)-np.floor(time1))*2+1)
  ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.85','w'],linewidth=0)
  ax.text(time1+0.1,110,'a',fontsize=12,fontname='Arial',fontweight='bold')
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  plt.text(time1+0.1,160,'Region M',fontsize=12,fontname='Arial',fontweight='bold')
  #plt.text(time1+0.1,360,'Region M',fontsize=12,fontname='Arial',fontweight='bold')
  #ax.text(time1+0.1,310,'a',fontsize=12,fontname='Arial',fontweight='bold')

  plt.subplot(gs[1])
  ax = plt.gca()
  plt.plot(times,tauds_M[:,0],'bo--',markersize=5)
  ax.set_ylim([163,175])
  plt.xlim([time1,time2])
  ax.set_yticks(np.arange(164,176,4))
  xTickPos = np.linspace(np.floor(time1)+0.25,np.ceil(time2)+0.25,(np.ceil(time2)-np.floor(time1))*2+1)
  ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.85','w'],linewidth=0)
  ax.set_xticklabels([])
  ax.text(time1+0.1,164,'b',fontsize=12,fontname='Arial',fontweight='bold')  
  ax.set_ylabel(r'$\tau_d$ (kPa)',fontsize=12,fontname='Arial')

  plt.subplot(gs[2])
  ax = plt.gca()
  plt.plot(times,taubs_M[:,0],'rs--',markersize=5)
  ax.set_ylim([161,173])
  ax.set_yticks(np.arange(162,174,4))
  plt.xticks(np.arange(2000,2018,1))
  labels=[]
  for i in range(2000,2017):
    labels.append('Jan \n'+str(i))
  plt.xlim([time1,time2])
  ax.set_xticklabels(labels,fontsize=12,fontname='Arial')
  ax.text(time1+0.1,162,'c',fontsize=12,fontname='Arial',fontweight='bold')
  ax.xaxis.set_minor_locator(AutoMinorLocator(2))
  xTickPos = np.linspace(np.floor(time1)+0.25,np.ceil(time2)+0.25,(np.ceil(time2)-np.floor(time1))*2+1)
  ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.85','w'],linewidth=0)
  ax.set_ylabel(r'$\tau_b$ (kPa)',fontsize=12,fontname='Arial')

  plt.tight_layout()
  plt.subplots_adjust(hspace=0.03,wspace=0.03,right=0.94,left=0.20,bottom=0.15,top=0.95)
  plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_inversion_timeseries_margin.pdf"),FORMAT='PDF',dpi=600)
  plt.close()


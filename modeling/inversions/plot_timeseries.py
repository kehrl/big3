import sys
import os
import matplotlib.pyplot as plt
from pylab import *
import argparse
import masklib, elmerreadlib, glaclib, colorlib, datelib
import scipy
from matplotlib.path import Path

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

# Input directory
DIR = os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/"+DIR+"/")

# Load flowline
x,y,zb,dists = glaclib.load_flowline(glacier,dx=100)

# Load results
DEMDIRs = os.listdir(DIR)
argsort = np.argsort(DEMDIRs)
nt = len(DEMDIRs)
n = 0
for k in argsort:
  DEMDIR = DEMDIRs[k]

  # Mesh boundaries
  extent = np.loadtxt(DIR+DEMDIR+"/inputs/mesh_extent.dat")
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
  betas[:,:,n] = np.reshape(scipy.interpolate.griddata((bed['x'],bed['y']),bed['beta'],(xb_flat,yb_flat)),(ny,nx))
  vsurf[:,:,n] = np.reshape(scipy.interpolate.griddata((bed['x'],bed['y']),surf['vsurfini'],(xb_flat,yb_flat)),(ny,nx))
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
taubs = np.ma.masked_array(taubs,masks)
tauds = np.ma.masked_array(tauds,masks)
betas = np.ma.masked_array(betas,masks)
tauds_flow = np.ma.masked_array(tauds_flow,masks)

# Extent
xmin = np.min(extent[:,0])
xmax = np.max(extent[:,0])
ymin = np.min(extent[:,1])
ymax = np.max(extent[:,1])

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

 # plt.figure(2)
 # ax2.plot(dists/1e3,flowbed['beta']**2,label=DEMDIR[3:])
  
 # plt.figure(3)
 # ax3.plot(dists/1e3,flowbed['taub']*1e3,label=DEMDIR[3:])

 # plt.figure(4)
 # alpha = np.diff(flowsurf['z'])/np.diff(dists)
 # taud = -1*917.0*9.8*(flowsurf['z'][1:]-flowbed['z'][1:])*alpha
 # ax4.plot(dists[1:]/1e3,taud/1e3,label=DEMDIR[3:])

 # plt.figure(5)
 # ax5.plot(dists[1:]/1e3,flowbed['taub'][1:]*1e3/(taud/1e3),label=DEMDIR[3:])

 # plt.figure(6)
 # plt.scatter(bed['x'],bed['y'],c=bed['taub']*1e3,vmin=0,vmax=500,lw=0)
 # plt.colorbar()
 # plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_inversion_"+DEMDIR[3:]+".pdf"),format='PDF')
 # plt.close()


#plt.figure(1)
#plt.legend(loc=2,ncol=2)
#plt.ylabel('Surface velocity (m/yr)')
#plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_inversion_surface_velocity.pdf"),format='PDF')

#plt.figure(2)
#plt.legend(ncol=2)
#plt.ylabel('Beta^2')
#plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_inversion_betas.pdf"),format='PDF')
#plt.close()

#plt.figure(3)
#plt.legend(ncol=2)
#plt.ylabel('tau_b (kPa)')
#plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_inversion_taub.pdf"),format='PDF')
#plt.close()

#plt.figure(4)
#plt.legend(ncol=2)
#plt.ylabel('tau_d (kPa)')
#plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_inversion_taud.pdf"),format='PDF')
#plt.close()

#plt.figure(5)
#plt.legend(ncol=2)
#plt.ylabel('tau_b/tau_d')
#plt.ylim([0,2])
#nonnan = np.where(~(np.isnan(flowbed['taub'])))[0]
#plt.plot([dists[nonnan[0]]/1e3,dists[nonnan[-1]]/1e3],[1,1],'k',lw=2)
#plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_inversion_taudb.pdf"),format='PDF')
#plt.close()

#plt.figure(1)
#jet = plt.get_cmap('RdBu_r') 
#cNorm  = matplotlib.colors.Normalize(vmin=-800, vmax=800)
#scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm, cmap=jet)
#ind = np.argmin(abs(dists--5.0e3))
#for i in range(0,len(times)):
#  colorval = scalarMap.to_rgba(vel_flow[ind,i]-np.mean(vel_flow[ind,:]))
#  plt.plot(dists,taubs_flow[:,i]*1e3,c=colorval[0])
#plt.plot(dists,np.mean(taubs_flow,axis=1)*1e3,'k')


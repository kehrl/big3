import glaclib, elmerreadlib
import matplotlib.pyplot as plt
import numpy as np
import argparse, os, sys
import matplotlib
import scipy.interpolate

##########
# Inputs #
##########

# Get inputs to file
parser = argparse.ArgumentParser()
parser.add_argument("-glacier",dest="glacier",required = True, 
        help = "Name of glacier (Kanger or Helheim)")
parser.add_argument("-mesh", dest="mesh", required = True,
        help = "Name of mesh directory") 
parser.add_argument("-runname",dest="runname",required = False, 
        default='iceshelf',
        help = "runname in mesh2d directory (default = 'iceshelf'.")
parser.add_argument("-t2",dest="t2",required = False, default=10000,
        help = "runname in mesh2d directory (default = 10000.")

args, _ = parser.parse_known_args(sys.argv)

glacier = args.glacier
meshname = args.mesh
runname = args.runname
t2 = int(args.t2)

DIRM = os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/"+meshname+"/")
DIRR = DIRM+'figures/'

if not(os.path.isdir(DIRM)):
  sys.exit("Mesh directory "+meshname+" does not exist for "+glacier+".")

# Create grid
dx=200
xgrid = np.arange(278000,313000,dx)
ygrid = np.arange(-2583000,-2552000,dx)

# Load flowline
xflow,yflow,zb,dist = glaclib.load_flowline(glacier)

# Get measurements
xdem,ydem,udem = elmerreadlib.input_file(DIRM+'inputs/udem.xy')
xdem,ydem,vdem = elmerreadlib.input_file(DIRM+'inputs/vdem.xy')
xdem,ydem,zsdem = elmerreadlib.input_file(DIRM+'inputs/zsdem.xy')

f = scipy.interpolate.RegularGridInterpolator((ydem,xdem),np.sqrt(udem**2+vdem**2))
xmesh,ymesh = np.meshgrid(xgrid,ygrid)
u_init = (f((ymesh.flatten(),xmesh.flatten()))).reshape(len(ygrid),len(xgrid))

f = scipy.interpolate.RegularGridInterpolator((ydem,xdem),zsdem)
zs_init = (f((ymesh.flatten(),xmesh.flatten()))).reshape(len(ygrid),len(xgrid))

# Points for sampling
dists_eul = np.array([-2,-5,-10,-15,-20,-30])*1e3

model_grid = elmerreadlib.pvtu_timeseries_grid(xgrid,ygrid,DIRM+'mesh2d/',runname,\
        ['velocity'],DIRM+'inputs/',layer='surface',debug=True,t2=t2,crop_mesh=False)

# Get model time
model_time = np.arange(0,np.shape(model_grid)[2])*(7/365.25)

# Get model results at sample points  
vel_model = np.zeros([len(model_time),len(dists_eul)])
zs_model = np.zeros([len(model_time),len(dists_eul)])
zs_meas = np.zeros([len(dists_eul),1])
vel_meas = np.zeros([len(dists_eul),1])
for i in range(0,len(dists_eul)):
  ind = np.argmin(abs(dists_eul[i]-dist))
  indx = np.argmin(abs(xflow[ind]-xgrid))
  indy = np.argmin(abs(yflow[ind]-ygrid))
  for j in range(0,len(model_time)):
    vel_model[j,i] = model_grid['velocity'][indy,indx,j]
    zs_model[j,i] = model_grid['z'][indy,indx,j]    
  zs_meas[i] = zs_init[indy,indx]
  vel_meas[i] = u_init[indy,indx]

###############
# Plot points #
###############

colors = ['r','b','g','limegreen','gold','k']
fig = plt.figure(figsize=(5,5))
gs = matplotlib.gridspec.GridSpec(2,1)

plt.subplot(gs[0, :])
ax = plt.gca()
for i in range(0,len(dists_eul)):
  plt.plot(model_time,(vel_model[:,i]-vel_meas[i])/1e3,color=colors[i],label='H{0:02d}'.format(int(-1*dists_eul[i]/1e3)))
ax.set_xticklabels([])
plt.legend(loc=0,borderpad=0.3,fontsize=8,numpoints=1,handlelength=0.2,handletextpad=0.5,labelspacing=0.1,ncol=2,columnspacing=0.8)
#plt.yticks(np.arange(-3,1,1),fontsize=8)
plt.ylabel('Velocity \n (km/yr)',fontsize=8,fontname='Arial')
plt.xlim([model_time[0],model_time[-1]])
plt.xticks(fontsize=8,fontname='Arial')
plt.yticks(fontsize=8,fontname='Arial')

plt.subplot(gs[1, :])
for i in range(0,len(dists_eul)):
  plt.plot(model_time,zs_model[:,i]-zs_meas[i],color=colors[i])
#plt.yticks(np.arange(-100,25,25),fontsize=8,fontname='Arial')
plt.ylabel('Elevation (m)',fontsize=8)
plt.xlim([model_time[0],model_time[-1]])
plt.xticks(fontsize=8,fontname='Arial')
plt.yticks(fontsize=8,fontname='Arial')

plt.tight_layout()
plt.subplots_adjust(hspace=0.03,wspace=0.02,top=0.96,right=0.95,left=0.15,bottom=0.07)
plt.savefig(DIRR+'model_relaxation_points.pdf')
plt.close()

##################################################
# Plot difference from observations through time #
##################################################

veldiff = np.zeros([len(model_time),1])
zsdiff = np.zeros([len(model_time),1])
veldiff_std = np.zeros([len(model_time),1])
zsdiff_std = np.zeros([len(model_time),1])
for i in range(0,len(model_time)):
  veldiff[i] = np.nanmean(abs(model_grid['velocity'][:,:,i]-u_init))
  zsdiff[i] = np.nanmean(abs(model_grid['z'][:,:,i]-zs_init))
  veldiff_std[i] = np.nanstd(abs(model_grid['velocity'][:,:,i]-u_init))
  zsdiff_std[i] = np.nanstd(abs(model_grid['z'][:,:,i]-zs_init))

fig = plt.figure(figsize=(5,5))
gs = matplotlib.gridspec.GridSpec(2,1)

plt.subplot(gs[0,:])
ax = plt.gca()
plt.plot(model_time,veldiff,'k')
#plt.plot(model_time,veldiff+veldiff_std,c='0.5')
ax.set_xticklabels([])
plt.ylabel('Velocity difference \n (m/yr)',fontsize=8,fontname='Arial')
plt.xlim([model_time[0],model_time[-1]])
plt.xticks(fontsize=8,fontname='Arial')
plt.yticks(fontsize=8,fontname='Arial')

plt.subplot(gs[1, :])
plt.plot(model_time,zsdiff,'k')
#plt.plot(model_time,zsdiff+zsdiff_std,c='0.5')
plt.ylabel('Elevation difference (m)',fontsize=8)
plt.xlim([model_time[0],model_time[-1]])
plt.xticks(fontsize=8,fontname='Arial')
plt.yticks(fontsize=8,fontname='Arial')

plt.tight_layout()
plt.subplots_adjust(hspace=0.03,wspace=0.02,top=0.96,right=0.95,left=0.15,bottom=0.07)
plt.savefig(DIRR+'model_relaxation_difference.pdf')
plt.close()

####################################################
# Plot modeled and measured velocity at end of run #
####################################################

fig = plt.figure(figsize=(6,4))
gs = matplotlib.gridspec.GridSpec(2,2)

plt.subplot(gs[0,0])
plt.imshow(model_grid['velocity'][:,:,-1],extent=[xgrid[0],xgrid[-1],ygrid[0],ygrid[-1]],origin='lower',clim=[0,8000])
plt.xticks([])
plt.yticks([])
cb=plt.colorbar(ticks=[0,4000,8000])
cb.set_label('Velocity (m/yr)',size=8,fontname='arial')
cb.ax.tick_params(labelsize=8)

plt.subplot(gs[0,1])
plt.imshow((model_grid['velocity'][:,:,-1]-u_init),extent=[xgrid[0],xgrid[-1],ygrid[0],ygrid[-1]],origin='lower',cmap='RdBu_r',clim=[-1000,1000])
plt.xticks([])
plt.yticks([])
cb = plt.colorbar(ticks=np.arange(-1000,1500,500))
cb.set_label('Relaxed-Measured (m/yr)',size=8,fontname='arial')
cb.ax.tick_params(labelsize=8)

plt.subplot(gs[1,0])
plt.imshow(model_grid['z'][:,:,-1],extent=[xgrid[0],xgrid[-1],ygrid[0],ygrid[-1]],origin='lower',clim=[0,1e3])
plt.xticks([])
plt.yticks([])
cb=plt.colorbar(ticks=[0,500,1000])
cb.set_label('Elevation (m asl)',size=8,fontname='arial')
cb.ax.tick_params(labelsize=8)

plt.subplot(gs[1,1])
plt.imshow(model_grid['z'][:,:,-1]-zs_init,extent=[xgrid[0],xgrid[-1],ygrid[0],ygrid[-1]],origin='lower',cmap='RdBu_r',clim=[-30,30])
plt.xticks([])
plt.yticks([])
cb=plt.colorbar(ticks=[-30,0,30])
cb.set_label('Relaxed-Measured (m/yr)',size=8,fontname='Arial')
cb.ax.tick_params(labelsize=8)

plt.tight_layout()
plt.subplots_adjust(hspace=0.08,wspace=0.15,top=0.95,right=0.93,left=0.01,bottom=0.01)
plt.savefig(DIRR+'model_relaxation_final.pdf')
plt.close()


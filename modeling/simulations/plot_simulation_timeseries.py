import elmerreadlib, glaclib, vellib, zslib, datelib, icefrontlib
import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import argparse, os, sys, cubehelix
import subprocess
from pandas import rolling_median


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
        default='terminusdriven',
        help = "runname in mesh2d directory (default = 'terminusdriven'.")
parser.add_argument("-t2",dest="t2",required = False, default='none',
        help = "runname in mesh2d directory (default = 'terminusdriven'.")

args, _ = parser.parse_known_args(sys.argv)

glacier = args.glacier
meshname = args.mesh
runname = args.runname
t2 = args.t2

DIRM=os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/"+meshname+"/")

if not(os.path.isdir(DIRM)):
  sys.exit("Mesh directory "+meshname+" does not exist for "+glacier+".")

# Load flowline
print "Loading flowline..."
xflow,yflow,zb,dist = glaclib.load_flowline(glacier)

# Create grid
dx=200
xgrid = np.arange(278000,313000,dx)
ygrid = np.arange(-2583000,-2552000,dx)

# Points for sampling
dists_eul = np.array([-2,-5,-10,-15,-20])*1e3

# Get terminus position along flowline
print "Loading terminus positions..."
terminus_val,terminus_time = icefrontlib.distance_along_flowline(xflow,yflow,dist,glacier,datatypes=['WV','Landsat8','TSX'])

# Get data for comparison
vel_val,vel_time,vel_term = vellib.velocity_along_flowline(xflow,yflow,dist,glacier,data='TSX')
#zs_val,zs_time = zslib.dem_along_flowline(xflow,yflow,glacier)

# Find total number of timesteps, if not specified
print "Loading model (grid and grounding lines).."
if t2 == 'none':
#  model_flow = elmerreadlib.pvtu_timeseries_flowline(xflow,yflow,DIRM+'mesh2d/',runname,\
#        ['velocity'],inputsdir=DIRM+'inputs/',layer='surface',debug=True)
  model_grid = elmerreadlib.pvtu_timeseries_grid(xgrid,ygrid,DIRM+'mesh2d/',runname,\
        ['velocity'],DIRM+'inputs/',layer='surface',debug=True)
else: 
#  model_flow = elmerreadlib.pvtu_timeseries_flowline(xflow,yflow,DIRM+'mesh2d/',runname,\
#        ['velocity'],DIRM+'inputs/',layer='surface',debug=True,t2=int(t2))
  model_grid = elmerreadlib.pvtu_timeseries_grid(xgrid,ygrid,DIRM+'mesh2d/',runname,\
        ['velocity'],DIRM+'inputs/',layer='surface',debug=True,t2=int(t2))

# Get mesh extent files
print "Getting mesh info..."
mesh_extent_x = np.loadtxt(DIRM+'inputs/mesh_timeseries_x.dat')
mesh_extent_y = np.loadtxt(DIRM+'inputs/mesh_timeseries_y.dat')
mesh_hole1 = np.loadtxt(DIRM+'inputs/mesh_hole1.dat')
mesh_hole2 = np.loadtxt(DIRM+'inputs/mesh_hole2.dat')

# Get model time
print "Getting model time..."
fid = open(DIRM+'mesh_info.txt','r')
lines = fid.readlines()
date1 = lines[1][7:-1]
fid.close()
model_time = datelib.date_to_fracyear(int(date1[0:4]),int(date1[4:6]),int(date1[6:8]))+((np.arange(0,np.shape(model_grid)[2])-1)/365.25)
    
print "Getting model results at sample points..."    
vel_model = np.zeros([len(model_time),len(dists_eul)])
zs_model = np.zeros([len(model_time),len(dists_eul)])
for i in range(0,len(dists_eul)):
  ind = np.argmin(abs(dists_eul[i]-dist))
  indx = np.argmin(abs(xflow[ind]-xgrid))
  indy = np.argmin(abs(yflow[ind]-ygrid))
  for j in range(0,len(model_time)):
    vel_model[j,i] = model_grid['velocity'][indy,indx,j]
    zs_model[j,i] = model_grid['z'][indy,indx,j]
  
print "Making plot"
cx = cubehelix.cmap(start=1.0,rot=-1.1,reverse=True,minLight=0.05,sat=2)
colors = ['r','b','g','limegreen','gold','k']

fig = plt.figure(figsize=(7.5,2))
gs = matplotlib.gridspec.GridSpec(3,3)
plt.subplot(gs[0:3,0])
ax = plt.gca()
p = plt.imshow(model_grid['velocity'][:,:,-1]/1e3,extent=[np.min(xgrid),np.max(xgrid),\
     np.min(ygrid),np.max(ygrid)],origin='lower',cmap=cx,clim=[0,8.2])
plt.plot(np.r_[mesh_hole1[:,0],mesh_hole1[0,0]],np.r_[mesh_hole1[:,1],mesh_hole1[0,1]],'k',linewidth=0.75,zorder=2)
plt.plot(np.r_[mesh_hole2[:,0],mesh_hole2[0,0]],np.r_[mesh_hole2[:,1],mesh_hole2[0,1]],'k',linewidth=0.75,zorder=2)
ind = np.where(mesh_extent_x[:,-1] != 0)[0]
plt.plot(mesh_extent_x[ind,-1],mesh_extent_y[ind,-1],'k',linewidth=0.75,zorder=2)
plt.xlim([xgrid[0],xgrid[-1]])
plt.ylim([ygrid[0],ygrid[-1]])
plt.xticks([])
plt.yticks([])
  
xmin,xmax = plt.xlim()
ymin,ymax = plt.ylim()  
path = matplotlib.path.Path([[0.55*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
  			[0.99*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
  			[0.99*(xmax-xmin)+xmin,0.72*(ymax-ymin)+ymin],
  			[0.55*(xmax-xmin)+xmin,0.72*(ymax-ymin)+ymin],
  			[0.55*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin]])
patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1,zorder=3)
#ax.add_patch(patch)
cbaxes = fig.add_axes([0.21, 0.9, 0.06, 0.03]) 
cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=[0,4,8]) 
ax.text(xmin+0.75*(xmax-xmin),ymin+0.75*(ymax-ymin),'Velocity',fontsize=8)
ax.text(xmin+0.76*(xmax-xmin),ymin+0.67*(ymax-ymin),'(km/yr)',fontsize=8)
cb.ax.tick_params(labelsize=8)
ax.plot([xmin+0.06*(xmax-xmin),xmin+0.06*(xmax-xmin)+2e3],[ymin+0.08*(ymax-ymin),ymin+0.08*(ymax-ymin)],'k',linewidth=1.5,zorder=4)
ax.plot([xmin+0.06*(xmax-xmin),xmin+0.06*(xmax-xmin)],[ymin+0.08*(ymax-ymin),ymin+0.06*(ymax-ymin)],'k',linewidth=1.5,zorder=4)
ax.plot([xmin+0.06*(xmax-xmin)+2e3,xmin+0.06*(xmax-xmin)+2e3],[ymin+0.08*(ymax-ymin),ymin+0.06*(ymax-ymin)],'k',linewidth=1.5,zorder=4)
ax.text(xmin+0.10*(xmax-xmin)+2e3,ymin+0.05*(ymax-ymin),'2 km',fontsize=8)
for i in range(0,len(dists_eul)):
  ind = np.argmin(abs(dists_eul[i]-dist))
  ax.plot(xflow[ind],yflow[ind],'o',color=colors[i],mec='k',mew=0.5,markersize=5)
  if i in [1,2,4]:
    ax.text(xflow[ind]-4.5e3,yflow[ind]-2e3,glacier[0]+'{0:02d}'.format(int(abs(dists_eul[i]/1e3))),fontsize=8,fontname='Arial')
  elif i == 0:
    ax.text(xflow[ind]-2e3,yflow[ind]-2.75e3,glacier[0]+'{0:02d}'.format(int(abs(dists_eul[i]/1e3))),fontsize=8,fontname='Arial')
  else:
    ax.text(xflow[ind],yflow[ind]+1e3,glacier[0]+'{0:02d}'.format(int(abs(dists_eul[i]/1e3))),fontsize=8,fontname='Arial')
ax.text(xmin+1e3,ymin+0.9*(ymax-ymin),'a',fontsize=8,fontweight='bold')
  
plt.subplot(gs[0:1,1:3])
ax = plt.gca()
plt.plot(terminus_time,terminus_val/1e3,'ko-',markersize=2)
plt.xticks(np.arange(2008,2017,1))
ax.set_xticklabels([])
plt.xlim([model_time[1],model_time[-1]])
plt.yticks(np.arange(-2,4,2),fontsize=8,fontname='Arial')
plt.ylim([-2.5,3])
plt.ylabel('Terminus \n (km)',fontsize=8,fontname='Arial')
#plt.text(2012.55,-2.2,'b',fontsize=8,fontname='Arial',fontweight='bold')
  
plt.subplot(gs[1:3,1:3])
ax = plt.gca()
for j in range(0,len(dists_eul)):
  ind = np.argmin(abs(dists_eul[j]-dist))
  plt.plot(vel_time,vel_val[ind,:]/1e3,'o',color=colors[j],mec='k',mew=0.5,markersize=3,label='H{0:02d}'.format(int(-1*dists_eul[j]/1e3)))
  plt.plot(model_time,vel_model[:,j]/1e3,color=colors[j])
plt.xticks(np.arange(2008,2017),fontsize=8,fontname='Arial')
plt.xlim([model_time[1],model_time[-1]])
plt.yticks(np.arange(4,30,2),fontsize=8,fontname='Arial')
plt.ylim([3,11.5])
plt.ylabel('Velocity \n (km/yr)',fontsize=8,fontname='Arial')
x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
ax.xaxis.set_major_formatter(x_formatter)
#plt.text(2012.55,3.3,'c',fontsize=8,fontname='Arial',fontweight='bold')
plt.legend(loc=2,borderpad=0.3,fontsize=8,numpoints=1,handlelength=0.2,handletextpad=0.5,labelspacing=0.1,ncol=3,columnspacing=0.8)

plt.tight_layout()
plt.subplots_adjust(hspace=0.2,wspace=0.25,top=0.97,right=0.98,left=0.01,bottom=0.1) 
  
plt.savefig(DIRM+'figures/model_'+runname+'.pdf',format='PDF',dpi=600)
plt.close()


import elmerreadlib, glaclib, vellib, zslib, datelib, icefrontlib, matplotlib
import matplotlib.pyplot as plt
import numpy as np
import argparse, os, matplotlib, sys, cubehelix
import subprocess


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
DIRO=os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_"+meshname+"_movie")

if not(os.path.isdir(DIRM)):
  sys.exit("Mesh directory "+meshname+" does not exist for "+glacier+".")

if not(os.path.isdir(DIRO)):
  os.makedirs(DIRO)


# Load flowline
print "Loading flowline..."
xflow,yflow,zb,dist = glaclib.load_flowline(glacier)

# Create grid
dx=100
xgrid = np.arange(278000,313000,dx)
ygrid = np.arange(-2583000,-2552000,dx)

# Points for sampling
dists_eul = np.array([-2,-5,-10,-15,-20])*1e3

# Get terminus position along flowline
print "Loading terminus positions..."
terminus_val,terminus_time = icefrontlib.distance_along_flowline(xflow,yflow,dist,glacier,datatypes=['WV','Landsat8','TSX'])

# Get data for comparison
vel_val,vel_time,vel_term = vellib.velocity_along_flowline(xflow,yflow,dist,glacier)
#zs_val,zs_time = zslib.dem_along_flowline(xflow,yflow,glacier)

# Find total number of timesteps, if not specified
print "Loading model (grid and grounding lines).."
if t2 == 'none':
#  model_flow = elmerreadlib.pvtu_timeseries_flowline(xflow,yflow,DIRM+'mesh2d/',runname,\
#        ['velocity'],inputsdir=DIRM+'inputs/',layer='surface',debug=True)
  model_grid = elmerreadlib.pvtu_timeseries_grid(xgrid,ygrid,DIRM+'mesh2d/',runname,\
        ['velocity'],DIRM+'inputs/',layer='surface',debug=True)
  model_gl = elmerreadlib.pvtu_timeseries_grounding_line(DIRM+'mesh2d/',runname,debug=True)
else: 
#  model_flow = elmerreadlib.pvtu_timeseries_flowline(xflow,yflow,DIRM+'mesh2d/',runname,\
#        ['velocity'],DIRM+'inputs/',layer='surface',debug=True,t2=int(t2))
  model_grid = elmerreadlib.pvtu_timeseries_grid(xgrid,ygrid,DIRM+'mesh2d/',runname,\
        ['velocity'],DIRM+'inputs/',layer='surface',debug=True,t2=int(t2))
  model_gl = elmerreadlib.pvtu_timeseries_grounding_line(DIRM+'mesh2d/',runname,debug=True,t2=int(t2))

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

print "Making plots..."
cx = cubehelix.cmap(start=1.0,rot=-1.1,reverse=True,minLight=0.05,sat=2)
colors = ['r','b','g','limegreen','gold','k']

for i in range(1,len(model_time)):
  fig = plt.figure(figsize=(5.5,4.2))
  gs = matplotlib.gridspec.GridSpec(6,2)
  year,month,day = datelib.fracyear_to_date(model_time[i])

  plt.subplot(gs[0:3,0])
  ax = plt.gca()
  p = plt.imshow(model_grid['velocity'][:,:,i]/365.25,extent=[np.min(xgrid),np.max(xgrid),\
     np.min(ygrid),np.max(ygrid)],origin='lower',cmap=cx,clim=[0,22])
  plt.contour(model_grid['velocity'][:,:,i]/365.25,extent=[np.min(xgrid),np.max(xgrid),\
     np.min(ygrid),np.max(ygrid)],origin='lower',levels=np.arange(0,25,5),cmap=cx,linewidths=1)
  plt.plot(model_gl['x'][:,i],model_gl['y'][:,i],'w.',markersize=0.5)
  plt.plot(np.r_[mesh_hole1[:,0],mesh_hole1[0,0]],np.r_[mesh_hole1[:,1],mesh_hole1[0,1]],'k',linewidth=0.75,zorder=2)
  plt.plot(np.r_[mesh_hole2[:,0],mesh_hole2[0,0]],np.r_[mesh_hole2[:,1],mesh_hole2[0,1]],'k',linewidth=0.75,zorder=2)
  ind = np.where(mesh_extent_x[:,i] != 0)[0]
  plt.plot(mesh_extent_x[ind,i],mesh_extent_y[ind,i],'k',linewidth=0.75,zorder=2)
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
  cbaxes = fig.add_axes([0.325, 0.945, 0.09, 0.02]) 
  cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=[0,10,20]) 
  ax.text(xmin+0.75*(xmax-xmin),ymin+0.75*(ymax-ymin),'Velocity',fontsize=8)
  ax.text(xmin+0.76*(xmax-xmin),ymin+0.67*(ymax-ymin),'(m d$^{-1}$)',fontsize=8)
  cb.ax.tick_params(labelsize=8)
  ax.plot([xmin+0.06*(xmax-xmin),xmin+0.06*(xmax-xmin)+2e3],[ymin+0.08*(ymax-ymin),ymin+0.08*(ymax-ymin)],'k',linewidth=1.5,zorder=4)
  ax.plot([xmin+0.06*(xmax-xmin),xmin+0.06*(xmax-xmin)],[ymin+0.08*(ymax-ymin),ymin+0.06*(ymax-ymin)],'k',linewidth=1.5,zorder=4)
  ax.plot([xmin+0.06*(xmax-xmin)+2e3,xmin+0.06*(xmax-xmin)+2e3],[ymin+0.08*(ymax-ymin),ymin+0.06*(ymax-ymin)],'k',linewidth=1.5,zorder=4)
  ax.text(xmin+0.10*(xmax-xmin)+2e3,ymin+0.05*(ymax-ymin),'2 km',fontsize=8)
  ax.text(xmin+0.02*(xmax-xmin),ymin+0.94*(ymax-ymin),'{0:02d} {1:02d} {2}'.format(int(month),int(np.round(day)),int(year)),fontname='Arial',fontsize=8)

  plt.subplot(gs[3:6,0])
  ax = plt.gca()
  if i > 0:
    diff = (model_grid['z'][:,:,i]-model_grid['z'][:,:,i-1])
    p = plt.imshow(diff*100,extent=[np.min(xgrid),np.max(xgrid),\
       np.min(ygrid),np.max(ygrid)],origin='lower',cmap='RdBu_r',clim=[-15,15])
  plt.plot(np.r_[mesh_hole1[:,0],mesh_hole1[0,0]],np.r_[mesh_hole1[:,1],mesh_hole1[0,1]],'k',linewidth=0.75)
  plt.plot(np.r_[mesh_hole2[:,0],mesh_hole2[0,0]],np.r_[mesh_hole2[:,1],mesh_hole2[0,1]],'k',linewidth=0.75)
  ind = np.where(mesh_extent_x[:,i] != 0)[0]
  plt.plot(mesh_extent_x[ind,i],mesh_extent_y[ind,i],'k',linewidth=0.75)
  for j in range(0,len(dists_eul)):
    ind = np.argmin(abs(dists_eul[j]-dist))
    plt.plot(xflow[ind],yflow[ind],'o',color=colors[j],mec='k',mew=0.5,markersize=4)
  plt.xlim([xgrid[0],xgrid[-1]])
  plt.ylim([ygrid[0],ygrid[-1]])
  plt.xticks([])
  plt.yticks([])
  cbaxes = fig.add_axes([0.325, 0.49, 0.09, 0.02]) 
  cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=[-10,0,10]) 
  ax.text(xmin+0.77*(xmax-xmin),ymin+0.76*(ymax-ymin),'dh/dt',fontsize=8)
  ax.text(xmin+0.74*(xmax-xmin),ymin+0.68*(ymax-ymin),'(cm d$^{-1}$)',fontsize=8)
  cb.ax.tick_params(labelsize=8)
  
  plt.subplot(gs[0:2,1])
  ax = plt.gca()
  plt.plot(terminus_time,terminus_val/1e3,'k')
  plt.xticks(np.arange(2008,2016,.5))
  ax.set_xticklabels([])
  plt.xlim([model_time[1],model_time[-1]])
  plt.ylim([-3,3])
  plt.ylabel('Terminus (km)',fontsize=10,fontname='Arial')
  
  plt.subplot(gs[2:4,1])
  ax = plt.gca()
  for j in range(0,len(dists_eul)):
    plt.plot(model_time[0:i+1],vel_model[0:i+1,j]/365.25,color=colors[j],label='H{0:02d}'.format(int(-1*dists_eul[j]/1e3)))
    #ind = np.where(vel_time < model_time[i+1])
    ind = np.argmin(abs(dists_eul[j]-dist))
    plt.plot(vel_time,vel_val[ind,:]/365.25,'o',color=colors[j],mec='k',mew=0.5)
  plt.xticks(np.arange(2008,2016,.5))
  plt.xlim([model_time[1],model_time[-1]])
  plt.yticks(np.arange(5,30,5))
  plt.ylim([10,24])
  plt.ylabel('Velocity (km/yr)',fontsize=10,fontname='Arial')
  x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
  ax.xaxis.set_major_formatter(x_formatter)
  
  plt.subplot(gs[4:,1])
  ax = plt.gca()
  for j in range(0,2):
    plt.plot(model_time[0:i+1],zs_model[0:i+1,j]-zs_model[0,j],color=colors[j],label='H{0:02d}'.format(int(-1*dists_eul[j]/1e3)))
    #ind = np.where(vel_time < model_time[i+1])
    #plt.plot(vel_time,vel_val[:,j]/365.25,'o',color=colors[j],mec='k',mew=0.5)
  plt.xticks(np.arange(2008,2016,.5))
  plt.xlim([model_time[1],model_time[-1]])
  plt.yticks(np.arange(-20,20,10))
  plt.ylim([-20,20])
  plt.ylabel('Elevation (m)',fontsize=10,fontname='Arial')
  x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
  ax.xaxis.set_major_formatter(x_formatter)

  plt.tight_layout()
  plt.subplots_adjust(hspace=0.1,wspace=0.25,top=0.98,right=0.98,left=0.02,bottom=0.08) 
  
  plt.savefig(DIRO+'/{0:04g}'.format(i)+'.png',format='PNG',dpi=400)
  plt.close()
  
fps=15
ffmpeg_in_opt = "-r %i" % fps
#ffmpeg_out_opt = "-r 5 -y -an -force_fps -s '%ix%i'" % newsize
#-an disables audio
#ffmpeg_out_opt = "-y -an -c:v libx264 -pix_fmt yuv420p"
ffmpeg_out_opt = "-an -vcodec mpeg4 -f mp4"
#scale='iw/4:-1'"

os.chdir(DIRO)
outmov = glacier+'_'+meshname+'.mp4'
#cmd = 'ffmpeg {0} -i %04d.jpg {1} {2}'.format(ffmpeg_in_opt, ffmpeg_out_opt, outmov)
#cmd = 'ffmpeg {0} -pattern_type glob -i *_clip.png {1} {2}'.format(ffmpeg_in_opt, ffmpeg_out_opt, outmov)
cmd = 'ffmpeg {0} -i %04d.png {1} {2}'.format(ffmpeg_in_opt, ffmpeg_out_opt, outmov)

print cmd
subprocess.call(cmd, shell=True)

#cmd = 'rm *.png'
#print cmd
#subprocess.call(cmd, shell=True)
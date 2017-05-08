import os
import sys
import numpy as np
import vellib, bedlib, zslib, floatlib, geotifflib, datelib, glaclib
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import AutoMinorLocator
import scipy.signal as signal
from scipy import stats
import scipy
import pylab, cubehelix
import matplotlib.cm as cmx
import matplotlib.colors as colors


##########
# Inputs #
##########

# Get arguments
args = sys.argv
glacier = args[1][:] # Options: Kanger, Helheim

time1 = 2011.
time2 = 2015.

# Map extent
if glacier == 'Kanger':
  xmin = 468000.
  xmax = 498000.
  ymin = -2299000.
  ymax = -2264000.
elif glacier == 'Helheim':
  xmin = 283000.
  xmax = 313000.
  ymin = -2587000.
  ymax = -2552000.

# Image for plotting
if glacier == "Helheim":
  imagetime = datelib.date_to_fracyear(2014,7,4)
  ximage,yimage,image = geotifflib.readrgb(os.path.join(os.getenv("DATA_HOME"),"Imagery/Landsat/Helheim/TIF/20140704140535_LC82330132014185LGN00.tif"))
elif glacier == "Kanger":
  imagetime = datelib.date_to_fracyear(2014,7,6)
  ximage,yimage,image = geotifflib.readrgb(os.path.join(os.getenv("DATA_HOME"),"Imagery/Landsat/Kanger/TIF/20140706135251_LC82310122014187LGN00.tif"))

# Load bed
xb,yb,zb = bedlib.morlighem_grid(xmin,xmax,ymin,ymax,verticaldatum='geoid')
#xb,yb,zb = bedlib.smith_grid(glacier,xmin,xmax,ymin,ymax,verticaldatum='geoid',model='aniso',smoothing=4,grid='structured')

# Load velocity variability
xvel,yvel,velall,veltrend,veldetrend,velrange,velcount,veltime = vellib.variability(glacier,time1,time2)

# Load elevation variability
xzs,yzs,zsall,zstrend,zsdetrend,zsrange,zscount,zstrend_error,zstime = zslib.variability(glacier,time1,time2)

# Find flotation conditions
xwv,ywv,zwv,timewv = zslib.dem_grid(glacier,xmin,xmax,ymin,ymax,years='all',resolution=32,verticaldatum='geoid')
xf,yf,zabovefloat = floatlib.extent(xwv,ywv,zwv,timewv,glacier,rho_i=917.0,rho_sw=1025.0,bedsource='cresis',verticaldatum='geoid')

cutoff = 5.
floatcond = np.zeros(len(xf))
floatcond[:] = float('nan')
for i in range(0,len(xf)):
  if np.nanmin(zabovefloat[i,:]) > cutoff:
    floatcond[i] = 1 #grounded
  elif np.nanmax(zabovefloat[i,:]) < -1*cutoff:
    floatcond[i] = -1 #floating
  elif (np.nanmax(zabovefloat[i,:]) > -1*cutoff) or (np.nanmin(zabovefloat[i,:]) < cutoff):
    floatcond[i] = 0 #changing basal conditions

# Get flowline to plot some points for reference to other figures
xflow,yflow,zbflow,distflow = glaclib.load_flowline(glacier,shapefilename='flowline_flightline',filt_len=2.0e3)

# Get points for plotting
dists_eul = -1.0*np.array([0,5.0,10.0,15.0,20.0]) # kilometers
coloptions=['w','b','g','limegreen','gold']

# Find where points are located along flowline
ind_eul=[]
for i in range(0,len(dists_eul)):
  ind_eul.append( (abs(distflow - dists_eul[i]*1e3)).argmin() )


############################
# Plot spatial variability #
############################

fig = plt.figure(figsize=(7.5,2.4))
matplotlib.rc('font',family='Arial')
gs = matplotlib.gridspec.GridSpec(1,4)
cx = cubehelix.cmap(start=1.2,rot=-1.1,reverse=True,minLight=0.1,sat=2)

plt.subplot(gs[0])
ax1 = plt.gca()
plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
p=plt.imshow(zb/1.0e3,extent=[np.min(xb),np.max(xb),np.min(yb),np.max(yb)],origin='lower',cmap='RdBu_r',clim=[-1.5,1.5])
plt.contour(zb,levels=[0],extent=[np.min(xb),np.max(xb),np.min(yb),np.max(yb)],origin='lower',colors='k',lw=2,zorder=2)
ax1.axes.set_xlim([xmin,xmax])
ax1.axes.set_ylim([ymin,ymax])
ax1.set_xticks([])
ax1.set_yticks([])
xmin1,xmax1 = plt.xlim()
ymin1,ymax1 = plt.ylim()
path = matplotlib.path.Path([[0.49*(xmax1-xmin1)+xmin,0.98*(ymax1-ymin1)+ymin1],
  			[0.98*(xmax1-xmin1)+xmin1,0.98*(ymax1-ymin1)+ymin1],
  			[0.98*(xmax1-xmin1)+xmin1,0.67*(ymax1-ymin1)+ymin1],
  			[0.49*(xmax1-xmin1)+xmin1,0.67*(ymax1-ymin1)+ymin1],
  			[0.49*(xmax1-xmin1)+xmin1,0.98*(ymax1-ymin1)+ymin1]])
patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1,zorder=3)
ax1.add_patch(patch)
cbaxes = fig.add_axes([0.14, 0.87, 0.085, 0.03]) 
cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=[-1,0,1]) 
cb.set_label('Bed elevation \n (km asl)',size=8,fontname='arial')
cb.ax.tick_params(labelsize=8)
for i in range(0,len(dists_eul)):
  ax1.plot(xflow[ind_eul[i]],yflow[ind_eul[i]],'o',color=coloptions[i],markersize=5,mec='k',mew=0.5)
  if glacier == 'Kanger':
    if i in [1,2,3,4]:
      ax1.text(xflow[ind_eul[i]]-3.5e3,yflow[ind_eul[i]]+1e3,glacier[0]+'{0:02d}'.format(int(abs(dists_eul[i]))),fontsize=8,fontname='Arial')
    else:
      ax1.text(xflow[ind_eul[i]]-4.5e3,yflow[ind_eul[i]]-500,glacier[0]+'{0:02d}'.format(int(abs(dists_eul[i]))),fontsize=8,fontname='Arial')
  else:
    if i in [3,4]:
      ax1.text(xflow[ind_eul[i]]-3.2e3,yflow[ind_eul[i]]+1e3,glacier[0]+'{0:02d}'.format(int(abs(dists_eul[i]))),fontsize=8,fontname='Arial')
    else:
      ax1.text(xflow[ind_eul[i]]-3.3e3,yflow[ind_eul[i]]-2.5e3,glacier[0]+'{0:02d}'.format(int(abs(dists_eul[i]))),fontsize=8,fontname='Arial')
#ax1.plot([xmin+0.61*(xmax-xmin),xmin+0.61*(xmax-xmin)+5e3],[ymin+0.73*(ymax-ymin),ymin+0.73*(ymax-ymin)],'k',linewidth=1.5)
#ax1.plot([xmin+0.61*(xmax-xmin),xmin+0.61*(xmax-xmin)],[ymin+0.73*(ymax-ymin),ymin+0.71*(ymax-ymin)],'k',linewidth=1.5)
#ax1.plot([xmin+0.61*(xmax-xmin)+5e3,xmin+0.61*(xmax-xmin)+5e3],[ymin+0.73*(ymax-ymin),ymin+0.71*(ymax-ymin)],'k',linewidth=1.5)
#ax1.text(xmin+0.64*(xmax-xmin)+5e3,ymin+0.7*(ymax-ymin),'5 km',fontsize=9)
ax1.text(0.05*(xmax1-xmin1)+xmin1,ymin1+0.05*(ymax1-ymin1),'a',fontweight='bold',fontname='arial',fontsize=9)

plt.subplot(gs[1])
ax2 = plt.gca()
ax2.axes.autoscale(False)
plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
norms = matplotlib.colors.BoundaryNorm(np.arange(0,2.1,0.2),cx.N)
p=plt.imshow(velrange/1.0e3,clim=[0,2],extent=[np.min(xvel),np.max(xvel),np.min(yvel),np.max(yvel)],origin='lower',cmap=cx)
#plt.contour(velstd_high/1.0e3,levels=[0.2],extent=[np.min(xvel),np.max(xvel),np.min(yvel),np.max(yvel)],origin='lower',colors='k',lw=2,zorder=2)
plt.contour(zb,levels=[0],extent=[np.min(xb),np.max(xb),np.min(yb),np.max(yb)],origin='lower',colors='k',lw=2,zorder=2)
ax2.set_xticks([])
ax2.set_yticks([])
ax2.axes.set_xlim([xmin,xmax])
ax2.axes.set_ylim([ymin,ymax])
xmin1,xmax1 = plt.xlim()
ymin1,ymax1 = plt.ylim()
path = matplotlib.path.Path([[0.49*(xmax1-xmin1)+xmin,0.98*(ymax1-ymin1)+ymin1],
  			[0.98*(xmax1-xmin1)+xmin1,0.98*(ymax1-ymin1)+ymin1],
  			[0.98*(xmax1-xmin1)+xmin1,0.67*(ymax1-ymin1)+ymin1],
  			[0.49*(xmax1-xmin1)+xmin1,0.67*(ymax1-ymin1)+ymin1],
  			[0.49*(xmax1-xmin1)+xmin1,0.98*(ymax1-ymin1)+ymin1]])
patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1,zorder=3)
ax2.add_patch(patch)
cbaxes = fig.add_axes([0.39, 0.87, 0.085, 0.03]) 
cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=[0,1,2]) 
cb.set_label(r"Velocity range"+"\n"+"(km yr$^{-1}$)",size=8,fontname='arial')
cb.ax.tick_params(labelsize=8)
minorticks = p.norm([0.5,1.5])
cb.ax.xaxis.set_ticks(minorticks, minor=True)
for i in range(0,len(dists_eul)):
  ax2.plot(xflow[ind_eul[i]],yflow[ind_eul[i]],'o',color=coloptions[i],markersize=5,mec='k',mew=0.5)
ax2.text(0.05*(xmax1-xmin1)+xmin1,ymin1+0.05*(ymax1-ymin1),'b',fontweight='bold',fontname='arial',fontsize=9)

plt.subplot(gs[2])
ax3 = plt.gca()
ax3.axes.autoscale(False)
plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
norms = matplotlib.colors.BoundaryNorm(np.arange(0,44,4),cx.N)
p=plt.imshow(zsrange,clim=[0,40],extent=[np.min(xzs),np.max(xzs),np.min(yzs),np.max(yzs)],origin='lower',cmap=cx)
#plt.contour(velstd_high/1.0e3,levels=[0.2],extent=[np.min(xvel),np.max(xvel),np.min(yvel),np.max(yvel)],origin='lower',colors='k',lw=2,zorder=2)
plt.contour(zb,levels=[0],extent=[np.min(xb),np.max(xb),np.min(yb),np.max(yb)],origin='lower',colors='k',lw=2,zorder=2)
#plt.contour(zsstd_high*1.96,levels=[5],extent=[np.min(xvel),np.max(xvel),np.min(yvel),np.max(yvel)],origin='lower',colors='k',lw=2)
ax3.set_xticks([])
ax3.set_yticks([])
ax3.axes.set_xlim([xmin,xmax])
ax3.axes.set_ylim([ymin,ymax])
xmin1,xmax1 = plt.xlim()
ymin1,ymax1 = plt.ylim()
path = matplotlib.path.Path([[0.49*(xmax1-xmin1)+xmin,0.98*(ymax1-ymin1)+ymin1],
  			[0.98*(xmax1-xmin1)+xmin1,0.98*(ymax1-ymin1)+ymin1],
  			[0.98*(xmax1-xmin1)+xmin1,0.67*(ymax1-ymin1)+ymin1],
  			[0.49*(xmax1-xmin1)+xmin1,0.67*(ymax1-ymin1)+ymin1],
  			[0.49*(xmax1-xmin1)+xmin1,0.98*(ymax1-ymin1)+ymin1]])
patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1,zorder=3)
ax3.add_patch(patch)
cbaxes = fig.add_axes([0.635, 0.87, 0.085, 0.03]) 
cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=[0,20,40]) 
cb.set_label(r'Elevation range'+'\n (m)',size=8,fontname='arial')
cb.ax.tick_params(labelsize=8)
minorticks = p.norm(np.arange(0,50,10))
cb.ax.xaxis.set_ticks(minorticks, minor=True)
for i in range(0,len(dists_eul)):
  ax3.plot(xflow[ind_eul[i]],yflow[ind_eul[i]],'o',color=coloptions[i],markersize=5,mec='k',mew=0.5)
ax3.text(0.05*(xmax1-xmin1)+xmin1,ymin1+0.05*(ymax1-ymin1),'c',fontweight='bold',fontname='arial',fontsize=9)

plt.subplot(gs[3])
ax3 = plt.gca()
ax3.axes.autoscale(False)
plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
plt.plot(0,0,'ro',label='Grounded',markersize=4)
plt.plot(0,0,'o',color='cyan',label='Floating',markersize=4)
plt.plot(0,0,'ko',label='Changing',markersize=4)
ind = np.where(floatcond == 1)
plt.plot(xf[ind],yf[ind],'r.',lw=0,markersize=2.75)
ind = np.where(floatcond == 0)
plt.plot(xf[ind],yf[ind],'k.',lw=0,markersize=2.75)
ind = np.where(floatcond == -1)
plt.plot(xf[ind],yf[ind],'.',color='cyan',lw=0,markersize=2.75)
#plt.contour(zb,levels=[0],extent=[np.min(xb),np.max(xb),np.min(yb),np.max(yb)],origin='lower',colors='k',lw=2,zorder=2)
#p=plt.imshow(zsstd_high*2,clim=[0,20],extent=[np.min(xzs),np.max(xzs),np.min(yzs),np.max(yzs)],origin='lower',cmap='jet')
#plt.contour(zsstd_high*1.96,levels=[5],extent=[np.min(xvel),np.max(xvel),np.min(yvel),np.max(yvel)],origin='lower',colors='k',lw=2)
ax3.set_xticks([])
ax3.set_yticks([])
ax3.axes.set_xlim([xmin,xmax])
ax3.axes.set_ylim([ymin,ymax])
xmin1,xmax1 = plt.xlim()
ymin1,ymax1 = plt.ylim()
path = matplotlib.path.Path([[0.54*(xmax1-xmin1)+xmin,0.98*(ymax1-ymin1)+ymin1],
  			[0.98*(xmax1-xmin1)+xmin1,0.98*(ymax1-ymin1)+ymin1],
  			[0.98*(xmax1-xmin1)+xmin1,0.67*(ymax1-ymin1)+ymin1],
  			[0.54*(xmax1-xmin1)+xmin1,0.67*(ymax1-ymin1)+ymin1],
  			[0.54*(xmax1-xmin1)+xmin1,0.98*(ymax1-ymin1)+ymin1]])
patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1,zorder=3)
ax3.add_patch(patch)
plt.legend(loc=1,numpoints=1,frameon=False,labelspacing=0.07,handletextpad=0.5,handlelength=0.2,fontsize=8)
ax3.plot([xmin1+0.6*(xmax1-xmin1),xmin+0.6*(xmax1-xmin1)+5e3],[ymin1+0.75*(ymax1-ymin1),ymin1+0.75*(ymax1-ymin1)],'k',linewidth=1.5,zorder=5)
ax3.plot([xmin1+0.6*(xmax1-xmin1),xmin1+0.6*(xmax1-xmin1)],[ymin1+0.75*(ymax1-ymin1),ymin1+0.73*(ymax1-ymin1)],'k',linewidth=1.5,zorder=5)
ax3.plot([xmin1+0.6*(xmax1-xmin1)+5e3,xmin1+0.6*(xmax1-xmin1)+5e3],[ymin1+0.75*(ymax1-ymin1),ymin1+0.73*(ymax1-ymin1)],'k',linewidth=1.5,zorder=5)
for i in range(0,len(dists_eul)):
  ax3.plot(xflow[ind_eul[i]],yflow[ind_eul[i]],'o',color=coloptions[i],markersize=5,mec='k',mew=0.5)
ax3.text(xmin1+0.64*(xmax1-xmin1)+5e3,ymin1+0.72*(ymax1-ymin1),'5 km',fontsize=8,fontname='arial')
ax3.text(0.05*(xmax1-xmin1)+xmin1,ymin1+0.05*(ymax1-ymin1),'d',fontweight='bold',fontname='arial',fontsize=9)


plt.tight_layout()
plt.subplots_adjust(hspace=0.05,wspace=0.05,right=0.99,left=0.01,bottom=0.02,top=0.98)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_variability.pdf"),FORMAT='PDF',dpi=600)
plt.close()


###############
# Plot trends #
###############

fig = plt.figure(figsize=(3.75,2.4))
matplotlib.rc('font',family='Arial')
gs = matplotlib.gridspec.GridSpec(1,2)
cx = cubehelix.cmap(start=1.2,rot=-1.1,reverse=False,minLight=0.1,sat=2,nlev=8)

plt.subplot(gs[0])
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
ax1 = plt.gca()
plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
norms = matplotlib.colors.BoundaryNorm(np.arange(-500,1,10),cx.N)
p=plt.imshow(veltrend,extent=[np.min(xvel),np.max(xvel),np.min(yvel),np.max(yvel)],origin='lower',clim=[-500,0],cmap=cx,norm=norms)
p.cmap.set_under('k')
p.cmap.set_over('w')
ax1.axes.set_xlim([xmin,xmax])
ax1.axes.set_ylim([ymin,ymax])
ax1.set_xticks([])
ax1.set_yticks([])
xmin1,xmax1 = plt.xlim()
ymin1,ymax1 = plt.ylim()
path = matplotlib.path.Path([[0.46*(xmax1-xmin1)+xmin,0.98*(ymax1-ymin1)+ymin1],
  			[0.98*(xmax1-xmin1)+xmin1,0.98*(ymax1-ymin1)+ymin1],
  			[0.98*(xmax1-xmin1)+xmin1,0.60*(ymax1-ymin1)+ymin1],
  			[0.46*(xmax1-xmin1)+xmin1,0.60*(ymax1-ymin1)+ymin1],
  			[0.46*(xmax1-xmin1)+xmin1,0.98*(ymax1-ymin1)+ymin1]])
patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1,zorder=3)
ax1.add_patch(patch)
cbaxes = fig.add_axes([0.28, 0.87, 0.17, 0.03]) 
cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=[-500,-250,0],extend='both') 
cb.set_label('Velocity trend \n (m yr$^{-2}$)',size=8,fontname='arial')
cb.ax.tick_params(labelsize=8,length=5)
for i in range(0,len(dists_eul)):
  ax1.plot(xflow[ind_eul[i]],yflow[ind_eul[i]],'o',color=coloptions[i],markersize=5,mew=0.5,mec='k')
for i in range(0,len(dists_eul)):
  ax1.plot(xflow[ind_eul[i]],yflow[ind_eul[i]],'o',color=coloptions[i],markersize=5,mew=0.5,mec='k')
  if glacier == 'Kanger':
    if i in [1,2,3,4]:
      ax1.text(xflow[ind_eul[i]]-3.5e3,yflow[ind_eul[i]]+1e3,glacier[0]+'{0:02d}'.format(int(abs(dists_eul[i]))),fontsize=8,fontname='Arial')
    else:
      ax1.text(xflow[ind_eul[i]]-2.5e3,yflow[ind_eul[i]]-2400,glacier[0]+'{0:02d}'.format(int(abs(dists_eul[i]))),fontsize=8,fontname='Arial')
  else:
    if i in [3,4]:
      ax1.text(xflow[ind_eul[i]]-3.2e3,yflow[ind_eul[i]]+1e3,glacier[0]+'{0:02d}'.format(int(abs(dists_eul[i]))),fontsize=8,fontname='Arial')
    else:
      ax1.text(xflow[ind_eul[i]]-3.3e3,yflow[ind_eul[i]]-2.5e3,glacier[0]+'{0:02d}'.format(int(abs(dists_eul[i]))),fontsize=8,fontname='Arial')
ax1.text(0.05*(xmax1-xmin1)+xmin1,0.05*(ymax1-ymin1)+ymin1,'a',weight='bold',fontsize=9)
ax1.plot([xmin1+0.54*(xmax1-xmin1),xmin+0.54*(xmax1-xmin1)+5e3],[ymin1+0.66*(ymax1-ymin1),ymin1+0.66*(ymax1-ymin1)],'k',linewidth=1.5,zorder=5)
ax1.plot([xmin1+0.54*(xmax1-xmin1),xmin1+0.54*(xmax1-xmin1)],[ymin1+0.66*(ymax1-ymin1),ymin1+0.64*(ymax1-ymin1)],'k',linewidth=1.5,zorder=5)
ax1.plot([xmin1+0.54*(xmax1-xmin1)+5e3,xmin1+0.54*(xmax1-xmin1)+5e3],[ymin1+0.66*(ymax1-ymin1),ymin1+0.64*(ymax1-ymin1)],'k',linewidth=1.5,zorder=5)
ax1.text(xmin1+0.58*(xmax1-xmin1)+5e3,ymin1+0.63*(ymax1-ymin1),'5 km',fontsize=8,fontname='arial')


plt.subplot(gs[1])
ax2 = plt.gca()
plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
norms = matplotlib.colors.BoundaryNorm(np.arange(-10,0,1),cx.N)
p=plt.imshow(zstrend,extent=[np.min(xzs),np.max(xzs),np.min(yzs),np.max(yzs)],origin='lower',cmap=cx,clim=[-10,0])
p.cmap.set_under('k')
p.cmap.set_over('w')
ax2.set_xticks([])
ax2.set_yticks([])
ax2.axes.set_xlim([xmin,xmax])
ax2.axes.set_ylim([ymin,ymax])
xmin1,xmax1 = plt.xlim()
ymin1,ymax1 = plt.ylim()
path = matplotlib.path.Path([[0.48*(xmax1-xmin1)+xmin,0.98*(ymax1-ymin1)+ymin1],
  			[0.98*(xmax1-xmin1)+xmin,0.98*(ymax1-ymin1)+ymin1],
  			[0.98*(xmax1-xmin1)+xmin,0.67*(ymax1-ymin1)+ymin1],
  			[0.48*(xmax1-xmin1)+xmin,0.67*(ymax1-ymin1)+ymin1],
  			[0.48*(xmax1-xmin1)+xmin,0.98*(ymax1-ymin1)+ymin1]])
patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1,zorder=3)
ax2.add_patch(patch)
cbaxes = fig.add_axes([0.77, 0.87, 0.17, 0.03]) 
cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=[-10,-5,0],extend='both') 
cb.set_label("Elevation trend \n (m yr$^{-1}$)",size=8,fontname='arial')
cb.ax.tick_params(labelsize=9)
ax2.text(0.05*(xmax1-xmin1)+xmin1,0.05*(ymax1-ymin1)+ymin1,'b',weight='bold',fontsize=9)
#cb.ax.xaxis.set_ticks(minorticks, minor=True)
cb.ax.tick_params(labelsize=8)
cb.ax.tick_params(which='both',length=5)
for i in range(0,len(dists_eul)):
  ax2.plot(xflow[ind_eul[i]],yflow[ind_eul[i]],'o',color=coloptions[i],markersize=5,mew=0.5,mec='k')

plt.tight_layout()
plt.subplots_adjust(hspace=0.05,wspace=0.05,right=0.98,left=0.02,bottom=0.02,top=0.98)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_trends.pdf"),FORMAT='PDF',dpi=600)
plt.close()

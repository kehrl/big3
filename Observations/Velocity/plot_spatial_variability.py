import os
import sys
import numpy as np
import vellib, bedlib, zslib, floatlib, geotifflib, datelib
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
xvel,yvel,velall,velmean,velstd,velcount,veltime = vellib.variability(glacier,time1,time2)

# Load elevation variability
xzs,yzs,zsall,zsmean,zsstd,zscount,zstime = zslib.variability(glacier,time1,time2)

# Find flotation conditions
xwv,ywv,zwv,timewv = zslib.dem_grid(glacier,xmin,xmax,ymin,ymax,years='all',resolution=32,verticaldatum='geoid')
xf,yf,zabovefloat = floatlib.extent(xwv,ywv,zwv,timewv,glacier,rho_i=917.0,rho_sw=1020.0,bedsource='cresis',verticaldatum='geoid')

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

velstd_high = np.array(velstd)
velstd_high[velcount < 0.5*len(veltime)] = float('nan')

fig = plt.figure(figsize=(7.5,2.4))
matplotlib.rc('font',family='Arial')
gs = matplotlib.gridspec.GridSpec(1,4)
cx = cubehelix.cmap(start=1.2,rot=-1.1,reverse=True,minLight=0.1,sat=2)

plt.subplot(gs[0])
ax1 = plt.gca()
plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
p=plt.imshow(zb/1.0e3,extent=[np.min(xb),np.max(xb),np.min(yb),np.max(yb)],origin='lower',cmap='RdBu_r',clim=[-1,1])
plt.contour(velstd_high*2/1.0e3,levels=[0.1],extent=[np.min(xvel),np.max(xvel),np.min(yvel),np.max(yvel)],origin='lower',colors='k',lw=2,zorder=2)
ax1.axes.set_xlim([xmin,xmax])
ax1.axes.set_ylim([ymin,ymax])
ax1.set_xticks([])
ax1.set_yticks([])
xmin1,xmax1 = plt.xlim()
ymin1,ymax1 = plt.ylim()
path = matplotlib.path.Path([[0.49*(xmax1-xmin1)+xmin,0.98*(ymax1-ymin1)+ymin1],
  			[0.98*(xmax1-xmin1)+xmin,0.98*(ymax1-ymin1)+ymin1],
  			[0.98*(xmax1-xmin1)+xmin,0.66*(ymax1-ymin1)+ymin1],
  			[0.49*(xmax1-xmin1)+xmin,0.66*(ymax1-ymin1)+ymin1],
  			[0.49*(xmax1-xmin1)+xmin,0.98*(ymax1-ymin1)+ymin1]])
patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1,zorder=3)
ax1.add_patch(patch)
cbaxes = fig.add_axes([0.15, 0.86, 0.085, 0.02]) 
cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=[-1,0,1]) 
cb.set_label('Bed elevation \n (km asl)',size=10,fontname='arial')
cb.ax.tick_params(labelsize=10)
#ax1.plot([xmin+0.61*(xmax-xmin),xmin+0.61*(xmax-xmin)+5e3],[ymin+0.73*(ymax-ymin),ymin+0.73*(ymax-ymin)],'k',linewidth=1.5)
#ax1.plot([xmin+0.61*(xmax-xmin),xmin+0.61*(xmax-xmin)],[ymin+0.73*(ymax-ymin),ymin+0.71*(ymax-ymin)],'k',linewidth=1.5)
#ax1.plot([xmin+0.61*(xmax-xmin)+5e3,xmin+0.61*(xmax-xmin)+5e3],[ymin+0.73*(ymax-ymin),ymin+0.71*(ymax-ymin)],'k',linewidth=1.5)
#ax1.text(xmin+0.64*(xmax-xmin)+5e3,ymin+0.7*(ymax-ymin),'5 km',fontsize=10)

plt.subplot(gs[1])
ax2 = plt.gca()
ax2.axes.autoscale(False)
plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
norms = matplotlib.colors.BoundaryNorm(np.arange(0,1.1,0.1),cx.N)
p=plt.imshow(velstd_high*2/1.0e3,clim=[0,1],extent=[np.min(xvel),np.max(xvel),np.min(yvel),np.max(yvel)],origin='lower',cmap=cx,norm=norms)
plt.contour(velstd_high*2/1.0e3,levels=[0.1],extent=[np.min(xvel),np.max(xvel),np.min(yvel),np.max(yvel)],origin='lower',colors='k',lw=2,zorder=2)
ax2.set_xticks([])
ax2.set_yticks([])
ax2.axes.set_xlim([xmin,xmax])
ax2.axes.set_ylim([ymin,ymax])
xmin1,xmax1 = plt.xlim()
ymin1,ymax1 = plt.ylim()
path = matplotlib.path.Path([[0.49*(xmax1-xmin1)+xmin,0.98*(ymax1-ymin1)+ymin1],
  			[0.98*(xmax1-xmin1)+xmin,0.98*(ymax1-ymin1)+ymin1],
  			[0.98*(xmax1-xmin1)+xmin,0.66*(ymax1-ymin1)+ymin1],
  			[0.49*(xmax1-xmin1)+xmin,0.66*(ymax1-ymin1)+ymin1],
  			[0.49*(xmax1-xmin1)+xmin,0.98*(ymax1-ymin1)+ymin1]])
patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1,zorder=3)
ax2.add_patch(patch)
cbaxes = fig.add_axes([0.39, 0.86, 0.085, 0.02]) 
cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=[0,0.5,1]) 
cb.set_label(r"Velocity 2-$\mathrm{\sigma}$"+"\n"+"(km yr$^{-1}$)",size=10,fontname='arial')
cb.ax.tick_params(labelsize=10)

plt.subplot(gs[2])
ax3 = plt.gca()
ax3.axes.autoscale(False)
plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
zsstd_high = np.array(zsstd)
zsstd_high[zscount < 1] = float('nan')
norms = matplotlib.colors.BoundaryNorm(np.arange(0,22,2),cx.N)
p=plt.imshow(zsstd_high*2,clim=[0,20],extent=[np.min(xzs),np.max(xzs),np.min(yzs),np.max(yzs)],origin='lower',cmap=cx,norm=norms)
plt.contour(velstd_high*2/1.0e3,levels=[0.1],extent=[np.min(xvel),np.max(xvel),np.min(yvel),np.max(yvel)],origin='lower',colors='k',lw=2,zorder=2)
#plt.contour(zsstd_high*1.96,levels=[5],extent=[np.min(xvel),np.max(xvel),np.min(yvel),np.max(yvel)],origin='lower',colors='k',lw=2)
ax3.set_xticks([])
ax3.set_yticks([])
ax3.axes.set_xlim([xmin,xmax])
ax3.axes.set_ylim([ymin,ymax])
xmin1,xmax1 = plt.xlim()
ymin1,ymax1 = plt.ylim()
path = matplotlib.path.Path([[0.49*(xmax1-xmin1)+xmin,0.98*(ymax1-ymin1)+ymin1],
  			[0.98*(xmax1-xmin1)+xmin,0.98*(ymax1-ymin1)+ymin1],
  			[0.98*(xmax1-xmin1)+xmin,0.66*(ymax1-ymin1)+ymin1],
  			[0.49*(xmax1-xmin1)+xmin,0.66*(ymax1-ymin1)+ymin1],
  			[0.49*(xmax1-xmin1)+xmin,0.98*(ymax1-ymin1)+ymin1]])
patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1,zorder=3)
ax3.add_patch(patch)
cbaxes = fig.add_axes([0.625, 0.86, 0.09, 0.02]) 
cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=[0,10,20]) 
cb.set_label(r'Elevation 2-$\sigma$'+'\n (m)',size=10,fontname='arial')
cb.ax.tick_params(labelsize=10)


plt.subplot(gs[3])
ax3 = plt.gca()
ax3.axes.autoscale(False)
plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
plt.plot(0,0,'r.',label='Grounded',markersize=5)
plt.plot(0,0,'b.',label='Floating',markersize=5)
plt.plot(0,0,'k.',label='Changing',markersize=5)
ind = np.where(floatcond == 1)
plt.plot(xf[ind],yf[ind],'r.',lw=0,markersize=2)
ind = np.where(floatcond == -1)
plt.plot(xf[ind],yf[ind],'b.',lw=0,markersize=2)
ind = np.where(floatcond == 0)
plt.plot(xf[ind],yf[ind],'k.',lw=0,markersize=2)
#p=plt.imshow(zsstd_high*2,clim=[0,20],extent=[np.min(xzs),np.max(xzs),np.min(yzs),np.max(yzs)],origin='lower',cmap='jet')
#plt.contour(zsstd_high*1.96,levels=[5],extent=[np.min(xvel),np.max(xvel),np.min(yvel),np.max(yvel)],origin='lower',colors='k',lw=2)
ax3.set_xticks([])
ax3.set_yticks([])
ax3.axes.set_xlim([xmin,xmax])
ax3.axes.set_ylim([ymin,ymax])
xmin1,xmax1 = plt.xlim()
ymin1,ymax1 = plt.ylim()
path = matplotlib.path.Path([[0.54*(xmax1-xmin1)+xmin,0.98*(ymax1-ymin1)+ymin1],
  			[0.98*(xmax1-xmin1)+xmin,0.98*(ymax1-ymin1)+ymin1],
  			[0.98*(xmax1-xmin1)+xmin,0.66*(ymax1-ymin1)+ymin1],
  			[0.54*(xmax1-xmin1)+xmin,0.66*(ymax1-ymin1)+ymin1],
  			[0.54*(xmax1-xmin1)+xmin,0.98*(ymax1-ymin1)+ymin1]])
patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1,zorder=3)
ax3.add_patch(patch)
plt.legend(loc=1,numpoints=1,frameon=False,labelspacing=0.07,handletextpad=0.5,handlelength=0.2,fontsize=10)
ax3.plot([xmin1+0.6*(xmax1-xmin1),xmin+0.6*(xmax1-xmin1)+5e3],[ymin1+0.73*(ymax1-ymin1),ymin1+0.73*(ymax1-ymin1)],'k',linewidth=1.5,zorder=5)
ax3.plot([xmin1+0.6*(xmax1-xmin1),xmin1+0.6*(xmax1-xmin1)],[ymin1+0.73*(ymax1-ymin1),ymin1+0.71*(ymax1-ymin1)],'k',linewidth=1.5,zorder=5)
ax3.plot([xmin1+0.6*(xmax1-xmin1)+5e3,xmin1+0.6*(xmax1-xmin1)+5e3],[ymin1+0.73*(ymax1-ymin1),ymin1+0.71*(ymax1-ymin1)],'k',linewidth=1.5,zorder=5)
ax3.text(xmin1+0.64*(xmax1-xmin1)+5e3,ymin1+0.7*(ymax1-ymin1),'5 km',fontsize=10,fontname='arial')

plt.tight_layout()
plt.subplots_adjust(hspace=0.05,wspace=0.05)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_variability.pdf"),FORMAT='PDF',dpi=400)
plt.close()

# Plot trends
# xzs,yzs,zsall,zsmean,zsstd,zscount,zstime = zslib.variability(glacier,time1,time2,data='TDM')

zstrend = np.zeros_like(zsstd)
zstrend_time1 = np.zeros_like(zsstd)
zstrend_time2 = np.zeros_like(zsstd)
zstrend_p = np.zeros_like(zsstd)
zstrend_count = np.zeros_like(zsstd)
zstrend_error = np.zeros_like(zsstd)
zstrend_r = np.zeros_like(zsstd)
zstrend[:,:] = float('nan')
zstrend_p[:,:] = float('nan')
zstrend_error[:,:] = float('nan')
zstrend_r[:,:] = float('nan')
ind = np.where((zstime > time1) & (zstime < time2))[0]
for j in range(0,len(yzs)):
  for i in range(0,len(xzs)):
    nonnan = np.where(~(np.isnan(zsall[j,i,ind])))[0]
    zstrend_count = len(nonnan)
    if len(nonnan) > 1:
      if (np.floor(np.min(zstime[ind[nonnan]]))==time1) and np.ceil(np.max(zstime[ind[nonnan]]))==time2:
        zstrend_time1[j,i] = np.min(zstime[ind[nonnan]])
        zstrend_time2[j,i] = np.max(zstime[ind[nonnan]])
        slope,intercept,r,p,std_err = stats.linregress(zstime[ind[nonnan]],zsall[j,i,ind[nonnan]])
        zstrend[j,i] = slope
        zstrend_error[j,i] = std_err
        zstrend_p[j,i] = p
        zstrend_r[j,i] = r
ind = np.where(zstrend_p > 0.05)
zstrend[ind] = float('nan')

veltrend = np.zeros_like(velstd)
veltrend_time1 = np.zeros_like(velstd)
veltrend_time2 = np.zeros_like(velstd)
veltrend_count = np.zeros_like(velstd)
veltrend_p = np.zeros_like(velstd)
veltrend_error = np.zeros_like(velstd)
veltrend_r = np.zeros_like(velstd)
veltrend_p[:,:] = float('nan')
veltrend[:,:] = float('nan')
veltrend_error[:,:] = float('nan')
veltrend_r[:,:] = float('nan')
ind = np.where((veltime > time1) & (veltime < time2))[0]
for j in range(0,len(yvel)):
  for i in range(0,len(xvel)):
    nonnan = np.where((~(np.isnan(velall[j,i,ind]))))[0]
    if len(nonnan) > 0.75*len(ind):
      if (np.floor(np.min(veltime[ind[nonnan]]))==time1) and np.ceil(np.max(veltime[ind[nonnan]]))==time2:
        slope,intercept,r,p,std_err = stats.linregress(veltime[ind[nonnan]],velall[j,i,ind[nonnan]])
        veltrend_count[j,i] = len(nonnan)
        veltrend[j,i] = slope
        veltrend_p[j,i] = p
        veltrend_error[j,i] = std_err
        veltrend_time1[j,i] = np.min(veltime[ind[nonnan]])
        veltrend_time2[j,i] = np.max(veltime[ind[nonnan]])
        veltrend_r[j,i] = r
ind = np.where(veltrend_p > 0.05)
veltrend[ind] = float('nan')

fig = plt.figure(figsize=(3.75,2.4))
matplotlib.rc('font',family='Arial')
gs = matplotlib.gridspec.GridSpec(1,2)
cx = cubehelix.cmap(start=1.2,rot=-1.1,reverse=False,minLight=0.1,sat=2,nlev=8)

plt.subplot(gs[0])
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
ax1 = plt.gca()
plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
norms = matplotlib.colors.BoundaryNorm(np.arange(-500,1,10),cx.N)
p=plt.imshow(veltrend,extent=[np.min(xvel),np.max(xvel),np.min(yvel),np.max(yvel)],origin='lower',clim=[-500,0],cmap=cx)
plt.contour(veltrend,0,extent=[np.min(xvel),np.max(xvel),np.min(yvel),np.max(yvel)],origin='lower',colors='k')
ax1.axes.set_xlim([xmin,xmax])
ax1.axes.set_ylim([ymin,ymax])
ax1.set_xticks([])
ax1.set_yticks([])
xmin1,xmax1 = plt.xlim()
ymin1,ymax1 = plt.ylim()
path = matplotlib.path.Path([[0.46*(xmax1-xmin1)+xmin,0.98*(ymax1-ymin1)+ymin1],
  			[0.98*(xmax1-xmin1)+xmin1,0.98*(ymax1-ymin1)+ymin1],
  			[0.98*(xmax1-xmin1)+xmin1,0.64*(ymax1-ymin1)+ymin1],
  			[0.46*(xmax1-xmin1)+xmin1,0.64*(ymax1-ymin1)+ymin1],
  			[0.46*(xmax1-xmin1)+xmin1,0.98*(ymax1-ymin1)+ymin1]])
patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1,zorder=3)
ax1.add_patch(patch)
cbaxes = fig.add_axes([0.3, 0.84, 0.17, 0.03]) 
cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=[-500,-250,0],extend='max') 
cb.set_label('du/dt \n (m yr$^{-2}$)',size=10,fontname='arial')
cb.ax.tick_params(labelsize=8.5,length=5)
ax1.text(0.05*(xmax1-xmin1)+xmin1,0.9*(ymax1-ymin1)+ymin1,'a',weight='bold',fontsize=10)

plt.subplot(gs[1])
ax2 = plt.gca()
plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
p=plt.imshow(zstrend,extent=[np.min(xzs),np.max(xzs),np.min(yzs),np.max(yzs)],origin='lower',cmap='RdBu_r',clim=[-10,10])
ax2.set_xticks([])
ax2.set_yticks([])
ax2.axes.set_xlim([xmin,xmax])
ax2.axes.set_ylim([ymin,ymax])
xmin1,xmax1 = plt.xlim()
ymin1,ymax1 = plt.ylim()
path = matplotlib.path.Path([[0.48*(xmax1-xmin1)+xmin,0.98*(ymax1-ymin1)+ymin1],
  			[0.98*(xmax1-xmin1)+xmin,0.98*(ymax1-ymin1)+ymin1],
  			[0.98*(xmax1-xmin1)+xmin,0.64*(ymax1-ymin1)+ymin1],
  			[0.48*(xmax1-xmin1)+xmin,0.64*(ymax1-ymin1)+ymin1],
  			[0.48*(xmax1-xmin1)+xmin,0.98*(ymax1-ymin1)+ymin1]])
patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1,zorder=3)
ax2.add_patch(patch)
cbaxes = fig.add_axes([0.755, 0.84, 0.16, 0.03]) 
cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=[-10,0,10]) 
cb.set_label("dH/dt \n (m yr$^{-1}$)",size=10,fontname='arial')
cb.ax.tick_params(labelsize=9)
ax2.text(0.05*(xmax1-xmin1)+xmin1,0.9*(ymax1-ymin1)+ymin1,'b',weight='bold',fontsize=10)
minorticks = p.norm(np.arange(-10, 15, 5))
cb.ax.xaxis.set_ticks(minorticks, minor=True)
cb.ax.tick_params(labelsize=9)
cb.ax.tick_params(which='both',length=5)
plt.tight_layout()
plt.subplots_adjust(hspace=0.03,wspace=0.03)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_trends.pdf"),FORMAT='PDF',dpi=300)
plt.close()

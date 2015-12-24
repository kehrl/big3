import os
import sys
import numpy as np
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"Util/Modules"))
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"BigThreeGlaciers/Tools"))
import elevation, geotiff, fracyear, icemask,flotation
import matplotlib.pyplot as plt
import matplotlib

# Get arguments
args = sys.argv
glacier = args[1][:] # Options: Kanger, Helheim


##########################
# Get surface elevations #
##########################

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

# DEMs
xdem,ydem,zdem,timedem,errordem = elevation.dem_grid(glacier,xmin,xmax,ymin,ymax,years='all',verticaldatum='geoid',return_error=True)

# Mask
xmask,ymask,zmask=icemask.load_grid(glacier,np.min(xdem),np.max(xdem),np.min(ydem),np.max(ydem),32)

######################
# Get Landsat images #
######################

# Image for plotting
if glacier == "Helheim":
  imagetime = fracyear.date_to_fracyear(2014,7,4)
  ximage,yimage,image = geotiff.readrgb(os.path.join(os.getenv("DATA_HOME"),"Imagery/Landsat/Helheim/TIF/20140704140535_LC82330132014185LGN00.tif"))
elif glacier == "Kanger":
  imagetime = fracyear.date_to_fracyear(2014,7,6)
  ximage,yimage,image = geotiff.readrgb(os.path.join(os.getenv("DATA_HOME"),"Imagery/Landsat/Kanger/TIF/20140706135251_LC82310122014187LGN00.tif"))

###############################
# Choose DEMs for subtraction #
###############################

if glacier == 'Helheim':
  time1 = [2013,2,9]
  ind1 = np.argmin(abs(fracyear.date_to_fracyear(time1[0],time1[1],time1[2])-timedem))
  time2 = [2013,10,31]
  ind2 = np.argmin(abs(fracyear.date_to_fracyear(time2[0],time2[1],time2[2])-timedem))
  time3 = [2012,3,16]
  ind3 = np.argmin(abs(fracyear.date_to_fracyear(time3[0],time3[1],time3[2])-timedem))
  time4 = [2012,12,5]
  ind4 = np.argmin(abs(fracyear.date_to_fracyear(time4[0],time4[1],time4[2])-timedem))
elif glacier == 'Kanger':
  time1 = [2012,5,22]
  time2 = [2012,12,17]
  time3 = [2012,5,22]
  time4 = [2012,12,17]
  ind1 = np.argmin(abs(fracyear.date_to_fracyear(time1[0],time1[1],time1[2])-timedem))
  ind2 = np.argmin(abs(fracyear.date_to_fracyear(time2[0],time2[1],time2[2])-timedem))
  ind3 = np.argmin(abs(fracyear.date_to_fracyear(time3[0],time3[1],time3[2])-timedem))
  ind4 = np.argmin(abs(fracyear.date_to_fracyear(time4[0],time4[1],time4[2])-timedem))

########################################################
# Find out if flotation condition changed between DEMs #
########################################################

xf,yf,zabovefloat = flotation.extent(xdem,ydem,zdem[:,:,[ind1,ind2]],timedem[[ind1,ind2]],glacier,rho_i=917.0,rho_sw=1020.0,bedsource='cresis',verticaldatum='geoid')

cutoff = 0.
floatcond = np.zeros(len(xf))
floatcond[:] = float('nan')

for i in range(0,len(xf)):
  if np.min(zabovefloat[i,:]) > cutoff:
    floatcond[i] = 1 #grounded
  elif np.max(zabovefloat[i,:]) < -1*cutoff:
    floatcond[i] = -1 #floating
  elif (np.max(zabovefloat[i,:]) > -1*cutoff) or (np.min(zabovefloat[i,:]) < cutoff):
    floatcond[i] = 0 #changing basal conditions



fig = plt.figure(figsize=(3.75,2.4))
matplotlib.rc('font',family='Arial')
gs = matplotlib.gridspec.GridSpec(1,2)
cx = matplotlib.cm.get_cmap('RdBu_r',8)

plt.subplot(gs[0])
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
ax1 = plt.gca()
plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
norms = matplotlib.colors.BoundaryNorm(np.arange(-20,20,5),cx.N)
p=plt.imshow(np.ma.masked_array(zdem[:,:,ind2]-zdem[:,:,ind1],zmask),extent=[np.min(xdem),np.max(xdem),np.min(ydem),np.max(ydem)],origin='lower',clim=[-20,20],cmap='RdBu_r')
#plt.contour(np.ma.masked_array(zdem[:,:,ind2]-zdem[:,:,ind1],zmask),0,colors='k',extent=[np.min(xdem),np.max(xdem),np.min(ydem),np.max(ydem)],origin='lower',fontsize=8)
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
cbaxes = fig.add_axes([0.3, 0.84, 0.16, 0.03]) 
cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=[-20,0,20]) 
cb.set_label('Difference \n (m)',size=10,fontname='arial')
minorticks = p.norm(np.arange(-20, 25, 10))
cb.ax.xaxis.set_ticks(minorticks, minor=True)
cb.ax.tick_params(labelsize=9)
cb.ax.tick_params(which='both',length=5)
ax1.text(0.05*(xmax1-xmin1)+xmin1,0.9*(ymax1-ymin1)+ymin1,'a',weight='bold',fontsize=10)

plt.subplot(gs[1])
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
ax1 = plt.gca()
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
plt.legend(loc=1,numpoints=1,frameon=False,labelspacing=0.07,handletextpad=0.5,handlelength=0.2,fontsize=10,borderpad=0)
ax1.plot([xmin1+0.56*(xmax1-xmin1),xmin+0.56*(xmax1-xmin1)+5e3],[ymin1+0.71*(ymax1-ymin1),ymin1+0.71*(ymax1-ymin1)],'k',linewidth=1.5,zorder=5)
ax1.plot([xmin1+0.56*(xmax1-xmin1),xmin1+0.56*(xmax1-xmin1)],[ymin1+0.71*(ymax1-ymin1),ymin1+0.69*(ymax1-ymin1)],'k',linewidth=1.5,zorder=5)
ax1.plot([xmin1+0.56*(xmax1-xmin1)+5e3,xmin1+0.56*(xmax1-xmin1)+5e3],[ymin1+0.71*(ymax1-ymin1),ymin1+0.69*(ymax1-ymin1)],'k',linewidth=1.5,zorder=5)
ax1.text(xmin1+0.60*(xmax1-xmin1)+5e3,ymin1+0.68*(ymax1-ymin1),'5 km',fontsize=10,fontname='arial')
ax1.text(0.05*(xmax1-xmin1)+xmin1,0.9*(ymax1-ymin1)+ymin1,'b',weight='bold',fontsize=10)

plt.tight_layout()
plt.subplots_adjust(hspace=0.03,wspace=0.03)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_seasonal.pdf"),FORMAT='PDF',dpi=400)
plt.close()

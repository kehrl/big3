# Calculate taud for Helheim for various time periods. 

import zslib, bedlib, vellib, datelib, geotifflib, glaclib, masklib
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
import scipy
import matplotlib
#from scipy.spatial import cKDTree

glacier = 'Kanger'

# Parameters
rho_i = 917.0
g = 9.8

if glacier == 'Kanger':
  xmin = 468000.
  xmax = 500800.
  ymin = -2300000.
  ymax = -2260200.
  ximage,yimage,image = geotifflib.readrgb(os.path.join(os.getenv("DATA_HOME"),"Imagery/Landsat/Helheim/TIF/20140704140535_LC82330132014185LGN00.tif"))
elif glacier == 'Helheim':
  xmin = 280200.
  xmax = 313000.
  ymin = -2585100.
  ymax = -2545300.
  ximage,yimage,image = geotifflib.readrgb(os.path.join(os.getenv("DATA_HOME"),"Imagery/Landsat/Kanger/TIF/20140706135251_LC82310122014187LGN00.tif"))

dx = 200.
x = np.arange(xmin-1e3,xmax+1e3,dx)
y = np.arange(ymin-1e3,ymax+1e3,dx)

xdem,ydem,zdem,timedem,errordem = zslib.dem_grid(glacier,xmin-2e3,xmax+2e3,ymin-2e3,ymax+2e3,years='all',verticaldatum='ellipsoid',return_error=True,data='TDM')
xbed,ybed,zbed = bedlib.morlighem_grid(xmin-2e3,xmax+2e3,ymin-2e3,ymax+2e3,verticaldatum='ellipsoid')

# Smooth surface DEMs
#zdem_blur = np.zeros_like(zdem)
#for k in range(0,len(timedem)):
#  zdem_blur[:,:,k] = scipy.ndimage.filters.gaussian_filter(zdem[:,:,k],sigma=2,truncate=4)
#  # "Blurring DEM over 17 pixels (roughly 500 m in each direction)..."

# Interpolate bed onto grid
x_grid,y_grid = np.meshgrid(x,y)
x_flat = x_grid.flatten()
y_flat = y_grid.flatten()
f = RegularGridInterpolator((ybed,xbed),zbed)
zbed_int_flat = f((y_flat,x_flat))
zbed_int = np.reshape(zbed_int_flat,[len(y),len(x)])

# Interpolate zs onto grid
zs = np.zeros([len(y),len(x),len(timedem)])
#for k in range(0,len(timedem)):
#  f = RegularGridInterpolator((ydem,xdem),zdem_blur[:,:,k])
#  zs_flat = f((y_flat,x_flat))
#  zs[:,:,k] = np.reshape(zs_flat,[len(y),len(x)])
xdem_grid,ydem_grid = np.meshgrid(xdem,ydem)
for j in range(0,len(y)):
  for i in range(0,len(x)):
    ind = np.where(np.sqrt((x_grid[j,i]-xdem_grid)**2+(y_grid[j,i]-ydem_grid)**2)<dx/2.)
    zs[j,i,:] = np.nanmean(zdem[ind],axis=0)

# Calculate grad(h) & H
dhdx = np.zeros([len(y),len(x),len(timedem)])
dhdy = np.zeros([len(y),len(x),len(timedem)])
H = np.zeros([len(y),len(x),len(timedem)])
taud = np.zeros([len(y),len(x),len(timedem)])
taud_u = np.zeros([len(y),len(x),len(timedem)])
H[:,:,:] = float('nan')
taud[:,:,:] = float('nan')
taud_u[:,:,] = float('nan')
dhdx[:,:,:] = float('nan')
dhdy[:,:,:] = float('nan')
u,v = vellib.inversion_3D(glacier,x,y,timedem[-1])
for k in range(0,len(timedem)):
  xmask,ymask,zmask = masklib.load_grid(glacier,x[0],x[-1],y[0],y[-1],y[1]-y[0],icefront_time=timedem[k],ice=1)

  for i in range(1,len(x)-1):
    for j in range(1,len(y)-1):
      if (zmask[j,i] == 1) and (zs[j,i,k] > zbed_int[j,i]):
        dhdx[j,i,k] = (zs[j,i+1,k]-zs[j,i-1,k])/(x[i+1]-x[i-1]) 
        dhdy[j,i,k] = (zs[j+1,i,k]-zs[j-1,i,k])/(y[j+1]-y[j-1]) 
        H[j,i,k] = zs[j,i,k]-zbed_int[j,i]
        taud[j,i,k] = rho_i*g*H[j,i,k]*np.sqrt((dhdx[j,i,k])**2+(dhdy[j,i,k])**2)
        taud_u[j,i,k] = -rho_i*g*H[j,i,k]*(dhdx[j,i,k]*u[j,i]+dhdy[j,i,k]*v[j,i])/np.sqrt(u[j,i]**2+v[j,i]**2)

mean_taud_u = np.nanmean(taud_u,axis=2)

for k in range(0,len(timedem)):
  fig = plt.figure(figsize=(3.4,4.0))
  ax = plt.gca()
  plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
  p = plt.imshow((taud_u[:,:,k])/1e3,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',vmin=-750,vmax=750,cmap='RdBu_r')
  plt.xticks([])
  plt.yticks([])
  ax.axis('equal')
  plt.xlim([xmin,xmax])
  plt.ylim([ymin,ymax])

  xmin,xmax = plt.xlim()
  ymin,ymax = plt.ylim()

  path = matplotlib.path.Path([[0.58*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
  			[0.98*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
  			[0.98*(xmax-xmin)+xmin,0.76*(ymax-ymin)+ymin],
  			[0.58*(xmax-xmin)+xmin,0.76*(ymax-ymin)+ymin],
  			[0.58*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin]])
  patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1)
  ax.add_patch(patch)
  cbaxes = fig.add_axes([0.61, 0.935, 0.28, 0.02]) 
  cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=np.arange(-500,1000,500)) 
  ax.text(xmin+0.60*(xmax-xmin),ymin+0.86*(ymax-ymin),r'Along-flow driving',fontsize=9)
  ax.text(xmin+0.68*(xmax-xmin),ymin+0.82*(ymax-ymin),'stress (kPa)',fontsize=9)
  #if label != 'none':
  #  ax.text(xmin+0.02*(xmax-xmin),ymin+0.93*(ymax-ymin),label,fontweight='bold',fontsize=10)
  cb.ax.tick_params(labelsize=10)
  ax.plot([xmin+0.61*(xmax-xmin),xmin+0.61*(xmax-xmin)+5e3],[ymin+0.80*(ymax-ymin),ymin+0.80*(ymax-ymin)],'k',linewidth=1.5)
  ax.plot([xmin+0.61*(xmax-xmin),xmin+0.61*(xmax-xmin)],[ymin+0.80*(ymax-ymin),ymin+0.78*(ymax-ymin)],'k',linewidth=1.5)
  ax.plot([xmin+0.61*(xmax-xmin)+5e3,xmin+0.61*(xmax-xmin)+5e3],[ymin+0.80*(ymax-ymin),ymin+0.78*(ymax-ymin)],'k',linewidth=1.5)
  ax.text(xmin+0.64*(xmax-xmin)+5e3,ymin+0.78*(ymax-ymin),'5 km',fontsize=10)

  year,month,day = datelib.fracyear_to_date(timedem[k])
  date = '{0}{1:02d}{2:02d}'.format(year,month,int(day))
  plt.tight_layout()
  plt.subplots_adjust(top=0.99,right=0.99,left=0.01,bottom=0.01) 
  plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_taud_"+date+"_flow.pdf"),format='PDF',transparent=True)
  plt.close()
  
  fig = plt.figure(figsize=(3.4,4.0))
  ax = plt.gca()
  plt.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
  cx = cubehelix.cmap(start=1.2,rot=-1.1,reverse=True,minLight=0.1,sat=2)
  p = plt.imshow((taud[:,:,k])/1e3,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',vmin=0,vmax=900,cmap=cx)
  plt.xticks([])
  plt.yticks([])
  ax.axis('equal')
  plt.xlim([xmin,xmax])
  plt.ylim([ymin,ymax])

  xmin,xmax = plt.xlim()
  ymin,ymax = plt.ylim()

  path = matplotlib.path.Path([[0.58*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
  			[0.98*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
  			[0.98*(xmax-xmin)+xmin,0.76*(ymax-ymin)+ymin],
  			[0.58*(xmax-xmin)+xmin,0.76*(ymax-ymin)+ymin],
  			[0.58*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin]])
  patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1)
  ax.add_patch(patch)
  cbaxes = fig.add_axes([0.61, 0.935, 0.28, 0.02]) 
  cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=np.arange(0,1200,400)) 
  ax.text(xmin+0.61*(xmax-xmin),ymin+0.86*(ymax-ymin),'Driving stress',fontsize=10)
  ax.text(xmin+0.72*(xmax-xmin),ymin+0.82*(ymax-ymin),'(kPa)',fontsize=10)
  #if label != 'none':
  #  ax.text(xmin+0.02*(xmax-xmin),ymin+0.93*(ymax-ymin),label,fontweight='bold',fontsize=10)
  cb.ax.tick_params(labelsize=10)
  ax.plot([xmin+0.61*(xmax-xmin),xmin+0.61*(xmax-xmin)+5e3],[ymin+0.80*(ymax-ymin),ymin+0.80*(ymax-ymin)],'k',linewidth=1.5)
  ax.plot([xmin+0.61*(xmax-xmin),xmin+0.61*(xmax-xmin)],[ymin+0.80*(ymax-ymin),ymin+0.78*(ymax-ymin)],'k',linewidth=1.5)
  ax.plot([xmin+0.61*(xmax-xmin)+5e3,xmin+0.61*(xmax-xmin)+5e3],[ymin+0.80*(ymax-ymin),ymin+0.78*(ymax-ymin)],'k',linewidth=1.5)
  ax.text(xmin+0.64*(xmax-xmin)+5e3,ymin+0.78*(ymax-ymin),'5 km',fontsize=10)

  year,month,day = datelib.fracyear_to_date(timedem[k])
  date = '{0}{1:02d}{2:02d}'.format(year,month,int(day))
  plt.tight_layout()
  plt.subplots_adjust(top=0.99,right=0.99,left=0.01,bottom=0.01) 
  plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_taud_"+date+"_mag.pdf"),format='PDF',transparent=True)
  plt.close()


xf,yf,zbf,dists = glaclib.load_flowline(glacier)

taud_u_flow = np.zeros([len(dists),len(timedem)])
for k in range(0,len(timedem)):
  f = RegularGridInterpolator((y,x),taud_u[:,:,k]-mean_taud_u,bounds_error=False,fill_value=float('nan'))
  taud_u_flow[:,k] = f((yf,xf))
  
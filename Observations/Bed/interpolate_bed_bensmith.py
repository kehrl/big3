import os
import numpy as np
import bedlib, crossoverlib, masklib
import matplotlib
import matplotlib.pyplot as plt

glacier = 'Helheim'

# Extent for interpolation
extent = np.loadtxt(os.path.join(os.getenv("HOME"),"Bigtmp/bed_bensmith_20150915/helheim_extent.dat"))

xmin = np.min(extent[:,0])-10.0e3
xmax = np.max(extent[:,0])+10.0e3
ymin = np.min(extent[:,1])-10.0e3
ymax = np.max(extent[:,1])+10.0e3

# Load velocities
x,y,v,vx,vy,ex,ey,time,int = geodat.readbinary(os.path.join(os.getenv("DATA_HOME"),"Velocity/Random/Greenland/AllGLVel/mosaicOffsets"))
xind1 = np.argmin(abs(x-xmin))
xind2 = np.argmin(abs(x-xmax))
yind1 = np.argmin(abs(y-ymin))
yind2 = np.argmin(abs(y-ymax))

x = x[xind1:xind2]
y = y[yind1:yind2]
v = v[yind1:yind2,xind1:xind2]
vx = vx[yind1:yind2,xind1:xind2]
vy = vy[yind1:yind2,xind1:xind2]

# Rescale vx,vy
vx_norm = np.array(vx/v)
vy_norm = np.array(vy/v)

# Load mask
xmask,ymask,icemask = masklib.load(glacier,np.min(x),np.max(x),np.min(y),np.max(y),500.0,ice=0)

# Fill in unknown azimuths
xgrid,ygrid = np.meshgrid(x,y)
indtodelete = np.where((xgrid < 325500.0) & (xgrid > 311500.0) & (ygrid > -2583050.0) & (ygrid < -2572930.0) & (v < 2000.0))
vx_norm[indtodelete] = float('NaN')
vy_norm[indtodelete] = float('NaN')
indtofill = np.where(np.isnan(vx_norm))
ind = np.where(~(np.isnan(vx_norm)) & (icemask==0))
xgood = xgrid[ind].ravel()
ygood = ygrid[ind].ravel()
vxgood = vx_norm[ind].ravel()
vygood = vy_norm[ind].ravel()
vx_norm[indtofill]=scipy.interpolate.griddata((xgood,ygood),vxgood,(x[indtofill[1]],y[indtofill[0]]),method='nearest')
vy_norm[indtofill]=scipy.interpolate.griddata((xgood,ygood),vygood,(x[indtofill[1]],y[indtofill[0]]),method='nearest')

vxmask = np.array(vx_norm)
vxmask[icemask==1] = 0.0
vymask = np.array(vy_norm)
vymask[icemask==1] = 0.0

# Load bed measurements
cresis_2001 = bedlib.cresis('2001',glacier,verticaldatum='ellipsoid')
cresis_composite = bedlib.cresis('all',glacier,verticaldatum='ellipsoid')

cresis = np.row_stack([cresis_2001,cresis_composite])

xcrossnew,ycrossnew,zcross1new,zcross2new,zcrossdiffnew = crossoverlib.find(cresis[:,0],cresis[:,1],cresis[:,2],maxdist=100)
xcrossnew,ycrossnew,tcross1new,tcross2new,tcrossdiffnew = crossoverlib.find(cresis[:,0],cresis[:,1],cresis[:,3],maxdist=100)


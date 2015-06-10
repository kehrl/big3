# Beginning of looking at surface mass balance for the models...

import os
import sys
import numpy as np
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Modules"))
import netCDF4
import coords

##########
# Inputs #
##########

# RACMO data
data = netCDF4.Dataset(os.path.join(os.getenv("HOME"),"Data/Climate/SurfaceBalance/Greenland/RACMO/daily_smb_2006-10.nc"))

# Glacier extent
DIRM = os.path.join(os.getenv("HOME"),"Models/Helheim/Meshes/3D/InverseClass/Inputs/")
extent = np.loadtxt(DIRM+"mesh_extent.dat")

#############
# Load data #
#############

# Load RACMO data
lat = np.array(data.variables['lat'][:])
lon = np.array(data.variables['lon'][:])
smb = np.array(data.variables['smb'][:])
time = np.array(data.variables['time'][:])
avesmb = np.mean(smb,0)

# Convert lat,lon to epsg 3413
xrac,yrac = coords.convert(lon,lat,4326,3413)

# Find domain of interest
xmin = np.min(extent[:,0]-1e4)
xmax = np.max(extent[:,0]+1e4)
ymin = np.min(extent[:,1]-1e4)
ymax = np.max(extent[:,1]+1e4)

xind = np.where((xrac >= xmin) & (xrac <= xmax))
xrac = xrac[xind]
yrac = yrac[xind]
avesmb = avesmb[xind]

yind = np.where((yrac >= ymin) & (yrac <= ymax))
xrac = xrac[yind]
yrac = yrac[yind]
avesmb = avesmb[yind]
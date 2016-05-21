# I'd like to look at the spatial distribution of strain rates through time at Helheim 
# and Kanger, so this is the beginning of a hack to do just that. 

import vellib
import numpy as np
import os
import geotifflib

glacier = 'Helheim'

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

file = os.path.join(os.getenv("DATA_HOME"),"Velocity/TSX/"+glacier+"/TIF/all-2008-2016")
x,y,u = geotifflib.read(file+'_vx.tif')
x,y,v = geotifflib.read(file+'_vy.tif')  

nx = len(x)
ny = len(y)

dudx = np.zeros_like(v)
dudy = np.zeros_like(v)
dvdx = np.zeros_like(v)
dvdy = np.zeros_like(v)

dudx[:,:] = float('nan')
dudy[:,:]  = float('nan')
dvdx[:,:]  = float('nan')
dvdy[:,:]  = float('nan')

for i in range(1,nx-1):
  for j in range(1,ny-1):
    dudx[j,i] = (u[j,i+1]-u[j,i-1])/(x[i+1]-x[i-1])
    dudy[j,i] = (u[j+1,i]-u[j-1,i])/(y[j+1]-y[j-1])
    dvdx[j,i] = (v[j,i+1]-v[j,i-1])/(x[i+1]-x[i-1])
    dvdy[j,i] = (v[j+1,i]-v[j-1,i])/(y[j+1]-y[j-1])


theta = np.arctan(v,u) # Flow direction
#theta = (1./2.)*np.arctan2((dudy+dvdx),(dudx-dvdy)) # Principal axis
e_xx = dudx*(np.cos(theta)**2)+dvdy*(np.sin(theta)**2)+(dvdx+dudy)*np.sin(theta)*np.cos(theta)
e_yy = dudx*(np.sin(theta)**2)+dvdy*(np.cos(theta)**2)-(dvdx+dudy)*np.sin(theta)*np.cos(theta)
#theta = (1./2.)*np.arctan2((dudy+dvdx),(dudx+dvdy))
#e_xx = dudx*(np.cos(theta)**2)+dvdy*(np.sin(theta)**2)

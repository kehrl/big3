# Look at spatial distribution of stresses from an inversion. Working on this for Kristin.

import elmerreadlib
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

DIRM = os.path.join(os.getenv("MODEL_HOME"),"Helheim/3D/DEM20120316/")
DIR = DIRM+"inversion_robin/lambda_1e10_20160525"
inputs = DIRM+"inputs/"
runname = "robin_beta"
bsurf = 5

# Get glacier mesh
extent = np.loadtxt(inputs+"mesh_extent.dat")
try:
  hole1 = np.loadtxt(inputs+"mesh_hole1.dat")
  hole2 = np.loadtxt(inputs+"mesh_hole2.dat")
  holes=[hole1,hole2]  
except:
  holes=[]

surf = elmerreadlib.saveline_boundary(DIR,runname,bsurf) 
x,y,sigxx = elmerreadlib.grid3d(surf,'stress1',holes,extent) 
x,y,sigyy = elmerreadlib.grid3d(surf,'stress2',holes,extent) 
x,y,sigxy = elmerreadlib.grid3d(surf,'stress4',holes,extent) 
x,y,uu = elmerreadlib.grid3d(surf,'vel1',holes,extent)
x,y,vv = elmerreadlib.grid3d(surf,'vel2',holes,extent)

# Flow direction
theta = np.zeros_like(uu)
theta[:,:] = float('nan')
sigflow = np.zeros_like(uu)
sigflow[:,:] = float('nan')

for i in range(0,len(x)):
  for j in range(0,len(y)):
    if uu[j,i] is not(np.ma.masked):
      theta[j,i] = np.arctan2(vv[j,i],uu[j,i])

      sigflow[j,i] = sigxx[j,i]*(np.cos(theta[j,i])**2)+sigyy[j,i]*(np.sin(theta[j,i])**2)+(2*sigxy[j,i])*np.sin(theta[j,i])*np.cos(theta[j,i])

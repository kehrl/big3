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

DIRM = os.path.join(os.getenv("MODEL_HOME"),"Helheim/3D/BASIN20120316/")
DIRR = DIRM+"mesh2d/"
DIR = DIRR+"inversion_adjoint/lambda_1e10_20160615_E3"
inputs = DIRM+"inputs/"
runname = "adjoint_beta"
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
x,y,sigxx = elmerreadlib.grid3d(surf,'stress1',holes,extent,dx=1.0e3) 
x,y,sigyy = elmerreadlib.grid3d(surf,'stress2',holes,extent,dx=1.0e3) 
x,y,sigxy = elmerreadlib.grid3d(surf,'stress4',holes,extent,dx=1.0e3) 
x,y,taub = elmerreadlib.grid3d(surf,'taub',holes,extent,dx=1.0e3) 
x,y,uu = elmerreadlib.grid3d(surf,'vel1',holes,extent,dx=1.0e3)
x,y,vv = elmerreadlib.grid3d(surf,'vel2',holes,extent,dx=1.0e3)

# Flow direction
theta = np.zeros_like(uu)
theta[:,:] = float('nan')
sigflow = np.zeros_like(uu)
sigflow[:,:] = float('nan')
taub_nomask = np.zeros_like(uu)
taub_nomask[:,:] = float('nan')

for i in range(0,len(x)):
  for j in range(0,len(y)):
    if uu[j,i] is not(np.ma.masked):
      theta[j,i] = np.arctan2(vv[j,i],uu[j,i])
      
      sigflow[j,i] = sigxx[j,i]*(np.cos(theta[j,i])**2)+sigyy[j,i]*(np.sin(theta[j,i])**2)+(2*sigxy[j,i])*np.sin(theta[j,i])*np.cos(theta[j,i])
      
      taub_nomask[j,i] = taub[j,i]

geotifflib.write_from_grid(x,y,np.flipud(sigflow),float('nan'),os.path.join(os.getenv("HOME"),"Bigtmp/sigma_alongflow.tif"))
geotifflib.write_from_grid(x,y,np.flipud(taub_nomask),float('nan'),os.path.join(os.getenv("HOME"),"Bigtmp/sigma_taub.tif"))

#fidsig = open(os.path.join(os.getenv("HOME"),"Bigtmp/sigma_alongflow.dat"),"w")
#fidsig.write("x y sigflow_MPa\n")
#fidtaub = open(os.path.join(os.getenv("HOME"),"Bigtmp/sigma_taub.dat"),"w")
#fidtaub.write("x y taub_MPa\n")
#for i in range(0,len(x)):
#  for j in range(0,len(y)):
#    fidsig.write('{0} {1} {2}\n'.format(x[i],y[j],sigflow[j,i]))
#    fidtaub.write('{0} {1} {2}\n'.format(x[i],y[j],taub[j,i]))
#fidsig.close()
#fidtaub.close()
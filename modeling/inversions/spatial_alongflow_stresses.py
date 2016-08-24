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
  xmin = 251000.
  xmax = 281000.
  ymin = -2595000.
  ymax = -2563000.
  #xmin = 283000.
  #xmax = 313000.
  #ymin = -2587000.
  #ymax = -2552000.

DIRM = os.path.join(os.getenv("MODEL_HOME"),"Helheim/3D/AquiferDS/")
DIRR = DIRM+"mesh2d/"
DIR = DIRR+"inversion_adjoint/lambda_1e12_20160822/"
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

dx=500.
surf = elmerreadlib.saveline_boundary(DIR,runname,bsurf,['stress','velocity','beta','eigenstress']) 
x,y,eigstress1 = elmerreadlib.grid3d(surf,'eigenstress 1',holes,extent,dx=dx) 
x,y,eigstress2 = elmerreadlib.grid3d(surf,'eigenstress 2',holes,extent,dx=dx) 
x,y,eigstress3 = elmerreadlib.grid3d(surf,'eigenstress 3',holes,extent,dx=dx) 
x,y,sigxx = elmerreadlib.grid3d(surf,'stress 1',holes,extent,dx=dx) 
x,y,sigyy = elmerreadlib.grid3d(surf,'stress 2',holes,extent,dx=dx) 
x,y,sigxy = elmerreadlib.grid3d(surf,'stress 4',holes,extent,dx=dx) 
x,y,taub = elmerreadlib.grid3d(surf,'taub',holes,extent,dx=dx) 
x,y,uu = elmerreadlib.grid3d(surf,'velocity 1',holes,extent,dx=dx)
x,y,vv = elmerreadlib.grid3d(surf,'velocity 2',holes,extent,dx=dx)

# Combine eigenstresses
eigstress_all = np.zeros([len(y),len(x),3])
eigstress_all[:,:,0] = eigstress1
eigstress_all[:,:,1] = eigstress2
eigstress_all[:,:,2] = eigstress3

# Flow direction
theta = np.zeros_like(uu)
theta[:,:] = float('nan')
sigflow = np.zeros_like(uu)
sigflow[:,:] = float('nan')
taub_nomask = np.zeros_like(uu)
taub_nomask[:,:] = float('nan')
eigstress = np.zeros_like(uu)
eigstress[:,:] = float('nan')

for i in range(0,len(x)):
  for j in range(0,len(y)):
    if uu[j,i] is not(np.ma.masked):
      # Along flow stress
      theta[j,i] = np.arctan2(vv[j,i],uu[j,i])
      # Principal stress
      #theta[j,i] = (1.0/2.0)*np.arctan2(2*sigxy[j,i],(sigxx[j,i]-sigyy[j,i]))
      
      sigflow[j,i] = sigxx[j,i]*(np.cos(theta[j,i])**2)+sigyy[j,i]*(np.sin(theta[j,i])**2)+(2*sigxy[j,i])*np.sin(theta[j,i])*np.cos(theta[j,i])
      
      taub_nomask[j,i] = taub[j,i]
      
      ind = np.argmax(abs(eigstress_all[j,i,:]))
      eigstress[j,i] = eigstress_all[j,i,ind]

geotifflib.write_from_grid(x,y,np.flipud(sigflow),float('nan'),os.path.join(os.getenv("HOME"),"Bigtmp/sigma_alongflow_DS.tif"))
geotifflib.write_from_grid(x,y,np.flipud(taub_nomask),float('nan'),os.path.join(os.getenv("HOME"),"Bigtmp/sigma_taub_DS.tif"))
geotifflib.write_from_grid(x,y,np.flipud(eigstress),float('nan'),os.path.join(os.getenv("HOME"),"Bigtmp/sigma_eigen_DS.tif"))

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
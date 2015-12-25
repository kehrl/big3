import os
import sys
from subprocess import call
import shutil
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
import elmer_read
import numpy as np
import matplotlib.pyplot as plt
import vtk
from vtk.util.numpy_support import vtk_to_numpy

DIRS=os.path.join(os.getenv("HOME"),"Code/Helheim/Modeling/SolverFiles/Flowline/")
DIRM=os.path.join(os.getenv("HOME"),"Models/Helheim/Meshes/Flowline/")
DIRR_rm=os.path.join(os.getenv("HOME"),"Models/Helheim/Results/Flowline/High/Retreat200m/")
flowline = np.loadtxt(os.path.join(os.getenv("HOME"),"Code/Helheim/Modeling/SolverFiles/Flowline/Inputs/flowline.dat"),skiprows=1)

runname1='flowline1'
runname4='flowline4_lag'

# Boundaries
bbed=1
bsurf=2
bcalve=3

#####################
# For 200 m retreat #
#####################

bed_rm = elmer_read.saveline_boundary(DIRR_rm+"Elmer/",runname4,bbed)
surf_rm = elmer_read.saveline_boundary(DIRR_rm+"Elmer/",runname4,bsurf)
surf1_rm = elmer_read.saveline_boundary(DIRR_rm+"Elmer/",runname1,bsurf)
terminus_rm = elmer_read.saveline_boundary(DIRR_rm+"Elmer/",runname4,bcalve)
terminus1_rm = elmer_read.saveline_boundary(DIRR_rm+"Elmer/",runname1,4)

# Initialize variables
dists = 82-np.array([0,5,10,20,28.9]) # kilometers
velocities_rm = np.zeros([len(np.unique(surf_rm['timestep']))+2,len(dists)])
time_rm = np.zeros([len(np.unique(surf_rm['timestep']))+2,1])
term_rm = np.zeros([len(np.unique(surf_rm['timestep']))+2,1])
sxx_rm = np.zeros([len(np.unique(surf_rm['timestep']))+2,1])
sxy_rm = np.zeros([len(np.unique(surf_rm['timestep']))+2,1])
termbed_rm = np.zeros([len(np.unique(surf_rm['timestep']))+2,1])

# Points before simulation
sortind=np.argsort(surf1_rm['coord1'])
velocities_rm[0,:] = np.interp(dists*1.0e3,surf1_rm['coord1'][sortind],surf1_rm['vel1'][sortind])
velocities_rm[1,:] = np.interp(dists*1.0e3,surf1_rm['coord1'][sortind],surf1_rm['vel1'][sortind])
time_rm[0] = [-1]
time_rm[1] = -0.01
term_rm[0] = terminus1_rm['coord1'][0]#-23
term_rm[1] = terminus1_rm['coord1'][0]
sxx_rm[0:2] = np.mean(terminus1_rm['sxx'])
termbed_rm[0:2] = np.interp(term_rm[0],flowline[:,0],flowline[:,3])

# Get velocity at select locations through time
timesteps=np.unique(surf_rm['timestep'])
for i in range(0,len(timesteps)):
  time_rm[i+2] = (i)*(1/2.0) 
  ind1 = np.where(surf_rm['timestep'] == timesteps[i])
  ind2 = np.where(terminus_rm['timestep'] == timesteps[i])
  coord1 = surf_rm['coord1'][ind1]
  vel1 = surf_rm['vel1'][ind1]
  bedind = (terminus_rm['coord2'][ind2]).argmin()
  
  sxy_rm[i+2] = terminus_rm['sxy'][ind2][bedind]
  sxx_rm[i+2] = np.mean(terminus_rm['sxx'][ind2])
  term_rm[i+2] = terminus_rm['coord1'][ind2][0]
  termbed_rm[i+2] = np.interp(term_rm[i+2],flowline[:,0],flowline[:,3])
  sortind=np.argsort(coord1)
  velocities_rm[i+2,:] = np.interp(dists*1e3,coord1[sortind],vel1[sortind])
  #plt.plot(coord1,vel1)
    
plt.figure(figsize=(8,9))
plt.subplot(311)
ax = plt.gca()
plt.plot([0,50],[0,0],'--',color='0.6',linewidth=3)
coloptions=['b','c','g','y','r']
for i in reversed(range(0,len(velocities_rm[0,:]))):
  plt.plot(time_rm,(velocities_rm[:,i]-velocities_rm[1,i])/1e3,'-',linewidth=2,color=coloptions[i],label=str(82-dists[i])+' km')
ax.set_xticklabels([])
plt.xticks(np.arange(0,50,10),fontsize=26)
plt.ylim([-0.5,2.1])
plt.yticks([0,1,2],fontsize=24)
plt.xlim([-1,50])
plt.ylabel('Relative velocity \n (km/yr)',fontsize=28)
plt.legend()

plt.subplot(312)
ax = plt.gca()
plt.plot([0,50],[0,0],'--',color='0.6',linewidth=2)
plt.plot(time_rm,(term_rm-term_rm[1])/1e3,'k-',linewidth=3)
plt.ylim([-1.1,0.4])
#plt.plot(time,termbed,'ko-',linewidth=2)
plt.xticks(np.arange(0,50,10),fontsize=24)
plt.xlim([-1,50])
plt.yticks([-1,-0.5,0],fontsize=26)
ax.set_xticklabels([])
plt.ylabel('Terminus \n (km)',fontsize=28)

plt.subplot(313)
plt.plot([-1,50],[termbed_rm[0],termbed_rm[0]],'--',color='0.6',linewidth=2)
plt.plot(time_rm,termbed_rm,'k-',linewidth=3)
plt.xticks(np.arange(0,50,10),fontsize=24)
plt.xlim([-1,50])
plt.ylim([-640,-570])
plt.yticks([-640,-610,-580],fontsize=26)
plt.xlabel('Days since calving event',fontsize=28)
plt.ylabel('Water depth \n (m asl)', fontsize=28)

plt.tight_layout()
plt.subplots_adjust(hspace=0.1,wspace=0)
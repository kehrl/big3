import os
import sys
from subprocess import call
import shutil
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Helheim/Observations/IceFronts"))
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Helheim/Observations/Velocity"))
import helheim_velocity, helheim_icefronts
import elmer_read
import numpy as np
import matplotlib.pyplot as plt
import vtk
from vtk.util.numpy_support import vtk_to_numpy

DIRS=os.path.join(os.getenv("HOME"),"Code/Helheim/Modeling/SolverFiles/Flowline/")
DIRM=os.path.join(os.getenv("HOME"),"Models/Helheim/Meshes/Flowline/")
DIRR_RKM=os.path.join(os.getenv("HOME"),"Models/Helheim/Results/Flowline/High/Retreat1000m/")
DIRR_rm=os.path.join(os.getenv("HOME"),"Models/Helheim/Results/Flowline/High/Retreat200m/")
flowline = np.loadtxt(os.path.join(os.getenv("HOME"),"Code/Helheim/Modeling/SolverFiles/Flowline/Inputs/flowline.dat"),skiprows=1)

runname1='flowline1'
runname4='flowline4_lag'

# Boundaries
bbed=1
bsurf=2
bcalve=3

#####################
# For 1000 m retreat #
#####################

bed_rkm = elmer_read.saveline_boundary(DIRR_RKM+"Elmer/",runname4,bbed)
surf_rkm = elmer_read.saveline_boundary(DIRR_RKM+"Elmer/",runname4,bsurf)
surf1_rkm = elmer_read.saveline_boundary(DIRR_RKM+"Elmer/",runname1,bsurf)
terminus_rkm = elmer_read.saveline_boundary(DIRR_RKM+"Elmer/",runname4,bcalve)
terminus1_rkm = elmer_read.saveline_boundary(DIRR_RKM+"Elmer/",runname1,4)

# Initialize variables
dists = 82-np.array([0,5,10,20,28.9]) # kilometers
velocities_rkm = np.zeros([len(np.unique(surf_rkm['timestep']))+2,len(dists)])
time_rkm = np.zeros([len(np.unique(surf_rkm['timestep']))+2,1])
term_rkm = np.zeros([len(np.unique(surf_rkm['timestep']))+2,1])
sxx_rkm = np.zeros([len(np.unique(surf_rkm['timestep']))+2,1])
sxy_rkm = np.zeros([len(np.unique(surf_rkm['timestep']))+2,1])
termbed_rkm = np.zeros([len(np.unique(surf_rkm['timestep']))+2,1])

# Points before simulation
sortind=np.argsort(surf1_rkm['coord1'])
velocities_rkm[0,:] = np.interp(dists*1.0e3,surf1_rkm['coord1'][sortind],surf1_rkm['vel1'][sortind])
velocities_rkm[1,:] = np.interp(dists*1.0e3,surf1_rkm['coord1'][sortind],surf1_rkm['vel1'][sortind])
time_rkm[0] = [-1]
time_rkm[1] = -0.01
term_rkm[0] = terminus1_rkm['coord1'][0]#-23
term_rkm[1] = terminus1_rkm['coord1'][0]
sxx_rkm[0:2] = np.mean(terminus1_rkm['sxx'])
termbed_rkm[0:2] = np.interp(term_rkm[0],flowline[:,0],flowline[:,3])

# Get velocity at select locations through time
timesteps=np.unique(surf_rkm['timestep'])
for i in range(0,len(timesteps)):
  time_rkm[i+2] = (i)*(1/2.0) 
  ind1 = np.where(surf_rkm['timestep'] == timesteps[i])
  ind2 = np.where(terminus_rkm['timestep'] == timesteps[i])
  coord1 = surf_rkm['coord1'][ind1]
  vel1 = surf_rkm['vel1'][ind1]
  bedind = (terminus_rkm['coord2'][ind2]).argmin()
  
  sxy_rkm[i+2] = terminus_rkm['sxy'][ind2][bedind]
  sxx_rkm[i+2] = np.mean(terminus_rkm['sxx'][ind2])
  term_rkm[i+2] = terminus_rkm['coord1'][ind2][0]
  termbed_rkm[i+2] = np.interp(term_rkm[i+2],flowline[:,0],flowline[:,3])
  sortind=np.argsort(coord1)
  velocities_rkm[i+2,:] = np.interp(dists*1e3,coord1[sortind],vel1[sortind])
  plt.plot(coord1,vel1,c='0.8',linewidth=2)

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

#Make figure    

plt.figure(figsize=(15,9))
gs = gridspec.GridSpec(3, 2, width_ratios=[2, 1]) 
plt.subplot(gs[0])
ax = plt.gca()
plt.plot([0,50],[0,0],'--',color='0.6',linewidth=3)
coloptions=['b','c','g','y','r']
for i in (range(0,len(velocities_rkm[0,:]))):
  plt.plot(time_rkm,(velocities_rkm[:,i]-velocities_rkm[1,i])/1e3,'-',linewidth=2,color=coloptions[i],label=str(82-dists[i])+' km')
ax.set_xticklabels([])
plt.xticks(np.arange(0,60,10),fontsize=26)
plt.ylim([-0.5,2.1])
plt.yticks([0,1,2],fontsize=24)
plt.xlim([-1,50])
plt.ylabel('Relative velocity \n (km/yr)',fontsize=28)
plt.legend(fontsize=14)

plt.subplot(gs[1])
ax = plt.gca()
plt.plot([0,50],[0,0],'--',color='0.6',linewidth=3)
coloptions=['b','c','g','y','r']
for i in (range(0,len(velocities_rm[0,:]))):
  plt.plot(time_rm,(velocities_rm[:,i]-velocities_rm[1,i])/1e3,'-',linewidth=2,color=coloptions[i],label=str(82-dists[i])+' km')
ax.set_xticklabels([])
plt.xticks(np.arange(0,60,10),fontsize=26)
plt.ylim([-0.5,2.1])
plt.yticks([0,1,2],fontsize=24)
plt.xlim([-1,25])
#plt.ylabel('Relative velocity \n (km/yr)',fontsize=28)
ax.set_yticklabels([])
plt.legend(fontsize=14)

plt.subplot(gs[2])
ax = plt.gca()
plt.plot([0,50],[term_rkm[0]/1e3-82,term_rkm[0]/1e3-82],'--',color='0.6',linewidth=2)
plt.plot(time_rkm,(term_rkm)/1e3-82,'k-',linewidth=3)
plt.ylim([1.9,3.4])
#plt.plot(time,termbed,'ko-',linewidth=2)
plt.xticks(np.arange(0,60,10),fontsize=24)
plt.xlim([-1,50])
plt.yticks([2,2.5,3.0],fontsize=26)
ax.set_xticklabels([])
plt.ylabel('Terminus \n (km)',fontsize=28)

plt.subplot(gs[3])
ax = plt.gca()
plt.plot([0,50],[term_rm[0]/1e3-82,term_rm[0]/1e3-82],'--',color='0.6',linewidth=2)
plt.plot(time_rm,(term_rm)/1e3-82,'k-',linewidth=3)
plt.ylim([1.9,3.4])
#plt.plot(time,termbed,'ko-',linewidth=2)
plt.xticks(np.arange(0,60,10),fontsize=24)
plt.xlim([-1,25])
plt.yticks([2,2.5,3.0],fontsize=26)
ax.set_xticklabels([])
ax.set_yticklabels([])
#plt.ylabel('Terminus \n (km)',fontsize=28)

plt.subplot(gs[4])
plt.plot([-1,50],[termbed_rkm[0],termbed_rkm[0]],'--',color='0.6',linewidth=2)
plt.plot(time_rkm,termbed_rkm,'k-',linewidth=3)
plt.xticks(np.arange(0,50,10),fontsize=24)
plt.xlim([-1,50])
plt.ylim([-640,-570])
plt.yticks([-640,-610,-580],fontsize=26)
plt.xlabel('                      Days since calving event',fontsize=28)
plt.ylabel('Bed \n (m asl)', fontsize=28)

plt.subplot(gs[5])
ax = plt.gca()
plt.plot([-1,50],[termbed_rm[0],termbed_rm[0]],'--',color='0.6',linewidth=2)
plt.plot(time_rm,termbed_rm,'k-',linewidth=3)
plt.xticks(np.arange(0,60,10),fontsize=24)
plt.xlim([-1,25])
plt.ylim([-640,-570])
plt.yticks([-640,-610,-580],fontsize=26)
ax.set_yticklabels([])
#plt.xlabel('Days since calving event',fontsize=28)
#plt.ylabel('Water depth \n (m asl)', fontsize=28)

plt.tight_layout()
plt.subplots_adjust(hspace=0.1,wspace=0.05)

# Compare modeled and measured
velmeasured,times=helheim_velocity.velocity_along_flowline(flowline[:,1],flowline[:,2],flowline[:,0])

plt.figure(figsize=(7,7))

plt.subplot(311)
ax=plt.gca()
meanvel=np.zeros([len(velflowlines[:,0]),1])
for i in range(0,len(velflowlines[:,0])):
  ind = np.where(~(np.isnan(velflowlines[i,:])))
  meanvel[i]=np.mean(velflowlines[i,ind])
  
plt.plot(flowline[:,0]/1e3-82,velflowlines[:,:]/1e3,c='b',alpha=0.3,linewidth=2)
plt.plot(flowline[:,0]/1e3-82,meanvel/1e3,c='b',linewidth=2)
plt.plot([-100],[-10],c='b',linewidth=2,label='Measured')
plt.plot([-100],[-10],c='r',linewidth=2,label='Simulated')
plt.xlim([-50,5])
plt.ylim([0,12])
plt.xticks([-40,-30,-20,-10,0],fontsize=24)
ax.set_xticklabels([])
plt.yticks([0,5,10],fontsize=24)
plt.legend(fontsize=20,handlelength=0,loc=2)

plt.subplot(312)
timesteps=np.unique(surf_rkm['timestep'])
ind1 = np.where(surf_rkm['timestep'] == timesteps[0])
coord1s = surf_rkm['coord1'][ind1]
vel1s = np.zeros([len(timesteps),len(coord1)])
for i in range(0,len(timesteps)):
  ind1 = np.where(surf_rkm['timestep'] == timesteps[i])
  coord1 = surf_rkm['coord1'][ind1]
  vel1 = surf_rkm['vel1'][ind1]
  vel1s[i,:]=np.interp(coord1s,coord1,vel1)
  plt.plot(coord1/1e3-82,vel1/1e3,c='r',alpha=0.4,linewidth=2)
plt.plot(coord1s/1e3-82,np.mean(vel1s,0)/1e3,c='r',linewidth=2)
plt.xlim([-50,5])
plt.xticks([-40,-30,-20,-10,0],fontsize=24)
plt.yticks([0,5,10],fontsize=24)
plt.ylabel('           Velocity (km/yr)',fontsize=28)


plt.subplot(313)
plt.plot(flowline[:,0]/1e3-82,flowline[:,3],'k',linewidth=2)
plt.xlim([-50,5])
plt.xticks([-40,-30,-20,-10,0],fontsize=24)
plt.ylim([-1000,800])
plt.yticks([-1000,-500,0,500],fontsize=24)
plt.ylabel('Bed (m asl)',fontsize=28)
plt.xlabel('Distance along flowline (km)',fontsize=28)

plt.tight_layout()
plt.subplots_adjust(hspace=0,wspace=0.05)

#reader = vtk.vtkXMLUnstructuredGridReader()
#reader.SetFileName(DIRR_RKM+"Retreat1000m/flowline4_lag.0011.vtu")
#reader.Update()

#sxx_vtk= reader.GetOutput().GetPointData().GetArray(16)
#sxx_all = vtk_to_numpy(sxx_vtk)

#nodes_vtk_array= reader.GetOutput().GetPoints().GetData()
#nodes_nummpy_array = vtk_to_numpy(nodes_vtk_array)
#x,y,z= nodes_nummpy_array[:,0] , nodes_nummpy_array[:,1] , nodes_nummpy_array[:,2]

#x_z=[]
#xlist = x[np.where(~np.isnan(x))]
#sxx_z = []
#satisfied = 0
#while not(satisfied):
#  ind1 = np.where(abs(x-np.min(xlist)) < 20)
#  ind2 = np.where(abs(xlist-np.min(xlist)) < 20)
#  x_z.append(np.mean(x[ind1]))
#  sxx_z.append(np.mean(sxx_all[ind1]))
#  xlist = np.delete(xlist,ind2[0])
#  if not(xlist.size):
#    satisfied = 1


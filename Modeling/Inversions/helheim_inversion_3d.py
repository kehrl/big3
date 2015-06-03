# This script let's you choose the velocity record for the basal inversion, and then runs 
# the Helheim inversion script.
#
# LMK, UW, 06/12/2014

import os
import shutil
import sys
import numpy as np
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
from geodat import *
from subprocess import call
import math
import glob
import numpy as np
import elmer_read 
import elmer_mesh as mesh
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

##########
# Inputs #
##########

# Model Resolution
RES = 'InverseClass'

# Directories
DIRS=os.path.join(os.getenv("HOME"),"Code/Helheim/Modeling/SolverFiles/3D/")
DIRM=os.path.join(os.getenv("HOME"),"Models/Helheim/Meshes/3D/"+RES+"/")
DIRR=os.path.join(os.getenv("HOME"),"Models/Helheim/Results/3D/"+RES+"/Inversion/")
DIRX=os.path.join(os.getenv("HOME"),"Data/ShapeFiles/Glaciers/3D/Helheim/")
Inputs=os.path.join(DIRM+"/Inputs/")

if not(os.path.exists(DIRR)):
  os.makedirs(DIRR)

# Copy model inputs into the solver file directory
input_files=os.listdir(DIRM+"Inputs/")
for file in input_files:
  shutil.copy(DIRM+"Inputs/"+file,DIRS+"Inputs/")
del input_files

# Boundary numbers 
bbed=3
bsurf=4
runname="beta"

# Regularization parameters (lambda)
regpars=['1e10','1e8','1e9','1e11','1e12','1e7','1e13','1e6','1e14']

# Mesh files for making pretty graphs
extent = np.genfromtxt(Inputs+"mesh_extent.dat")
hole1 = np.genfromtxt(Inputs+"mesh_hole1.dat")
hole2 = np.genfromtxt(Inputs+"mesh_hole2.dat")
holes = []
holes.append({'xy': np.array(hole1[:,0:2])})
holes.append({'xy': np.array(hole2[:,0:2])})
del hole1, hole2

############################################################
# Run inversion solver file for different values of lambda #
############################################################

print "\n## Running elmer inversion code ##\n"

fid = open(DIRS+"ELMERSOLVER_STARTINFO","w")
fid.write('robin_beta_temp.sif')
fid.close()

fid_info = open(DIRR+"summary.dat","a")
fid_info.write('Lambda Nsim Cost Norm RelPrec_G \n')
fid_info.close()
n=-1
for regpar in regpars:
  for filename in glob.glob(DIRM+"elmer/robin_beta*"):
    os.remove(filename)
  n=n+1
  os.chdir(DIRS)
  fid1 = open('robin_beta_temperature.sif', 'r')
  fid2 = open('robin_beta_temp.sif', 'w')
  lines=fid1.readlines()
  for line in lines:
    line=line.replace('$Lambda=1.0e10', '$Lambda={}'.format(regpar))
    line=line.replace('Low','{}'.format(RES))
    fid2.write(line)
  fid1.close()
  fid2.close()
  del fid1, fid2
  call(["mpirun","-quiet","-n","4","elmersolver_mpi"])
  os.system('rm robin_beta_temp.sif')
  
  #####################################
  # Write cost values to summary file #
  ##################################### 
  fid = open(DIRS+"cost_robin_beta.dat","r")
  lines = fid.readlines()
  line=lines[-1]
  p=line.split()
  nsim = float(p[0])
  cost1 = float(p[1])
  cost2 = float(p[2])
  norm = float(p[3]) 
  fid.close()
  fid_info = open(DIRR+"summary".dat","a")
  fid_info.write('{} {} {} {} {}\n'.format(regpar,nsim,cost1,cost2,norm))
  fid_info.close()
  del fid
  
  #######################################
  # Combine elmer results into one file #
  #######################################
  bed = elmer_read.saveline_boundary(DIRM+"/elmer/",runname,bbed)
  surf = elmer_read.saveline_boundary(DIRM+"/elmer/",runname,bsurf)
  
  os.rename(DIRM+"/elmer/"+runname+".dat",DIRR+runname+"_"+regpar+".dat")
  os.rename(DIRM+"/elmer/"+runname+".dat.names",DIRR+runname+"_"+regpar+".dat.names")
  os.rename(DIRS+"M1QN3_robin_beta.out",DIRR+"M1QN3_"+regpar+".out")
  os.rename(DIRS+"gradientnormadjoint_robin_beta.dat",DIRR+"gradient_"+regpar+".dat")
  os.rename(DIRS+"cost_robin_beta.dat",DIRR+"cost_"+regpar+".dat")
  
  names = os.listdir(DIRM+"/elmer")
  os.chdir(DIRM+"/elmer")
  if not os.path.exists(DIRR+"lambda_"+regpar):
    os.makedirs(DIRR+"lambda_"+regpar)
  for name in names:
    if name.endswith('vtu') and name.startswith('robin'):
      os.rename(name,DIRR+"lambda_"+regpar+"/"+name)
  
  
################################
# Output friction coefficients #
################################
#Linear Beta square
fid = open(Inputs+"beta_linear.xyz",'w')
fid.write('{0}\n'.format(len(bed['node'])))
for i in range(0,len(bed['node'])):
  fid.write('{0} {1} {2:.4f} {3}\n'.format(bed['coord1'][i],bed['coord2'][i],bed['coord3'][i],bed['beta'][i]**2))
fid.close()

#Weertman coefficient
fid = open(Inputs+"beta_weertman.xyz",'w')
fid.write('{0}\n'.format(len(bed['node'])))
for i in range(0,len(bed['node'])):
  coeff=(bed['beta'][i]**2)*(bed['vel'][i]**(2.0/3.0))
  fid.write('{0} {1} {2:.4f} {3}\n'.format(bed['coord1'][i],bed['coord2'][i],bed['coord3'][i],coeff))
fid.close()

###########################################################
# Make figures of basal shear stress and velocity anomaly #
###########################################################
# Plot basal shear stress
taub_3D=elmer_read.grid3d(bed,'taub',holes,extent)
velmod_3D=elmer_read.grid3d(surf,'vel',holes,extent)
velmes_3D=elmer_read.grid3d(surf,'velmes',holes,extent)

limits=[min(taub_3D[0]),max(taub_3D[0]),max(taub_3D[1]),min(taub_3D[1])]
plotnames=['basin','terminus']
plotlimits=np.array([[min(taub_3D[0]),max(taub_3D[0]),min(taub_3D[1]),max(taub_3D[1])],
            [287000,310000,-2585000,-2558000]])

for i in range(0,len(plotlimits)):
  plt.clf()
  plt.figure(figsize=(12,12))
  ax1 = plt.subplot(221)
  plt.imshow(velmes_3D[2]/1000,extent=limits)
  plt.gca().invert_yaxis()
  plt.clim([0,8])
  cf = plt.colorbar(orientation='horizontal',fraction=0.3,pad=0.03)#,ax=ax1,shrink=0.4)
  plt.setp( ax1.get_xticklabels(), visible=False)
  plt.setp( ax1.get_yticklabels(), visible=False)
  plt.title('Measured velocity')
  cf.set_ticks([0,2,4,6])
  cf.set_label('km yr$^{-1}$')
  plt.ylim([plotlimits[i,2],plotlimits[i,3]])
  plt.xlim(plotlimits[i,0],plotlimits[i,1])

  ax2 = plt.subplot(222)
  plt.imshow(velmod_3D[2]/1000,extent=limits)
  plt.gca().invert_yaxis()
  plt.clim([0,8])
  cf = plt.colorbar(orientation='horizontal',fraction=0.3,pad=0.03)#,ax=ax1,shrink=0.4)
  plt.setp( ax2.get_xticklabels(), visible=False)
  plt.setp( ax2.get_yticklabels(), visible=False)
  plt.title('Modelled velocity')
  cf.set_ticks([-1,0,1])
  cf.set_label('km yr$^{-1}$')
  #cf.ax.tick_params(labelsize=20)
  plt.xlim(plotlimits[i,0],plotlimits[i,1])
  plt.ylim([plotlimits[i,2],plotlimits[i,3]]) 

  ax3 = plt.subplot(223)
  plt.imshow((velmod_3D[2]-velmes_3D[2])/1000,cmap='RdBu_r',extent=limits)
  plt.gca().invert_yaxis()
  plt.clim([-1,1])
  cf = plt.colorbar(orientation='horizontal',fraction=0.3,pad=0.03)#,ax=ax1,shrink=0.4)
  plt.setp( ax3.get_xticklabels(), visible=False)
  plt.setp( ax3.get_yticklabels(), visible=False)
  plt.title('Velocity anomaly')
  cf.set_label('km yr$^{-1}$')
  cf.set_ticks([-1,0,1])
  plt.xlim(plotlimits[i,0],plotlimits[i,1])
  plt.ylim([plotlimits[i,2],plotlimits[i,3]])

  ax4 = plt.subplot(224)
  plt.imshow(taub_3D[2]*10**6,norm=matplotlib.colors.LogNorm(),extent=limits)
  plt.gca().invert_yaxis()
  plt.clim([0.1,10**6])
  cf = plt.colorbar(orientation='horizontal',fraction=0.3,pad=0.03)#,ax=ax1,shrink=0.4)
  plt.setp( ax4.get_xticklabels(), visible=False)
  plt.setp( ax4.get_yticklabels(), visible=False)
  plt.title('Basal shear stress')
  cf.set_ticks([1e1,1e3,1e6])
  cf.set_label('Pa')
  plt.xlim(plotlimits[i,0],plotlimits[i,1])
  plt.ylim([plotlimits[i,2],plotlimits[i,3]])
  plt.savefig(DIRR+"results_"+plotnames[i]+".jpg")
   

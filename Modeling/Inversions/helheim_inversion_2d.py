# This script let's you choose the velocity record for the basal inversion, and then runs 
# the Helheim inversion script.
#
# LMK, UW, 06/12/2014

import os
import sys
import shutil
import numpy as np
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"Util/Modules"))
from subprocess import call
import math
import scipy.interpolate as interpolate
import glob
import elmer_read

##########
# Inputs #
##########

#MESHNAME='DEM20120624Low'
MESHNAME='DEM20120624Low'

glacier = 'Helheim'

# Regularization parameters (lambda)

#regpars=['1e7','1e8','5e8','1e9','5e9','1e10','5e10','1e11','5e11','1e12','5e12','1e13','5e13','1e14','5e14','1e15','5e15','1e16','5e16','1e17']

#regpars=["1e10"]
regpars=['1e14']

# Directories
DIRS=os.path.join(os.getenv("CODE_HOME"),"BigThreeGlaciers/Modeling/SolverFiles/Flowline/")
DIRM=os.path.join(os.getenv("MODEL_HOME"),glacier+"/Meshes/Flowline/"+MESHNAME+"/")
DIRR=os.path.join(os.getenv("MODEL_HOME"),glacier+"/Results/Flowline/"+MESHNAME+"/Inversion/")
DIRELMERLIB = os.path.join(os.getenv("CODE_HOME"),"BigThreeGlaciers/Modeling/Elmerlib/")

# Flowline
file_flowline_in=DIRM+"Inputs/flowline.dat"

# Elmer results name
runname = 'beta'
bbed=1
bsurf=2

#####################
# Set up model runs #
#####################

# Copy model inputs into the solver file directory
input_files=os.listdir(DIRM+"Inputs/")
for file in input_files:
  shutil.copy(DIRM+"Inputs/"+file,DIRS+"Inputs/")
del input_files

# Compile fortran functions for flowline model
readfiles = glob.glob(DIRELMERLIB+"Flowline/"+"*.f90")

with open(DIRELMERLIB+"Flowline.f90", "wb") as outfile:
  for f in readfiles:
    with open(f, "rb") as infile:
      outfile.write(infile.read())
      outfile.write("\n")
outfile.close()

os.chdir(DIRELMERLIB)
call(["elmerf90","-o","Flowline.so","Flowline.f90"])
del DIRELMERLIB

############################################################
# Run inversion solver file for different values of lambda #
############################################################

if not(os.path.isdir(DIRR)):
  os.makedirs(DIRR)

fid = open(DIRS+"ELMERSOLVER_STARTINFO","w")
fid.write('temp.sif')
fid.close()

fid_info = open(DIRR+"summary.dat","a")
fid_info.write('Lambda Nsim Cost1 Cost2 Norm \n')
fid_info.close()
n=-1
for regpar in regpars:
  n=n+1
  os.chdir(DIRS)
  fid1 = open('robin_beta.sif', 'r')
  fid2 = open('temp.sif', 'w')
  lines=fid1.readlines()
  for line in lines:
    line=line.replace('$Lambda=1e10', '$Lambda={0}'.format(regpar))
    line=line.replace('Mesh_Input','{0}'.format("../../../../../Models/Helheim/Meshes/Flowline/"+MESHNAME))
    fid2.write(line)
  fid1.close()
  fid2.close()
  del fid1, fid2
  call(["mpirun","-n","4","elmersolver_mpi"])
  #os.system('rm robin_beta_temp.sif')
  
  ###############################################################################
  # Get number of simulations, cost value, and mean velocity difference between #
  # modelled and measured velocities for each regularization parameter          #
  ###############################################################################
  fid = open(DIRS+"cost_robin_beta.dat","r")
  lines = fid.readlines()
  line=lines[-1]
  p=line.split()
  nsim = float(p[0])
  cost1 = float(p[1])
  cost2 = float(p[2])
  norm = float(p[3]) 
  fid.close()
  fid_info = open(DIRR+"summary.dat","a")
  fid_info.write('{0} {1} {2} {3} {4}\n'.format(regpar,nsim,cost1,cost2,norm))
  fid_info.close()
  del fid
  
  #####################################################
  # Save results into different files and directories #
  #####################################################
  
  surf_bed = elmer_read.saveline_boundary(DIRM+"elmer/",runname,bbed)
  surf_sur = elmer_read.saveline_boundary(DIRM+"elmer/",runname,bsurf)
  
  os.rename(DIRM+"elmer/"+runname+".dat",DIRR+runname+"_"+regpar+".dat")
  os.rename(DIRM+"elmer/"+runname+".dat.names",DIRR+runname+"_"+regpar+".dat.names")
  os.rename(DIRS+"M1QN3_robin_beta.out",DIRR+"M1QN3_"+regpar+".out")
  os.rename(DIRS+"cost_robin_beta.dat",DIRR+"cost_"+regpar+".dat")
  os.rename(DIRS+"gradientnormadjoint_robin_beta.dat",DIRR+"gradient_"+regpar+".dat")
  
  names = os.listdir(DIRM+"elmer")
  os.chdir(DIRM+"elmer")
  if not os.path.exists(DIRR+"lambda_"+regpar):
    os.makedirs(DIRR+"lambda_"+regpar)
  for name in names:
    if name.endswith('vtu') and name.startswith('robin'):
      os.rename(name,DIRR+"lambda_"+regpar+"/"+name)
  
fid_info.close()
del fid_info

##############################################
# Export friction coefficient for future simulations #
##############################################
#Linear Beta square
fid = open(DIRM+"/Inputs/beta_linear.xy",'w')
fid.write('{0}\n'.format(len(surf_bed['node'])))
ind=np.argsort(surf_bed['coord1'])
for i in ind:
  fid.write('{0} {1} {2:.4g}\n'.format(surf_bed['coord1'][i],surf_bed['coord2'][i],surf_bed['beta'][i]**2))
fid.close()

#Weertman coefficient
fid = open(DIRM+"/Inputs/beta_weertman.dat",'w')
coeff=(surf_bed['beta'][ind[0]]**2)*(surf_bed['vel1'][ind[0]]**(2.0/3.0))
fid.write('{0} {1:.4g}\n'.format(surf_bed['coord1'][ind[0]]-50,coeff))
for i in ind:
  coeff=(surf_bed['beta'][i]**2)*(surf_bed['vel1'][i]**(2.0/3.0))
  fid.write('{0} {1:.4g}\n'.format(surf_bed['coord1'][i],coeff))
fid.write('{0} {1:.4g}\n'.format(surf_bed['coord1'][i]+5000,coeff))
fid.close()


print "\n## Done running helheim_inversion_2d.py ##\n"

# This code runs a steady state Helheim simulation.

import os
import sys
from subprocess import call
import shutil
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
import elmer_read
import numpy as np
import matplotlib.pyplot as plt

##########
# Inputs #
##########

MESHNAME='High' 

# Directories
DIRS=os.path.join(os.getenv("HOME"),"Code/Helheim/Modeling/SolverFiles/Flowline/")
DIRM=os.path.join(os.getenv("HOME"),"Models/Helheim/Meshes/Flowline/"+MESHNAME+"/")
DIRR=os.path.join(os.getenv("HOME"),"Models/Helheim/Results/Flowline/"+MESHNAME+"/")

# Solver file
solverfile='flowline2.sif'

#solverfile = 'flowline2.sif'
#runname='flowline2'
#solverfile='flowline4_lagrangian_nofree.sif'
runname='flowline2'

# Boundaries
bbed=1
bsurf=2
bcalve=3

############################################
# Check that fortran routines are compiled #
############################################
os.chdir(os.path.join(os.getenv("HOME"),"Code/Helheim/Modeling/Elmerlib"))
routines = ['Flowline','USF_Calving','USF_Zs']
#for routine in routines:
  #if not(os.path.exists(DIRS+"../../lib/"+routine+".so")):
  #call(["elmerf90","-o",routine+".so",routine+".f90"])
  
####################################################
# Copy model inputs into the solver file directory #
####################################################
input_files=os.listdir(DIRM+"Inputs/")
for file in input_files:
  shutil.copy(DIRM+"Inputs/"+file,DIRS+"Inputs/")
del input_files

###################################
# Change mesh file in solver file #
###################################
os.chdir(DIRS)
fid1 = open(solverfile, 'r')
fid2 = open('temp.sif', 'w')
lines=fid1.readlines()
for line in lines:
  line=line.replace('Mesh_Input','{}'.format("../../../../../Models/Helheim/Meshes/Flowline/"+MESHNAME))
  fid2.write(line)
fid1.close()
fid2.close()


###################
# Run elmersolver #
###################
os.chdir(DIRS)
fid = open(DIRS+"ELMERSOLVER_STARTINFO","w")
fid.write('{}'.format('temp.sif'))
fid.close()

call(["mpirun","-n","4","elmersolver_mpi"])

##############################
# Load results from saveline #
##############################
bed = elmer_read.saveline_boundary(DIRM+"Elmer",runname,bbed)
surf = elmer_read.saveline_boundary(DIRM+"Elmer",runname,bsurf)
surf1 = elmer_read.saveline_boundary(DIRM+"Elmer",'flowline1',bsurf)
terminus = elmer_read.saveline_boundary(DIRM+"Elmer",runname,bcalve)

coord1=np.sort(surf['coord1'][np.where(surf['timestep'] == 1)])
elevation=np.zeros([len(coord1),np.max(surf['timestep'])])
timesteps=np.unique(surf['timestep'])
for i in range(0,len(timesteps)):
  ind1 = np.where(surf['timestep'] == timesteps[i])
  x = surf['coord1'][ind1]
  y = surf['coord2'][ind1]
  sortind=np.argsort(x)
  elevation[:,i]=y[sortind]

  
# This code runs a steady state Helheim simulation.

import os
import sys
from subprocess import call
import shutil, glob
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"Util/Modules"))
import elmer_read
import numpy as np
import matplotlib.pyplot as plt

##########
# Inputs #
##########

MESHNAME='DEM20120624Low' 

# Directories
DIRS=os.path.join(os.getenv("CODE_HOME"),"BigThreeGlaciers/Modeling/SolverFiles/Flowline/")
DIRM=os.path.join(os.getenv("MODEL_HOME"),"Helheim/Meshes/Flowline/"+MESHNAME+"/")
DIRR=os.path.join(os.getenv("MODEL_HOME"),"Helheim/Results/Flowline/"+MESHNAME+"/")
DIRELMERLIB = os.path.join(os.getenv("CODE_HOME"),"BigThreeGlaciers/Modeling/Elmerlib/")

# Solver file
solverfile='flowline4.sif'
runname=solverfile[0:-4]

# Variables
dt = 1/365.25
totaldt = 10

# Boundaries
bbed=1
bsurf=2
bcalve=3

############################################
# Check that fortran routines are compiled #
############################################
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
  line=line.replace('totaldt','{}'.format(totaldt))
  line=line.replace('$(dt)','{}'.format(dt))
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

call(["mpiexec","-n","4","elmersolver_mpi"])

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

  

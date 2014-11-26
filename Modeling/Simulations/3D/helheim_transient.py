# This code puts runs a 3D Helheim simulation, using one of the solver files and a mesh given
# by "RES". It then saves the results in a directory in "Results/Models/"+"RES".
#
#
# LMK, UW, 10/16/2014

import os
from subprocess import call
import shutil

##########
# Inputs #
##########

# Mesh
RES="Low_Normal"

# Solver file
solverfile='model6.sif'

# Directories
DIRS=os.path.join(os.getenv("HOME"),"Code/Modeling/SolverFiles/3D/Helheim/")
DIRM=os.path.join(os.getenv("HOME"),"Models/Meshes/3D/Helheim/"+RES+"/")
DIRR=os.path.join(os.getenv("HOME"),"Models/Results/3D/Helheim/"+RES+"/Inversion/")
DIRX=os.path.join(os.getenv("HOME"),"Data/ShapeFiles/Glaciers/3D/Helheim/")

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
  line=line.replace('Mesh_Input','{}'.format("../../../../../Models/Meshes/3D/Helheim/"+RES))
  fid2.write(line)
fid1.close()
fid2.close()

###################
# Run Elmersolver #
###################
os.chdir(DIRS)
fid = open(DIRS+"ELMERSOLVER_STARTINFO","w")
fid.write('{}'.format('temp.sif'))
fid.close()

call(["mpirun","-n","4","elmersolver_mpi"])
os.system('rm temp.sif')
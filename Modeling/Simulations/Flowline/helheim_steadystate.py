# This code runs a steady state Helheim simulation.

import os
from subprocess import call
import shutil
import glob

##########
# Inputs #
##########

MESHNAME='High' 

# Directories
DIRS=os.path.join(os.getenv("HOME"),"Code/Helheim/Modeling/SolverFiles/Flowline/")
DIRM=os.path.join(os.getenv("HOME"),"Models/Helheim/Meshes/Flowline/"+MESHNAME+"/")
DIRR=os.path.join(os.getenv("HOME"),"Models/Helheim/Results//Flowline/")
DIRELMERLIB = os.path.join(os.getenv("HOME"),"Code/Helheim/Modeling/Elmerlib/")

# Solver file
solverfile='flowline1.sif'

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



call(["mpirun","-quiet","-n","4","elmersolver_mpi"])
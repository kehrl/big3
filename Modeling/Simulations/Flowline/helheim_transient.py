# This code runs a steady state Helheim simulation.

import os
from subprocess import call
import shutil

##########
# Inputs #
##########

MESHNAME='Worldview_Advance' 

# Directories
DIRS=os.path.join(os.getenv("HOME"),"Code/Helheim/Modeling/SolverFiles/Flowline/")
DIRM=os.path.join(os.getenv("HOME"),"Models/Helheim/Meshes/Flowline/"+MESHNAME+"/")
DIRR=os.path.join(os.getenv("HOME"),"Models/Helheim/Results/Flowline/"+MESHNAME+"/")

# Solver file
#solverfile='flowline1.sif'
solverfile='flowline4_lagrangian_nofree.sif'

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
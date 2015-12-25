# This script let's you choose the velocity record for the basal inversion, and then runs 
# the Helheim inversion script.
#
# LMK, UW, 06/12/2014

import os
import shutil
import sys
import numpy as np
from subprocess import call
import math
import glob
import numpy as np
import elmerreadlib, meshlib
import numpy as np
import argparse

##########
# Inputs #
##########

# Get inputs to file
parser = argparse.ArgumentParser()
parser.add_argument("-mesh", dest="meshname", required = True,
        help = "Name of meshlib")
parser.add_argument("-n", dest="npartitions", required = True,
        help = "Number of partitions.")
parser.add_argument("-regpar", dest="regpar", required = False,
		default='1e10',help = "Regularization parameter.")
parser.add_argument("-method", dest="method", required = False,
		default='adjoint',help = "adjoint or robin.")

args, _ = parser.parse_known_args(sys.argv)
RES = args.meshname
partitions = args.npartitions
regpar = args.regpar
method = args.method

# Model Resolution
glacier = 'Helheim'

# Directories
DIRS=os.path.join(os.getenv("CODE_HOME"),"big3/modeling/solverFiles/3D/")
DIRM=os.path.join(os.getenv("MODEL_HOME"),glacier+"/Meshes/3D/"+RES+"/")
DIRR=os.path.join(os.getenv("MODEL_HOME"),glacier+"/Results/3D/"+RES+"/Inversion/")
DIRX=os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/Glaciers/3D/"+glacier)
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
runname=method+"_beta"

############################################################
# Run inversion solver file for different values of lambda #
############################################################

print "\n## Running elmer inversion code ##\n"

fid = open(DIRS+"ELMERSOLVER_STARTINFO","w")
fid.write('temp.sif')
fid.close()

fid_info = open(DIRR+"summary.dat","a")
fid_info.write('Lambda Nsim Cost Norm RelPrec_G \n')
fid_info.close()

#for filename in glob.glob(DIRM+"elmer/robin_beta*"):
#  os.remove(filename)
os.chdir(DIRS)
fid1 = open(method+'_beta.sif', 'r')
fid2 = open('temp.sif', 'w')
lines=fid1.readlines()
for line in lines:
  line=line.replace('$Lambda=1.0e10', '$Lambda={}'.format(regpar))
  line=line.replace('Mesh_Input','{}'.format("../../../../../Models/Helheim/Meshes/3D/"+RES))
  fid2.write(line)
fid1.close()
fid2.close()
del fid1, fid2
call(["mpiexec","-np",partitions,"elmersolver_mpi"])
os.system('rm temp.sif')
  
#####################################
# Write cost values to summary file #
##################################### 

fid = open(DIRS+"cost_"+method+"_beta.dat","r")
lines = fid.readlines()
line=lines[-1]
p=line.split()
nsim = float(p[0])
cost1 = float(p[1])
cost2 = float(p[2])
norm = float(p[3]) 
fid.close()
fid_info = open(DIRR+"summary.dat","a")
fid_info.write('{} {} {} {} {}\n'.format(regpar,nsim,cost1,cost2,norm))
fid_info.close()
del fid

#######################################
# Combine elmer results into one file #
#######################################

bed = elmerreadlib.saveline_boundary(DIRM+"/elmer/",runname,bbed)
surf = elmerreadlib.saveline_boundary(DIRM+"/elmer/",runname,bsurf)

os.rename(DIRM+"/elmer/"+runname+".dat",DIRR+runname+method+"_"+regpar+"_beta.dat")
os.rename(DIRM+"/elmer/"+runname+".dat.names",DIRR+runname+method+"_"+regpar+"_beta.dat.names")
os.rename(DIRS+"M1QN3_"+method+"_beta.out",DIRR+"M1QN3_"+method+"_"+regpar+"_beta.out")
os.rename(DIRS+"gradientnormadjoint_"+method+"_beta.dat",DIRR+"gradient_"+method+"_"+regpar+"_beta.dat")
os.rename(DIRS+"cost_"+method+"_beta.dat",DIRR+"cost_"+method+"_"+regpar+"_beta.dat")

names = os.listdir(DIRM+"/elmer")
os.chdir(DIRM+"/elmer")
if not os.path.exists(DIRR+"lambda_"+regpar):
  os.makedirs(DIRR+"lambda_"+regpar)
for name in names:
  if name.endswith('vtu') and name.startswith(method):
    os.rename(name,DIRR+"lambda_"+regpar+"/"+name)
  
  
################################
# Output friction coefficients #
################################

#Linear Beta square
fid = open(Inputs+"beta_linear.xyz",'w')
fid.write('{0}\n'.format(len(bed['node'])))
for i in range(0,len(bed['node'])):
  fid.write('{0} {1} {2:.4f} {3}\n'.format(bed['coord1'][i],bed['coord2'][i],\
  		bed['coord3'][i],bed['beta'][i]**2))
fid.close()

#Weertman coefficient
fid = open(Inputs+"beta_weertman.xyz",'w')
fid.write('{0}\n'.format(len(bed['node'])))
for i in range(0,len(bed['node'])):
  coeff=(bed['beta'][i]**2)*(bed['vel'][i]**(2.0/3.0))
  fid.write('{0} {1} {2:.4f} {3}\n'.format(bed['coord1'][i],bed['coord2'][i],bed['coord3'][i],coeff))
fid.close() 


# Remove files in solver Input directory
for file in os.listdir(DIRS+"Inputs/"):
  file_path = os.path.join(DIRS+"Inputs/", file)
  try:
   if os.path.isfile(file_path):
     os.unlink(file_path)  
  except:
    pass

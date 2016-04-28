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
import datetime
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

def get_arguments():

  # Get inputs to file
  parser = argparse.ArgumentParser()
  parser.add_argument("-mesh", dest="mesh", required = True,
        help = "Name of meshlib") 
  parser.add_argument("-n", dest="n", required = True,
        help = "Number of partitions.")
  parser.add_argument("-regpar", dest="regpar", required = False,
		default='1e10',help = "Regularization parameter.")
  parser.add_argument("-method", dest="method", required = False,
		default='adjoint',help = "adjoint or robin.")
  parser.add_argument("-extrude", dest="extrude", type=int,required = False,\
		default=5,\
		help = "Number of extrusion levels.")

  args, _ = parser.parse_known_args(sys.argv)


  return args

def main():

  ##########
  # inputs #
  ##########

  args = get_arguments()

  RES = args.mesh
  partitions = str(args.n)
  regpar = str(args.regpar)
  method = args.method
  extrude = str(args.extrude)

  # Get current date
  now = datetime.datetime.now()
  date = '{0}{1:02.0f}{2:02.0f}'.format((now.year),(now.month),(now.day))

  # Model Resolution
  glacier = 'Helheim'

  # Directories
  DIRS=os.path.join(os.getenv("CODE_HOME"),"big3/modeling/solverfiles/3D/")
  DIRM=os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/"+RES+"/")
  DIRR=os.path.join(DIRM+'inversion_'+method+'/')
  DIRX=os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/Glaciers/3D/"+glacier)
  inputs=os.path.join(DIRM+"/inputs/")

  extent = np.loadtxt(inputs+"mesh_extent.dat")
  hole1 = np.loadtxt(inputs+"mesh_hole1.dat")
  hole2 = np.loadtxt(inputs+"mesh_hole2.dat")
  holes=[hole1,hole2]  

  if not(os.path.exists(DIRR)):
    os.makedirs(DIRR)

  # Boundary numbers 
  bbed=4
  bsurf=5
  runname=method+"_beta"

  #############################
  # Run inversion solver file #
  #############################

  print "\n## Running elmer inversion code ##\n"

  if not(os.path.exists(DIRR+"summary.dat")):
    fid_info = open(DIRR+"summary.dat","w")
    fid_info.write('Lambda Nsim Cost Norm RelPrec_G \n')

  #for filename in glob.glob(DIRM+"elmer/robin_beta*"):
  #  os.remove(filename)
  os.chdir(DIRM)
  fid1 = open(DIRS+method+'_beta.sif', 'r')
  fid2 = open(DIRM+method+'_beta_'+date+'.sif', 'w')

  lines=fid1.readlines()
  for line in lines:
    line=line.replace('Extruded Mesh Levels=Integer', 'Extruded Mesh Levels=Integer {}'.format(extrude))
    line=line.replace('$Lambda=1.0e10', '$Lambda={}'.format(regpar))
    fid2.write(line)
  fid1.close()
  fid2.close()
  del fid1, fid2

  fid3 = open(DIRM+'ELMERSOLVER_STARTINFO','w')
  fid3.write(method+'_beta_'+date+'.sif')
  fid3.close()

  call(["mpiexec","ElmerSolver_mpi"])
	
  #####################################
  # Write cost values to summary file #
  ##################################### 

  fid = open(DIRM+"cost_"+method+"_beta.dat","r")
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

  DIRR_lambda = DIRR+"lambda_"+regpar+"_"+date+"/"

  names = os.listdir(DIRM+"/mesh2d")
  os.chdir(DIRM+"/mesh2d")
  if not os.path.exists(DIRR_lambda):
    os.makedirs(DIRR_lambda)
  for name in names:
    if name.endswith('vtu') and name.startswith(method):
      os.rename(name,DIRR_lambda+name)

  bed = elmerreadlib.saveline_boundary(DIRM+"/mesh2d/",runname,bbed)
  surf = elmerreadlib.saveline_boundary(DIRM+"/mesh2d/",runname,bsurf)

  os.rename(DIRM+"/mesh2d/"+runname+".dat",DIRR_lambda+runname+".dat")
  os.rename(DIRM+"/mesh2d/"+runname+".dat.names",DIRR_lambda+runname+".dat.names")
  os.rename(DIRM+"M1QN3_"+method+"_beta.out",DIRR_lambda+"M1QN3_"+method+"_beta.out")
  os.rename(DIRM+"gradientnormadjoint_"+method+"_beta.dat",DIRR_lambda+"gradient_"+runname+".dat")
  os.rename(DIRM+"cost_"+method+"_beta.dat",DIRR_lambda+"cost_"+runname+".dat")

	################################
	# Output friction coefficients #
	################################

  # Linear Beta square
  fid = open(inputs+"beta_linear.xyz",'w')
  fid.write('{0}\n'.format(len(bed['node'])))
  for i in range(0,len(bed['node'])):
    fid.write('{0} {1} {2:.4f} {3}\n'.format(bed['coord1'][i],bed['coord2'][i],\
				bed['coord3'][i],bed['beta'][i]**2))
  fid.close()

  # Weertman coefficient
  fid = open(inputs+"beta_weertman.xyz",'w')
  fid.write('{0}\n'.format(len(bed['node'])))
  for i in range(0,len(bed['node'])):
    coeff=(bed['beta'][i]**2)*(bed['vel'][i]**(2.0/3.0))
    fid.write('{0} {1} {2:.4f} {3}\n'.format(bed['coord1'][i],bed['coord2'][i],bed['coord3'][i],coeff))
  fid.close() 

  xgrid,ygrid,taubgrid = elmerreadlib.grid3d(bed,'taub',holes,extent)
  plt.imshow(taubgrid*1e4,origin='lower',clim=[0,500])
  plt.xticks([])
  plt.yticks([])
  plt.colorbar()
  cb.set_label('Basal shear stress (kPa)')
  plt.savefig(DIRR+'lambda_'+regpar+'_'+date+'.pdf',format='PDF',dpi=400)
  plt.close()
	
if __name__ == "__main__":
  main()

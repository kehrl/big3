#!/usr/bin/python

'''
This script let's you choose the velocity record for the basal inversion, and then runs 
the inversion solver file for either the adjoint or robin methods.

LMK, UW, 06/12/2014
'''

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
import matplotlib, elmerrunlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

def get_arguments():

  # Get inputs to file
  parser = argparse.ArgumentParser()
  parser.add_argument("-glacier",dest="glacier",required = True, 
        help = "Name of glacier (Kanger or Helheim)")
  parser.add_argument("-mesh", dest="mesh", required = True,
        help = "Name of mesh directory") 
  #parser.add_argument("-front", dest="frontbc", required = True,
  #      help = "Calving front boundary condition (velocity or pressure).") 
  parser.add_argument("-n", dest="n", required = True,
        help = "Number of partitions.")
  parser.add_argument("-regpar", dest="regpar", required = False,
       default='1e10',help = "Regularization parameter.")
  parser.add_argument("-restartsolverfile",dest="restartfile",required = False,\
       default="none",help = "Name of restart solver file.")
  parser.add_argument("-restartposition",dest="restartposition",required = False,\
       default=0,type=int,help = "Restart position in results file (if applicable).")
  parser.add_argument("-temperature",dest="temperature",required  = False,\
       default=-10.0,help = "Use modeled or constant temperature.") 
  parser.add_argument("-itmax",dest="itmax",required = False,\
       default=1000,help = "Maximum number of steady state iterations.")


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
  #frontbc = str(args.frontbc)
  glacier = args.glacier
  restartfile = args.restartfile
  restartposition = args.restartposition
  temperature = args.temperature
  itmax = args.itmax

  # Directories
  DIRS=os.path.join(os.getenv("CODE_HOME"),"big3/modeling/solverfiles/ssa/")
  DIRM=os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/"+RES+"/")
  DIRR=os.path.join(DIRM+'mesh2d/inversion_adjoint/')
  inputs=os.path.join(DIRM+"/inputs/")

  extent = np.loadtxt(inputs+"mesh_extent.dat")
  try:
    hole1 = np.loadtxt(inputs+"mesh_hole1.dat")
    hole2 = np.loadtxt(inputs+"mesh_hole2.dat")
    holes=[hole1,hole2]  
  except:
    holes=[]
  
  if not(os.path.exists(DIRR)):
    os.makedirs(DIRR)

  # Boundary numbers 
  bbed=4
  bsurf=5
  runname="adjoint_beta_ssa"

#   if temperature == 'model':
#     temperature_text="""
#   Constant Temperature = Variable Coordinate 1, Coordinate 2
#     Real Procedure "USF_Init.so" "ModelTemperature" """
#   else:
#     try:
#       float(temperature)
#       temperature_text="""
#   Constant Temperature = Real """+str(temperature)
#     except:
#       sys.exit("Unknown temperature of "+temperature)

  ################################################
  # Compile fortran libraries for adjoint solver #
  ################################################
  
  # Check that fortran files are compiled
  DIRELMERLIB = os.path.join(os.getenv("CODE_HOME"),"big3/modeling/elmerlib/")
  readfiles = os.listdir(DIRELMERLIB+"AdjointSSA/")
  outputfile = DIRELMERLIB+"AdjointSSASolvers.so"

  allstring=''
  for f in readfiles:
    if f.endswith('.F90'):
      allstring = allstring+f+' '

  os.chdir(DIRELMERLIB+"AdjointSSA/")
  call(["elmerf90","-o",outputfile,allstring])
  del DIRELMERLIB


  #############################
  # Run inversion solver file #
  #############################

  if not(os.path.exists(DIRR+"summary.dat")):
    fid_info = open(DIRR+"summary.dat","w")
    fid_info.write('Lambda Nsim Cost Norm RelPrec_G \n')

  if restartfile !='none':
    solverfile = restartfile
    date = restartfile[-12:-4]
    
    # Move previous timesteps to output directory
    DIRR_lambda = DIRR+"lambda_"+regpar+"_"+date+"/"
    names = os.listdir(DIRM+"/mesh2d")
    if not os.path.exists(DIRR_lambda):
      os.makedirs(DIRR_lambda)
    for name in names:
      if name.endswith('vtu') and name.startswith('adjoint'):
        os.rename(DIRM+"/mesh2d/"+name,DIRR_lambda+name)
#     try:
#       os.rename(DIRM+"M1QN3_"+method+"_beta.out",DIRR_lambda+"M1QN3_"+method+"_beta_beforerestart.out")
#       os.rename(DIRM+"gradientnormadjoint_"+method+"_beta.dat",DIRR_lambda+"gradient_"+runname+"_beforerestart.dat")
#       os.rename(DIRM+"cost_"+method+"_beta.dat",DIRR_lambda+"cost_"+runname+"_beforerestart.dat")
#     except:
#       try:
#         os.rename(DIRR_lambda+"M1QN3_"+method+"_beta.out",DIRR_lambda+"M1QN3_"+method+"_beta_beforerestart.out")
#         os.rename(DIRR_lambda+"gradientnormadjoint_"+method+"_beta.dat",DIRR_lambda+"gradient_"+runname+"_beforerestart.dat")
#         os.rename(DIRR_lambda+"cost_"+method+"_beta.dat",DIRR_lambda+"cost_"+runname+"_beforerestart.dat")
#       except:
#         pass

  else: 
    # Get current date
    now = datetime.datetime.now()
    date = '{0}{1:02.0f}{2:02.0f}'.format((now.year),(now.month),(now.day))
    restartposition = 0
 
    solverfile = 'adjoint_beta_'+regpar+'_'+date+'.sif'
    os.chdir(DIRM)
    fid1 = open(DIRS+'adjoint_beta_ssa.sif', 'r')
    fid2 = open(DIRM+solverfile, 'w')

    lines=fid1.read()
    lines=lines.replace('{Lambda}', '{0}'.format(regpar))
    lines=lines.replace('{ItMax}', '{0}'.format(int(itmax)))
    #lines=lines.replace('{Temperature}', '{0}'.format(temperature_text))
    fid2.write(lines)
    fid1.close()
    fid2.close()
    del fid1, fid2

  returncode = elmerrunlib.run_elmer(DIRM+solverfile,n=partitions)
	
  #####################################
  # Write cost values to summary file #
  ##################################### 

#   fid = open(DIRM+"cost_"+method+"_beta.dat","r")
#   lines = fid.readlines()
#   line=lines[-1]
#   p=line.split()
#   nsim = float(p[0])
#   cost1 = float(p[1])
#   cost2 = float(p[2])
#   norm = float(p[3]) 
#   fid.close()
#   fid_info = open(DIRR+"summary.dat","a")
#   fid_info.write('{} {} {} {} {}\n'.format(regpar,nsim+restartposition,cost1,cost2,norm))
#   fid_info.close()
#   del fid

	#######################################
	# Combine elmer results into one file #
	#######################################

#   DIRR_lambda = DIRR+"lambda_"+regpar+"_"+date+"/"
# 
#   names = os.listdir(DIRM+"/mesh2d")
#   if not os.path.exists(DIRR_lambda):
#     os.makedirs(DIRR_lambda)
#   for name in names:
#     if name.endswith('pvtu') and name.startswith(method):
#       os.rename(DIRM+"/mesh2d/"+name,DIRR_lambda+'{0}{1:04d}{2}'.format(name[0:-9],int(name[-9:-5])+restartposition,'.pvtu'))
#     elif name.endswith('vtu') and name.startswith(method):
#       os.rename(DIRM+"/mesh2d/"+name,DIRR_lambda+'{0}{1:04d}{2}'.format(name[0:-8],int(name[-8:-4])+restartposition,'.vtu'))
#     elif name.startswith(method) and 'result' in name:
#       os.rename(DIRM+"/mesh2d/"+name,DIRR_lambda+name)

#   bed = elmerreadlib.saveline_boundary(DIRM+"/mesh2d/",runname,bbed,['velocity','beta'])
#   surf = elmerreadlib.saveline_boundary(DIRM+"/mesh2d/",runname,bsurf,['vsurfini','velocity'])
# 
#   # Move saveline results
#   files = os.listdir(DIRM+"/mesh2d/")
#   for file in files:
#     if file.startswith(runname) and not file.endswith('names') and ('.dat' in file):
#       os.rename(DIRM+"/mesh2d/"+file,DIRR_lambda+file)
#     if file.startswith(runname) and file.endswith('names'):
#       os.rename(DIRM+"/mesh2d/"+file,DIRR_lambda+file)
#   
#   # Move outputs for optimization
#   os.rename(DIRM+"M1QN3_"+method+"_beta.out",DIRR_lambda+"M1QN3_"+method+"_beta.out")
#   os.rename(DIRM+"gradientnormadjoint_"+method+"_beta.dat",DIRR_lambda+"gradient_"+runname+".dat")
#   os.rename(DIRM+"cost_"+method+"_beta.dat",DIRR_lambda+"cost_"+runname+".dat")
  
	################################
	# Output friction coefficients #
	################################

  # Linear Beta square
#   fid = open(inputs+"beta_linear.xyz",'w')
#   fid.write('{0}\n'.format(len(bed['Node Number'])))
#   for i in range(0,len(bed['Node Number'])):
#     fid.write('{0} {1} {2:.4f} {3}\n'.format(bed['x'][i],bed['y'][i],\
# 				bed['z'][i],bed['beta'][i]**2))
#   fid.close()
# 
#   # Weertman coefficient
#   fid = open(inputs+"beta_weertman.xyz",'w')
#   fid.write('{0}\n'.format(len(bed['Node Number'])))
#   for i in range(0,len(bed['Node Number'])):
#     coeff=(bed['beta'][i]**2)*(bed['velocity'][i]**(2.0/3.0))
#     fid.write('{0} {1} {2:.4f} {3}\n'.format(bed['x'][i],bed['y'][i],bed['z'][i],coeff))
#   fid.close() 
# 
#   xgrid,ygrid,taubgrid = elmerreadlib.grid3d(bed,'taub',holes,extent)
#   xgrid,ygrid,vmodgrid = elmerreadlib.grid3d(surf,'velocity',holes,extent)
#   xgrid,ygrid,vmesgrid = elmerreadlib.grid3d(surf,'vsurfini',holes,extent)
#   
#   plt.figure(figsize=(6.5,3))
#   plt.subplot(121)
#   plt.imshow(taubgrid*1e3,origin='lower',clim=[0,500])
#   plt.xticks([])
#   plt.yticks([])
#   cb = plt.colorbar(ticks=np.arange(0,600,100))
#   cb.ax.tick_params(labelsize=10)
#   cb.set_label('Basal shear stress (kPa)')
#   
#   plt.subplot(122)
#   plt.imshow(vmodgrid-vmesgrid,origin='lower',clim=[-500,500],cmap='RdBu_r')
#   plt.xticks([])
#   plt.yticks([])
#   cb = plt.colorbar(ticks=np.arange(-500,600,100))
#   cb.ax.tick_params(labelsize=10)
#   cb.set_label('Modeled-Measured (m/yr)')
#   
#   plt.tight_layout()
#   
#   plt.savefig(DIRR+'lambda_'+regpar+'_'+date+'.pdf',format='PDF',dpi=400)
#   plt.close()
# 	
if __name__ == "__main__":
  main()

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
import elmerreadlib, meshlib, flowparameterlib
import numpy as np
import argparse
import datetime
import scipy.interpolate
import matplotlib, elmerrunlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

def get_arguments():

  # Get inputs to file
  parser = argparse.ArgumentParser()
  parser.add_argument("-glacier",dest="glacier",required = True, 
        help = "Name of glacier (Kanger or Helheim).")
  parser.add_argument("-mesh", dest="mesh", required = True,
        help = "Name of mesh directory.") 
  parser.add_argument("-front", dest="frontbc", required = True,
        help = "Calving front boundary condition (velocity or pressure).") 
  parser.add_argument("-n", dest="n", required = True,
        help = "Number of partitions.")
  parser.add_argument("-regpar", dest="regpar", required = False,
       default='1e10',help = "Regularization parameter. Default is 1e10.")
  parser.add_argument("-restartsolverfile",dest="restartfile",required = False,\
       default="none",help = "Name of restart solver file. Default is none.")
  parser.add_argument("-restartposition",dest="restartposition",required = False,\
       default=0,type=int,help = "Restart position in results file (if applicable).")
  parser.add_argument("-temperature",dest="temperature",required  = False,\
       default=-10.0,help = "Ice temperature in deg C (or 'model'). Default is -10.") 
  parser.add_argument("-itmax",dest="itmax",required = False,\
       default=500,help = "Maximum number of steady state iterations. Default is 500.")


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
  frontbc = str(args.frontbc)
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

  if temperature == 'model':
     temperature_text="""
  mu = Variable Coordinate 1, Coordinate 2\r
    Real Procedure "USF_Init.so" "SSAViscosity"""""
  else:
    try:
      print float(temperature)
      A = flowparameterlib.arrhenius(273.15+float(temperature))
      E = 3
      temperature_text="""
  mu = Real $(("""+'2*'+str(E)+'*'+str(A[0])+'*'+'yearinsec)^(-1.0/3.0)*1.0e-6)'
    except:
      sys.exit("Unknown temperature of "+temperature)

  if frontbc == 'velocity' or frontbc == 'dirichlet':
    frontbc_text = """
  SSAVelocity 1= Equals vsurfini 1\r
  SSAVelocity 2= Equals vsurfini 2\r
  Adjoint 1 = Real 0.0\r
  Adjoint 2 = Real 0.0"""
  elif frontbc == 'pressure' or frontbc == 'neumann':
    frontbc_text = """
  Calving Front = Logical True"""

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
    lines=lines.replace('{Temperature}', '{0}'.format(temperature_text))
    lines=lines.replace('{FrontBC}', '{0}'.format(frontbc_text))
    fid2.write(lines)
    fid1.close()
    fid2.close()
    del fid1, fid2

  returncode = elmerrunlib.run_elmer(DIRM+solverfile,n=partitions)
	
  #####################################
  # Write cost values to summary file #
  ##################################### 

  fidcost = open(DIRM+"cost.dat","r")
  lines = fidcost.readlines()
  line=lines[-1]
  p=line.split()
  nsim = float(p[0])
  costsur = float(p[1])
  fidcost.close()
  
  fidcostreg = open(DIRM+"costreg.dat")
  lines = fidcostreg.readlines()
  line=lines[-1]
  p=line.split()
  nsim = float(p[0])
  costbed = float(p[1])#/float(regpar)
  fidcostreg.close()

  costtot = costsur+float(regpar)*costbed  

  fid_info = open(DIRR+"summary.dat","a")
  fid_info.write('{} {} {} {} {}\n'.format(regpar,nsim+restartposition,costtot,costsur,costbed))
  fid_info.close()
  del fidcost, fidcostreg

  #######################################
  # Combine elmer results into one file #
  #######################################

  DIRR_lambda = DIRR+"lambda_"+regpar+"_"+date+"/"
 
  names = os.listdir(DIRM+"/mesh2d")
  if not os.path.exists(DIRR_lambda):
    os.makedirs(DIRR_lambda)
  for name in names:
    if name.endswith('pvtu') and name.startswith('adjoint'):
      os.rename(DIRM+"/mesh2d/"+name,DIRR_lambda+'{0}{1:04d}{2}'.format(name[0:-9],int(name[-9:-5])+restartposition,'.pvtu'))
    elif name.endswith('vtu') and name.startswith('adjoint'):
      os.rename(DIRM+"/mesh2d/"+name,DIRR_lambda+'{0}{1:04d}{2}'.format(name[0:-8],int(name[-8:-4])+restartposition,'.vtu'))
    elif name.startswith('adjoint') and 'result' in name:
      os.rename(DIRM+"/mesh2d/"+name,DIRR_lambda+name)
  
  bed = elmerreadlib.pvtu_file(DIRR_lambda+'adjoint_beta_ssa0001.pvtu',['ssavelocity','beta','vsurfini'])
   
  # Move outputs for optimization
  os.rename(DIRM+"M1QN3_adjoint_beta_ssa.out",DIRR_lambda+"M1QN3_adjoint_beta_ssa.out")
  os.rename(DIRM+"cost.dat",DIRR_lambda+"cost.dat")
  os.rename(DIRM+"costreg.dat",DIRR_lambda+"costreg.dat")
  os.rename(DIRM+"gradientnormadjoint_adjoint_beta_ssa.dat",DIRR_lambda+"gradientnormadjoint_adjoint_beta_ssa.dat")

  ################################
  # Output friction coefficients #
  ################################

  # Gridded linear beta square
  x,y,u = elmerreadlib.input_file(inputs+"udem.xy", dim=2)
  xx,yy = np.meshgrid(x,y)
  beta_linear = scipy.interpolate.griddata((bed['x'],bed['y']),bed['beta']**2,(xx,yy),\
      method='nearest')
  beta_linear_lin = scipy.interpolate.griddata((bed['x'],bed['y']),bed['beta']**2,(xx,yy),\
      method='linear')
  beta_weertman_lin = scipy.interpolate.griddata((bed['x'],bed['y']),\
      (bed['beta']**2)/(bed['ssavelocity']**(-2.0/3.0)), (xx,yy), method='linear')
  beta_weertman = scipy.interpolate.griddata((bed['x'],bed['y']),\
      (bed['beta']**2)*(bed['ssavelocity']**(-2.0/3.0)),(xx,yy),method='nearest')
  ind = np.where(~(np.isnan(beta_linear_lin)))
  beta_linear[ind] = beta_linear_lin[ind]
  ind = np.where(~(np.isnan(beta_weertman_lin)))
  beta_weertman[ind] = beta_weertman_lin[ind]  
  del beta_linear_lin, beta_weertman_lin
 
  # Output original results
  fidr = open(inputs+"beta_linear.dat",'w')
  fidw = open(inputs+"beta_weertman.dat",'w')
  fidr.write('{}\n'.format(len(bed['x'])))
  fidw.write('{}\n'.format(len(bed['x'])))
  for i in range(0,len(bed['x'])):
    fidr.write('{0} {1} {2:.6f}\n'.format(bed['x'][i],bed['y'][i],bed['beta'][i]**2))
    fidw.write('{0} {1} {2:.6f}\n'.format(bed['x'][i],bed['y'][i],(bed['beta'][i]**2)/(bed['ssavelocity'][i]**(-2.0/3.0))))
  fidr.close()
  fidw.close()
  del fidr, fidw
  
  fidl = open(inputs+"beta_linear.xy",'w')
  fidw = open(inputs+"beta_weertman.xy",'w')
  fidl.write('{}\n{}\n'.format(len(x),len(y)))
  fidw.write('{}\n{}\n'.format(len(x),len(y)))
  for i in range(0,len(x)):
    for j in range(0,len(y)):
      fidl.write('{0} {1} {2}\n'.format(x[i],y[j],beta_linear[j,i]))
      fidw.write('{0} {1} {2}\n'.format(x[i],y[j],beta_weertman[j,i]))
  fidl.close()
  fidw.close()
  del beta_linear, beta_weertman, u, x, y, fidl, fidw, xx, yy
  
  xgrid,ygrid,taubgrid = elmerreadlib.grid3d(bed,'taub',holes,extent)
  xgrid,ygrid,vmodgrid = elmerreadlib.grid3d(bed,'ssavelocity',holes,extent)
  xgrid,ygrid,vmesgrid = elmerreadlib.grid3d(bed,'vsurfini',holes,extent)
   
  plt.figure(figsize=(6.5,3))
  plt.subplot(121)
  plt.imshow(taubgrid*1e3,origin='lower',clim=[0,500])
  plt.xticks([])
  plt.yticks([])
  cb = plt.colorbar(ticks=np.arange(0,600,100))
  cb.ax.tick_params(labelsize=10)
  cb.set_label('Basal shear stress (kPa)')
   
  plt.subplot(122)
  plt.imshow(vmodgrid-vmesgrid,origin='lower',clim=[-200,200],cmap='RdBu_r')
  plt.xticks([])
  plt.yticks([])
  cb = plt.colorbar(ticks=np.arange(-500,600,100))
  cb.ax.tick_params(labelsize=10)
  cb.set_label('Modeled-Measured (m/yr)')
   
  plt.tight_layout()
   
  plt.savefig(DIRR+'lambda_'+regpar+'_'+date+'.pdf',format='PDF',dpi=400)
  plt.close()
 	
if __name__ == "__main__":
  main()

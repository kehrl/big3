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
  parser.add_argument("-glacier",dest="glacier",required = True,\
        help = "Name of glacier (Kanger or Helheim).")
  parser.add_argument("-mesh", dest="mesh", required = True,\
        help = "Name of mesh directory.") 
  parser.add_argument("-n", dest="n", required = True,\
        help = "Number of partitions.")
  parser.add_argument("-front", dest="frontbc", required = True,\
        help = "Calving front boundary condition (velocity or pressure).")
  parser.add_argument("-regpar", dest="regpar", required = False,
       default='1e10',help = "Regularization parameter. Default is 1e10.")
  parser.add_argument("-method", dest="method", required = False,
       default='adjoint',help = "Adjoint or robin. Default is adjoint.")
  parser.add_argument("-extrude", dest="extrude", type=int,required = False,\
       default=10,help = "Number of extrusion levels. Default is 10.")
  parser.add_argument("-sidewall", dest="sidewallbc", required = False,\
        default='velocity',help = "Sidewall boundary condition (velocity or friction). Default is velocity.") 
  parser.add_argument("-slipcoefficient", dest="slipcoefficient", required = False,
        default='1.0E-3',help = "Sidewall boundary slip coefficient. Default is 1.0E-3.")   
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
  method = args.method
  extrude = str(args.extrude)
  frontbc = str(args.frontbc)
  sidewallbc = str(args.sidewallbc)
  slipcoefficient = args.slipcoefficient
  glacier = args.glacier
  restartfile = args.restartfile
  restartposition = args.restartposition
  temperature = args.temperature
  itmax = args.itmax

  # Directories
  DIRS=os.path.join(os.getenv("CODE_HOME"),"big3/modeling/solverfiles/3D/")
  DIRM=os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/"+RES+"/")
  DIRR=os.path.join(DIRM+'mesh2d/inversion_'+method+'/')
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
  runname=method+"_beta"

  # Grab boundary condition for front -- it will be different if we are using an actual
  # terminus position (neumann, pressure BC) or an outflow boundary (dirichlet, velocity BC)
  if method=='adjoint':
    if (frontbc == 'neumann') or (frontbc == 'pressure'):
      frontbc_text ="""
  Adjoint Force BC = Logical True
  Flow Force BC = Logical True
  External Pressure = Variable Coordinate 3 !we are in MPa units
  Real MATC "-1.0*waterpressure(tx)*1.0E-06" """
    elif (frontbc == 'dirichlet') or (frontbc == 'velocity'):
      frontbc_text = """
  Velocity 1 = Variable Coordinate 1
    Real procedure "USF_Init.so" "UIni"
  Velocity 2 = Variable Coordinate 1
    Real procedure "USF_Init.so" "VIni"
  ! Dirichlet BC => Dirichlet = 0 for Adjoint
  Adjoint 1 = Real 0.0
  Adjoint 2 = Real 0.0"""
    else:
      sys.exit("Unknown BC for front of glacier.")
  elif method == 'robin':
    if (frontbc == 'neumann') or (frontbc == 'pressure'):
      frontbc_text = """
  Flow Force BC = Logical True
  External Pressure = Variable Coordinate 3 !we are in MPa units
  Real MATC "-1.0*waterpressure(tx)*1.0E-06" """
    elif (frontbc == 'dirichlet') or (frontbc == 'velocity'):
      frontbc_text ="""
  ! Dirichlet BCs
  Velocity 1 = Variable Coordinate 1, Coordinate 2
    Real procedure "USF_Init.so" "UIni"
  Velocity 2 = Variable Coordinate 1, Coordinate 2
    Real procedure "USF_Init.so" "VIni"
        
  ! Dirichlet BC => Same Dirichlet
  VeloD 1 = Variable Coordinate 1, Coordinate 2
    Real procedure "USF_Init.so" "UIni"
  VeloD 2 = Variable Coordinate 1, Coordinate 2
    Real procedure "USF_Init.so" "VIni" """

  # Set boundary condition for sidewalls (either measured velocity or slip coefficient).
  if method == 'adjoint':
    if sidewallbc == 'friction':
      sidewallbc_text = """
  Normal-Tangential Velocity = Logical True
  Normal-Tangential Adjoint = Logical True
  
  Adjoint Force BC = Logical True
  
  Velocity 1 = Real 0.0e0
  Adjoint 1 = Real 0.0e0
  
  Slip Coefficient 2 = Real """+slipcoefficient+"""
  Slip Coefficient 3 = Real """+slipcoefficient
    elif sidewallbc == 'freeslip':
      sidewallbc_text = """
  Normal-Tangential Velocity = Logical True
  Normal-Tangential Adjoint = Logical True

  Adjoint Force BC = Logical True

  Velocity 1 = Real 0.0e0
  Adjoint 1 = Real 0.0e0"""
    elif sidewallbc == 'velocity':
      sidewallbc_text = """
  !! Dirichlet BC 
  Velocity 1 = Variable Coordinate 1
    Real procedure "USF_Init.so" "UWa"
  Velocity 2 = Variable Coordinate 1
    Real procedure "USF_Init.so" "VWa"
  
  ! Dirichlet BC => Dirichlet = 0 for Adjoint
  Adjoint 1 = Real 0.0
  Adjoint 2 = Real 0.0"""
    else:
      sys.exit("Unknown sidewall BC of "+sidewallbc)
  if method == 'robin':
    if sidewallbc == 'friction':
      sidewallbc_text = """
  Normal-Tangential Velocity = Logical True
  Normal-Tangential VeloD = Logical True

  Flow Force BC = Logical True

  Velocity 1 = Real 0.0
  VeloD 1 = Real 0.0e0
  Slip Coefficient 2 = Real """+slipcoefficient+"""
  Slip Coefficient 3 = Real """+slipcoefficient
    elif sidewallbc == 'freeslip':
      sidewallbc_text = """
  Normal-Tangential Velocity = Logical True
  Normal-Tangential VeloD = Logical True

  Flow Force BC = Logical True

  Velocity 1 = Real 0.0
  VeloD 1 = Real 0.0e0"""
    elif sidewallbc == 'velocity':
      sidewallbc_text = """
  ! Dirichlet BCs
  Velocity 1 = Variable Coordinate 1
    Real procedure "USF_Init.so" "UWa"
  Velocity 2 = Variable Coordinate 1
    Real procedure "USF_Init.so" "VWa"

  ! Dirichlet BC => Same Dirichlet
  VeloD 1 = Variable Coordinate 1
    Real procedure "USF_Init.so" "UWa"
  VeloD 2 = Variable Coordinate 1
    Real procedure "USF_Init.so" "VWa" """
    else:
      sys.exit("Unknown sidewall BC of "+sidewallbc)    
      
  if temperature == 'model':
    temperature_text="""
  Viscosity = Variable Coordinate 1, Coordinate 2
    Real Procedure "USF_Init.so" "ModelViscosity"
  Constant Temperature = Variable Coordinate 1, Coordinate 2
    Real Procedure "USF_Init.so" "ModelTemperature" """
  else:
    try:
      A = flowparameterlib.arrhenius(273.15+float(temperature))
      E = 3
      temperature_text="""
  Viscosity = Real $(("""+'2*'+str(E)+'*'+str(A[0])+'*'+'yearinsec)^(-1.0/3.0)*1.0e-6)\r'+"""
  Constant Temperature = Real """+str(temperature)
    except:
      sys.exit("Unknown temperature of "+temperature)

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
      if name.endswith('vtu') and name.startswith(method):
        os.rename(DIRM+"/mesh2d/"+name,DIRR_lambda+name)
    try:
      os.rename(DIRM+"M1QN3_"+method+"_beta.out",DIRR_lambda+"M1QN3_"+method+"_beta_beforerestart.out")
      os.rename(DIRM+"gradientnormadjoint_"+method+"_beta.dat",DIRR_lambda+"gradient_"+runname+"_beforerestart.dat")
      os.rename(DIRM+"cost_"+method+"_beta.dat",DIRR_lambda+"cost_"+runname+"_beforerestart.dat")
    except:
      try:
        os.rename(DIRR_lambda+"M1QN3_"+method+"_beta.out",DIRR_lambda+"M1QN3_"+method+"_beta_beforerestart.out")
        os.rename(DIRR_lambda+"gradientnormadjoint_"+method+"_beta.dat",DIRR_lambda+"gradient_"+runname+"_beforerestart.dat")
        os.rename(DIRR_lambda+"cost_"+method+"_beta.dat",DIRR_lambda+"cost_"+runname+"_beforerestart.dat")
      except:
        pass

  else: 
    # Get current date
    now = datetime.datetime.now()
    date = '{0}{1:02.0f}{2:02.0f}'.format((now.year),(now.month),(now.day))
    restartposition = 0
 
    solverfile = method+'_beta_'+regpar+'_'+date+'.sif'
    os.chdir(DIRM)
    fid1 = open(DIRS+method+'_beta_grounded.sif', 'r')
    fid2 = open(DIRM+solverfile, 'w')

    lines=fid1.read()
    lines=lines.replace('{Extrude}', '{0}'.format(extrude))
    lines=lines.replace('{Lambda}', '{0}'.format(regpar))
    lines=lines.replace('{FrontBC}', '{0}'.format(frontbc_text))
    lines=lines.replace('{SidewallBC}', '{0}'.format(sidewallbc_text))
    lines=lines.replace('{ItMax}', '{0}'.format(int(itmax)))
    lines=lines.replace('{Temperature}', '{0}'.format(temperature_text))
    fid2.write(lines + os.linesep)
    fid1.close()
    fid2.close()
    del fid1, fid2

  returncode = elmerrunlib.run_elmer(DIRM+solverfile,n=partitions)
	
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
  fid_info.write('{} {} {} {} {}\n'.format(regpar,nsim+restartposition,cost1,cost2,norm))
  fid_info.close()
  del fid

  #######################################
  # Move results to one directory #
  #######################################

  DIRR_lambda = DIRR+"lambda_"+regpar+"_"+date+"/"

  # Something wrong with saveline boundary, so moving to reading pvtu_files
  #bed = elmerreadlib.saveline_boundary(DIRM+"/mesh2d/",runname,bbed,['velocity 1','velocity 2','velocity 3','beta'])
  #surf = elmerreadlib.saveline_boundary(DIRM+"/mesh2d/",runname,bsurf,['vsurfini 1','vsurfini 2','velocity 1','velocity 2','velocity 3'])
  data = elmerreadlib.pvtu_file(DIRM+"/mesh2d/"+runname+"0001.pvtu",['beta','velocity','vsurfini 1','vsurfini 2'])
  surf = elmerreadlib.values_in_layer(data,'surf')
  bed = elmerreadlib.values_in_layer(data,'bed')

  names = os.listdir(DIRM+"/mesh2d")
  if not os.path.exists(DIRR_lambda):
    os.makedirs(DIRR_lambda)
  for name in names:
    if name.endswith('pvtu') and name.startswith(method):
      os.rename(DIRM+"/mesh2d/"+name,DIRR_lambda+'{0}{1:04d}{2}'.format(name[0:-9],int(name[-9:-5])+restartposition,'.pvtu'))
    elif name.endswith('vtu') and name.startswith(method):
      os.rename(DIRM+"/mesh2d/"+name,DIRR_lambda+'{0}{1:04d}{2}'.format(name[0:-8],int(name[-8:-4])+restartposition,'.vtu'))
    elif name.startswith(method) and 'result' in name:
      os.rename(DIRM+"/mesh2d/"+name,DIRR_lambda+name)

  # Move saveline results
  files = os.listdir(DIRM+"/mesh2d/")
  for file in files:
    if file.startswith(runname) and ('.dat' in file):
      os.rename(DIRM+"/mesh2d/"+file,DIRR_lambda+file)
  
  # Move outputs for optimization
  os.rename(DIRM+"M1QN3_"+method+"_beta.out",DIRR_lambda+"M1QN3_"+method+"_beta.out")
  os.rename(DIRM+"gradientnormadjoint_"+method+"_beta.dat",DIRR_lambda+"gradient_"+runname+".dat")
  os.rename(DIRM+"cost_"+method+"_beta.dat",DIRR_lambda+"cost_"+runname+".dat")
  
  ################################
  # Output friction coefficients #
  ################################

  # Gridded linear beta square
  x,y,u = elmerreadlib.input_file(inputs+"udem.xy", dim=2)
  xx,yy = np.meshgrid(x,y)
  beta_linear = scipy.interpolate.griddata((bed['x'],bed['y']),\
    bed['beta']**2, (xx,yy), method='nearest')
  beta_linear_lin = scipy.interpolate.griddata((bed['x'],bed['y']),\
    bed['beta']**2, (xx,yy), method='linear')
  beta_weertman = scipy.interpolate.griddata((bed['x'],bed['y']),\
      (bed['beta']**2)/(bed['velocity']**(-2.0/3.0)), (xx,yy), method='nearest')
  beta_weertman_lin = scipy.interpolate.griddata((bed['x'],bed['y']),\
      (bed['beta']**2)/(bed['velocity']**(-2.0/3.0)), (xx,yy), method='linear')
  #if glacier == 'Helheim':
  #  # To get rid of edge effects near terminus, we take an average beta from farther upstream
  #  # and use that near the terminus
  #  ind_ave = np.where((xx > 302000) & (xx < 306000) & (yy < -2576000) & (yy > -2578000))
  #  ind = np.where((xx > 306000) & (yy < -2572000) & (yy > -2583000))
  #  beta_linear_lin[ind] = np.mean(beta_linear_lin[ind_ave])
  #  beta_weertman_lin[ind] = np.mean(beta_weertman_lin[ind_ave])
  ind = np.where(~(np.isnan(beta_linear_lin)))
  beta_linear[ind] = beta_linear_lin[ind]
  ind = np.where(~(np.isnan(beta_weertman_lin)))
  beta_weertman[ind] = beta_weertman_lin[ind]
  del beta_linear_lin, beta_weertman_lin

  # Output original results
  fidr = open(inputs+"beta_linear_"+regpar+"_FS_"+RES+".dat",'w')
  fidw = open(inputs+"beta_weertman_"+regpar+"_FS_"+RES+".dat",'w')  
  fidt = open(inputs+"taub_"+regpar+"_FS_"+RES+".dat",'w')
  fidr.write('{}\n'.format(len(bed['x'])))
  fidw.write('{}\n'.format(len(bed['x'])))
  fidt.write('{}\n'.format(len(bed['x'])))
  for i in range(0,len(bed['x'])):
    fidr.write('{0} {1} {2:.16f}\n'.format(bed['x'][i],bed['y'][i],bed['beta'][i]**2))  
    fidw.write('{0} {1} {2:.16f}\n'.format(bed['x'][i],bed['y'][i],(bed['beta'][i]**2)/(bed['velocity'][i]**(-2.0/3.0))))
    fidt.write('{0} {1} {2:.16f}\n'.format(bed['x'][i],bed['y'][i],bed['taub'][i]*1e3))
  fidr.close()
  fidw.close()
  fidt.close()
  del fidr, fidw, fidt

  # Output gridded results
  fidl = open(inputs+"beta_linear.xy",'w')
  fidw = open(inputs+"beta_weertman.xy",'w')
  fidl.write('{}\n{}\n'.format(len(x),len(y)))
  fidw.write('{}\n{}\n'.format(len(x),len(y)))
  for i in range(0,len(x)):
    for j in range(0,len(y)):
      fidl.write('{0} {1} {2:.6f}\n'.format(x[i],y[j],beta_linear[j,i]))
      fidw.write('{0} {1} {2:.6f}\n'.format(x[i],y[j],beta_weertman[j,i]))
  fidl.close()
  fidw.close()
  del beta_linear, beta_weertman, u, x, y, fidl, fidw, xx, yy

  xgrid,ygrid,taubgrid = elmerreadlib.grid3d(bed,'taub',holes,extent)
  xgrid,ygrid,vmodgrid = elmerreadlib.grid3d(surf,'velocity',holes,extent)
  xgrid,ygrid,vmesgrid1 = elmerreadlib.grid3d(surf,'vsurfini 1',holes,extent)
  xgrid,ygrid,vmesgrid2 = elmerreadlib.grid3d(surf,'vsurfini 2',holes,extent)

  plt.figure(figsize=(6.5,3))
  plt.subplot(121)
  plt.imshow(taubgrid*1e3,origin='lower',clim=[0,500])
  plt.xticks([])
  plt.yticks([])
  cb = plt.colorbar(ticks=np.arange(0,600,100))
  cb.ax.tick_params(labelsize=10)
  cb.set_label('Basal shear stress (kPa)')
  
  plt.subplot(122)
  plt.imshow(vmodgrid-np.sqrt(vmesgrid1**2+vmesgrid2**2),origin='lower',clim=[-200,200],cmap='RdBu_r')
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

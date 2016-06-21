import os
import shutil
import sys
import numpy as np
from subprocess import call
import math
import glob
import numpy as np
import meshlib
import numpy as np
import argparse
import datetime
import matplotlib, elmerrunlib, elmerreadlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

def get_arguments():

  # Get inputs to file
  parser = argparse.ArgumentParser()
  parser.add_argument("-glacier",dest="glacier",required = True, 
        help = "Name of glacier (Kanger or Helheim)")
  parser.add_argument("-mesh", dest="mesh", required = True,
        help = "Name of mesh") 
  parser.add_argument("-front", dest="frontbc", required = True,
        help = "Calving front boundary condition (velocity or pressure).") 
  parser.add_argument("-n", dest="n", required = True,
        help = "Number of partitions.")
  parser.add_argument("-s", dest="simulation", required = False,
        help = "Simulation type (steady or transient).",default='steady')
  parser.add_argument("-t", dest='iterations', required = False,
        help = "Number of steady state or transient iterations.", default=20)
  parser.add_argument("-dt", dest='timestepsize', required = False,
        help = "Time step size in units of year.", default=1)        
  parser.add_argument("-extrude", dest="extrude", type=int,required = False,\
       default=5,help = "Number of extrusion levels.")
  parser.add_argument("-restartfile",dest="restartfile",required = False,\
       default="none",help = "Name of restart file.")
  parser.add_argument("-restartposition",dest="restartposition",required = False,\
       default=0,type=int,\
       help = "Restart position in results file (if applicable.")  

  args, _ = parser.parse_known_args(sys.argv)


  return args

def main():

  ##########
  # inputs #
  ##########

  args = get_arguments()

  RES = args.mesh
  partitions = str(args.n)
  extrude = str(args.extrude)
  frontbc = str(args.frontbc)
  glacier = args.glacier
  simulation = args.simulation
  iterations = args.iterations
  dt = args.timestepsize
  restartfile = args.restartfile
  restartposition = args.restartposition

  # Directories
  DIRS=os.path.join(os.getenv("CODE_HOME"),"big3/modeling/solverfiles/3D/")
  DIRM=os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/"+RES+"/")
  DIRR=os.path.join(DIRM+'mesh2d/temperature/')
  DIRX=os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/Glaciers/3D/"+glacier)
  inputs=os.path.join(DIRM+"/inputs/")
  
  if not(os.path.exists(DIRR)):
    os.makedirs(DIRR)

  # Boundary numbers 
  bbed=4
  bsurf=5
  runname="temperature"

  # Grab boundary condition for front -- it will be different if we are using an actual
  # terminus position (neumann, pressure BC) or an outflow boundary (dirichlet, velocity BC)
  if (frontbc == 'neumann') or (frontbc == 'pressure'):
    frontbc_text ="""
  Flow Force BC = Logical True
  External Pressure = Variable Coordinate 3 !we are in MPa units
  Real MATC "-1.0*waterpressure(tx)*1.0E-06" """
  elif (frontbc == 'dirichlet') or (frontbc == 'velocity'):
    frontbc_text = """
  Velocity 1 = Variable Coordinate 1
    Real procedure "USF_Init.so" "UWa"
  Velocity 2 = Variable Coordinate 1
    Real procedure "USF_Init.so" "VWa" """
  else:
    sys.exit("Unknown BC for front of glacier.")

  # Set up text for simulation type (steady or transient) and the number of 
  # timesteps or iterations we want the simulation to run for.
  if simulation == 'steady':
    simulation_text = """
  Simulation Type = Steady State
  Steady State Max Iterations = {0}
  Steady State Min Iterations = 1""".format(iterations)
  elif simulation == 'transient':
    simulation_text = """
  Simulation Type = Transient
  Timestepping Method = "BDF"
  BDF Order = 1
  Timestep Intervals = {0}     
  Timestep Sizes = {1}""".format(iterations,dt)
  
  # Add in restart file and restart position into simulation text if we are using it
  if restartfile != 'none':
    simulation_text = simulation_text+"""
  Restart File = {0}
  Restart Position = {1}""".format(restartfile,restartposition)
  
  ###############################
  # Run temperature solver file #
  ###############################

  # Get current date
  now = datetime.datetime.now()
  date = '{0}{1:02.0f}{2:02.0f}'.format((now.year),(now.month),(now.day))
 
  solverfile = runname+'_'+date+'.sif'
  os.chdir(DIRM)
  fid1 = open(DIRS+runname+'.sif', 'r')
  fid2 = open(DIRM+solverfile, 'w')

  lines=fid1.read()
  lines=lines.replace('{Simulation Type}', '{0}'.format(simulation_text))
  lines=lines.replace('{Extrude}', '{0}'.format(extrude))
  lines=lines.replace('{FrontBC}', '{0}'.format(frontbc_text))
  fid2.write(lines)
  fid1.close()
  fid2.close()
  del fid1, fid2

  returncode = elmerrunlib.run_elmer(DIRM+solverfile,n=partitions)
  
  #####################
  # Move output files #
  #####################
  
  DIRR_output = DIRR+"temperature_"+date+"/"
  
  names = os.listdir(DIRM+"/mesh2d")
  if not os.path.exists(DIRR_output):
    os.makedirs(DIRR_output)
  for name in names:
    if name.endswith('pvtu') and name.startswith("temperature"):
      os.rename(DIRM+"/mesh2d/"+name,DIRR_output+name)
    elif name.endswith('vtu') and name.startswith("temperature"):
      os.rename(DIRM+"/mesh2d/"+name,DIRR_output+name)
    elif name.startswith("temperature") and 'result' in name:
      os.rename(DIRM+"/mesh2d/"+name,DIRR_output+name)

  bed = elmerreadlib.saveline_boundary(DIRM+"/mesh2d/",runname,bbed)
  surf = elmerreadlib.saveline_boundary(DIRM+"/mesh2d/",runname,bsurf)

  os.rename(DIRM+"/mesh2d/"+runname+".dat",DIRR_output+runname+".dat")
  os.rename(DIRM+"/mesh2d/"+runname+".dat.names",DIRR_output+runname+".dat.names")
  
if __name__ == "__main__":
  main()

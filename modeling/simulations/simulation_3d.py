# This code runs a 3D Helheim simulation, using one of the solver files.
#
# LMK, UW, 10/16/2014

import os, sys, datetime
import argparse
import elmerrunlib

##########
# Inputs #
##########

# Get inputs to file
parser = argparse.ArgumentParser()
parser.add_argument("-glacier",dest="glacier",required = True, 
        help = "Name of glacier (Kanger or Helheim)")
parser.add_argument("-mesh", dest="mesh", required = True,
        help = "Name of mesh directory") 
parser.add_argument("-sif",dest="solverfile",required = True, 
        help = "SIF file to run")
parser.add_argument("-front", dest="frontbc", required = False,default='pressure',
        help = "Calving front boundary condition (velocity or pressure).") 
parser.add_argument("-n", dest="n", required = True,
        help = "Number of partitions.")
parser.add_argument("-extrude", dest="extrude", type=int,required = False,\
       default=10,help = "Number of extrusion levels.")
parser.add_argument("-restartfile",dest="restartfile",required = False,\
       default="none",help = "Name of restart solver file (if applicable).")
parser.add_argument("-restartposition",dest="restartposition",required = False,\
       default=0,type=int,help = "Restart position in results file (if applicable).")
parser.add_argument("-temperature",dest="temperature",required  = False,\
       default=-10.0,help = "Use modeled or constant temperature (-10.0 deg C).") 
parser.add_argument("-nt",dest="nt",required  = False,type=str,\
       default='10',help = "Number of timesteps (10).") 
parser.add_argument("-dt",dest="dt",required  = False,type=str,\
       default='1/365.25',help = "Timestep size (1/365.25, i.e., 1 day).") 
parser.add_argument("-slipcoefficient",dest="slipcoefficient",required  = False,type=str,\
       default='1.0E-3',help = "Sidewall slip coefficient.")
parser.add_argument("-slidinglaw",dest="slidinglaw",required  = False,type=str,\
       default='linear',help = "Sliding law (linear or weertman).")     

args, _ = parser.parse_known_args(sys.argv)

RES = args.mesh
partitions = int(args.n)
extrude = str(args.extrude)
frontbc = str(args.frontbc)
solverfile_in = args.solverfile
glacier = args.glacier
restartfile = args.restartfile
restartposition = args.restartposition
temperature = args.temperature
dt = args.dt
nt = args.nt
slipcoefficient = args.slipcoefficient
slidinglaw = args.slidinglaw

# Directories
DIRS=os.path.join(os.getenv("CODE_HOME"),"big3/modeling/solverfiles/3D/")
DIRM=os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/"+RES+"/")
DIRR=os.path.join(DIRM+"mesh2d/simulation/")
inputs=os.path.join(DIRM+"inputs/")

# Check that solver file exists and remove the '.sif' at the end, if it's present
if restartfile !='none':
  solverfile_in = restartfile
  date = restartfile[-12:-4]
  if solverfile_in.endswith('.sif'):
    solverfile_in = solverfile_in[0:-4]
else:
  if solverfile_in.endswith('.sif'):
    solverfile_in = solverfile_in[0:-4]
  if not(os.path.exists(DIRS+solverfile_in+'.sif')):
    sys.exit('No solverfile with name '+solverfile_in+'.sif')

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
    Real procedure "USF_Init.so" "VWa"""
else:
  sys.exit("Unknown BC for front of glacier.")

if temperature == 'model':
  temperature_text="""
  Constant Temperature = Variable Coordinate 1, Coordinate 2
    Real Procedure "USF_Init.so" "ModelTemperature" """
else:
  try:
    float(temperature)
    temperature_text="""
  Constant Temperature = Real """+str(temperature)
  except:
    sys.exit("Unknown temperature of "+temperature)

if slidinglaw == 'linear':
  beta_data_file = 'inputs/beta_linear.xy'
  sliding_exponent='1.0/1.0'
elif slidinglaw == 'weertman':
  beta_data_file = 'inputs/beta_weertman.xy'
  sliding_exponent='1.0/3.0'
else:
  sys.exit("Unknown sliding law of "+slidinglaw)
  
if restartposition > 0:
  print "To restart the model you will also need to:\n"+\
    "(1) move the results file to the desired mesh directory\n"+\
    "(2) set Exec Solver = Never in StructuredMeshMapper, and \n"+\
    "(3) set Restart Position = 1 and Restart File in the Simulation section\n"+\
    "Remember that we can only restart from the last timestep of a run because previous results files are overwritten"

###################################
# Change mesh file in solver file #
###################################

if restartfile !='none':
  solverfile_out = solverfile_in
else:
  # Get current date
  now = datetime.datetime.now()
  date = '{0}{1:02.0f}{2:02.0f}'.format((now.year),(now.month),(now.day))
 
  os.chdir(DIRM)
  fid1 = open(DIRS+solverfile_in+'.sif', 'r')
  solverfile_out = solverfile_in+'_'+date
  fid2 = open(DIRM+solverfile_out+'.sif', 'w')

  lines=fid1.read()
  lines=lines.replace('{Extrude}', '{0}'.format(extrude))
  lines=lines.replace('{FrontBC}', '{0}'.format(frontbc_text))
  lines=lines.replace('{Temperature}', '{0}'.format(temperature_text))
  lines=lines.replace('{TimeStepSize}', '$({0})'.format(dt))
  lines=lines.replace('{TimeSteps}', '{0}'.format(nt))
  lines=lines.replace('{SlipCoefficient}', '{0}'.format(slipcoefficient))
  lines=lines.replace('{BetaDataFile}', '{0}'.format(beta_data_file))
  lines=lines.replace('{SlidingExponent}', '{0}'.format(sliding_exponent))
  lines=lines.replace('{nrestart}','{0}'.format(restartposition))

  fid2.write(lines)
  fid1.close()
  fid2.close()
  del fid1, fid2

#####################################
# Start zip/unzip process if wanted #
#####################################


if solverfile_in.startswith('terminusdriven') or solverfile_in == 'checkmeshes':
  os.chdir(DIRM)
  
  for i in range(max(1,restartposition),25+restartposition):
    if not(os.path.isdir('mesh{0:04d}'.format(i))):
      if os.path.isfile('mesh{0:04d}.tar.gz'.format(i)):
        os.system('tar -xzf mesh{0:04d}.tar.gz'.format(i))
    if not(os.path.isfile(inputs+'smb{0:04d}.xy'.format(i))):
      os.system('tar -xf '+inputs+'smb.tar -C '+inputs+' smb{0:04d}.xy'.format(i))

  import logging
  logging.basicConfig()
  
  from apscheduler.schedulers.background import BackgroundScheduler
  sched = BackgroundScheduler()

  def check_files():
    import os
    itmax = int(1)
    files = os.listdir(DIRM+"mesh2d")
    for file in files:
      if (file.endswith('.pvtu')) and (int(file[-9:-5]) > itmax):
        itmax = int(file[-9:-5])
    # print "itmax = "+str(itmax)
      
    # First check to make sure we have adequate mesh and smb files
    
    for i in range(itmax,itmax+25):
      if not(os.path.isdir('mesh{0:04d}'.format(i))):
        if os.path.isfile('mesh{0:04d}.tar.gz'.format(i)):
          os.system('tar -xzf mesh{0:04d}.tar.gz'.format(i))
      if not(os.path.isfile(inputs+'smb{0:04d}.xy'.format(i))):
        os.system('tar -xf '+inputs+'smb.tar -C '+inputs+' smb{0:04d}.xy'.format(i))
     
    # Now remove old mesh files and zip vtu files
    for i in range(max(1,restartposition),itmax-1):
      if os.path.isdir('mesh{0:04d}'.format(i)):
        os.system('rm -r '+'mesh{0:04d}'.format(i))
      if os.path.isfile(inputs+'smb{0:04d}.xy'.format(i)):
        os.system('rm '+inputs+'smb{0:04d}.xy'.format(i))
      if os.path.isfile('mesh2d/terminusdriven{0:04d}.pvtu'.format(i)) and not(os.path.isfile('mesh2d/terminusdriven{0:04d}.pvtu.tar.gz'.format(i))):  
        os.chdir(DIRM+"mesh2d")
        os.system('tar -czf terminusdriven{0:04d}.pvtu.tar.gz '.format(i)+\
                'terminusdriven*{0:04d}.'.format(i)+'*vtu')
        os.system('rm terminusdriven*{0:04d}.'.format(i)+'*vtu')
        os.chdir(DIRM)
      
  job = sched.add_job(check_files,'interval',minutes=5)
  sched.start()

###################
# Run Elmersolver #
###################

returncode = elmerrunlib.run_elmer(DIRM+solverfile_out+'.sif',n=partitions,email=True)

if solverfile_in.startswith('terminusdriven') or (solverfile_in == 'checkmeshes'):

  job.remove()
  print "Stopped unzipping/zipping mesh files and vtu."
  
  # Remove all remaining mesh and tar remaining vtu files
  os.chdir(DIRM)
  for i in range(max(1,restartposition),int(nt)+25):
    if os.path.isdir('mesh{0:04d}'.format(i)):
      os.system('rm -r '+'mesh{0:04d}'.format(i))
    if os.path.isfile('mesh2d/terminusdriven{0:04d}.pvtu'.format(i)):
      os.chdir(DIRM+"mesh2d")
      os.system('tar -czf '+'terminusdriven{0:04d}.pvtu.tar.gz '.format(i)+\
                'terminusdriven*{0:04d}.'.format(i)+'*vtu')
      os.system('rm '+'terminusdriven*{0:04d}.'.format(i)+'*vtu')
      os.chdir(DIRM)

#!/usr/bin/python
 
import subprocess
import time
import os

#glacier = 'Kanger'
glacier = 'Helheim'

# Mesh geometry

meshshp = 'glacier_extent_basin_front'
extrude = 10
bname = 'smith'
bmodel = 'aniso'
bsmooth = '4'
#lc = '300 500 750 2500'
lc = '400 800 1500 5000'
dx = 200

if glacier == 'Helheim':
  date = '20140127'
elif glacier == 'Kanger':
  date = '20110213' 

# Inversion options:
method = 'adjoint'
if method == 'adjoint':
  regpar = '5e11'
elif method == 'robin':
  regpar = '1e10' 

# Temperature SIF options
temp_simulation = 'steady'
temp_iterations = 60
temp_timestep = 1
temperature = -10.0 # Constant temperature for initial inversion

# Options for PBS submission
queue = 'long'
model = 'ivy'
nparts = 160
ncpus = 20
runtime = '48:00:00'

if meshshp.endswith('nofront'):
  frontBC = 'pressure'
else:
  frontBC = 'velocity'

#################
# Generate mesh #
#################

# Output mesh name
meshname = 'BASIN'+date+'_'+regpar

# Create mesh
command = "python /u/lkehrl/Code/big3/modeling/meshing/"+\
          "mesh_3d.py -glacier {0} -mesh {1} -d {2} -bname {3} -bmodel {4} -bsmooth {5} -lc {6} -n {7} -output {8} -dx {9}".format(glacier,meshshp,date,bname,bmodel,bsmooth,lc,nparts,meshname,dx)
print command
os.system(command)

###########################
# Submit inversion to PBS #
###########################
    
# Customize job options
job_name = glacier+"_"+date+"_basin"
walltime = runtime
processors = "select={0}:ncpus={1}:mpiprocs={2}:model={3}".format(nparts/ncpus,ncpus,ncpus,model)
command = "python /u/lkehrl/Code/big3/modeling/inversions/inversion_3d.py"+\
          " -glacier {0} -method {1} -regpar {2} -mesh {3} -extrude {4} -front {5} -n {6} -temperature {7}".format(glacier,method,regpar,meshname,extrude,frontBC,nparts,temperature)
dir = "/nobackupp8/lkehrl/Models/"+glacier+"/3D/"+meshname+"/"

job_string = """
#PBS -S /bin/bash
#PBS -M kehrl@uw.edu
#PBS -m abe
#PBS -N %s
#PBS -l walltime=%s
#PBS -l %s
#PBS -o %s
#PBS -e %s
source /u/lkehrl/.profile
source /u/dlilien/sw/elmer/.bashrc_pleiades
cd %s
%s""" % (job_name, walltime, processors, dir+"PBS_"+method+"_"+regpar+".out",dir+"PBS_"+method+"_"+regpar+".err",dir,command)
     
print command
os.chdir(dir)
fid = open("PBS_"+method+"_"+regpar+".pbs","w")
fid.write(job_string)
fid.close()
try:
  subprocess.call(['qsub','-q',queue,'PBS_'+method+'_'+regpar+'.pbs'])
except:
  print "Couldn't submit job for %s for lambda=%s" % (meshname,regpar)
  
########################################################
# Set up temperature run for after inversion completes #
########################################################

# Customize job options
job_name = glacier+"_"+date+"_basin_temperature"
walltime = runtime
processors = "select={0}:ncpus={1}:mpiprocs={2}:model={3}".format(nparts/ncpus,ncpus,ncpus,model)
command = "python /u/lkehrl/Code/big3/modeling/inversions/temperature_3d.py"+\
          " -glacier {0} -s {1} -t {2} -dt {3} -mesh {4} -extrude {5} -front {6} -n {7}".format(glacier,temp_simulation,temp_iterations,temp_timestep,meshname,extrude,frontBC,nparts)
dir = "/nobackupp8/lkehrl/Models/"+glacier+"/3D/"+meshname+"/"

job_string = """
#PBS -S /bin/bash
#PBS -M kehrl@uw.edu
#PBS -m abe
#PBS -N %s
#PBS -l walltime=%s
#PBS -l %s
#PBS -o %s
#PBS -e %s
source /u/lkehrl/.profile
source /u/dlilien/sw/elmer/.bashrc_pleiades
cd %s
%s""" % (job_name, walltime, processors, dir+"PBS_temperature.out",dir+"PBS_temperature.err",dir,command)
     
print command
os.chdir(dir)
fid = open("PBS_temperature.pbs","w")
fid.write(job_string)
fid.close()

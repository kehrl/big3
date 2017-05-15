import subprocess
import time
import os

glacier = 'Helheim'

# Mesh characteristics and geometry
meshshp = 'glacier_extent_inversion_nofront'
extrude = 12
bname = 'smith'
bmodel = 'aniso'
bsmooth = '4'
bottomsurface = 'iceshelf'
temperature = -10.0
#lc = '500 1000 1000 2000'
lc = '150 250 500 700'

# Info terminus driven model
date1 = '20080613'
date2 = '20090613'
dt = '1/365.25'
timeseries = 'True'

runname = 'TD_'+date1+'_'+date2

# Inversion options
method = 'adjoint'
regpar = '5e11'
frontBC = 'pressure'
slipcoefficient = '1.0E-3'
sidewallBC = 'friction'

# Options for PBS submission
queue = 'long'
model = 'has'
nparts = 24
ncpus = 24
runtime = '30:00:00'

# Mesh directory
dir = "/nobackupp8/lkehrl/Models/"+glacier+"/3D/"+runname+"/"

#################
# Create meshes #
#################
command = "python /u/lkehrl/Code/big3/modeling/meshing/"+\
          "mesh_3d.py "+"-glacier {0} -mesh {1}".format(glacier,meshshp)+\
          " -bname {0} -bmodel {1} -bsmooth {2}".format(bname,bmodel,bsmooth)+\
          " -zb {0} -temperature {1} -lc {2}".format(bottomsurface,temperature,lc)+\
          " -timeseries {0} -dt {1} -d {2}".format(timeseries,dt,date1)+\
          " -d2 {0} -n {1} -output {2} ".format(date2,nparts,runname)            
print command
os.system(command)

###############################
# Create inversion PBS script #
###############################
command = "python /u/lkehrl/Code/big3/modeling/inversions/"+\
          "inversion_3d.py -glacier {0} -mesh {1}".format(glacier,runname)+\
          " -method {0} -regpar {1} -extrude {2}".format(method,regpar,extrude)+\
          " -front {0} -n {1} -sidewall {2}".format(frontBC,nparts,sidewallBC)+\
          " -slipcoefficient {0} -temperature {1}".format(slipcoefficient, temperature)

job_name = glacier+"_"+runname+"_inversion"
walltime = runtime
processors =     processors = "select={0}:ncpus={1}:mpiprocs={2}:model={3}".format(nparts/ncpus,ncpus,ncpus,model)

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
source /u/dlilien/sw/elmer/.bashrc_pleiades_haswell
cd %s
%s""" % (job_name, walltime, processors, dir+"PBS_inversion.out",dir+"PBS_inversion.err",dir,command)

print command  
os.chdir(dir)
fid = open("PBS_inversion.pbs","w")
fid.write(job_string)
fid.close()

####################################
# Create simulation PBS job script #
####################################

# Get number of timesteps from mesh files
nt = 0
files = os.listdir(dir)
for file in files:
  if file.startswith('mesh0') and not(file.endswith('.msh')):
    if int(file[-4:]) > nt:
      nt = int(file[-4:])
   
command = "python /u/lkehrl/Code/big3/modeling/simulations/"+\
          "simulation_3d.py -glacier {0} -sif terminusdriven.sif".format(glacier)+\
          " -mesh {0} -extrude {1} -n {2}".format(runname,extrude,nparts)+\
          " -dt {0} -nt {1} -temperature {2}".format(dt,nt,temperature) 

job_name = glacier+"_"+runname+"_terminusdriven"
walltime = runtime
processors = "select={0}:ncpus={1}:mpiprocs={2}:model={3}".format(nparts/ncpus,ncpus,ncpus,model)

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
source /u/dlilien/sw/elmer/.bashrc_pleiades_haswell
cd %s
%s""" % (job_name, walltime, processors, dir+"PBS_terminusdriven.out",dir+"PBS_terminusdriven.err",dir,command)

print command
os.chdir(dir)
fid = open("PBS_terminusdriven.pbs","w")
fid.write(job_string)
fid.close() 

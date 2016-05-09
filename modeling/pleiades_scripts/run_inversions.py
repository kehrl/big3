#!/usr/bin/python
 
from popen2 import popen2
import time

# Date for mesh
date = '20120624'

# Mesh geometry
glacier = 'Helheim'
meshshp = 'glacier_extent_inversion_front.shp'
extrude = 20
date = '20120624'
bname = 'smith'
bmodel = 'aniso'
bsmooth = '4'
lc = '2000 2000 2000 2000'

# Output mesh name
meshname = 'TEST'

# Inversion options
method = 'robin'
regpars = ['1e10'] 


# Options for 
model = 'ivy'
nparts = 60
ncpus = 20

if meshshp.endswith('nofront'):
  frontBC = 'pressure'
else:
  frontBC = 'velocity'

###############
# Create mesh #
###############

command = "python /u/lkehrl/Code/big3/modeling/meshing/"+\
  "helheim_mesh_3d.py -mesh {0} -d {1} -bname {2} -bmodel {3} -bsmooth {4} -lc {5} -n {6} -output {7}".format(meshshp,date,bname,bmodel,bsmooth,lc,nparts,meshname)

############################
# Submit inversions to PBS #
############################

# Loop over PBS jobs
for regpar in regpars:
 
    # Open a pipe to the qsub command.
    #output, input = popen2('qsub')
     
    # Customize your options here
    job_name = "lambda_%s" % regpar
    walltime = "6:00:00"
    processors = "select={0}:ncpus={1}:mpiprocs={2}:model={3}".format(nparts/ncpus,ncpus,ncpus,model)
    command = "python /u/lkehrl/Code/big3/modeling/inversions/helheim_inversion_3d.py"+\
              " -method {0} -regpar {1} -mesh {2} -extrude {3} -front {4}".format(method,regpar,meshname,extrude,frontBC)
    dir = "/nobackupp8/lkehrl/Models/"+glacier+"/3D/"+meshname
 
    job_string = """
    #PBS -S /bin/bash
    #PBS -M kehrl@uw.edu
    #PBS -m ab
    #PBS -N %s
    #PBS -l walltime=%s
    #PBS -l %s
    #PBS -o ./output/%s.out
    #PBS -e ./error/%s.err
    cd %s
    %s""" % (job_name, walltime, processors, job_name, job_name, command,dir)
     
    # Send job_string to qsub
    #input.write(job_string)
    #input.close()
     
    # Print your job and the system response to the screen as it's submitted
    #print job_string
    #print output.read()
     
    time.sleep(0.1)
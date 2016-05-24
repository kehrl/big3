#!/usr/bin/python
 
import subprocess
import time
import os

# Mesh geometry
#glacier = 'Kanger'
glacier = 'Helheim'
meshshp = 'glacier_extent_inversion_front'
extrude = 18
#date= '20120213' # Kanger
date = '20120316' # Helheim
bname = 'morlighem'
bmodel = 'aniso'
bsmooth = '4'
lc = '300 300 500 500'

# Output mesh name
meshname = 'DEM'+date+'_lowerres'

# Inversion options
method = 'robin'
regpars = ['1e10'] 


# Options for PBS submission
queue = 'long'
model = 'ivy'
nparts = 80
ncpus = 20
runtime = '12:00:00'

if meshshp.endswith('nofront'):
  frontBC = 'pressure'
else:
  frontBC = 'velocity'

###############
# Create mesh #
###############

command = "python /u/lkehrl/Code/big3/modeling/meshing/"+\
          "mesh_3d.py -glacier {0} -mesh {1} -d {2} -bname {3} -bmodel {4} -bsmooth {5} -lc {6} -n {7} -output {8}".format(glacier,meshshp,date,bname,bmodel,bsmooth,lc,nparts,meshname)

print command
os.system(command)

############################
# Submit inversions to PBS #
############################

# Loop over PBS jobs
for regpar in regpars:
 
    # Open a pipe to the qsub command.
    #output, input = popen2('qsub')
     
    # Customize your options here
    job_name = glacier+"_"+date+"_lambda"+regpar
    walltime = runtime
    processors = "select={0}:ncpus={1}:mpiprocs={2}:model={3}".format(nparts/ncpus,ncpus,ncpus,model)
    command = "python /u/lkehrl/Code/big3/modeling/inversions/inversion_3d.py"+\
              " -glacier {0} -method {1} -regpar {2} -mesh {3} -extrude {4} -front {5} -n {6}".format(glacier,method,regpar,meshname,extrude,frontBC,nparts)
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
      print "Couldn't submit job for regularization %s" % regpar 

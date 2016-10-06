#!/usr/bin/python
 
import subprocess
import time
import os

#glacier = 'Kanger'
glacier = 'Helheim'

# Mesh geometry

meshshp = 'glacier_extent_inversion_front'
extrude = 10
#bname = 'smith'
bname = 'morlighem'
bmodel = 'aniso'
bsmooth = '4'
bottomsurface = 'bed' # or 'iceshelf'
temperature = 'model'
lc = '300 300 300 300'
#lc = '1000 1000 4000 5000'

if glacier == 'Helheim':
  dates = ['20120316']
  #dates = ['20110319','20110615','20110828','20111116',\
  #         '20120316','20120624','20120908','20121205',\
  #         '20130209','20130508','20130804','20131031',\
  #         '20140127','20140509','20140731','20141016']
elif glacier == 'Kanger':
  dates = ['20120213']
  #dates = ['20110708','20110826','20111106',\
  #         '20120213','20120522','20121012','20121217',\
  #         '20130210','20130714','20131004','20131204',\
  #         '20140213']  

# Inversion options
method = 'adjoint'
regpars = ['5e11']
#regpars = ['1e8','1e9','1e10','1e11','5e11','1e12','1e13','1e14'] 


# Options for PBS submission
queue = 'normal'
model = 'ivy'
nparts = 60
ncpus = 20
runtime = '8:00:00'

if meshshp.endswith('nofront'):
  frontBC = 'pressure'
else:
  frontBC = 'velocity'

############################
# Submit inversions to PBS #
############################

# Loop over PBS jobs
for date in dates:

  # Output mesh name
  meshname = 'DEM'+date+'_variableA'

  # Create mesh
  command = "python /u/lkehrl/Code/big3/modeling/meshing/"+\
          "mesh_3d.py -glacier {0} -mesh {1} -d {2} -bname {3} -bmodel {4} -bsmooth {5} -lc {6} -n {7} -output {8} -zb {9} -temperature {10}".format(glacier,meshshp,date,bname,bmodel,bsmooth,lc,nparts,meshname,bottomsurface,temperature)
  print command
  os.system(command)

 
  for regpar in regpars:
    
    # Customize job options
    job_name = glacier+"_"+date+"_lambda"+regpar
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
    #try:
    #  subprocess.call(['qsub','-q',queue,'PBS_'+method+'_'+regpar+'.pbs'])
    #except:
    #  print "Couldn't submit job for %s for lambda=%s" % (meshname,regpar)

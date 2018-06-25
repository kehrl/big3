#!/usr/bin/python
 
import subprocess, shutil, os
import time

glacier = 'Kanger'
temperature = '-10.0'#'-10.0'#'model'
sliding_exponent = 1
beta_date = 'average'
regpar = '1e12'

# Other options
extrude = 12
DIRM = os.path.join(os.getenv("MODEL_HOME"),glacier+'/3D/')

# Get dates
if glacier == 'Helheim':
  #dates = ['20120316']
  dates = ['20040803','20050829','20060825','20070912',\
           '20080814','20100624',\
           '20110319','20110615','20110828','20111116',\
           '20120316','20120624','20120908','20121205',\
           '20130209','20130508','20130804','20131031',\
           '20140127','20140509','20140731']
elif glacier == 'Kanger':
  #dates = ['20010712','20030801','20050621','20050831']
  dates = ['20010712','20030801','20050621','20050831',\
           '20060505','20060708','20060919',\
           '20070728','20090817','20100603',\
           '20110308','20110708','20110826','20111106',\
           '20120213','20120522','20121012','20121217',\
           '20130210','20130714','20131004','20131204',\
           '20140213','20150808']  

# Options for PBS submission
queue = 'normal'
model = 'has'
nparts = 24
ncpus = 24
runtime = '0:20:00'

# Get sliding law string
if sliding_exponent == 1:
  slidinglaw = 'linear'
elif sliding_exponent == 3:
  slidinglaw = 'weertman'
else:
  slidinglaw = 'm'+str(int(sliding_exponent))

# Get temperature text
if temperature == 'model':
  temperature_text = 'modelT'
else:
  temperature_text = 'constantT'

# Set up beta suffix
if beta_date == 'average':
  beta_suffix = regpar+'_FS_average_2011_2016_'+temperature_text
else:
  beta_suffix = regpar+'_FS_DEM'+beta_date+'_'+temperature_text
# Get beta file name in case it needs to be copied from another directory
beta_file = DIRM+'DEM'+beta_date+'_'+temperature_text+'/inputs/beta_'+\
	slidinglaw+'_'+beta_suffix+'.dat'

############################
# Submit jobs to PBS #
############################

#Loop over PBS jobs
for date in dates:

  # Output mesh name
  meshname = 'DEM'+date+'_'+temperature_text

  if beta_date == '':
    beta_suffix = regpar+'_FS_DEM'+date+'_'+temperature_text
  else:
    if not(os.path.isfile(DIRM+meshname+'/inputs/beta_'+slidinglaw+'_'+beta_suffix+'.dat')):
      shutil.copy(beta_file,DIRM+meshname+'/inputs/')

  # Customize job options
  job_name = glacier+"_"+date+"_"+beta_suffix+"_"+slidinglaw
  walltime = runtime
  processors = "select={0}:ncpus={1}:mpiprocs={2}:model={3}".format(nparts/ncpus,ncpus,ncpus,model)
  command = "python /u/lkehrl/Code/big3/modeling/simulations/simulation_3d.py"+\
              " -glacier {0} -sif steady.sif -mesh {1} -extrude {2}".format(glacier,meshname,extrude)+\
	      " -front velocity -n {0} -temperature {1} -sidewall velocity".format(nparts,temperature)+\
	      " -beta {0} -slidinglaw {1}".format(beta_suffix,slidinglaw)
  dir = "/nobackupp8/lkehrl/Models/"+glacier+"/3D/"+meshname+"/"

  filename = "PBS_steady_"+beta_suffix+"_"+slidinglaw
  job_string = """
  #PBS -S /bin/bash
  #PBS -M kehrl@uw.edu
  #PBS -W group_list=s1877
  #PBS -m a
  #PBS -N %s
  #PBS -l walltime=%s
  #PBS -l %s
  #PBS -o %s
  #PBS -e %s
  source /u/lkehrl/.profile
  source /u/lkehrl/.bashrc_pleiades_haswell_sles12
  cd %s
  %s""" % (job_name, walltime, processors, dir+filename+".out",dir+filename+".err",dir,command)
     
  print command
  os.chdir(dir)
  fid = open(filename+".pbs","w")
  fid.write(job_string)
  fid.close()
  try:
    subprocess.call(['qsub','-q',queue,filename+".pbs"])
  except:
    print "Couldn't submit job for %s" % (filename)

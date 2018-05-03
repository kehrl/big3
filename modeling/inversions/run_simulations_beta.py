# This script will loop through models and run simulations for a given beta field and sliding exponent.
#
# Laura Kehrl, University of Washington, 10 April 2018

import subprocess
import time
import os
import shutil

#glacier = 'Kanger'
glacier = 'Helheim'
temperature = '-10.0'
sliding_exponent = 1
beta_date = '20120316'#'20120316'
regpar = '1e13'
nparts = 4

DIRM = os.path.join(os.getenv("MODEL_HOME"),glacier+'/3D/')

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
    beta_suffix = regpar+'_SSA_average_2000_2016_'+temperature_text
else:
    beta_suffix = regpar+'_SSA_DEM'+beta_date+'_'+temperature_text
# Get beta file name in case it needs to be copied from another directory
beta_file = DIRM+'DEM'+beta_date+'_'+temperature_text+'/inputs/beta_'+\
            slidinglaw+'_'+beta_suffix+'.dat'

# Dates for running the simulation
if glacier == 'Helheim':
    dates = ['20040803','20050829','20060825','20070912',\
           '20080814','20100624',\
           '20110319','20110615','20110828','20111116',\
           '20120316','20120624','20120908','20121205',\
           '20130209','20130508','20130804','20131031',\
           '20140127','20140509','20140731','20141016']
elif glacier == 'Kanger':
    dates = ['20010712','20030801','20050621','20050831',\
           '20060505','20060708','20060919',\
           '20070728','20090817','20100603',\
           '20110308','20110708','20110826','20111106',\
           '20120213','20120522','20121012','20121217',\
           '20130210','20130714','20131004','20131204',\
           '20140213','20150808']  

# Set up a file for saving model cost
fid_info = open(DIRM+"cost_files/cost_"+beta_suffix+"_"+slidinglaw+".dat","w")
fid_info.write('Date Cost\n')

# Loop over jobs
for date in dates:

    meshname = 'DEM'+date+'_'+temperature_text
  
    if beta_date == '':
        beta_suffix = regpar+'_SSA_DEM'+date+'_'+temperature_text
    else:
        if not(os.path.isfile(DIRM+meshname+'/inputs/beta_'+slidinglaw+'_'+beta_suffix+'.dat')):
            shutil.copy(beta_file,DIRM+meshname+'/inputs/')

    # Customize job options
    command = "python "+os.path.join(os.getenv("CODE_HOME"),"big3/modeling/simulations/simulation_ssa.py")+\
              " -glacier {0} -mesh {1} -n {2} -sif steady_ssa.sif".format(glacier,meshname,nparts)+\
              " -temperature {0} -slidinglaw {1} -beta {2}".format(temperature,slidinglaw,beta_suffix)

    print command
    os.system(command)

    fidf = open(DIRM+meshname+'/cost.dat','r')
    lines = fidf.readlines()
    p = lines[-1].split()
    fid_info.write('{0} {1}\n'.format(date,p[1]))

fid_info.close()

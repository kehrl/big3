import subprocess
import time
import os
import shutil

#glacier = 'Kanger'
glacier = 'Helheim'
temperature = 'model'
slidinglaw = 'weertman'
beta_date = ''#'20120316'
regpar = '1e13'
nparts = 4

DIRM = os.path.join(os.getenv("MODEL_HOME"),glacier+'/3D/')
if temperature == 'model':
    beta_suffix = regpar+'_SSA_DEM'+beta_date+'_modelT'
    beta_file = DIRM+'DEM'+beta_date+'_modelT/inputs/beta_'+\
            slidinglaw+'_'+beta_suffix+'.dat'
else:
    beta_suffix = regpar+'_SSA_DEM'+beta_date+'_constantT'
    beta_file = DIRM+'DEM'+beta_date+'_constantT/inputs/beta_'+\
            slidinglaw+'_'+beta_suffix+'.dat'

if glacier == 'Helheim':
    dates = ['20040803','20050829','20060825','20070912',\
           '20080814','20100624',\
           '20110319','20110615','20110828','20111116',\
           '20120316','20120624','20120908','20121205',\
           '20130209','20130508','20130804','20131031',\
           '20140127','20140509','20140731']
elif glacier == 'Kanger':
    dates = ['20010712','20030801','20050621','20050831',\
           '20060505','20060708','20060919',\
           '20070728','20090817','20100603',\
           '20110308','20110708','20110826','20111106',\
           '20120213','20120522','20121012','20121217',\
           '20130210','20130714','20131004','20131204',\
           '20140213','20150808']  

fid_info = open(DIRM+"cost_files/cost_"+beta_suffix+"_"+slidinglaw+".dat","w")
fid_info.write('Date Cost\n')

#Loop over jobs
for date in dates:

    if temperature == 'model':
        meshname = 'DEM'+date+'_modelT'
    else:
        meshname = 'DEM'+date+'_constantT'
  
    if beta_date == '':
        if temperature == 'model':
            beta_suffix = regpar+'_SSA_DEM'+date+'_modelT'
        else:
            beta_suffix = regpar+'_SSA_DEM'+date+'_constantT'
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

    if not(date == beta_date) and not(beta_date == ''):
        try:
            os.system("rm "+DIRM+meshname+'/inputs/beta_'+slidinglaw+'_'+beta_suffix+'.dat')
        except:
            pass

fid_info.close()

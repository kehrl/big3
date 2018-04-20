# Calculate an average beta for a given sliding law. The script "calculate_beta_sliding_exponent.py" may need to
# be run first to calculate the beta for each inversion for a particular sliding law exponent m.
#
# Laura Kehrl, University of Washington, 11 April 2018

import os
import numpy as np
import elmerreadlib
import datelib

glacier = 'Kanger'
regpar = '1e12'
temperature = '-10.0'
time1 = 2000.
time2 = 2016.
slidinglaw = 'linear'
SSA = True

if SSA:
    modelname = 'SSA'
else:
    modelname = 'FS'

if temperature == 'model':
    temperature_text = 'modelT'
else:
    temperature_text = 'constantT'

DIRM = os.path.join(os.getenv("MODEL_HOME"),glacier+'/3D/')

dirs = os.listdir(DIRM)
n = 0
for dir in dirs:
    if temperature_text in dir and dir.endswith('T'):
        date = dir[3:11]
        time = datelib.date_to_fracyear(int(date[0:4]),int(date[4:6]),(date[6:]))
        if (time >= time1) and (time <= time2):
            data = np.loadtxt(DIRM+dir+'/inputs/beta_'+slidinglaw+'_'+regpar+'_'+modelname+'_DEM'+\
                    date+'_'+temperature_text+'.dat',skiprows=1)
            if n == 0:
                x = data[:,0]
                y = data[:,1]
                betas = data[:,2]
            elif x[0] == data[0,0]:
                betas = np.column_stack([betas,data[:,2]])
            else:
                print "We have a problem"
            n = n+1

ave_beta = np.mean(betas,axis=1)

for dir in dirs:
    if temperature_text in dir and dir.endswith('T'):
        fid = open(DIRM+dir+'/inputs/beta_'+slidinglaw+'_'+regpar+\
                '_'+modelname+'_average_'+str(int(time1))+'_'+str(int(time2))+\
                '_'+temperature_text+'.dat','w')
        fid.write('{0}\n'.format(len(x)))
        for i in range(0,len(x)):
            fid.write('{0} {1} {2:.16f}\n'.format(x[i],y[i],ave_beta[i]))
        fid.close()

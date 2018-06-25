# Cycle through model outputs and export a beta friction coefficient for a given sliding law exponent.
#
# Laura Kehrl, 12 April 2018

import os
import numpy as np
import elmerreadlib


regpar = '1e13'
glacier = 'Kanger'
SSA = True
temperature = 'constantT'
m = 10.0

DIRM = os.path.join(os.getenv("MODEL_HOME"),glacier+'/3D/')
dirs = os.listdir(DIRM)

if SSA:
    vname = 'ssavelocity'
    modelname = 'SSA'
else:
    vname = 'velocity'
    modelname = 'FS'

for dir in dirs:
    if temperature in dir and not('Lcurve' in dir):
        
        # Get date of file
        date = dir[3:11]
        dirs2 = os.listdir(DIRM+dir+'/mesh2d/inversion_adjoint/')
        for dir2 in dirs2:
            if (regpar in dir2) and not(dir2.endswith('.pdf')):
                lambda_dir = dir2

        # Load data
        if SSA:
            bed = elmerreadlib.pvtu_file(DIRM+dir+'/mesh2d/inversion_adjoint/'+lambda_dir+\
                '/adjoint_beta_ssa0001.pvtu',['beta',vname])
        else:
            data = elmerreadlib.pvtu_file(DIRM+dir+'/mesh2d/inversion_adjoint/'+lambda_dir+\
                '/adjoint_beta0001.pvtu',['beta',vname])
            bed = elmerreadlib.values_in_layer(data,'bed')
            
        # Set up file name
        if m == 1 or m == 'linear':
            fid = open(DIRM+dir+'/inputs/beta_linear_'+regpar+'_'+modelname+'_DEM'+date+'_'+temperature+'.dat','w')
            m = 1.0
        elif m == 3 or m == 'weertman':
            fid = open(DIRM+dir+'/inputs/beta_weertman_'+regpar+'_'+modelname+'_DEM'+date+'_'+temperature+'.dat','w')
            m = 3.0
        else:
            fid = open(DIRM+dir+'/inputs/beta_m'+str(int(m))+'_'+regpar+'_'+modelname+'_DEM'+date+'_'+temperature+'.dat','w')

        fid.write('{0}\n'.format(len(bed['x'])))
        for i in range(0,len(bed)):
            fid.write('{0} {1} {2:.16f}\n'.format(bed['x'][i],bed['y'][i],\
                (bed['betasquared'][i]/(bed[vname][i]**(-float(m-1)/float(m))))))
        fid.close()
            
        del bed
        print dir

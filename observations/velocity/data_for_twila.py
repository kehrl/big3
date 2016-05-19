import os
import sys
import icefrontlib, fluxlib, datelib
import numpy as np

DIR = os.path.join(os.getenv("HOME"),"Bigtmp/Twila/")

glacier = 'Helheim'

######################
# Terminus positions #
######################

terminus_val,terminus_time,terminus_source = icefrontlib.box_method(glacier,imagesource = True,time1=2008.,time2=2016)

fid = open(DIR+"terminus_2008_2016_box_05182016.dat","w")

fid.write("#Time Position_m Satellite \n")
for i in range(0,len(terminus_time)):
  fid.write('{0:.6f} {1:.6f} {2}\n'.format(terminus_time[i],terminus_val[i],terminus_source[i])) 
fid.close()

############
# Ice flux #
############

time_Q_u, Q_u, Hbar_u, ubar_u, error_u, across_u= fluxlib.fluxgate(glacier,fluxgate_filename='fluxgate_for_twila',bedsource='cresis',dl=10.0,timing='velocity')
time_Q_zs, Q_zs, Hbar_zs, ubar_zs, error_zs, across_zs = fluxlib.fluxgate(glacier,fluxgate_filename='fluxgate_for_twila',bedsource='cresis',dl=10.0,timing='elevation')

fid = open(DIR+"iceflux_velocity_05182016.dat","w")

fid.write("#Time IceFlux_m3_yr-1 Across_m3 \n")

ind = np.where(~(np.isnan(Q_u)))[0]
for i in ind:
  fid.write('{0:.6f} {1:.2f} {2:.2f}\n'.format(time_Q_u[i],Q_u[i],across_u[i])) 
fid.close()

fid = open(DIR+"iceflux_elevation_05182016.dat","w")

fid.write("#Time IceFlux_m3_yr-1 Across_m3\n")

ind = np.where(~(np.isnan(Q_zs)))[0]
for i in ind:
  fid.write('{0:.6f} {1:.2f} {2:.2f}\n'.format(time_Q_zs[i],Q_zs[i],across_zs[i])) 
fid.close()
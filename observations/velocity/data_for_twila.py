import os
import sys
import icefrontlib, fluxlib, datelib
import numpy as np

DIR = os.path.join(os.getenv("HOME"),"Bigtmp/Twila/")

######################
# Terminus positions #
######################

#terminus_val,terminus_time,terminus_source = icefrontlib.box_method(glacier,imagesource = True,time1=2008.,time2=2016)
#
#fid = open(DIR+"terminus_2008_2016_box_05182016.dat","w")
#
#fid.write("#Time Position_m Satellite \n")
#for i in range(0,len(terminus_time)):
#  fid.write('{0:.6f} {1:.6f} {2}\n'.format(terminus_time[i],terminus_val[i],terminus_source[i])) 
#fid.close()

########################
# Ice flux for Helheim #
########################

time_Q_u, Q_u, Hbar_u, ubar_u, error_u, across_u= fluxlib.fluxgate('Helheim',fluxgate_filename='fluxgate_for_twila',bedsource='cresis',dl=10.0,timing='velocity')
#time_Q_zs, Q_zs, Hbar_zs, ubar_zs, error_zs, across_zs = fluxlib.fluxgate(glacier,fluxgate_filename='fluxgate_for_twila',bedsource='cresis',dl=10.0,timing='elevation')

fid = open(DIR+"iceflux_helheim_11062016.dat","w")

fid.write("#Time IceFlux_m3_yr-1 Across_m3 \n")

ind = np.where(~(np.isnan(Q_u)))[0]
for i in ind:
  fid.write('{0:.6f} {1:.2f} {2:.2f}\n'.format(time_Q_u[i],Q_u[i],across_u[i])) 
fid.close()


#####################
# Flux for Midgaard #
#####################

# Note that I originally calculated ice fluxes for both the west and east branch, but 
# average ice thickness is ~30 m near the east branch terminus, so it's not really contributing
# anything to total flux. I found that the errors for the east branch were roughly the same
# magnitude as the estimated flux and ~five orders of magnitude smaller than the flux from the
# west branch.

time_mid,Q_mid,Hbar_mid,ubar_mid,error_mid,across_mid = fluxlib.fluxgate('Midgaard',fluxgate_filename='fluxgate_midgaard_west.shp',bedsource='morlighem',dl=10.0,timing='velocity')

fid = open(DIR+"iceflux_midgaard_west_11062016.dat","w")

fid.write("#Time IceFlux_m3_yr-1 Across_m3 \n")

ind = np.where(~(np.isnan(Q_mid)))[0]
for i in ind:
  fid.write('{0:.6f} {1:.2f} {2:.2f}\n'.format(time_mid[i],Q_mid[i],across_mid[i])) 
fid.close()

###################
# Flux for Fenris #
###################

# OK. So we don't have any TSX velocities for Fenris (at least I think). So let's just
# play around with Howat's velocities and see what happens.

time_fen,Q_fen,Hbar_fen,ubar_fen,error_fen,across_fen = fluxlib.fluxgate('Fenris','fluxgate_fenris.shp',bedsource='morlighem',dl=10.0,timing='velocity')

fid = open(DIR+"iceflux_fenris_11062016.dat","w")

fid.write("#Time IceFlux_m3_yr-1 Across_m3\n")

ind = np.where(~(np.isnan(Q_fen)))[0]
for i in ind:
  fid.write('{0:.6f} {1:.2f} {2:.2f}\n'.format(time_fen[i],Q_fen[i],across_fen[i])) 
fid.close()


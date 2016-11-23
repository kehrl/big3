import numpy as np
import matplotlib.pyplot as plt
import vellib, bedlib, zslib, glaclib, icefrontlib, floatlib, datelib

##########
# Inputs #
##########

# Get Glacier name
glacier = 'Helheim'

# Start and end time
time1 = 2011.0 #start time for plot
time2 = 2016.0 # end time for plot

# Sample points 
dists_eul = -1.0*np.array([2,5.0,10.0,15.0,20.0]) # kilometers

#############
# Constants #
#############

g = 9.81       # gravity
rho_i = 917.   # ice density
rho_sw = 1020. # seawater density

#############
# Load data #
#############

# Load flowline
x,y,zb,dists = glaclib.load_flowline(glacier,shapefilename='flowline_flightline',filt_len=2.0e3,bedsource='cresis',verticaldatum='geoid')

# Get velocities at points
#ind_eul=[]
#for i in range(0,len(dists_eul)):
#  ind_eul.append( (abs(dists - dists_eul[i]*1e3)).argmin() )
#vel_val,vel_time,vel_error = vellib.velocity_at_eulpoints(x[ind_eul],y[ind_eul],glacier)

# Get elevations at points
# zpt_dem,zpterror_dem,zpt_time_dem = zslib.dem_at_pts(x[ind_eul],y[ind_eul],glacier,years='all',verticaldatum='geoid',cutoff='terminus',method='average',radius=200.)

# Get velocities along flowline
vel_dem,vel_time,vel_term = vellib.velocity_along_flowline(x,y,dists,glacier,cutoff='terminus',data='TSX')
ind = np.where((vel_time > time1) & (vel_time < time2))[0]
vel_dem = vel_dem[:,ind]
vel_time = vel_time[ind]

# Terminus position
#terminus_val, terminus_time = icefrontlib.distance_along_flowline(x,y,dists,glacier,type='icefront',time1=time1,time2=time2)
#time = terminus_time
terminus_time = vel_time
terminus_val = vel_term

# Get elevations along flowline
zs_dem,time_dem = zslib.dem_along_flowline(x,y,glacier,years='all',cutoff='terminus',verticaldatum='geoid',filt_len=2.0e3,data='TDM')
ind = np.where((time_dem > time1) & (time_dem < time2))[0]
zs_dem = zs_dem[ind,:]
time_dem = time_dem[ind]

time = time_dem[:,0]
#time = vel_time

#########################################
# Calculate terminus force for flowline #
#########################################

# This is the force at any given point along the flowline, if the glacier were to retreat to
# that location. Note that the assumed surface height will greatly affect the result.

# Assuming the ice is at flotation
h_float = floatlib.height(zb)
H_float = h_float-zb
Ft_float = (g*rho_i/2)*( (1-rho_sw/rho_i)*H_float**2 + (rho_sw/rho_i)*h_float*(2*H_float-h_float) )

# Using average surface heights
h_ave = np.nanmean(zs_dem,axis=0)
H_ave = h_ave-zb
Ft_ave = (g*rho_i/2)*( (1-rho_sw/rho_i)*H_ave**2 + (rho_sw/rho_i)*h_ave*(2*H_ave-h_ave) )

####################################################################################
# Calculate terminus force & its contributing to the "driving stress" through time #
####################################################################################

#h = np.zeros(len(time_dem))
#for i in range(0,len(time_dem)):
#  term = np.interp(time_dem[i,0],terminus_time,terminus_val)
#  ind = np.argmin(abs(term-dists))
#  h[i] = np.nanmean(zs_dem[i,ind-20:ind+1])

# For the moment, we will assume the surface height is always 60 m.
h = 70

Ftot_time = np.zeros(len(time))
for i in range(0,len(time)):
  ind = np.argmin(abs(dists-terminus_val[i]))
  H = h-zb[ind]
  Ftot_time[i] = (g*rho_i/2)*( (1-rho_sw/rho_i)*H**2 + (rho_sw/rho_i)*h*(2*H-h) )

# Length over which the terminus force acts
L_o = 20.e3

# Contributing of terminus force to the driving stress
tau_f = np.zeros([len(dists),len(time)])
for i in range(0,len(time)):
  ind1 = np.argmin(abs(dists-terminus_val[i]))
  ind2 = np.argmin(abs(dists--L_o))
  L = np.linspace(0,1.,ind1-ind2)
  tau_f[ind2:ind1,i] = L*Ftot_time[i]/np.sum(L)
  
############################################################  
# Calculate driving stress along the flowline through time #  
############################################################

tau_d = np.zeros([len(dists),len(time)])
for i in range(0,len(time)):
  ind = np.argmin(abs(time_dem[:,0]-time[i]))
  zs = zs_dem[ind,:]
  dhdx = np.zeros(len(dists))
  dhdx[0:-1] = abs(np.diff(zs)/np.diff(dists))
  for j in range(50,len(dists)-50):
    tau_d[j,i] = rho_i*g*(zs[j]-zb[j])*np.nanmean(dhdx[j-50:j+50])

###################################
# Reference velocity and stresses #
###################################

time_ref = datelib.date_to_fracyear(2012,03,16)

ind = np.argmin(abs(time-time_ref))
tau_f_ref = tau_f[:,ind]
tau_d_ref = tau_d[:,ind]

ind = np.argmin(abs(time_ref-vel_time))
vel_ref = vel_dem[:,ind]

################################
# Calculate modeled velocities #
################################

fudge=0
ratios = np.zeros([len(dists),len(time)])
velocities = np.zeros([len(dists),len(time)])
velocities_measured = np.zeros([len(dists),len(time)])
for i in range(0,len(time)):
  for j in range(0,len(dists)):
    ratios[j,i] = ((np.nanmean(tau_d[j,:])+tau_f[j,i]+fudge)/(np.nanmean(tau_d[j,:])+fudge+tau_f_ref[j]))**3
  ind = np.where((ratios[:,i] > 5) | (ratios[:,i] < 0.5))
  ratios[ind,i] = float('NaN')
  velocities[:,i] = ratios[:,i]*vel_ref
  ind = np.argmin(abs(vel_time-time[i]))
  velocities_measured[:,i] = vel_dem[:,ind].T
  
data_ratios = np.zeros([len(dists),len(time)])
for i in range(0,len(time)):
  data_ratios[:,i] = velocities_measured[:,i]/vel_ref
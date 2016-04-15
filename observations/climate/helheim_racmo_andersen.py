# This plot looks at SMB for Helheim from RACMO and Andersen observations. We need to know how
# much RACMO underestimates ablation rates along the trunk of Helheim. What should
# ablation rates be there?
#
# LMK, UW, 9/16/2015

import os
import numpy as np
import climlib, glaclib,icefrontlib, vellib, datelib, coordlib
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import AutoMinorLocator

#################################################################################
# Get coordinates for distance D along the flowline for both Kanger and Helheim #
#################################################################################
D = [10.0e3] # m inland from average terminus position

# Helheim:
xH,yH,zbH,distsH = glaclib.load_flowline('Helheim')
ind = []
for i in range(0,len(D)):
  ind.append(np.argmin(abs(distsH+D[i])))
xptH = xH[ind]
yptH = yH[ind]
terminus_valH,terminus_timeH = icefrontlib.distance_along_flowline(xH,yH,distsH,'Helheim',type='icefront')
vel_valH,vel_timeH,vel_errorH = vellib.velocity_at_eulpoints(xptH,yptH,'Helheim',data='TSX')

del ind

###################
# Load RACMO data #
###################

filt_len = 'none' # filter length in days for timeseries

# Helheim
xracH,yracH,runoffH,timeH = climlib.racmo_at_pts(xptH,yptH,'runoff',filt_len=filt_len)
xracH,yracH,zsH,timeHz = climlib.racmo_at_pts(xptH,yptH,'zs',filt_len='none')
xracH,yracH,t2mH,timeH = climlib.racmo_at_pts(xptH,yptH,'t2m',filt_len=filt_len)
xracH,yracH,smbH,timeH = climlib.racmo_at_pts(xptH,yptH,'smb',filt_len=filt_len)
xracH,yracH,precipH,timeH = climlib.racmo_at_pts(xptH,yptH,'precip',filt_len=filt_len)

##############################
# Numbers from Andersen 2010 #
##############################

# Their measurements were from an AWS station installed at xA,yA

rho_i=900.

xA,yA = coordlib.convert(-38.44,66.42,4326,3413)

xrac_A,yrac_A,zs_A,time_z_A = climlib.racmo_at_pts(xA,yA,'zs',filt_len=filt_len)
xrac_A,yrac_A,smb_A,time_A = climlib.racmo_at_pts(xA,yA,'smb',filt_len=filt_len)
xrac_A,yrac_A,runoff_A,time_A = climlib.racmo_at_pts(xA,yA,'runoff',filt_len=filt_len)

# 2007
day1_2007 = datelib.doy_to_fracyear(2007,208)
day2_2007 = datelib.doy_to_fracyear(2007,235)
meltrate_obs_2007 = -0.89/(day2_2007-day1_2007)
meltrate_obs_2007_downstream = -0.0444*365.25

ind1 = np.argmin(abs(time_z_A-day1_2007))
ind2 = np.argmin(abs(time_z_A-day2_2007))
meltrate_rac_zs_2007 = (zs_A[ind2]-zs_A[ind1])/(time_z_A[ind2]-time_z_A[ind1])

ind1 = np.argmin(abs(time_A-day1_2007))
ind2 = np.argmin(abs(time_A-day2_2007))
meltrate_rac_smb_2007 = np.mean((smb_A[ind1:ind2+1]))*365.25/rho_i
meltrate_rac_runoff_2007 = -1*np.mean((runoff_A[ind1:ind2+1]))*365.25/rho_i

# 2008
day1_2008 = datelib.doy_to_fracyear(2008,182)
day2_2008 = datelib.doy_to_fracyear(2008,232)
meltrate_obs_2008 = -1.59/(day2_2008-day1_2008)
meltrate_obs_2008_downstream = -0.0436*365.25

ind1 = np.argmin(abs(time_z_A-day1_2008))
ind2 = np.argmin(abs(time_z_A-day2_2008))
meltrate_rac_zs_2008 = (zs_A[ind2]-zs_A[ind1])/(time_z_A[ind2]-time_z_A[ind1])

ind1 = np.argmin(abs(time_A-day1_2008))
ind2 = np.argmin(abs(time_A-day2_2008))
meltrate_rac_smb_2008 = np.mean((smb_A[ind1:ind2+1]))*365.25/rho_i
meltrate_rac_runoff_2008 = -1*np.mean((runoff_A[ind1:ind2+1]))*365.25/rho_i


###########################################################
# Quantify beginning and end of melt season for each year #
###########################################################

yH_t2m,day1H_t2m,day2H_t2m,meltlengthH_t2m,totalH_t2m = climlib.seasonlength(timeH,t2mH,'t2m')
yH_runoff,day1H_runoff,day2H_runoff,meltlengthH_runoff,totalH_runoff = climlib.seasonlength(timeH,runoffH,'runoff')
yH_smb,day1H_smb,day2_smb,meltlengthH_smb,totalH_smb = climlib.seasonlength(timeH,smbH,'smb')

###########################################
# Terminus position, velocity, and runoff #
###########################################

# Helheim
plt.figure(figsize=(6.5,4))
gs = matplotlib.gridspec.GridSpec(3,1)

plt.subplot(gs[0,:])
ax = plt.gca()
plt.plot(vel_timeH,vel_valH/1.0e3,'go',markersize=3)
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1H_runoff[i],5.5],[day1H_runoff[i],7.5],[day2H_runoff[i],7.5],[day2H_runoff[i],5.5]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.75',edgecolor='none',lw=0)
  ax.add_patch(patch)
plt.xlim([time1,time2])
plt.ylim([5.5,7.5])
plt.yticks([5.5,6.0,6.5,7.0,7.5],fontsize=10)
plt.ylabel('Velocity \n (km/yr)',fontsize=10)
ax.set_xticklabels([])

plt.subplot(gs[1,:])
ax = plt.gca()
plt.plot(terminus_timeH,terminus_valH/1.0e3,'k.')
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1H_runoff[i],-2.5],[day1H_runoff[i],2.5],[day2H_runoff[i],2.5],[day2H_runoff[i],-2.5]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.75',edgecolor='none',lw=0)
  ax.add_patch(patch)
plt.xlim([time1,time2])
plt.ylim([-2,2.5])
plt.yticks([-2.0,-1.0,0,1.0,2.0],fontsize=10)
plt.ylabel('Terminus \n (km)',fontsize=10)
ax.set_xticklabels([])

plt.subplot(gs[2,:])
ax = plt.gca()
ax2 = ax.twinx()
ax.set_zorder(ax2.get_zorder()+1)
ax.set_frame_on(False)
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1H_runoff[i],0],[day1H_runoff[i],5000],[day2H_runoff[i],5000],[day2H_runoff[i],0]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.75',edgecolor='none',lw=0)
  ax2.add_patch(patch)
  ax2.plot([day1H_runoff[i],day1H_runoff[i],day2H_runoff[i],day2H_runoff[i]],[0,totalH_runoff[i],totalH_runoff[i],0],'r',lw=1.5)
  if yH_runoff[i] >= time1 and yH_runoff[i] < time2:
    ax.text(day1H_runoff[i]+0.03,58,str(int(meltlengthH_runoff[i]))+' d',fontsize=10)
plt.xlim([time1,time2])
ax.plot(timeH,runoffH,'k',linewidth=1.5)
x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
ax.xaxis.set_major_formatter(x_formatter)
ax.set_ylabel('Runoff \n (kg m$^{-2}$ d$^{-1}$)',fontsize=10)
ax2.set_ylabel('Total Runoff \n (kg m$^{-2}$ d$^{-1}$)',fontsize=10,color='r')
ax.set_yticks([0,20,40,60])
ax.set_ylim([0,70])
ax2.set_ylim([0,5000])
ax.tick_params(axis='both',labelsize=10)
ax2.tick_params(axis='both',labelsize=10)
ax2.set_yticks(np.arange(0,5000,1000))
for tl in ax2.get_yticklabels():
    tl.set_color('r')

plt.subplots_adjust(hspace=0.05,wspace=0) 

plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Helheim_racmo.pdf"),FORMAT='PDF',dpi=800)
plt.close()

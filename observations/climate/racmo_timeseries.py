# This file makes a plot of the RACMO timeseries for both Helheim and Kanger from 2008-2015.
#
# LMK, UW, 9/16/2015

import os
import numpy as np
import climlib, glaclib,icefrontlib, vellib, datelib
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import AutoMinorLocator

#################################################################################
# Get coordinates for distance D along the flowline for both Kanger and Helheim #
#################################################################################
D = 10.0e3 # m inland from average terminus position

# Helheim:
xH,yH,zbH,distsH = glaclib.load_flowline('Helheim')
ind = np.argmin(abs(distsH+D))
xptH = xH[ind]
yptH = yH[ind]
terminus_valH,terminus_timeH = icefrontlib.distance_along_flowline(xH,yH,distsH,'Helheim',type='icefront')
vel_valH,vel_timeH,vel_errorH = vellib.velocity_at_eulpoints(xptH,yptH,'Helheim',data='TSX')

del xH,yH,zbH,distsH,ind

# Kanger:
xK,yK,zbK,distsK = glaclib.load_flowline('Kanger')
ind = np.argmin(abs(distsK+D))
xptK = xK[ind]
yptK = yK[ind]
terminus_valK,terminus_timeK = icefrontlib.distance_along_flowline(xK,yK,distsK,'Kanger',type='icefront')
vel_valK,vel_timeK,vel_errorK = vellib.velocity_at_eulpoints(xptK,yptK,'Kanger',data='TSX')

del xK,yK,zbK,distsK,ind

###################
# Load RACMO data #
###################

filt_len = 31.0 # filter length in days for timeseries

# Kanger
xracK,yracK,runoffK,timeK = climlib.racmo_at_pts(xptK,yptK,'runoff',filt_len=filt_len)
xracK,yracK,zsK,timeKz = climlib.racmo_at_pts(xptK,yptK,'zs',filt_len='none')
xracK,yracK,t2mK,timeK = climlib.racmo_at_pts(xptK,yptK,'t2m',filt_len=filt_len)
xracK,yracK,smbK,timeK = climlib.racmo_at_pts(xptK,yptK,'smb',filt_len=filt_len)
xracK,yracK,precipK,timeK = climlib.racmo_at_pts(xptK,yptK,'precip',filt_len=filt_len)

# Helheim
xracH,yracH,runoffH,timeH = climlib.racmo_at_pts(xptH,yptH,'runoff',filt_len=filt_len)
xracH,yracH,zsH,timeHz = climlib.racmo_at_pts(xptH,yptH,'zs',filt_len='none')
xracH,yracH,t2mH,timeH = climlib.racmo_at_pts(xptH,yptH,'t2m',filt_len=filt_len)
xracH,yracH,smbH,timeH = climlib.racmo_at_pts(xptH,yptH,'smb',filt_len=filt_len)
xracH,yracH,precipH,timeH = climlib.racmo_at_pts(xptH,yptH,'precip',filt_len=filt_len)

#############
# Make plot #
#############

time1=2008.0
time2=2015.0
plt.figure(figsize=(4,4))
gs = matplotlib.gridspec.GridSpec(3,1)

plt.subplot(gs[0,:])
ax = plt.gca()
plt.plot([time1,time2],[0,0],'k')
plt.plot(timeH,t2mH,'r',linewidth=1.5,label='Helheim')
plt.plot(timeK,t2mK,'b',linewidth=1.5,label='Kanger')
plt.ylabel('Temperature ($^o$C)')
plt.xlim([time1,time2])
ax.set_xticklabels([])

plt.subplot(gs[1,:])
ax = plt.gca()
plt.plot([time1,time2],[0,0],'k')
plt.plot(timeH,smbH*365.25/917.0,'r',linewidth=1.5,label='Helheim')
plt.plot(timeK,smbK*365.25/917.0,'b',linewidth=1.5,label='Kanger')
plt.legend(loc=2,fontsize=10,numpoints=1,handlelength=0.4,labelspacing=0.05,ncol=2,columnspacing=0.7,handletextpad=0.2)
plt.ylabel('SMB (m a$^{-1}$)')
plt.xlim([time1,time2])
ax.set_xticklabels([])

plt.subplot(gs[2,:])
ax = plt.gca()
plt.plot(timeH,runoffH*365.25/917.0,'r',linewidth=1.5)
plt.plot(timeK,runoffK*365.25/917.0,'b',linewidth=1.5)
plt.ylabel('Runoff (m a$^{-1}$)')
plt.xlim([time1,time2])
x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
ax.xaxis.set_major_formatter(x_formatter)

plt.subplots_adjust(hspace=0.05,wspace=0)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/racmo_timeseries.pdf"),FORMAT='PDF',dpi=800)
plt.close()

###########################################################
# Quantify beginning and end of melt season for each year #
###########################################################

yH_t2m,day1H_t2m,day2H_t2m,meltlengthH_t2m,totalH_t2m = climlib.seasonlength(timeH,t2mH,'t2m')
yH_runoff,day1H_runoff,day2H_runoff,meltlengthH_runoff,totalH_runoff = climlib.seasonlength(timeH,runoffH,'runoff')
yH_smb,day1H_smb,day2_smb,meltlengthH_smb,totalH_smb = climlib.seasonlength(timeH,smbH,'smb')

yK_t2m,day1K_t2m,day2K_t2m,meltlengthK_t2m,totalK_t2m = climlib.seasonlength(timeK,t2mK,'t2m')
yK_runoff,day1K_runoff,day2K_runoff,meltlengthK_runoff,totalK_runoff = climlib.seasonlength(timeK,runoffK,'runoff')
yK_smb,day1K_smb,day2_smb,meltlengthK_smb,totalK_smb = climlib.seasonlength(timeK,smbH,'smb')

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

# Kanger
plt.figure(figsize=(6.5,4))
gs = matplotlib.gridspec.GridSpec(3,1)

plt.subplot(gs[0,:])
ax = plt.gca()
plt.plot(vel_timeK,vel_valK/1.0e3,'go',markersize=3)
for i in range(0,len(yK_runoff)):
  path = matplotlib.path.Path([[day1K_runoff[i],5.5],[day1K_runoff[i],7.5],[day2K_runoff[i],7.5],[day2K_runoff[i],5.5]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.75',edgecolor='none',lw=0)
  ax.add_patch(patch)
plt.xlim([time1,time2])
plt.ylim([5.5,7.5])
plt.yticks([5.5,6.0,6.5,7.0,7.5],fontsize=10)
plt.ylabel('Velocity \n (km/yr)',fontsize=10)
ax.set_xticklabels([])

plt.subplot(gs[1,:])
ax = plt.gca()
plt.plot(terminus_timeK,terminus_valK/1.0e3,'k.')
for i in range(0,len(yK_runoff)):
  path = matplotlib.path.Path([[day1K_runoff[i],-3.5],[day1K_runoff[i],3.5],[day2K_runoff[i],3.5],[day2K_runoff[i],-3.5]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.75',edgecolor='none',lw=0)
  ax.add_patch(patch)
plt.xlim([time1,time2])
plt.ylim([-3.5,3])
plt.yticks([-3.0,-2.0,-1.0,0,1.0,2.0,3.0],fontsize=10)
plt.ylabel('Terminus \n (km)',fontsize=10)
ax.set_xticklabels([])

plt.subplot(gs[2,:])
ax = plt.gca()
ax2 = ax.twinx()
ax.set_zorder(ax2.get_zorder()+1)
ax.set_frame_on(False)
for i in range(0,len(yK_runoff)):
  path = matplotlib.path.Path([[day1K_runoff[i],0],[day1K_runoff[i],5000],[day2K_runoff[i],5000],[day2K_runoff[i],0]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.75',edgecolor='none',lw=0)
  ax2.add_patch(patch)
  ax2.plot([day1K_runoff[i],day1K_runoff[i],day2K_runoff[i],day2K_runoff[i]],[0,totalK_runoff[i],totalK_runoff[i],0],'r',lw=1.5)
  if yK_runoff[i] >= time1 and yK_runoff[i] < time2:
    ax.text(day1K_runoff[i]+0.03,55,str(int(meltlengthK_runoff[i]))+' d',fontsize=10)
plt.xlim([time1,time2])
ax.plot(timeK,runoffK,'k',linewidth=1.5)
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

plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Kanger_racmo.pdf"),FORMAT='PDF',dpi=800)
plt.close()
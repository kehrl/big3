# This script compares elevation change from the DEMs to anticipated thinning rates from 
# the fluxgate method.

# LMK, UW, 8/31/2015

import os
import sys
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"Util/Modules/"))
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"BigThreeGlaciers/Tools/"))
import glacier_flowline, elevation, icefronts, fluxgate, climate
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib
import geotiff
 

###########
# WV DEMs #
###########

# Load DEMs	
xwv_H,ywv_H,zwv_H,timewv_H,errorwv_H = elevation.dem_grid('Helheim',285000.0,320000.0,-2588000.0,-2566000.0,years='all',verticaldatum='ellipsoid',method='nearest',return_error=True)
xwv_K,ywv_K,zwv_K,timewv_K,errorwv_K = elevation.dem_grid('Kanger',449800.0,503000.0,-2302000.0,-2266000.0,years='all',verticaldatum='ellipsoid',method='nearest',return_error=True)

###########################
# Load terminus positions #
###########################

x_H,y_H,zb_H,dists_H = glacier_flowline.load('Helheim',shapefilename='flowline_flightline')
x_K,y_K,zb_K,dists_K = glacier_flowline.load('Kanger')

terminus_val_H, terminus_time_H = icefronts.distance_along_flowline(x_H,y_H,dists_H,'Helheim',type='icefront',time1=2008.,time2=2016.)
terminus_val_K, terminus_time_K = icefronts.distance_along_flowline(x_K,y_K,dists_K,'Kanger',type='icefront',time1=2008.,time2=2016.)

##################
# Load fluxgates #
##################

xgate_H,ygate_H = fluxgate.fluxbox_geometry('Helheim',"fluxgate3")
xgate_K,ygate_K = fluxgate.fluxbox_geometry('Kanger',"fluxgate3")

##################
# SMB from RACMO #
##################

xrac_H,yrac_H,smbrac_H,timerac_H = climate.racmo_at_pts(np.mean(xgate_H),np.mean(ygate_H),'smb',filt_len='none')
xrac_K,yrac_K,smbrac_K,timerac_K = climate.racmo_at_pts(np.mean(xgate_K),np.mean(ygate_K),'smb',filt_len='none')


###############################################
# dH/dt from the fluxgate, WV, and comparison #
###############################################

# Helheim
flux_time_H,flux_dH_H = fluxgate.fluxgate_thinning('Helheim',"fluxgate3",bedsource='smith')
wv_time_H,wv_dH_H = fluxgate.dem_thinning('Helheim',xwv_H,ywv_H,zwv_H,timewv_H,errorwv_H,"fluxgate3")
dH_time_H,dH_flux_H,dH_dem_H,dH_smb_H= fluxgate.compare_thinning_rates(wv_time_H,wv_dH_H,flux_time_H,flux_dH_H,timerac_H,smbrac_H,rho_i=900.0)


# Kanger
flux_time_K,flux_dH_K = fluxgate.fluxgate_thinning('Kanger',"fluxgate3",bedsource='cresis')
wv_time_K,wv_dH_K = fluxgate.dem_thinning('Kanger',xwv_K,ywv_K,zwv_K,timewv_K,errorwv_K,"fluxgate3")
dH_time_K,dH_flux_K,dH_dem_K,dH_smb_K = fluxgate.compare_thinning_rates(wv_time_K,wv_dH_K,flux_time_K,flux_dH_K,timerac_K,smbrac_K,rho_i=900.0)

plt.figure(figsize=(3,3))
matplotlib.rc('font',family='Arial')
plt.plot([-100,100],[-100,100],c='k',lw=1)
plt.errorbar(dH_dem_K[:,0],dH_flux_K[:,0],yerr=dH_flux_K[:,1],xerr=dH_dem_K[:,1],fmt='o',color='0.7',markersize=3,zorder=3,label='Kanger')
plt.errorbar(dH_dem_H[:,0],dH_flux_H[:,0],yerr=dH_flux_H[:,1],xerr=dH_dem_H[:,1],fmt='ko',markersize=3,zorder=3,label='Helheim')
plt.xticks(np.arange(-100,120,20),fontsize=8)
plt.yticks(np.arange(-100,120,20),fontsize=8)
plt.ylim([-100,80])
plt.xlim([-100,80])
plt.ylabel('Flux thinning rate (m/yr)',fontsize=8)
plt.xlabel('DEM thinning rate (m/yr)',fontsize=8)
plt.tight_layout()
plt.legend(loc=4,numpoints=1,fontsize=8)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+"Thinning_comparison.pdf"),FORMAT='PDF')
plt.close()

plt.figure(figsize=(3,3))
matplotlib.rc('font',family='Arial')
plt.plot([0,0],[-150,150],'k')
plt.plot([-4,4],[0,0],'k')
interped = np.interp(flux_time_H,terminus_time_H,terminus_val_H)
plt.errorbar(interped/1e3,flux_dH_H[:,0],fmt='o',yerr=flux_dH_H[:,1],markersize=3,c='k',label='Flux')
interped = np.interp(wv_time_H[:,0],terminus_time_H,terminus_val_H)
plt.errorbar(interped/1e3,wv_dH_H[:,0],fmt='o',yerr=wv_dH_H[:,1],markersize=3,c='r',label='DEM')
plt.xlabel('Terminus position (km)',fontsize=8)
plt.ylabel('Thinning rate (m/yr)',fontsize=8)
plt.yticks(np.arange(-150,200,75),fontsize=8)
plt.xticks(np.arange(-4,5,1),fontsize=8)
plt.xlim([-2,2])
plt.ylim([-150,150])
plt.legend(loc=2,fontsize=8,numpoints=1)
plt.tight_layout()
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Helheim_terminus_thinning.pdf"),FORMAT='PDF')
plt.close()

plt.figure(figsize=(3,3))
matplotlib.rc('font',family='Arial')
plt.plot([0,0],[-150,150],'k')
plt.plot([-4,4],[0,0],'k')
interped = np.interp(flux_time_K,terminus_time_K,terminus_val_K)
plt.errorbar(interped/1e3,flux_dH_K[:,0],fmt='o',yerr=flux_dH_K[:,1],markersize=3,c='k',label='Flux')
interped = np.interp(wv_time_K[:,0],terminus_time_K,terminus_val_K)
plt.errorbar(interped/1e3,wv_dH_K[:,0],fmt='o',yerr=wv_dH_K[:,1],markersize=3,c='r',label='DEM')
plt.xlabel('Terminus position (km)',fontsize=8)
plt.ylabel('Thinning rate (m/yr)',fontsize=8)
plt.yticks(np.arange(-150,200,75),fontsize=8)
plt.xticks(np.arange(-4,5,1),fontsize=8)
plt.xlim([-2,2])
plt.ylim([-150,150])
plt.legend(loc=2,fontsize=8,numpoints=1)
plt.tight_layout()
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Kanger_terminus_thinning.pdf"),FORMAT='PDF')
plt.close()


####################################################################
# OK. Let's try to really compare DEM and inferred thinning rates. #
####################################################################

time1 = 2008.5
time2 = 2015.5
plt.figure(figsize=(6.5,6.5))

ax = plt.gca()
nonnan = np.where(~(np.isnan(flux_dH_gate3[:,0])))[0]
plt.errorbar(wv_time_gate3[:,0],wv_dH_gate3[:,0],xerr=wv_time_gate3[:,1],yerr=wv_dH_gate3[:,1],fmt='o',c='r',markersize=3)
plt.errorbar(flux_time_gate3[nonnan],flux_dH_gate3[nonnan,0],yerr=flux_dH_gate3[nonnan,1],fmt='o-',color='k',markersize=3)
#plt.plot(flux_time_gate3_cresis,flux_dH_gate3_cresis[:,0],'g+',markersize=3)
x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
ax.xaxis.set_major_formatter(x_formatter)

nonnan = np.where(~(np.isnan(flux_dH_gate3[:,0])))[0]
plt.errorbar(wv_time_gate1[:,0],wv_dH_gate1[:,0],xerr=wv_time_gate1[:,1],yerr=wv_dH_gate1[:,1],fmt='o',c='r',markersize=3)
plt.errorbar(flux_time_gate1[nonnan],flux_dH_gate1[nonnan,0],flux_dH_gate1[nonnan,1],fmt='o',color='m',markersize=3)
plt.xlim([time1,time2])

plt.figure()
plt.plot([-150,80],[-150,80],'k--',lw=1.5)
plt.errorbar(dH_dem_gate1[:,0],dH_flux_gate1[:,0],xerr=dH_dem_gate1[:,1],yerr=dH_flux_gate1[:,1],fmt='go')
plt.errorbar(dH_dem_gate3[:,0]*10,dH_flux_gate3[:,0],xerr=dH_dem_gate3[:,1],yerr=dH_flux_gate3[:,1],fmt='ro')
plt.gca().set_aspect('equal', adjustable='box')
plt.xlim([-120,80])
plt.ylim([-120,80])
plt.xlabel('DEM thinning rate (m/yr)',fontsize=10)
plt.ylabel('Velocity thinning rate (m/yr)',fontsize=10)
plt.xticks(np.arange(-100,120,40),fontsize=10)
plt.yticks(np.arange(-100,120,40),fontsize=10)
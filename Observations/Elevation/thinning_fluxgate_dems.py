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

plt.figure(figsize=(2.7,2.7))
matplotlib.rc('font',family='Arial')
plt.plot([-100,100],[-100,100],c='k',lw=1)
plt.errorbar(dH_dem_K[:,0],dH_flux_K[:,0],yerr=dH_flux_K[:,1],xerr=dH_dem_K[:,1],fmt='o',color='0.7',markersize=3,zorder=3,label='Kanger')
plt.errorbar(dH_dem_H[:,0],dH_flux_H[:,0],yerr=dH_flux_H[:,1],xerr=dH_dem_H[:,1],fmt='ko',markersize=3,zorder=3,label='Helheim')
plt.xticks(np.arange(-120,120,40),fontsize=10)
plt.yticks(np.arange(-120,120,40),fontsize=10)
plt.ylim([-100,80])
plt.xlim([-100,80])
plt.ylabel('Flux dH/dt (m yr$^{-1}$)',fontsize=11)
plt.xlabel('DEM dH/dt (m yr$^{-1}$)',fontsize=11)
plt.tight_layout()
plt.legend(loc=4,borderpad=0.3,fontsize=10,numpoints=1,handlelength=0.7,labelspacing=0.05,columnspacing=0.7,handletextpad=0.5)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+"Thinning_comparison.pdf"),FORMAT='PDF',dpi=600)
plt.close()

plt.figure(figsize=(3,3))
matplotlib.rc('font',family='Arial')
plt.plot([0,0],[-150,150],'k')
plt.plot([-4,4],[0,0],'k')
plt.errorbar(-200,-200,fmt='ko',markersize=3,label='Flux')
plt.errorbar(-200,-200,fmt='k^',markersize=3,label='DEM')
years = range(2008,2016)
colors=['b','skyblue','g','limegreen','gold','orange','r','maroon']
n=0
for year in years:
  ind = np.where(np.floor(flux_time_H)==year)[0]
  interped = np.interp(flux_time_H[ind],terminus_time_H,terminus_val_H)
  plt.errorbar(interped/1e3,flux_dH_H[ind,0],fmt='o',yerr=flux_dH_H[ind,1],markersize=3,c=colors[n])
  ind = np.where(np.floor(wv_time_H)==year)[0]
  for i in ind:
    time = np.arange(wv_time_H[i,0]-wv_time_H[i,1],wv_time_H[i,0]+wv_time_H[i,1],1/365.25)
    interped = np.interp(time,terminus_time_H,terminus_val_H)
    plt.errorbar(np.mean(interped)/1e3,wv_dH_H[i,0],fmt='^',yerr=wv_dH_H[i,1],markersize=3,c=colors[n])
  plt.plot(-200,-200,'s',markersize=3,color=colors[n],label=year)
  n=n+1
n=0
for year in years:
  ind = np.where(np.floor(flux_time_H)==year)[0]
  interped = np.interp(flux_time_H[ind],terminus_time_H,terminus_val_H)
  plt.plot(interped/1e3,flux_dH_H[ind,0],'o',markersize=3,c=colors[n])
  ind = np.where(np.floor(wv_time_H)==year)[0]
  for i in ind:
    time = np.arange(wv_time_H[i,0]-wv_time_H[i,1],wv_time_H[i,0]+wv_time_H[i,1],1/365.25)
    interped = np.interp(time,terminus_time_H,terminus_val_H)
    plt.plot(np.mean(interped)/1e3,wv_dH_H[i,0],'^',markersize=3,c=colors[n])
  n=n+1
plt.errorbar(-200,-200,fmt='w.',label=' ')
plt.errorbar(-200,-200,fmt='w.',label=' ')
plt.xlabel('Terminus position (km)',fontsize=11)
plt.ylabel('dH/dt (m yr$^{-1}$)',fontsize=11)
plt.yticks(np.arange(-120,200,60),fontsize=11)
plt.xticks(np.arange(-4,5,1),fontsize=11)
plt.xlim([-2,2])
plt.ylim([-130,130])
plt.legend(loc=2,fontsize=11,numpoints=1,ncol=3,handlelength=0.2,handletextpad=0.5,labelspacing=0.2,columnspacing=0.8)
plt.tight_layout()
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Helheim_terminus_thinning.pdf"),FORMAT='PDF',dpi=600)
plt.close()

plt.figure(figsize=(3,3))
matplotlib.rc('font',family='Arial')
plt.plot([0,0],[-150,150],'k')
plt.plot([-4,4],[0,0],'k')
plt.errorbar(-200,-200,fmt='ko',markersize=3,label='Flux')
plt.errorbar(-200,-200,fmt='k^',markersize=3,label='DEM')
years = range(2008,2016)
colors=['b','skyblue','g','limegreen','gold','orange','r','maroon']
n=0
for year in years:
  ind = np.where(np.floor(flux_time_K)==year)[0]
  interped = np.interp(flux_time_K[ind],terminus_time_K,terminus_val_K)
  plt.errorbar(interped/1e3,flux_dH_K[ind,0],fmt='o',yerr=flux_dH_K[ind,1],markersize=3,c=colors[n])
  ind = np.where(np.floor(wv_time_K)==year)[0]
  for i in ind:
    time = np.arange(wv_time_K[i,0]-wv_time_K[i,1],wv_time_K[i,0]+wv_time_K[i,1],1/365.25)
    interped = np.interp(time,terminus_time_K,terminus_val_K)
    plt.errorbar(np.mean(interped)/1e3,wv_dH_K[i,0],fmt='^',yerr=wv_dH_K[i,1],markersize=3,c=colors[n])
  plt.plot(-200,-200,'s',markersize=3,color=colors[n],label=year)
  n=n+1
plt.errorbar(-200,-200,fmt='w.',label=' ')
plt.errorbar(-200,-200,fmt='w.',label=' ')
plt.xlabel('Terminus position (km)',fontsize=11)
plt.ylabel('dH/dt (m yr$^{-1}$)',fontsize=11)
plt.yticks(np.arange(-80,200,40),fontsize=11)
plt.xticks(np.arange(-4,5,1),fontsize=11)
plt.xlim([-4,4])
plt.ylim([-80,80])
plt.legend(loc=2,fontsize=11,numpoints=1,ncol=3,handlelength=0.2,handletextpad=0.5,labelspacing=0.2,columnspacing=0.8)
plt.tight_layout()
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Kanger_terminus_thinning.pdf"),FORMAT='PDF',dpi=600)
plt.close()



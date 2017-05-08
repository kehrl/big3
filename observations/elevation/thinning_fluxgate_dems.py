# This script compares elevation change from the DEMs to anticipated thinning rates from 
# the fluxgate method.

# LMK, UW, 8/31/2015

import os
import sys
import glaclib, zslib, icefrontlib, fluxlib, climlib, vellib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib

###########
# dem DEMs #
###########

# Load DEMs	
xdem_H,ydem_H,zdem_H,timedem_H,errordem_H = zslib.dem_grid('Helheim',285000.0,320000.0,-2588000.0,-2566000.0,years='all',verticaldatum='ellipsoid',method='nearest',return_error=True)
xdem_K,ydem_K,zdem_K,timedem_K,errordem_K = zslib.dem_grid('Kanger',449800.0,503000.0,-2302000.0,-2266000.0,years='all',verticaldatum='ellipsoid',method='nearest',return_error=True)

#############################################################
# Load terminus positions, velocities, & surface elevations #
#############################################################

# Load flowlines
x_H,y_H,zb_H,dists_H = glaclib.load_flowline('Helheim',shapefilename='flowline_flightline')
x_K,y_K,zb_K,dists_K = glaclib.load_flowline('Kanger',shapefilename='flowline_flightline')

# Load terminus positions
terminus_val_H, terminus_time_H = icefrontlib.distance_along_flowline(x_H,y_H,dists_H,'Helheim',type='icefront',time1=2008.,time2=2016.)
terminus_val_K, terminus_time_K = icefrontlib.distance_along_flowline(x_K,y_K,dists_K,'Kanger',type='icefront',time1=2008.,time2=2016.)

# Get indices along flowlines for velocities, elevations
dists_eul = -1*np.array([2.,5.,7,10.,15.])

ind_eul_H=[]
for i in range(0,len(dists_eul)):
  ind_eul_H.append( (abs(dists_H - dists_eul[i]*1e3)).argmin() )
ind_eul_K=[]
for i in range(0,len(dists_eul)):
  ind_eul_K.append( (abs(dists_K - dists_eul[i]*1e3)).argmin() )

# Load velocities
vel_val_H,vel_time_H,vel_error_H = vellib.velocity_at_eulpoints(x_H[ind_eul_H],y_H[ind_eul_H],'Helheim',data='TSX')
vel_val_K,vel_time_K,vel_error_K = vellib.velocity_at_eulpoints(x_K[ind_eul_K],y_K[ind_eul_K],'Kanger',data='TSX')

# Get longitudinal strain rates for flowline in flux box
indbox_H = np.where((x_H > 303940.) & (x_H < 307012.))[0]
vel_val_H_sec,vel_time_H_sec,vel_error_H_sec = vellib.velocity_at_eulpoints(x_H[indbox_H],y_H[indbox_H],'Helheim',data='TSX')
strain_H = np.zeros([len(vel_time_H_sec),2])
for i in range(0,len(vel_time_H_sec)):
  strain_H[i,0] = np.mean(np.diff(vel_val_H_sec[i,:])/(dists_H[1]-dists_H[0]))
  strain_H[i,1] = np.mean(vel_val_H_sec)*0.02/(dists_H[indbox_H[-1]]-dists_H[indbox_H[0]])
del vel_val_H_sec,vel_time_H_sec,vel_error_H_sec

indbox_K = np.where((x_K > 483686.) & (x_K < 485974.))[0]
#indbox_K = np.where((x_K > 486758) & (x_K < 489415))[0]
vel_val_K_sec,vel_time_K_sec,vel_error_K_sec = vellib.velocity_at_eulpoints(x_K[indbox_K],y_K[indbox_K],'Kanger',data='TSX')
strain_K = np.zeros([len(vel_time_K_sec),2])
for i in range(0,len(vel_time_K_sec)):
  strain_K[i,0] = np.mean(np.diff(vel_val_K_sec[i,:])/(dists_K[1]-dists_K[0]))
  strain_K[i,1] = np.std(np.diff(vel_val_K_sec[i,:])/(dists_K[1]-dists_K[0])) 

del vel_val_K_sec,vel_time_K_sec,vel_error_K_sec

# Load elevations
#zpt_atm_H,zptstd_atm_H,time_atm_H = zslib.atm_at_pts(x_H[ind_eul_H],y_H[ind_eul_H],'Helheim',years='all',maxdist=500.,verticaldatum='geoid',method='average',cutoff='terminus')
zpt_dem_H,zpterror_dem_H,time_dem_H = zslib.dem_at_pts(x_H[ind_eul_H],y_H[ind_eul_H],'Helheim',years='all',verticaldatum='geoid',cutoff='terminus',method='average',radius=500.)

# Calculate average surface slope for flowline in flux box
zpt_dem_H_sec,zpterror_dem_H_sec,time_dem_H_sec = zslib.dem_at_pts(x_H[indbox_H],y_H[indbox_H],'Helheim',years='all',verticaldatum='ellipsoid',cutoff='terminus',method='average',radius=500.)
slope_H = np.zeros([len(time_dem_H_sec),3])
slope_H[:,0] = time_dem_H_sec
slope_H[:,1] = (zpt_dem_H_sec[:,0]-zpt_dem_H_sec[:,-1])/(dists_H[indbox_H[-1]]-dists_H[indbox_H[0]])
slope_H[:,2] = 2*zpterror_dem_H_sec/(dists_H[indbox_H[-1]]-dists_H[indbox_H[0]])
del zpt_dem_H_sec,time_dem_H_sec,zpterror_dem_H_sec

#zpt_atm_K,zptstd_atm_K,time_atm_K = zslib.atm_at_pts(x_K[ind_eul_K],y_K[ind_eul_K],'Kanger',years='all',maxdist=500.,verticaldatum='geoid',method='average',cutoff='terminus')
zpt_dem_K,zpterror_dem_K,time_dem_K = zslib.dem_at_pts(x_K[ind_eul_K],y_K[ind_eul_K],'Kanger',years='all',verticaldatum='geoid',cutoff='terminus',method='average',radius=500.)

# Calculate average surface slope for flowline in flux box
zpt_dem_K_sec,zpterror_dem_K_sec,time_dem_K_sec = zslib.dem_at_pts(x_K[indbox_K],y_K[indbox_K],'Kanger',years='all',verticaldatum='ellipsoid',cutoff='terminus',method='average',radius=500.)
slope_K = np.zeros([len(time_dem_K_sec),3])
slope_K[:,0] = time_dem_K_sec
slope_K[:,1] = (zpt_dem_K_sec[:,0]-zpt_dem_K_sec[:,-1])/(dists_K[indbox_K[-1]]-dists_K[indbox_K[0]])
slope_K[:,2] = 2*zpterror_dem_K_sec/(dists_K[indbox_K[-1]]-dists_K[indbox_K[0]])
del zpt_dem_K_sec,time_dem_K_sec,zpterror_dem_K_sec


##################
# Load fluxgates #
##################

xgate_H,ygate_H = fluxlib.fluxbox_geometry('Helheim',"fluxgate1")
xgate_K,ygate_K = fluxlib.fluxbox_geometry('Kanger',"fluxgate1")

#####################
# Climate variables #
#####################

filt_len = 11.

xrac_H,yrac_H,smbrac_H,timerac_H = climlib.racmo_at_pts(np.mean(xgate_H),np.mean(ygate_H),'smb',filt_len=filt_len)
xrac_K,yrac_K,smbrac_K,timerac_K = climlib.racmo_at_pts(np.mean(xgate_K),np.mean(ygate_K),'smb',filt_len=filt_len)

xrac_H,yrac_H,runrac_H,timerac_H = climlib.racmo_at_pts(np.mean(xgate_H),np.mean(ygate_H),'runoff',filt_len=filt_len)
xrac_K,yrac_K,runrac_K,timerac_K = climlib.racmo_at_pts(np.mean(xgate_K),np.mean(ygate_K),'runoff',filt_len=filt_len)

xrac_H,yrac_H,zsrac_H,timezsrac_H = climlib.racmo_at_pts(np.mean(xgate_H)-10e3,np.mean(ygate_H)+10e3,'zs',filt_len='none')
xrac_K,yrac_K,zsrac_K,timezsrac_K = climlib.racmo_at_pts(np.mean(xgate_K),np.mean(ygate_K),'zs',filt_len='none')

yH_runoff,day1H_runoff,day2H_runoff,meltlengthH_runoff,totalH_runoff = climlib.seasonlength(timerac_H,runrac_H,'runoff')
yH_smb,day1H_smb,day2_smb,meltlengthH_smb,totalH_smb = climlib.seasonlength(timerac_H,smbrac_H,'smb')

yK_runoff,day1K_runoff,day2K_runoff,meltlengthK_runoff,totalK_runoff = climlib.seasonlength(timerac_K,runrac_K,'runoff')
yK_smb,day1K_smb,day2_smb,meltlengthK_smb,totalK_smb = climlib.seasonlength(timerac_K,smbrac_K,'smb')


###############################################
# dH/dt from the fluxgate, dem, and comparison #
###############################################

# Note comparison doesn't work right now because fluxgate thinning rates have no 
# uncertainty estimates

# Helheim
time_H,dH_H,Q_H,hbar_H,ubar_H,smb_H,width_H,area_H = fluxlib.fluxgate_thinning('Helheim','fluxgate1','cresis',10.,timing='velocity')
time_H_zs,dH_H_zs,Q_H_zs,hbar_H_zs,ubar_H_zs,smb_H_zs,width_H_zs,area_H_zs = fluxlib.fluxgate_thinning('Helheim','fluxgate1','cresis',10.,timing='elevation')
dem_time_H,dem_dH_H = fluxlib.dem_thinning('Helheim',xdem_H,ydem_H,zdem_H,timedem_H,errordem_H,"fluxgate1")

# Kanger
time_K,dH_K,Q_K,hbar_K,ubar_K,smb_K,width_K,area_K = fluxlib.fluxgate_thinning('Kanger','fluxgate1','cresis',10.,timing='velocity')
time_K_zs,dH_K_zs,Q_K_zs,hbar_K_zs,ubar_K_zs,smb_K_zs,width_K_zs,area_K_zs = fluxlib.fluxgate_thinning('Kanger','fluxgate1','cresis',10.,timing='elevation')
dem_time_K,dem_dH_K = fluxlib.dem_thinning('Kanger',xdem_K,ydem_K,zdem_K,timedem_K,errordem_K,"fluxgate1")

# Compare thinning rates
dH_time_H,dH_flux_H,dH_dem_H = fluxlib.compare_thinning_rates(dem_time_H,dem_dH_H,time_H,dH_H,smb_H)
dH_time_K,dH_flux_K,dH_dem_K = fluxlib.compare_thinning_rates(dem_time_K,dem_dH_K,time_K,dH_K,smb_K)

plt.figure(figsize=(4,4))
matplotlib.rc('font',family='Arial')
plt.plot([-110,100],[-110,100],c='k',lw=1)
plt.errorbar(dH_dem_K[:,0]/365.25*100,dH_flux_K[:,0]/365.25*100,yerr=dH_flux_K[:,1]/365.25*100,xerr=dH_dem_K[:,1]/365.25*100,fmt='o',color='0.7',markersize=3,zorder=3,label='Kangerdlugssuaq')
plt.errorbar(dH_dem_H[:,0]/365.25*100,dH_flux_H[:,0]/365.25*100,yerr=dH_flux_H[:,1]/365.25*100,xerr=dH_dem_H[:,1]/365.25*100,fmt='ro',markersize=3,zorder=3,label='Helheim')
plt.xticks(np.arange(-30,30,10),fontsize=8)
plt.yticks(np.arange(-30,30,10),fontsize=8)
plt.ylim([-35,22])
plt.xlim([-35,22])
plt.ylabel('Flux-gate $dh/dt$ (cm d$^{-1}$)',fontsize=9)
plt.xlabel('DEM $dh/dt$ (cm d$^{-1}$)',fontsize=9)
plt.tight_layout()
plt.legend(loc=4,borderpad=0.3,fontsize=8,numpoints=1,handlelength=0.7,labelspacing=0.05,columnspacing=0.7,handletextpad=0.5)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+"Thinning_comparison.pdf"),FORMAT='PDF',dpi=600)
plt.close()

########################################
# Plot components of fluxgate analysis #
########################################

plt.figure(figsize=(3.75,4.5))
gs = matplotlib.gridspec.GridSpec(5,1)
matplotlib.rc('font',family='Arial')

time1 = 2008.
time2 = 2016.

# Helheim

plt.subplot(gs[0:2])
ax = plt.gca()
plt.plot([time1,time2],[0,0],'k')
plt.plot(timerac_H,smbrac_H/900.*100,'0.4',lw=1.2,label='SMB')
#plt.plot(timezsrac_H[1:],np.convolve(np.diff(zsrac_H)/np.diff(timezsrac_H),np.ones((11,))/11.,mode='same'),c='0.4',lw=1.2,label='SMB')
#plt.plot(time_H,smb_H,c='0.5',lw=0.75,label='SMB')
plt.errorbar(dem_time_H[:,0],dem_dH_H[:,0]/365.25*100,xerr=dem_time_H[:,1],yerr=dem_dH_H[:,1]/365.25*100,fmt='bo',mec='k',mew=0.5,markersize=3,capsize=1,lw=0.5,label='DEM')
#plt.errorbar(time_H_2,dH_H_2,label='Helheim',fmt='r.',markersize=4,capsize=1,lw=0.5)
plt.errorbar(time_H,(dH_H[:,0]+smb_H)/365.25*100,yerr=3.,fmt='rs',mec='k',mew=0.5,markersize=3,capsize=1,lw=0.5,label='Flux-gate')
#time_H_2,dH_H_2,Q_H_2,hbar_H_2,ubar_H_2 = fluxlib.fluxgate_thinning('Helheim','fluxgate1','cresis',10.,timing='velocity')
#plt.plot(timezsrac_K[1:],np.diff(zsrac_K)/np.diff(timezsrac_K),c='r',lw=0.5)
plt.xticks(range(2008,2017),fontsize=8,fontname='Arial')
plt.xlim([time1,time2])
plt.ylabel(r'$dh/dt$ (cm d$^{-1}$)',fontname='Arial',fontsize=9)
plt.yticks(np.arange(-40,40,10),fontsize=8,fontname='Arial')
plt.ylim([-43,32])
ax.set_xticklabels([])
#xTickPos = np.linspace(np.floor(time1)-0.25,np.ceil(time2)-0.25,(np.ceil(time2)-np.floor(time1))*2+1)
#ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.85','w'],linewidth=0)
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1H_runoff[i],-160],[day1H_runoff[i],160],[day2H_runoff[i],160],[day2H_runoff[i],-160]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
  ax.add_patch(patch)
plt.text(time1+0.2,25,'a',fontsize=9,fontname='Arial',fontweight='bold')  
plt.legend(loc=1,numpoints=1,ncol=3,handletextpad=0.2,fontsize=8,columnspacing=0.05,markerscale=1,handlelength=1.2,borderpad=0.2)

plt.subplot(gs[4])
ax = plt.gca()
#plt.plot(time_H_2,Q_H_2[:,0]*917./1e12,'k.-',markersize=4,label = 'Flux in')
#plt.plot([2008,2017],np.array([30,30])*1e9*917/1e12,'k',lw=1.5,label='Balance')
ind = np.where(~(np.isnan(Q_H[:,0])))[0]
plt.plot(time_H[ind],Q_H[ind,0]/1e9,'b.',markersize=5,label = 'Inflow')
ind = np.where(~(np.isnan(Q_H[:,1])))[0]
plt.plot(time_H[ind],Q_H[ind,1]/1e9,'r.',markersize=5, label = 'Outflow')
plt.ylabel('Ice flux \n'+r'(km$^3$ yr$^{-1}$)',fontsize=9,fontname='Arial')
plt.xticks(range(2008,2017))
plt.xlim([time1,time2])
plt.yticks(np.arange(20,42,4),fontsize=8,fontname='Arial')
plt.ylim([25,42])
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1H_runoff[i],-100],[day1H_runoff[i],100],[day2H_runoff[i],100],[day2H_runoff[i],-100]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
  ax.add_patch(patch)
plt.text(time1+0.2,39,'d',fontsize=9,fontname='Arial',fontweight='bold')
plt.legend(loc=1,numpoints=1,ncol=3,handletextpad=0.2,fontsize=8,columnspacing=0.4,handlelength=0.3,borderpad=0.2)
labels=[]
for i in range(2000,2017):
  labels.append('Jan \n'+str(i))
plt.xticks(range(2000,2017))
plt.xlim([time1,time2])
ax.set_xticklabels(labels,fontsize=8,fontname='Arial')

plt.subplot(gs[2])
ax = plt.gca()
#plt.errorbar(time_H,(ubar_H[:,1]-ubar_H[:,0])/1e3,yerr=(ubar_H[:,1]+ubar_H[:,0])/1e3*0.03/np.sqrt(2000./100),fmt='ko',markersize=2.5,capsize=1,lw=0.5)
plt.errorbar(time_H,(ubar_H[:,1]-ubar_H[:,0])/1e3,yerr=0,fmt='ko',markersize=2.5,capsize=1,lw=0.5)
#plt.ylabel(r'$\bar{u}_{in}-\bar{u}_{out}$'+'\n'+r'(km yr$^{-1}$)',fontsize=9,fontname='Arial')
plt.ylabel(r'$\bar{u}_{o}-\bar{u}_{i}$'+'\n'+r'(km yr$^{-1}$)',fontsize=9,fontname='Arial')
plt.yticks(np.arange(0,4,0.5),fontsize=8,fontname='Arial')
#plt.ylim([1.1,2.8])
plt.ylim([0.4,1.6])
plt.xticks(range(2008,2017))
plt.xlim([time1,time2])
ax.set_xticklabels([])
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1H_runoff[i],-100],[day1H_runoff[i],100],[day2H_runoff[i],100],[day2H_runoff[i],-100]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
  ax.add_patch(patch)
plt.text(time1+0.2,1.4,'b',fontsize=9,fontname='Arial',fontweight='bold')
#plt.legend(loc=1,numpoints=1,ncol=2,handletextpad=0.1,fontsize=8,columnspacing=0.05,handlelength=1)

plt.subplot(gs[3])
ax = plt.gca()
#plt.errorbar(time_H_zs,(hbar_H_zs[:,0]-hbar_H_zs[:,1]),yerr=hbar_H_zs[:,2]+hbar_H_zs[:,3],fmt='ko',markersize=2.5,capsize=1,lw=0.5)
plt.errorbar(time_H_zs,(hbar_H_zs[:,0]-hbar_H_zs[:,1]),yerr=0,fmt='ko',markersize=2.5,capsize=1,lw=0.5)
#plt.ylabel(r'$\bar{H}_in-\bar{H}_out$'+'\n (m)',fontsize=9,fontname='Arial')
plt.ylabel(r'$\bar{H}_{i}-\bar{H}_{o}$'+'\n (m)',fontsize=9,fontname='Arial')
plt.xticks(range(2008,2017))
#plt.yticks(np.arange(15,45,5),fontsize=8,fontname='Arial')
plt.yticks(np.arange(0,120,5),fontsize=8,fontname='Arial')
plt.ylim([80,93])
#plt.ylim([17,32])
plt.xlim([time1,time2])
ax.set_xticklabels([])
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1H_runoff[i],-100],[day1H_runoff[i],200],[day2H_runoff[i],200],[day2H_runoff[i],-100]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
  ax.add_patch(patch)
plt.text(time1+0.2,90.5,'c',fontsize=9,fontname='Arial',fontweight='bold')
ax.set_xticklabels([])

plt.tight_layout()
plt.subplots_adjust(hspace=0.05,wspace=-0.2,top=0.98,right=0.95,left=0.16,bottom=0.075) 
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Helheim_fluxgate_components.pdf"),FORMAT='PDF',dpi=600)
plt.close()

# Kanger

plt.figure(figsize=(3.75,4.5))
gs = matplotlib.gridspec.GridSpec(5,1)
matplotlib.rc('font',family='Arial')

time1 = 2008.
time2 = 2016.

plt.subplot(gs[0:2])
ax = plt.gca()
plt.plot([time1,time2],[0,0],'k')
plt.plot(timerac_K,smbrac_K/900.*100,'0.4',lw=1.2,label='SMB')
plt.errorbar(dem_time_K[:,0],dem_dH_K[:,0]/365.25*100,xerr=dem_time_K[:,1],yerr=dem_dH_K[:,1]/365.25*100,fmt='bo',mec='k',mew=0.5,markersize=3,capsize=1,lw=0.5,label='DEM')
#plt.errorbar(time_K_2,dH_K_2,label='Helheim',fmt='r.',markersize=4,capsize=1,lw=0.5)
plt.errorbar(time_K,(dH_K[:,0]+smb_K)/365.25*100,yerr=5.,fmt='rs',mec='k',mew=0.5,markersize=3,capsize=1,lw=0.5,label='Flux')
#plt.errorbar(time_K_zs,dH_K_zs[:,0]+smb_K_zs,fmt='bo',markersize=2.5,capsize=1,lw=0.5,label='Flux')
#time_K_2,dH_K_2,Q_K_2,hbar_K_2,ubar_K_2 = fluxlib.fluxgate_thinning('Helheim','fluxgate1','cresis',10.,timing='velocity')
#plt.plot(timezsrac_K[1:],np.diff(zsrac_K)/np.diff(timezsrac_K),c='r',lw=0.5)
plt.xticks(range(2008,2017),fontsize=8,fontname='Arial')
plt.xlim([time1,time2])
plt.ylabel(r'$dh/dt$ (cm d$^{-1}$)',fontname='Arial',fontsize=9)
plt.yticks(np.arange(-40,40,10),fontsize=8,fontname='Arial')
plt.ylim([-37,32])
ax.set_xticklabels([])
#xTickPos = np.linspace(np.floor(time1)-0.25,np.ceil(time2)-0.25,(np.ceil(time2)-np.floor(time1))*2+1)
#ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.85','w'],linewidth=0)
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1K_runoff[i],-160],[day1K_runoff[i],160],[day2K_runoff[i],160],[day2K_runoff[i],-160]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
  ax.add_patch(patch)
plt.text(time1+0.2,25,'a',fontsize=9,fontname='Arial',fontweight='bold')  
plt.legend(loc=1,numpoints=1,ncol=3,handletextpad=0.2,fontsize=8,columnspacing=0.05,markerscale=1,handlelength=1.2,borderpad=0.2)

plt.subplot(gs[4])
ax = plt.gca()
#plt.plot(time_K_2,Q_K_2[:,0]*917./1e12,'k.-',markersize=4,label = 'Flux in')
#plt.plot([2008,2017],np.array([22.9,22.9])*1e9*917/1e12,'k',lw=1.5,label='Balance')
ind = np.where(~(np.isnan(Q_K[:,0])))[0]
plt.plot(time_K[ind],Q_K[ind,0]/1e9,'b.',markersize=5,label = 'Inflow')
ind = np.where(~(np.isnan(Q_K[:,1])))[0]
plt.plot(time_K[ind],Q_K[ind,1]/1e9,'r.',markersize=5, label = 'Outflow')
plt.ylabel('Ice flux \n'+r'(km$^{3}$ yr$^{-1}$)',fontsize=9,fontname='Arial')
plt.xticks(range(2008,2017))
plt.xlim([time1,time2])
plt.yticks(np.arange(20,40,4),fontsize=8,fontname='Arial')
plt.ylim([23,33.5])
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1K_runoff[i],-100],[day1K_runoff[i],100],[day2K_runoff[i],100],[day2K_runoff[i],-100]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
  ax.add_patch(patch)
plt.text(time1+0.2,31.5,'d',fontsize=9,fontname='Arial',fontweight='bold')
plt.legend(loc=1,numpoints=1,ncol=3,handletextpad=0.2,fontsize=8,columnspacing=0.4,handlelength=0.3,borderpad=0.2)
labels=[]
for i in range(2000,2017):
  labels.append('Jan \n'+str(i))
plt.xticks(range(2000,2017))
plt.xlim([time1,time2])
ax.set_xticklabels(labels,fontsize=8,fontname='Arial')

plt.subplot(gs[2])
ax = plt.gca()
#plt.errorbar(time_K,(ubar_K[:,1]-ubar_K[:,0])/1e3,yerr=(ubar_K[:,1]+ubar_K[:,0])/1e3*0.03/np.sqrt(2000./100),fmt='ko',markersize=2.5,capsize=1,lw=0.5)
plt.errorbar(time_K,(ubar_K[:,1]-ubar_K[:,0])/1e3,yerr=0,fmt='ko',markersize=2.5,capsize=1,lw=0.5)
plt.ylabel(r'$\bar{u}_{o}-\bar{u}_{i}$'+'\n'+r'(km yr$^{-1}$)',fontsize=9,fontname='Arial')
plt.yticks(np.arange(0,4,0.2),fontsize=8,fontname='Arial')
#plt.ylim([0.1,0.6])
plt.ylim([1.3,2.0])
plt.xticks(range(2008,2017))
plt.xlim([time1,time2])
ax.set_xticklabels([])
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1K_runoff[i],-100],[day1K_runoff[i],100],[day2K_runoff[i],100],[day2K_runoff[i],-100]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
  ax.add_patch(patch)
plt.text(time1+0.2,1.87,'b',fontsize=9,fontname='Arial',fontweight='bold')
#plt.legend(loc=1,numpoints=1,ncol=2,handletextpad=0.1,fontsize=8,columnspacing=0.05,handlelength=1)

plt.subplot(gs[3])
ax = plt.gca()
#plt.errorbar(time_K_zs,(hbar_K_zs[:,0]-hbar_K_zs[:,1]),yerr=hbar_K_zs[:,2]+hbar_K_zs[:,3],fmt='ko',markersize=2.5,capsize=1,lw=0.5)
plt.errorbar(time_K_zs,(hbar_K_zs[:,0]-hbar_K_zs[:,1]),yerr=0,fmt='ko',markersize=2.5,capsize=1,lw=0.5)
plt.ylabel(r'$\bar{H}_{i}-\bar{H}_{o}$'+'\n (m)',fontsize=9,fontname='Arial')
plt.xticks(range(2008,2017))
#plt.yticks(np.arange(0,1,0.001),fontsize=8,fontname='Arial')
#plt.ylim([0.0185,0.0235])
plt.yticks(np.arange(70,310,5),fontsize=8,fontname='Arial')
plt.ylim([75,90])
plt.xlim([time1,time2])
ax.set_xticklabels([])
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1K_runoff[i],-100],[day1K_runoff[i],100],[day2K_runoff[i],100],[day2K_runoff[i],-100]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
  ax.add_patch(patch)
plt.text(time1+0.2,87,'c',fontsize=9,fontname='Arial',fontweight='bold')
ax.set_xticklabels([])

plt.tight_layout()
plt.subplots_adjust(hspace=0.05,wspace=-0.2,top=0.98,right=0.95,left=0.16,bottom=0.075) 
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Kanger_fluxgate_components.pdf"),FORMAT='PDF',dpi=600)
plt.close()


###########################################################
# Plot terminus vs. thinning rate for individual glaciers #
###########################################################

fig=plt.figure(figsize=(3,2.75))
matplotlib.rc('font',family='Arial')
plt.plot([0,0],[-50,30],'k',zorder=0)
plt.plot([-4,4],[0,0],'k',zorder=0)
#plt.errorbar(-200,-200,fmt='ko',markersize=3,label='Flux')
#plt.errorbar(-200,-200,fmt='k^',markersize=3,label='DEM')

ind = np.where((time_H > time1) & (time_H < time2))[0]
interped = np.interp(time_H[ind],terminus_time_H,terminus_val_H)
plt.errorbar(interped/1e3,dH_H[ind,0]/365.25*100,3.,fmt='k.',marker=None,mew=0,markersize=4,zorder=0)
plt.scatter(interped/1e3,dH_H[ind,0]/365.25*100,c=(time_H[ind]-np.floor(time_H[ind]))*365,marker='s',s=15,vmin=0,vmax=365,cmap='viridis')
#interped = np.interp(dem_time_H[:,0],terminus_time_H,terminus_val_H)
#plt.plot(interped/1e3,dH_dem_H[:,0]/365.25*100,'ko',markerfacecolor='b',markersize=4)

plt.xlabel('Terminus position (km)',fontsize=8)
plt.ylabel(r'Dynamic $dh/dt$ (cm d$^{-1}$)',fontsize=8)
plt.yticks(np.arange(-40,30,10),fontsize=8)
plt.xticks(np.arange(-4,5,1),fontsize=8)
plt.xlim([-4,4])
plt.ylim([-42,15])
#plt.legend(loc=2,fontsize=11,numpoints=1,ncol=3,handlelength=0.2,handletextpad=0.5,labelspacing=0.2,columnspacing=0.8)

plt.text(-3.6,10,'a',fontsize=9,fontname='Arial',fontweight='bold')
plt.text(0.9,-32,'Day of Year',fontsize=8)
cbaxes = fig.add_axes([0.6, 0.24, 0.34, 0.03]) 
cb = plt.colorbar(cax=cbaxes, ticks=[0,100,200,300], orientation='horizontal')
cb.ax.tick_params(labelsize=8)

plt.tight_layout()
plt.subplots_adjust(hspace=0.0,wspace=0,top=0.98,right=0.97,left=0.17,bottom=0.15) 
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Helheim_terminus_thinning.pdf"),FORMAT='PDF',dpi=600)
plt.close()

fig = plt.figure(figsize=(3,2.75))
matplotlib.rc('font',family='Arial')
plt.plot([0,0],[-50,30],'k',zorder=0)
plt.plot([-4,4],[0,0],'k',zorder=0)
#plt.errorbar(-200,-200,fmt='ko',markersize=3,label='Flux')
#plt.errorbar(-200,-200,fmt='k^',markersize=3,label='DEM')

ind = np.where((time_K > time1) & (time_K < time2))[0]
interped = np.interp(time_K[ind],terminus_time_K,terminus_val_K)
plt.errorbar(interped/1e3,dH_K[ind,0]/365.25*100,3.,fmt='k.',marker=None,mew=0,markersize=4,zorder=0)
plt.scatter(interped/1e3,dH_K[ind,0]/365.25*100,c=(time_K[ind]-np.floor(time_K[ind]))*365,marker='s',s=15,vmin=0,vmax=365,cmap='viridis')
#interped = np.interp(dem_time_K[:,0],terminus_time_K,terminus_val_K)
#plt.plot(interped/1e3,dH_dem_K[:,0]/365.25*100,'ko',markerfacecolor='b',markersize=4)

plt.text(-3.6,10,'b',fontsize=9,fontname='Arial',fontweight='bold')
plt.xlabel('Terminus position (km)',fontsize=8)
plt.ylabel(r'Dynamic $dh/dt$ (cm d$^{-1}$)',fontsize=8)
plt.yticks(np.arange(-40,30,10),fontsize=8)
plt.xticks(np.arange(-4,5,1),fontsize=8)
plt.xlim([-4,4])
plt.ylim([-42,15])

plt.text(0.9,-32,'Day of Year',fontsize=8)
cbaxes = fig.add_axes([0.6, 0.24, 0.34, 0.03]) 
cb = plt.colorbar(cax=cbaxes, ticks=[0,100,200,300], orientation='horizontal')
cb.ax.tick_params(labelsize=8)

#plt.legend(loc=2,fontsize=11,numpoints=1,ncol=3,handlelength=0.2,handletextpad=0.5,labelspacing=0.2,columnspacing=0.8)
plt.tight_layout()
plt.subplots_adjust(hspace=0.0,wspace=0,top=0.98,right=0.97,left=0.17,bottom=0.15) 
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Kanger_terminus_thinning.pdf"),FORMAT='PDF',dpi=600)
plt.close()

##############################
# Plot both glaciers at once #
##############################

plt.figure(figsize=(3.75,5))
gs = matplotlib.gridspec.GridSpec(6,1)

time1 = 2010.
time2 = 2016.

plt.subplot(gs[0,:])
ax = plt.gca()
nonnan = np.where(~(np.isnan(vel_val_H)))[0]
ind = np.argmin(abs(vel_time_H[nonnan]-time1))
plt.plot(vel_time_H,(vel_val_H[:,0]-vel_val_H[nonnan[ind],0])/1e3,'k.',label='Helheim',markersize=3)
nonnan = np.where(~(np.isnan(vel_val_K)))[0]
ind = np.argmin(abs(vel_time_K[nonnan]-time1))
plt.plot(vel_time_K,(vel_val_K[:,0]-vel_val_K[ind,0])/1e3,'r.',label='Kanger',markersize=3)
plt.ylabel('Speed \n (km yr$^{-1}$)',fontsize=9,fontname='Arial')
plt.yticks(range(-2,4),fontsize=8,fontname='Arial')
plt.ylim([-1.5,3])
plt.xticks(range(2008,2017))
plt.xlim([time1,time2])
ax.set_xticklabels([])
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1H_runoff[i],-100],[day1H_runoff[i],100],[day2H_runoff[i],100],[day2H_runoff[i],-100]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
  ax.add_patch(patch)
plt.legend(loc=1,numpoints=1,ncol=2,handletextpad=0.02,fontsize=8,columnspacing=0.05,markerscale=3,handlelength=1)
plt.text(time1+0.2,1.6,'a',fontsize=9,fontname='Arial',fontweight='bold')
  
plt.subplot(gs[1,:])
ax = plt.gca()
#plt.plot(vel_time_H,(vel_val_H[:,0]-vel_val_H[:,1])/(((dists_eul[0]-dists_eul[1])*1e3),'k.',label='Helheim',markersize=3)
#plt.plot(vel_time_K,(vel_val_K[:,0]-vel_val_K[:,1])/(dists_eul[0]-dists_eul[1])*1e3),'r.-',label='Kanger',markersize=3)
plt.plot(vel_time_H,strain_H[:,0],'k.',label='Helheim',markersize=4)
plt.plot(vel_time_K,strain_K[:,0],'r.',label='Kanger',markersize=4)
plt.ylabel('Strain rate \n (yr$^{-1}$)',fontsize=9,fontname='Arial')
plt.ylim([0.1,0.7])
plt.yticks(np.arange(0,0.8,0.2),fontsize=8,fontname='Arial')
plt.xticks(range(2008,2017))
plt.xlim([time1,time2])
ax.set_xticklabels([])
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1H_runoff[i],-100],[day1H_runoff[i],100],[day2H_runoff[i],100],[day2H_runoff[i],-100]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
  ax.add_patch(patch)
plt.text(time1+0.2,0.55,'b',fontsize=9,fontname='Arial',fontweight='bold')

plt.subplot(gs[2,:])
ax = plt.gca()
#nonnan = np.where(~(np.isnan(time_atm_H)))[0]
#ind = np.argmin(abs(time_atm_H[nonnan]-time1))
#ind2 = np.argmin(timezsrac_H-time_atm_H[nonnan[ind]])
#plt.plot(timezsrac_H,zsrac_H-zsrac_H[ind2],c='k',lw=0.5)
#plt.errorbar(time_atm_H,(zpt_atm_H[:,0]-zpt_atm_H[nonnan[ind],0]),fmt='k.',label='Helheim',markersize=3,capsize=1)
#plt.errorbar(time_dem_H,(zpt_dem_H[:,0]-zpt_atm_H[nonnan[ind],0]),yerr=zpterror_dem_H,fmt='k.',label='Helheim',markersize=3,capsize=1)
#nonnan = np.where(~(np.isnan(time_atm_K)))[0]
#ind = np.argmin(abs(time_atm_K[nonnan]-time1))
#ind2 = np.argmin(timezsrac_K-time_atm_K[nonnan[ind]])
#plt.plot(timezsrac_H,zsrac_K-zsrac_K[ind2],c='r',lw=0.5)
#plt.errorbar(time_atm_K,(zpt_atm_K[:,0]-zpt_atm_K[nonnan[ind],0]),fmt='r.',label='Kanger',markersize=3,capsize=1)
#plt.errorbar(time_dem_K,(zpt_dem_K[:,0]-zpt_atm_K[nonnan[ind],0]),yerr=zpterror_dem_K,fmt='r.',label='Kanger',markersize=3,capsize=1)
plt.ylabel('Elevation \n (m)',fontsize=9,fontname='Arial')
plt.xticks(range(2008,2017))
plt.xlim([time1,time2])
plt.yticks(np.arange(-40,40,20),fontsize=8,fontname='Arial')
plt.ylim([-60,25])
ax.set_xticklabels([])
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1H_runoff[i],-100],[day1H_runoff[i],100],[day2H_runoff[i],100],[day2H_runoff[i],-100]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
  ax.add_patch(patch)
plt.text(time1+0.2,10,'c',fontsize=9,fontname='Arial',fontweight='bold')


plt.subplot(gs[-3,:])
ax = plt.gca()
plt.plot([time1,time2],[0,0],'k')
plt.plot(timezsrac_H[1:],np.diff(zsrac_H)/np.diff(timezsrac_H),c='0.5',lw=0.75)
plt.errorbar(dem_time_H[:,0],dem_dH_H[:,0],xerr=dem_time_H[:,1],yerr=dem_dH_H[:,1],label='Helheim',fmt='k.',markersize=3,capsize=1,lw=0.5)
plt.errorbar(dem_time_K[:,0],dem_dH_K[:,0],xerr=dem_time_K[:,1],yerr=dem_dH_K[:,1],label='Helheim',fmt='r.',markersize=3,capsize=1,lw=0.5)
#plt.plot(timezsrac_K[1:],np.diff(zsrac_K)/np.diff(timezsrac_K),c='r',lw=0.5)
plt.xticks(range(2008,2017),fontsize=8,fontname='Arial')
plt.xlim([time1,time2])
plt.ylabel('dH/dt \n (m/yr)',fontname='Arial',fontsize=9)
plt.yticks(np.arange(-80,120,40),fontsize=8,fontname='Arial')
plt.ylim([-80,80])
ax.set_xticklabels([])
#xTickPos = np.linspace(np.floor(time1)-0.25,np.ceil(time2)-0.25,(np.ceil(time2)-np.floor(time1))*2+1)
#ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.85','w'],linewidth=0)
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1H_runoff[i],-100],[day1H_runoff[i],100],[day2H_runoff[i],100],[day2H_runoff[i],-100]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
  ax.add_patch(patch)
plt.text(time1+0.2,50,'d',fontsize=9,fontname='Arial',fontweight='bold')

plt.subplot(gs[-2,:])
ax = plt.gca()
plt.plot(terminus_time_H,terminus_val_H/1e3,'k.',lw=1.5,label='Helheim',markersize=3)
plt.plot(terminus_time_K,terminus_val_K/1e3,'r.',lw=1.5,label='Kanger',markersize=3)
plt.ylabel('Terminus \n (km)',fontsize=9,fontname='Arial')
plt.yticks(np.arange(-2,6,2),fontsize=8,fontname='Arial')
plt.ylim([-4,4])
plt.xticks(range(2008,2017),fontsize=8,fontname='Arial')
plt.xlim([time1,time2])
ax.set_xticklabels([])
#xTickPos = np.linspace(np.floor(time1)-0.25,np.ceil(time2)-0.25,(np.ceil(time2)-np.floor(time1))*2+1)
#ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.85','w'],linewidth=0)
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1H_runoff[i],-4],[day1H_runoff[i],4],[day2H_runoff[i],4],[day2H_runoff[i],-4]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
  ax.add_patch(patch)
plt.text(time1+0.2,2.5,'e',fontsize=9,fontname='Arial',fontweight='bold')


plt.subplot(gs[-1,:])
ax = plt.gca()
plt.plot(timerac_H,runrac_H,'k',lw=1,label='Helheim')
plt.plot(timerac_K,runrac_K,'r',lw=1,label='Kanger')
labels=[]
for i in range(2000,2017):
  labels.append('Jan \n'+str(i))
plt.xticks(range(2000,2017))
ax.set_xticklabels(labels,fontsize=8,fontname='Arial')
plt.xlim([time1,time2])
plt.ylim([0,60])
plt.yticks(np.arange(0,80,20),fontsize=8,fontname='Arial')
plt.ylabel('Runoff \n (kg m$^{-2}$ d$^{-1}$)',fontsize=9,fontname='Arial')
#xTickPos = np.linspace(np.floor(time1)-0.25,np.ceil(time2)-0.25,(np.ceil(time2)-np.floor(time1))*2+1)
#ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.85','w'],linewidth=0)
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1H_runoff[i],0],[day1H_runoff[i],60],[day2H_runoff[i],60],[day2H_runoff[i],0]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
  ax.add_patch(patch)
plt.text(time1+0.2,45,'f',fontsize=9,fontname='Arial',fontweight='bold')


plt.tight_layout()
plt.subplots_adjust(hspace=0.1,wspace=0) 
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/both_glaciers.pdf"),FORMAT='PDF')
plt.close()

#################### 
# New plot for all #
####################

plt.figure(figsize=(7.5,5))
gs = matplotlib.gridspec.GridSpec(6,2)

time1 = 2010.
time2 = 2016.

plt.subplot(gs[1,0])
ax = plt.gca()
plt.plot(vel_time_H,vel_val_H[:,0]/1e3,'b.',label=str(int((dists_eul[0])))+' km',markersize=3)
plt.plot(vel_time_H,vel_val_H[:,1]/1e3,'k.',label=str(int((dists_eul[1])))+' km',markersize=3)
plt.plot(vel_time_H,vel_val_H[:,2]/1e3,'r.',label=str(int((dists_eul[2])))+' km',markersize=3)
plt.ylabel('Speed \n (km yr$^{-1}$)',fontsize=9,fontname='Arial')
plt.yticks(range(5,11),fontsize=8,fontname='Arial')
plt.ylim([4.5,10])
plt.xticks(range(2008,2017))
plt.xlim([time1,time2])
ax.set_xticklabels([])
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1H_runoff[i],-100],[day1H_runoff[i],100],[day2H_runoff[i],100],[day2H_runoff[i],-100]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
  ax.add_patch(patch)
plt.legend(loc=1,numpoints=1,ncol=3,handletextpad=0.02,fontsize=8,columnspacing=0.05,markerscale=2,handlelength=1,borderpad=0.2)
#plt.text(time1+0.2,9.5,'a',fontsize=9,fontname='Arial',fontweight='bold')

plt.subplot(gs[1,1])
ax = plt.gca()
plt.plot(vel_time_K,vel_val_K[:,0]/1e3,'b.',label=str(int((dists_eul[0])))+' km',markersize=3)
plt.plot(vel_time_K,vel_val_K[:,1]/1e3,'k.',label=str(int((dists_eul[1])))+' km',markersize=3)
plt.plot(vel_time_K,vel_val_K[:,2]/1e3,'r.',label=str(int((dists_eul[2])))+' km',markersize=3)
plt.yticks(range(5,12),fontsize=8,fontname='Arial')
plt.ylim([5.2,10.7])
plt.xticks(range(2008,2017))
plt.xlim([time1,time2])
ax.set_xticklabels([])
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1K_runoff[i],-100],[day1K_runoff[i],100],[day2K_runoff[i],100],[day2K_runoff[i],-100]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
  ax.add_patch(patch)
#plt.legend(loc=1,numpoints=1,ncol=2,handletextpad=0.02,fontsize=8,columnspacing=0.05,markerscale=3,handlelength=1)
#plt.text(time1+0.2,9.5,'a',fontsize=9,fontname='Arial',fontweight='bold')

  
plt.subplot(gs[2,0])
ax = plt.gca()
plt.plot(vel_time_H,(vel_val_H[:,0]-vel_val_H[:,1])/((dists_eul[0]-dists_eul[1])*1e3),'b.',label='5-7 km',markersize=3)
plt.plot(vel_time_H,(vel_val_H[:,1]-vel_val_H[:,2])/((dists_eul[1]-dists_eul[2])*1e3),'r.',label='7-10 km',markersize=3)
#plt.plot(vel_time_H,strain_H[:,0],'k.',label='Helheim',markersize=4)
plt.ylabel('Strain rate \n (yr$^{-1}$)',fontsize=9,fontname='Arial')
plt.yticks(np.arange(0,0.8,0.2),fontsize=8,fontname='Arial')
plt.ylim([-0.05,0.75])
plt.xticks(range(2008,2017))
plt.xlim([time1,time2])
ax.set_xticklabels([])
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1H_runoff[i],-100],[day1H_runoff[i],100],[day2H_runoff[i],100],[day2H_runoff[i],-100]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
  ax.add_patch(patch)
plt.legend(loc=1,numpoints=1,ncol=3,handletextpad=0.02,fontsize=8,columnspacing=0.05,markerscale=2,handlelength=1,borderpad=0.2)
#plt.text(time1+0.2,0.55,'b',fontsize=9,fontname='Arial',fontweight='bold')

plt.subplot(gs[2,1])
ax = plt.gca()
#plt.plot(vel_time_H,(vel_val_H[:,0]-vel_val_H[:,1])/(((dists_eul[0]-dists_eul[1])*1e3),'k.',label='Helheim',markersize=3)
plt.plot(vel_time_K,(vel_val_K[:,0]-vel_val_K[:,1])/((dists_eul[0]-dists_eul[1])*1e3),'b.',label='Kanger',markersize=3)
plt.plot(vel_time_K,(vel_val_K[:,1]-vel_val_K[:,2])/((dists_eul[1]-dists_eul[2])*1e3),'r.',label='Kanger',markersize=3)
#plt.plot(vel_time_K,strain_K[:,0],'k.',label='Kanger',markersize=4)
plt.yticks(np.arange(0,1.2,0.2),fontsize=8,fontname='Arial')
plt.ylim([0.2,1])
plt.xticks(range(2008,2017))
plt.xlim([time1,time2])
ax.set_xticklabels([])
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1K_runoff[i],-100],[day1K_runoff[i],100],[day2K_runoff[i],100],[day2K_runoff[i],-100]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
  ax.add_patch(patch)
#plt.text(time1+0.2,0.55,'b',fontsize=9,fontname='Arial',fontweight='bold')

plt.subplot(gs[3,0])
ax = plt.gca()
#nonnan = np.where(~(np.isnan(time_atm_H)))[0]
#ind = np.argmin(abs(time_atm_H[nonnan]-time1))
#ind2 = np.argmin(timezsrac_H-time_atm_H[nonnan[ind]])
#plt.plot(timezsrac_H,zsrac_H-zsrac_H[ind2],c='k',lw=0.5)
#plt.errorbar(time_atm_H,(zpt_atm_H[:,0]-zpt_atm_H[nonnan[ind],0]),fmt='b.',label='Helheim',markersize=3,capsize=1)
#plt.errorbar(time_dem_H,(zpt_dem_H[:,0]-zpt_atm_H[nonnan[ind],0]),yerr=zpterror_dem_H,fmt='b.',label='Helheim',markersize=3,capsize=1)
plt.ylabel('Elevation \n (m)',fontsize=9,fontname='Arial')
plt.xticks(range(2008,2017))
plt.xlim([time1,time2])
plt.yticks(np.arange(-40,40,20),fontsize=8,fontname='Arial')
plt.ylim([-30,30])
ax.set_xticklabels([])
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1H_runoff[i],-100],[day1H_runoff[i],100],[day2H_runoff[i],100],[day2H_runoff[i],-100]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
  ax.add_patch(patch)
#plt.text(time1+0.2,10,'c',fontsize=9,fontname='Arial',fontweight='bold')

plt.subplot(gs[3,1])
ax = plt.gca()
#nonnan = np.where(~(np.isnan(time_atm_K)))[0]
#ind = np.argmin(abs(time_atm_K[nonnan]-time1))
#ind2 = np.argmin(timezsrac_K-time_atm_K[nonnan[ind]])
#plt.plot(timezsrac_H,zsrac_K-zsrac_K[ind2],c='r',lw=0.5)
#plt.errorbar(time_atm_K,(zpt_atm_K[:,0]-zpt_atm_K[nonnan[ind],0]),fmt='b.',label='Kanger',markersize=3,capsize=1)
#plt.errorbar(time_dem_K,(zpt_dem_K[:,0]-zpt_atm_K[nonnan[ind],0]),yerr=zpterror_dem_K,fmt='b.',label='Kanger',markersize=3,capsize=1)
plt.xticks(range(2008,2017))
plt.xlim([time1,time2])
plt.yticks(np.arange(-40,40,20),fontsize=8,fontname='Arial')
plt.ylim([-55,5])
ax.set_xticklabels([])
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1K_runoff[i],-100],[day1K_runoff[i],100],[day2K_runoff[i],100],[day2K_runoff[i],-100]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
  ax.add_patch(patch)
#plt.text(time1+0.2,10,'c',fontsize=9,fontname='Arial',fontweight='bold')


plt.subplot(gs[-2,0])
ax = plt.gca()
plt.plot([time1,time2],[0,0],'k')
plt.plot(timezsrac_H[1:],np.diff(zsrac_H)/np.diff(timezsrac_H),c='0.5',label='SMB',lw=0.75)
plt.errorbar(dem_time_H[:,0],dem_dH_H[:,0],xerr=dem_time_H[:,1],yerr=dem_dH_H[:,1],label='Obs',fmt='k.',markersize=3,capsize=1,lw=0.5)
plt.xticks(range(2008,2017),fontsize=8,fontname='Arial')
plt.xlim([time1,time2])
plt.ylabel('dH/dt \n (m/yr)',fontname='Arial',fontsize=9)
plt.yticks(np.arange(-80,120,40),fontsize=8,fontname='Arial')
plt.ylim([-80,80])
ax.set_xticklabels([])
#xTickPos = np.linspace(np.floor(time1)-0.25,np.ceil(time2)-0.25,(np.ceil(time2)-np.floor(time1))*2+1)
#ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.85','w'],linewidth=0)
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1H_runoff[i],-100],[day1H_runoff[i],100],[day2H_runoff[i],100],[day2H_runoff[i],-100]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
  ax.add_patch(patch)
#plt.text(time1+0.2,50,'d',fontsize=9,fontname='Arial',fontweight='bold')
plt.legend(loc=1,numpoints=1,ncol=1,handletextpad=0.02,fontsize=8,columnspacing=0.05,markerscale=2,handlelength=1,borderpad=0.2)

plt.subplot(gs[-2,1])
ax = plt.gca()
plt.plot([time1,time2],[0,0],'k')
plt.plot(timezsrac_K[1:],np.diff(zsrac_K)/np.diff(timezsrac_K),c='0.5',lw=0.75)
plt.errorbar(dem_time_K[:,0],dem_dH_K[:,0],xerr=dem_time_K[:,1],yerr=dem_dH_K[:,1],fmt='k.',markersize=3,capsize=1,lw=0.5)
plt.xticks(range(2008,2017),fontsize=8,fontname='Arial')
plt.xlim([time1,time2])
plt.yticks(np.arange(-80,120,40),fontsize=8,fontname='Arial')
plt.ylim([-80,80])
ax.set_xticklabels([])
#xTickPos = np.linspace(np.floor(time1)-0.25,np.ceil(time2)-0.25,(np.ceil(time2)-np.floor(time1))*2+1)
#ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.85','w'],linewidth=0)
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1K_runoff[i],-100],[day1K_runoff[i],100],[day2K_runoff[i],100],[day2K_runoff[i],-100]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
  ax.add_patch(patch)
#plt.text(time1+0.2,50,'d',fontsize=9,fontname='Arial',fontweight='bold')

plt.subplot(gs[0,0])
ax = plt.gca()
plt.plot(terminus_time_H,terminus_val_H/1e3,'k.',lw=1.5,label='Helheim',markersize=3)
plt.ylabel('Terminus \n (km)',fontsize=9,fontname='Arial')
plt.yticks(np.arange(-2,6,2),fontsize=8,fontname='Arial')
plt.ylim([-4,4])
plt.xticks(range(2008,2017),fontsize=8,fontname='Arial')
plt.xlim([time1,time2])
ax.set_xticklabels([])
#xTickPos = np.linspace(np.floor(time1)-0.25,np.ceil(time2)-0.25,(np.ceil(time2)-np.floor(time1))*2+1)
#ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.85','w'],linewidth=0)
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1H_runoff[i],-4],[day1H_runoff[i],4],[day2H_runoff[i],4],[day2H_runoff[i],-4]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
  ax.add_patch(patch)
#plt.text(time1+0.2,2.5,'e',fontsize=9,fontname='Arial',fontweight='bold')

plt.subplot(gs[0,1])
ax = plt.gca()
plt.plot(terminus_time_K,terminus_val_K/1e3,'k.',lw=1.5,label='Kanger',markersize=3)
plt.yticks(np.arange(-2,6,2),fontsize=8,fontname='Arial')
plt.ylim([-4,4])
plt.xticks(range(2008,2017),fontsize=8,fontname='Arial')
plt.xlim([time1,time2])
ax.set_xticklabels([])
#xTickPos = np.linspace(np.floor(time1)-0.25,np.ceil(time2)-0.25,(np.ceil(time2)-np.floor(time1))*2+1)
#ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.85','w'],linewidth=0)
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1K_runoff[i],-4],[day1K_runoff[i],4],[day2K_runoff[i],4],[day2K_runoff[i],-4]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
  ax.add_patch(patch)
#plt.text(time1+0.2,2.5,'e',fontsize=9,fontname='Arial',fontweight='bold')

plt.subplot(gs[-1,0])
ax = plt.gca()
plt.plot(timerac_H,runrac_H,'k',lw=1,label='Helheim')
labels=[]
for i in range(2000,2017):
  labels.append('Jan \n'+str(i))
plt.xticks(range(2000,2017))
ax.set_xticklabels(labels,fontsize=8,fontname='Arial')
plt.xlim([time1,time2])
plt.ylim([0,60])
plt.yticks(np.arange(0,80,20),fontsize=8,fontname='Arial')
plt.ylabel('Runoff \n (kg m$^{-2}$ d$^{-1}$)',fontsize=9,fontname='Arial')
#xTickPos = np.linspace(np.floor(time1)-0.25,np.ceil(time2)-0.25,(np.ceil(time2)-np.floor(time1))*2+1)
#ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.85','w'],linewidth=0)
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1H_runoff[i],0],[day1H_runoff[i],60],[day2H_runoff[i],60],[day2H_runoff[i],0]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
  ax.add_patch(patch)
#plt.text(time1+0.2,45,'f',fontsize=9,fontname='Arial',fontweight='bold')

plt.subplot(gs[-1,1])
ax = plt.gca()
plt.plot(timerac_K,runrac_K,'k',lw=1,label='Kanger')
labels=[]
for i in range(2000,2017):
  labels.append('Jan \n'+str(i))
plt.xticks(range(2000,2017))
ax.set_xticklabels(labels,fontsize=8,fontname='Arial')
plt.xlim([time1,time2])
plt.ylim([0,60])
plt.yticks(np.arange(0,80,20),fontsize=8,fontname='Arial')
#xTickPos = np.linspace(np.floor(time1)-0.25,np.ceil(time2)-0.25,(np.ceil(time2)-np.floor(time1))*2+1)
#ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.85','w'],linewidth=0)
for i in range(0,len(yH_runoff)):
  path = matplotlib.path.Path([[day1K_runoff[i],0],[day1K_runoff[i],60],[day2K_runoff[i],60],[day2K_runoff[i],0]])
  patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
  ax.add_patch(patch)
#plt.text(time1+0.2,45,'f',fontsize=9,fontname='Arial',fontweight='bold')


plt.tight_layout()
plt.subplots_adjust(hspace=0.1,wspace=0.2) 
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/both_glaciers2.pdf"),FORMAT='PDF')
plt.close()

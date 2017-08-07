import os
import sys
import numpy as np
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"Modules/demtools"))
import glaclib, zslib, floatlib, datelib, demshadelib, bedlib, icefrontlib, vellib, fluxlib, climlib
import matplotlib.pyplot as plt
import matplotlib, scipy
from matplotlib.ticker import AutoMinorLocator
import scipy.signal as signal
import gmtColormap, subprocess
import cubehelix
cpt_rainbow = gmtColormap.get_rainbow()
plt.register_cmap(cmap=cpt_rainbow)


# Get arguments
args = sys.argv
glacier = args[1][:] # Options: Kanger, Helheim

# What bed to use
if glacier == 'Helheim':
  bedsource = 'cresis'
elif glacier == 'Kanger':
  bedsource = 'cresis'

if glacier == 'Helheim':
  xmin = 304700.0
  xmax = 313600.0
  ymin = -2582000.0
  ymax = -2573000.0
elif glacier == 'Kanger':
  xmin = 483700.0
  xmax = 497200.0
  ymin = -2297500.0
  ymax = -2283500.0

dists_eul = [-2.0e3,-5.0e3,-10.0e3,-15.0e3,-20.e3]

time1 = 2008.
time2 = 2016.

############ 
# Flowline #
############

x,y,zb,dists = glaclib.load_flowline(glacier,shapefilename='flowline_flightline',filt_len=2.0e3,verticaldatum='geoid',bedsource='cresis')

############################################
# Get ice front position and calving style #
############################################

terminus_val, terminus_time = icefrontlib.distance_along_flowline(x,y,dists,glacier,type='icefront',time1=time1,time2=time2)
calvingstyle = icefrontlib.calving(glacier)

############################################
# Get velocity and surface elevation at pt #
############################################

# Find index
ind_eul=[]
ind_eul=[]
for i in range(0,len(dists_eul)):
  ind_eul.append( (abs(dists - dists_eul[i])).argmin() )


# Get velocities
vel_val,vel_time,vel_error = vellib.velocity_at_eulpoints(x[ind_eul],y[ind_eul],glacier)

# Get surface elevations
zpt_atm,zptstd_atm,time_atm = zslib.atm_at_pts(x[ind_eul],y[ind_eul],glacier,years='all',maxdist=250.,verticaldatum='geoid',method='average',cutoff='terminus')
zpt_dem,zpterror_dem,time_dem = zslib.dem_at_pts(x[ind_eul],y[ind_eul],glacier,years='all',verticaldatum='geoid',cutoff='terminus',method='average',radius=250.)
zpt_wv,zpterror_wv,time_wv = zslib.dem_at_pts(x[ind_eul],y[ind_eul],glacier,years='all',verticaldatum='geoid',cutoff='terminus',method='average',radius=250.,data='WV')
zpt_tdm,zpterror_tdm,time_tdm = zslib.dem_at_pts(x[ind_eul],y[ind_eul],glacier,years='all',verticaldatum='geoid',cutoff='terminus',method='average',radius=250.,data='TDM')
zpt_spirit,zpterror_spirit,time_spirit = zslib.dem_at_pts(x[ind_eul],y[ind_eul],glacier,years='all',verticaldatum='geoid',cutoff='terminus',method='average',radius=250.,data='SPIRIT')

######################
# Get thinning rates #
######################

xzdem,yzdem,zsdem,zstimedem,errordem = zslib.dem_grid(glacier,xmin,xmax,ymin,ymax,
	years='all',verticaldatum='geoid',return_error=True)
dem_time,dem_dH = fluxlib.dem_thinning(glacier,xzdem,yzdem,zsdem,zstimedem,errordem,"fluxgate1",type='rate')
flux_time,flux_dH,Q,hbar,ubar,smb_flux,width_flux,area_flux = fluxlib.fluxgate_thinning(glacier,"fluxgate1",bedsource=bedsource,timing='velocity')
xflux,yflux = fluxlib.fluxbox_geometry(glacier,"fluxgate1")

xrac,yrac,smbrac,timerac = climlib.racmo_at_pts(np.mean(xflux),np.mean(yflux),'smb',filt_len=14.0)
xrac,yrac,runrac,timerac = climlib.racmo_at_pts(np.mean(xflux),np.mean(yflux),'runoff',filt_len=14.0)
#xrac,yrac,zsrac,timeraczs = climlib.racmo_at_pts(np.mean(xflux),np.mean(yflux),'zs',filt_len=14.0)
#xsif,ysif,sif,timesif = climlib.SIF_at_pts(np.mean(xflux),np.mean(yflux),filt_len=14.)

year_runoff,day1_runoff,day2_runoff,meltlength_runoff,total_runoff = climlib.seasonlength(timerac,runrac,'runoff')

#####################
# Get velocity DEMs #
#####################

xvdem,yvdem,vdem,vtrend,vdetrend,vrange,vcount,vtimedem = vellib.variability(glacier,time1,time2)

##############
# Get images #
##############

images,imagetimes,imagesats = glaclib.load_satimages(glacier,xmin,xmax,ymin,ymax,time1,time2)

####################################
# Get bed elevations near flowline #
####################################

# Get radar thicknesses close to flightline
cresis = bedlib.cresis('all',glacier)
if glacier == 'Helheim':
	cresis2001 = bedlib.cresis('2001',glacier)
	cresis = np.row_stack([cresis,cresis2001])

cutoff = 200.
dcresis = []
zcresis = []
tcresis = []
for i in range(0,len(cresis[:,0])):
	mindist = np.min(np.sqrt((cresis[i,0]-x)**2+(cresis[i,1]-(y))**2))
	if mindist < cutoff:
		minind = np.argmin(np.sqrt((cresis[i,0]-x)**2+(cresis[i,1]-(y))**2))
		dcresis.append(dists[minind])
		zcresis.append(cresis[i,2])
		tcresis.append(cresis[i,4])
dcresis = np.array(dcresis)
zcresis = np.array(zcresis)

#############################
# Make plots for each image #
#############################

if glacier == 'Helheim':
  ptind=0
elif glacier == 'Kanger':
  ptind=1

cx = cubehelix.cmap(start=1.2,rot=-1.1,reverse=True,minLight=0.1,sat=2,maxLight=0.9)

n=0  
nonnan = np.where(~(np.isnan(zpt_dem[:,ptind])))[0]
for time in np.unique(time_dem[nonnan]):

  zsind = np.where(zstimedem == time)[0]

  fig = plt.figure(figsize=(11,9))
  gs = matplotlib.gridspec.GridSpec(12,3)
  matplotlib.rc('font',family='Arial')
  
  plt.subplot(gs[0:4,0])
  imageind = np.argmin(abs(imagetimes-zstimedem[zsind[0]]))
  ax=plt.gca()
  plt.imshow(images[imageind][2],extent=[np.min(images[imageind][0]),np.max(images[imageind][0]),np.min(images[imageind][1]),np.max(images[imageind][1])],origin='lower',cmap='Greys_r')
  plt.xticks([])
  plt.yticks([])
  plt.xlim([xmin,xmax])
  plt.ylim([ymin,ymax])
  year,month,day = datelib.fracyear_to_date(imagetimes[imageind])
  path = matplotlib.path.Path([[0.02*(xmax-xmin)+xmin,0.82*(ymax-ymin)+ymin],
  			[0.02*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
  			[0.38*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
  			[0.38*(xmax-xmin)+xmin,0.82*(ymax-ymin)+ymin],
  			[0.02*(xmax-xmin)+xmin,0.82*(ymax-ymin)+ymin]])
  patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1,zorder=3)
  ax.add_patch(patch)
  plt.text(xmin+0.04*(xmax-xmin),0.92*(ymax-ymin)+ymin,'{0:4d}-{1:02d}-{2:02d}'.format(year,month,int(np.round(day))),fontsize=12,fontname='Arial')
  plt.text(xmin+0.04*(xmax-xmin),0.845*(ymax-ymin)+ymin,imagesats[imageind],fontsize=12,fontname='Arial')
  plt.text(xmin+0.04*(xmax-xmin),0.03*(ymax-ymin)+ymin,'A',fontweight='bold',color='w',fontsize=12,fontname='Arial')


  plt.subplot(gs[4:8,0])
  ax = plt.gca()
  velind = np.argmin(abs(vtimedem-zstimedem[zsind[0]]))
  p=plt.imshow(vdem[:,:,velind]/1e3,extent=[xvdem[0],xvdem[-1],yvdem[0],yvdem[-1]],origin='lower',clim=[4,10],cmap=cx)
  plt.plot(x,y,'k',lw=2)
  if glacier == 'Helheim':
    plt.plot(x[ind_eul[0]],y[ind_eul[0]],'o',color=[0.1,0.1,200./255],mec='k',mew=0.5,markersize=10)
    plt.text(x[ind_eul[0]]-900,y[ind_eul[0]]-800,'H02',fontsize=12,fontname='Arial')
  elif glacier == 'Kanger':
    plt.plot(x[ind_eul[1]],y[ind_eul[1]],'o',color=[182./255,219./255,1],mec='k',mew=0.5,markersize=10)
    plt.text(x[ind_eul[1]]-1000,y[ind_eul[1]]-1200,'K05',fontsize=12,fontname='Arial')
  plt.xticks([])
  plt.yticks([])
  plt.xlim([xmin,xmax])
  plt.ylim([ymin,ymax])
  year,month,day = datelib.fracyear_to_date(vtimedem[velind])
  path = matplotlib.path.Path([[0.02*(xmax-xmin)+xmin,0.9*(ymax-ymin)+ymin],
  			[0.02*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
  			[0.38*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
  			[0.38*(xmax-xmin)+xmin,0.9*(ymax-ymin)+ymin],
  			[0.02*(xmax-xmin)+xmin,0.9*(ymax-ymin)+ymin]])
  patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1,zorder=3)
  ax.add_patch(patch)
  plt.text(xmin+0.04*(xmax-xmin),0.92*(ymax-ymin)+ymin,'{0:4d}-{1:02d}-{2:02d}'.format(year,month,int(np.round(day))),fontsize=12,fontname='Arial')
  plt.text(xmin+0.04*(xmax-xmin),0.03*(ymax-ymin)+ymin,'B',fontweight='bold',color='k',fontsize=12,fontname='Arial')
  path = matplotlib.path.Path([[0.46*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
  			[0.98*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
  			[0.98*(xmax-xmin)+xmin,0.73*(ymax-ymin)+ymin],
  			[0.46*(xmax-xmin)+xmin,0.73*(ymax-ymin)+ymin],
  			[0.46*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin]])
  patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1,zorder=3)
  ax.add_patch(patch)
  cbaxes = fig.add_axes([0.157, 0.645, 0.085, 0.01]) 
  cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=[4,6,8,10]) 
  cb.set_label('Velocity'+r' (km yr$^{-1}$)',fontsize=12,fontname='Arial')

  
  plt.subplot(gs[8:12,0])
  ax = plt.gca()
  #zsind[0] = np.argmin(abs(zstimedem-imagetimes[imageind]))
  if len(zsind) > 1:
    plt.imshow(demshadelib.hillshade(np.nanmean(zsdem[:,:,zsind],axis=2)),extent=[xmin,xmax,ymin,ymax],origin='lower',cmap='Greys_r')
    p = plt.imshow(np.nanmean(zsdem[:,:,zsind],axis=2),extent=[xmin,xmax,ymin,ymax],origin='lower',clim=[0,400],cmap='cpt_rainbow',alpha=0.6)  
  else:
    plt.imshow(demshadelib.hillshade(zsdem[:,:,zsind[0]]),extent=[xmin,xmax,ymin,ymax],origin='lower',cmap='Greys_r')
    p = plt.imshow(zsdem[:,:,zsind[0]],extent=[xmin,xmax,ymin,ymax],origin='lower',clim=[0,400],cmap='cpt_rainbow',alpha=0.6)
  plt.plot(x,y,'k',lw=2)
  if glacier == 'Helheim':
    plt.plot(x[ind_eul[0]],y[ind_eul[ptind]],'o',color=[0.1,0.1,200./255],mec='k',mew=0.5,markersize=10)
    plt.text(x[ind_eul[0]]-900,y[ind_eul[0]]-800,'H02',fontsize=12,fontname='Arial')
  elif glacier == 'Kanger':
    plt.plot(x[ind_eul[1]],y[ind_eul[1]],'o',color=[182./255,219./255,1],mec='k',mew=0.5,markersize=10)
    plt.text(x[ind_eul[1]]-1000,y[ind_eul[1]]-1200,'K05',fontsize=12,fontname='Arial')
  plt.xlim([xmin,xmax])
  plt.ylim([ymin,ymax])
  plt.xticks([])
  plt.yticks([])
  path = matplotlib.path.Path([[0.02*(xmax-xmin)+xmin,0.9*(ymax-ymin)+ymin],
  			[0.02*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
  			[0.38*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
  			[0.38*(xmax-xmin)+xmin,0.9*(ymax-ymin)+ymin],
  			[0.02*(xmax-xmin)+xmin,0.9*(ymax-ymin)+ymin]])
  patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1,zorder=3)
  ax.add_patch(patch)
  year,month,day = datelib.fracyear_to_date(zstimedem[zsind[0]])
  plt.text(xmin+0.04*(xmax-xmin),0.92*(ymax-ymin)+ymin,'{0:4d}-{1:02d}-{2:02d}'.format(year,month,int(np.round(day))),fontsize=12,fontname='Arial')
  plt.text(xmin+0.04*(xmax-xmin),0.03*(ymax-ymin)+ymin,'C',fontweight='bold',color='w',fontsize=12,fontname='Arial')
  path = matplotlib.path.Path([[0.49*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
  			[0.98*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
  			[0.98*(xmax-xmin)+xmin,0.73*(ymax-ymin)+ymin],
  			[0.49*(xmax-xmin)+xmin,0.73*(ymax-ymin)+ymin],
  			[0.49*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin]])
  patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1,zorder=3)
  ax.add_patch(patch)
  cbaxes = fig.add_axes([0.16, 0.33, 0.085, 0.01]) 
  cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=[0,200,400]) 
  cb.set_label(r'Elevation'+' (m asl)',fontsize=12,fontname='Arial')
  
  plt.subplot(gs[0:2,1:])
  ax = plt.gca()
  ax.tick_params(axis='x',which='both',direction='in')
  ax.xaxis.set_ticks_position('both')
  for i in range(0,len(year_runoff)):
    path = matplotlib.path.Path([[day1_runoff[i],-5],[day1_runoff[i],5],[day2_runoff[i],5],[day2_runoff[i],-5]])
    patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
    ax.add_patch(patch)
  ind = np.where(calvingstyle[:,1] == 'Tabular')[0]
  for i in ind:
    if i == ind[1]:
      plt.bar(float(calvingstyle[i,0])-0.025,10.0,width=0.03,color=[0.6,0.6,1],edgecolor='none',label='Tabular',bottom=-5)
    elif i !=0:
      plt.bar(float(calvingstyle[i,0])-0.025,10.0,width=0.03,color=[0.6,0.6,1],edgecolor='none',bottom=-5)
  ind = np.where((calvingstyle[:,1] == 'Mixed'))[0]
  for i in ind:
    if i == ind[1]:
      plt.bar(float(calvingstyle[i,0])-0.025,10.0,width=0.03,color='gold',edgecolor='none',label='Mixed',bottom=-5)
    elif i != 0:
      plt.bar(float(calvingstyle[i,0])-0.025,10.0,width=0.03,color='gold',edgecolor='none',bottom=-5)
  ind = np.where((calvingstyle[:,1] == 'Domino'))[0]
  for i in ind:
    if i == ind[1]:
      plt.bar(float(calvingstyle[i,0])-0.025,10.0,width=0.03,color=[1,0.6,0.6],edgecolor='none',label='Non-tabular',bottom=-5)
    elif i != 0:
      plt.bar(float(calvingstyle[i,0])-0.025,10.0,width=0.03,color=[1,0.6,0.6],edgecolor='none',bottom=-5)
  plt.plot(terminus_time,terminus_val/1e3,'k.')
  plt.xticks(range(2000,2017))
  plt.xlim([time1,time2])
  ax.set_xticklabels([])
  if glacier == 'Helheim':
    plt.ylim([-2.5,2.5])
    plt.text(2008.2,-2.1,'D',fontweight='bold',fontsize=12,fontname='Arial')
  elif glacier == 'Kanger':
    plt.ylim([-4.5,4.5])
    plt.yticks(np.arange(-4,6,2))
    plt.text(2008.2,-4,'D',fontweight='bold',fontsize=12,fontname='Arial')
  plt.ylabel('Terminus \n (km)',fontsize=12,fontname='Arial')
  plt.legend(loc=2,fontsize=10,numpoints=1,handlelength=0.3,labelspacing=0.05,ncol=3,columnspacing=0.7,handletextpad=0.2,borderpad=0.25)
  plt.plot([zstimedem[zsind[0]],zstimedem[zsind[0]]],[-5.5,5.5],'k--',lw=2)
  
  
  plt.subplot(gs[2:4,1:])
  ax = plt.gca()
  ax.tick_params(axis='x',which='both',direction='in')
  ax.xaxis.set_ticks_position('both')
  for i in range(0,len(year_runoff)):
    path = matplotlib.path.Path([[day1_runoff[i],5],[day1_runoff[i],12],[day2_runoff[i],12],[day2_runoff[i],5]])
    patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
    ax.add_patch(patch)
  if glacier == 'Helheim':
    plt.plot(vel_time,vel_val[:,ptind]/1e3,'o',color=[0.1,0.1,200./255],mec='k',mew=0.5,markersize=5)
  elif glacier == 'Kanger':
    plt.plot(vel_time,vel_val[:,ptind]/1e3,'o',color=[182./255,219./255,1],mec='k',mew=0.5,markersize=5)
  plt.xticks(range(2000,2017))
  plt.xlim([time1,time2])
  ax.set_xticklabels([])
  if glacier == 'Helheim':
    plt.yticks(np.arange(6,12,2))
    plt.text(2008.2,6.4,'E',fontweight='bold',fontsize=12,fontname='Arial')
    plt.ylim([6,11])
  if glacier == 'Kanger':
    plt.yticks(np.arange(7,12,1))
    plt.text(2008.2,7.5,'E',fontweight='bold',fontsize=12,fontname='Arial')
    plt.ylim([7,11])
  plt.ylabel('Glacier velocity \n'+r'(km yr$^{-1}$)',fontsize=12,fontname='Arial')
  plt.plot([zstimedem[zsind[0]],zstimedem[zsind[0]]],[0,12],'k--',lw=2)
  
  plt.subplot(gs[4:6,1:])
  if glacier == 'Helheim':
    yerr = 3.
  elif glacier == 'Kanger':
    yerr = 5.
  ax = plt.gca()
  ax.tick_params(axis='x',which='both',direction='in')
  ax.xaxis.set_ticks_position('both')
  plt.plot([time1,time2],[0,0],'k')
  plt.plot(timerac,smbrac/900.*100,'0.4',lw=1.5,label='SMB')
  #plt.plot(timeraczs[1:],np.convolve(np.diff(zsrac)/np.diff(timeraczs),np.ones((11,))/11.,mode='same'),c='0.4',lw=1.2,label='SMB')
  plt.errorbar(dem_time[:,0],dem_dH[:,0]/365.25*100,xerr=dem_time[:,1],yerr=dem_dH[:,1]/365.25*100,fmt='o',color=[182./255,219./255,1],ecolor='k',mew=0.5,mec='k',markersize=5,capsize=1,lw=1,label='DEM')
  plt.errorbar(flux_time,(flux_dH[:,0]+smb_flux)/365.25*100,yerr=yerr,xerr=5.5/365.25,fmt='s',color=[200./255,0,0],markersize=5,mew=0.5,mec='k',capsize=1,lw=1,ecolor='k',label='Flux-gate')
  #time_H_2,dH_H_2,Q_H_2,hbar_H_2,ubar_H_2 = fluxlib.fluxgate_thinning('Helheim','fluxgate1','cresis',10.,timing='velocity')
  #plt.plot(timezsrac_K[1:],np.diff(zsrac_K)/np.diff(timezsrac_K),c='r',lw=0.5)
  plt.xticks(range(2008,2017),fontsize=8,fontname='Arial')
  plt.xlim([time1,time2])
  plt.ylabel('$dh/dt$ \n'+r'(cm d$^{-1}$)',fontsize=12,fontname='Arial')
  plt.yticks(np.arange(-40,40,10),fontname='Arial')
  plt.yticks(np.arange(-40,40,20))
  if glacier == 'Kanger':
    plt.ylim([-37,29])
    plt.text(2008.2,-32,'F',fontweight='bold',fontsize=12,fontname='Arial')
  elif glacier == 'Helheim':
    plt.ylim([-42,29])
    plt.text(2008.2,-35,'F',fontweight='bold',fontsize=12,fontname='Arial')
  ax.set_xticklabels([])
  #xTickPos = np.linspace(np.floor(time1)-0.25,np.ceil(time2)-0.25,(np.ceil(time2)-np.floor(time1))*2+1)
  #ax.bar(xTickPos, [max(plt.ylim())-min(plt.ylim())] * len(xTickPos), (xTickPos[1]-xTickPos[0]), bottom=min(plt.ylim()), color=['0.85','w'],linewidth=0)
  for i in range(0,len(year_runoff)):
    path = matplotlib.path.Path([[day1_runoff[i],-180],[day1_runoff[i],180],[day2_runoff[i],180],[day2_runoff[i],-180]])
    patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
    ax.add_patch(patch)
  plt.legend(loc=2,numpoints=1,ncol=3,handletextpad=0.2,fontsize=10,columnspacing=0.05,markerscale=1,handlelength=1.2,borderpad=0.2)
  plt.plot([zstimedem[zsind[0]],zstimedem[zsind[0]]],[-110,110],'k--',lw=2)

  plt.subplot(gs[6:8,1:])
  ax = plt.gca()
  ax.tick_params(axis='x',which='both',direction='in')
  ax.xaxis.set_ticks_position('both')
  for i in range(0,len(year_runoff)):
    path = matplotlib.path.Path([[day1_runoff[i],-120],[day1_runoff[i],100],[day2_runoff[i],100],[day2_runoff[i],-120]])
    patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
    ax.add_patch(patch)
  if glacier == 'Helheim':
    plt.fill_between([time1,time2],np.ones(2)*floatlib.height(zb[ind_eul[ptind]]-50)-floatlib.height(zb[ind_eul[ptind]]),np.ones(2)*floatlib.height(zb[ind_eul[ptind]]+50)-floatlib.height(zb[ind_eul[ptind]]),\
          alpha=0.15,facecolor=[0.1,0.1,200./255],edgecolor=[0.1,0.1,200./255],antialiased=True,zorder=2)
    plt.plot([time1,time2],[0,0],':',color=[0.1,0.1,200./255],mec='k',mew=0.5,linewidth=1.5)
    plt.errorbar(time_wv,zpt_wv[:,ptind]-floatlib.height(zb[ind_eul[ptind]]),capsize=1,yerr=zpterror_wv,fmt='o',color=[0.1,0.1,200./255],markersize=5,mec='k',mew=0.5,label='WV')
    plt.errorbar(time_tdm,zpt_tdm[:,0]-floatlib.height(zb[ind_eul[0]]),capsize=1,yerr=zpterror_tdm,fmt='^',color=[0.1,0.1,200./255],markersize=5,mec='k',mew=0.5,label='TDM')
    plt.plot(time_atm,zpt_atm[:,0]-floatlib.height(zb[ind_eul[0]]),'d',color=[0.1,0.1,200./255],markersize=5,mec='k',mew=0.5,label='ATM')
    plt.yticks(np.arange(-10,40,10))
    plt.ylim([-15,25])
    plt.text(2008.2,-5,'G',fontweight='bold',fontsize=12,fontname='Arial')
  elif glacier == 'Kanger':
    plt.fill_between([time1,time2],np.ones(2)*floatlib.height(zb[ind_eul[1]]-50)-floatlib.height(zb[ind_eul[1]]),np.ones(2)*floatlib.height(zb[ind_eul[1]]+50)-floatlib.height(zb[ind_eul[1]]),\
    	alpha=0.25,facecolor=[182./255,219./255,1],edgecolor=[182./255,219./255,1],antialiased=True,zorder=2)
    plt.plot([time1,time2],[0,0],':',color=[182./255,219./255,1],mec='k',mew=0.5,linewidth=1.5)
    plt.errorbar(time_wv,zpt_wv[:,1]-floatlib.height(zb[ind_eul[1]]),capsize=1,yerr=zpterror_wv,fmt='o',color=[182./255,219./255,1],mec='k',mew=0.5,markersize=5,label='WV')
    plt.errorbar(time_tdm,zpt_tdm[:,1]-floatlib.height(zb[ind_eul[1]]),capsize=1,yerr=zpterror_tdm,fmt='^',color=[182./255,219./255,1],mec='k',mew=0.5,markersize=5,label='TDM')
    plt.errorbar(time_spirit,zpt_spirit[:,1]-floatlib.height(zb[ind_eul[1]]),capsize=1,yerr=zpterror_spirit,fmt='v',color=[182./255,219./255,1],mec='k',mew=0.5,markersize=4,label='SPIRIT')
    plt.plot(time_atm,zpt_atm[:,1]-floatlib.height(zb[ind_eul[1]]),'d',color=[182./255,219./255,1],markersize=5,mec='k',mew=0.5,label='ATM')
    plt.ylim([-25,40])
    plt.text(2008.2,-7,'G',fontweight='bold',fontsize=12,fontname='Arial')
  plt.legend(loc=3,borderpad=0.3,handleheight=0.1,fontsize=10,numpoints=1,handlelength=0.4,labelspacing=0.05,ncol=4,columnspacing=0.7,handletextpad=0.5)
  plt.xlim([time1,time2])
  plt.xticks(np.arange(2008,2017,1),fontsize=12,fontname='Arial')
  plt.ylabel('Height above \n flotation (m)',fontsize=12,fontname='Arial')
  x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
  ax.xaxis.set_major_formatter(x_formatter)
  plt.plot([zstimedem[zsind[0]],zstimedem[zsind[0]]],[-110,110],'k--',lw=2)

  plt.subplot(gs[9:,1:])
  ind = np.where(dists < -4.8e3)[0]
  plt.plot(dcresis/1e3,zcresis,'.',color='0.7',label='Measured bed')
  if glacier == 'Helheim':
    plt.plot(dists/1e3,zb,'k',lw=1,label='Smoothed bed')
  elif glacier == 'Kanger':
    ind = np.where(dists < -4.8e3)[0]
    plt.plot(dists[ind]/1e3,zb[ind],color='k',linewidth=1.5,label='Smoothed bed')
  #plt.plot([-10,5],[0,0],'0.8',lw=0.5)
  #if i != 0:
  #   plt.plot(dists/1e3,f((y,x)),c='0.7')
  #   plt.plot(dists/1e3,floatlib.shelfbase(f(((y,x)))),c='0.7')
  if len(zsind) > 1:
    f = scipy.interpolate.RegularGridInterpolator((yzdem,xzdem),np.nanmean(zsdem[:,:,zsind],axis=2),bounds_error = False,method='linear',fill_value=float('nan'))
  else:
    f = scipy.interpolate.RegularGridInterpolator((yzdem,xzdem),zsdem[:,:,zsind[0]],bounds_error = False,method='linear',fill_value=float('nan'))
  plt.plot(dists/1e3,floatlib.height(zb),'k:',lw=2,label='Flotation threshold')
  plt.plot(dists/1e3,f((y,x)),color=[0.1,0.1,200./255],lw=2,label='Ice surface')
  plt.plot(dists/1e3,floatlib.shelfbase(f(((y,x)))),'--',color=[0.1,0.1,200./255],lw=2,dashes=[2,1],label='Ice bottom')
  plt.legend(loc=4,borderpad=0.3,handleheight=0.1,fontsize=10,framealpha=1,numpoints=1,handlelength=1.0,labelspacing=0.05,columnspacing=0.7,handletextpad=0.5)
  if glacier == 'Helheim':
    plt.xticks(range(-10,5,2))
    plt.yticks(np.arange(-1200,400,200))
    plt.ylim([-950,220])
    plt.xlim([-5,4])
    plt.text(-4.8,-50,'H',fontweight='bold',fontsize=12,fontname='Arial')
  elif glacier == 'Kanger':
    plt.ylim([-1200,220])
    plt.xticks(np.arange(-10,5,2))
    plt.yticks(np.arange(-1200,400,200))
    plt.xlim([-7,4])
    plt.text(-6.7,-50,'H',fontweight='bold',fontsize=12,fontname='Arial')
  #plt.title(str(year)+'-'+str(month)+'-'+str(int(day)))
  plt.ylabel('Elevation (m asl)',fontsize=12,fontname='Arial')
  plt.xlabel('Distance along flowline (km)',fontsize=12,fontname='Arial')


  plt.tight_layout()
  plt.subplots_adjust(hspace=0.1,wspace=0.4,top=0.98,right=0.98,left=0.02,bottom=0.05) 
  
  plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/movie/"+'{0:03g}'.format(n)+'.png'),format='PNG',dpi=200)
  plt.close()
  n=n+1

##############
# Make movie #
##############

fps=1
ffmpeg_in_opt = "-r %i" % fps
#ffmpeg_out_opt = "-r 5 -y -an -force_fps -s '%ix%i'" % newsize
#-an disables audio
ffmpeg_out_opt = "-y -an -c:v libx264 -pix_fmt yuv420p"
#scale='iw/4:-1'"

os.chdir(os.path.join(os.getenv("HOME"),"Bigtmp/movie/"))
outmov = glacier+'.mp4'
#cmd = 'ffmpeg {0} -i %04d.jpg {1} {2}'.format(ffmpeg_in_opt, ffmpeg_out_opt, outmov)
#cmd = 'ffmpeg {0} -pattern_type glob -i *_clip.png {1} {2}'.format(ffmpeg_in_opt, ffmpeg_out_opt, outmov)
cmd = 'ffmpeg {0} -i %03d.png {1} {2}'.format(ffmpeg_in_opt, ffmpeg_out_opt, outmov)

print cmd
subprocess.call(cmd, shell=True)

#cmd = 'rm *.png'
#print cmd
#subprocess.call(cmd, shell=True)


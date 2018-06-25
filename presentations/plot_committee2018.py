# Make some Helheim/Kanger plots for committee meeting, winter 2018

import datelib, meshlib, geotifflib, glaclib, zslib, fluxlib, vellib, masklib, icefrontlib, climlib, floatlib, zslib, bedlib
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib
import cubehelix
from mpl_toolkits.basemap import Basemap
from matplotlib.ticker import AutoMinorLocator
from matplotlib.path import Path
from matplotlib.patches import PathPatch

#############
# Load data #
#############

# Flowlines
x_H,y_H,zb_H,dists_H = glaclib.load_flowline('Helheim',shapefilename='flowline_flightline',filt_len=2.0e3,bedsource='cresis')
x_K,y_K,zb_K,dists_K = glaclib.load_flowline('Kanger',shapefilename='flowline_flightline',filt_len=2.0e3,bedsource='cresis')

# Get inversion dates
times_inv_K = []; times_inv_H = []
DIR = os.path.join(os.getenv("MODEL_HOME"),"Helheim/Results/INV_SSA_ModelT")
files = os.listdir(DIR)
for file in files:
    if file.endswith('taub.tif'):
        times_inv_H.append(datelib.date_to_fracyear(int(file[3:7]),int(file[7:9]),int(file[9:11])))
DIR = os.path.join(os.getenv("MODEL_HOME"),"Kanger/Results/INV_SSA_ModelT")
files = os.listdir(DIR)
for file in files:
    if file.endswith('taub.tif'):
        times_inv_K.append(datelib.date_to_fracyear(int(file[3:7]),int(file[7:9]),int(file[9:11])))

# Ice-front positions
terminus_val_H, terminus_time_H = icefrontlib.distance_along_flowline(x_H,y_H,dists_H,'Helheim',type='icefront',time1=2000.,time2=2016.5)
terminus_val_K, terminus_time_K = icefrontlib.distance_along_flowline(x_K,y_K,dists_K,'Kanger',type='icefront',time1=2000.,time2=2016.5)

# Locations where we want velocities and surface elevations
dists_eul = -1*np.array([2.,5.,10.,15.,20.])
ind_eul_H=[]
for i in range(0,len(dists_eul)):
    ind_eul_H.append( (abs(dists_H - dists_eul[i]*1e3)).argmin() )
ind_eul_K=[]
for i in range(0,len(dists_eul)):
    ind_eul_K.append( (abs(dists_K - dists_eul[i]*1e3)).argmin() )

# Load velocities
vel_val_H,vel_time_H,vel_error_H = vellib.velocity_at_eulpoints(x_H[ind_eul_H],y_H[ind_eul_H],'Helheim',data='all')
vel_val_K,vel_time_K,vel_error_K = vellib.velocity_at_eulpoints(x_K[ind_eul_K],y_K[ind_eul_K],'Kanger',data='all')
vel_val_H_op,vel_time_H_op,vel_error_H_op = vellib.howat_optical_at_pts(x_H[ind_eul_H],y_H[ind_eul_H],'Helheim')
vel_val_K_op,vel_time_K_op,vel_error_K_op = vellib.howat_optical_at_pts(x_K[ind_eul_K],y_K[ind_eul_K],'Kanger')

# Load elevations
zpt_dem_H,zpterror_dem_H,time_dem_H = zslib.dem_at_pts(x_H[ind_eul_H],y_H[ind_eul_H],'Helheim',years='all',verticaldatum='geoid',cutoff='none',method='average',radius=200.)
zpt_dem_K,zpterror_dem_K,time_dem_K = zslib.dem_at_pts(x_K[ind_eul_K],y_K[ind_eul_K],'Kanger',years='all',verticaldatum='geoid',cutoff='terminus',method='average',radius=200.)

# Color options for plotting
coloptions=['r','b','g','limegreen','gold']



###############
# Main figure #
###############
plot_inv = 0
for glacier in ['Kanger','Helheim']:

    vbars = 'none'
    plot_calving = 0

    if glacier == 'Helheim':
        x = x_H; y = y_H; zb = zb_H; dists = dists_H
        vel_val = vel_val_H; vel_time = vel_time_H
        vel_val_op = vel_val_H_op; vel_time_op = vel_time_H_op
        terminus_val = terminus_val_H; terminus_time = terminus_time_H
        time_dem = time_dem_H; zpt_dem = zpt_dem_H; zpterror_dem = zpterror_dem_H
        ind_eul = ind_eul_H
        times_inv = times_inv_H
    elif glacier == 'Kanger':  
        x = x_K; y = y_K; zb = zb_K; dists = dists_K
        vel_val = vel_val_K; vel_time = vel_time_K
        vel_val_op = vel_val_K_op; vel_time_op = vel_time_K_op
        terminus_val = terminus_val_K; terminus_time = terminus_time_K
        time_dem = time_dem_K; zpt_dem = zpt_dem_K; zpterror_dem = zpterror_dem_K
        ind_eul = ind_eul_K
        times_inv = times_inv_K
    
    plt.figure(figsize=(8,6))
    time1 = 2000.0; time2 = 2016.75; 
    years = np.arange(np.floor(time1),np.ceil(time2)+1)
 
    gs = matplotlib.gridspec.GridSpec(6,1)
    matplotlib.rc('font',family='Arial')

    # Plot terminus
    plt.subplot(gs[0, :])
    ax = plt.gca()

    #for i in range(2000,2017):
        #plt.plot([i,i],[-5,120],color='0.5',lw=0.5)    
    if vbars == 'runoff':
        for i in range(0,len(year_runoff)):
            path = matplotlib.path.Path([[day1_runoff[i],-5],[day1_runoff[i],5],[day2_runoff[i],5],\
                    [day2_runoff[i],-5]])
            patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
            ax.add_patch(patch)
    elif vbars == 'seasonal':
        for i in range(0,len(years)):
            path = matplotlib.path.Path([[years[i]+0.25,-6],[years[i]+0.25,12],[years[i]+0.75,12],[years[i]+0.75,-6]])
            patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
            ax.add_patch(patch)  

    nonnan = np.where(~(np.isnan(terminus_val)))[0]
    plt.plot(terminus_time[nonnan],terminus_val[nonnan]/1e3,'ko',linewidth=1,markersize=3.5)
    x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax.xaxis.set_major_formatter(x_formatter)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    plt.xticks(range(2000,2017),fontsize=16,fontname="Arial")
    #ax.xaxis.tick_top()
    labels=[]
    for i in range(2000,2017):
        labels.append('Jan \n'+str(i))
    plt.xticks(np.arange(2000,2017,2))
    ax.set_xticklabels([])
    #ax.set_xticklabels(labels,fontsize=11,fontname='Arial')
    plt.xlim([time1,time2])
    plt.yticks(np.arange(-6,8,3),fontsize=16,fontname="Arial")
    plt.ylabel('Terminus \n (km)',fontsize=16,fontname="Arial")
    ax.tick_params('both', length=6, width=1.25, which='major')
    ax.tick_params('x', length=6, width=1.25, which='minor')
    if plot_inv:
        for time in times_inv:
            plt.plot([time,time],[-10,10],'b',lw=2)
    plt.ylim([-5,6])

    # Plot velocities
    ax = plt.subplot(gs[1:4, :]) 
    coloptions=['r','b','g','limegreen','gold']
    markoptions=['o','o','o','o','o','o']

    #for i in range(2000,2017):
    #    plt.plot([i,i],[-5,120],color='0.5',lw=0.5)
    x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax.xaxis.set_major_formatter(x_formatter)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    matplotlib.rc('font',family="Arial",)
    ax.tick_params('both', length=6, width=1.25, which='major')
    ax.tick_params('both', length=6, width=1.25, which='minor')
    for i in range(0,len(dists_eul)):
        nonnan = np.where(~(np.isnan(vel_val[:,i])))[0]
        plt.plot(vel_time[nonnan],(vel_val[nonnan,i])/1e3,markoptions[i],color=coloptions[i],\
            label=glacier[0]+'{0:02d}'.format(int(abs(dists_eul[i]))),markersize=5,mec='k')
        nonnan = np.where(vel_time_op[:,0] < 2008)[0]
        plt.plot(vel_time_op[nonnan,0],vel_val_op[nonnan,i]/1e3,'^',color=coloptions[i],markersize=5,mec='k')
    plt.legend(loc=2,borderpad=0.3,fontsize=14,numpoints=1,handlelength=0.2,handletextpad=0.5,labelspacing=0.1,\
            ncol=2,columnspacing=0.8)
    plt.yticks(np.arange(2,14,2),fontsize=16,fontname="Arial")
    plt.ylabel('Glacier speed \n (km/yr)',fontsize=16,fontname="Arial")
    labels=[]
    plt.xticks(np.arange(2000,2017,2),fontsize=16,fontname="Arial")
    ax.set_xticklabels([])
    if vbars == 'runoff':
        for i in range(0,len(year_runoff)):
            path = matplotlib.path.Path([[day1_runoff[i],-3],[day1_runoff[i],12],[day2_runoff[i],12],[day2_runoff[i],-3]])
            patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
            ax.add_patch(patch)  
    elif vbars == 'seasonal':
        for i in range(0,len(years)):
            path = matplotlib.path.Path([[years[i]+0.25,-3],[years[i]+0.25,12],[years[i]+0.75,12],[years[i]+0.75,-3]])
            patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
            ax.add_patch(patch)  
    plt.xlim([time1,time2])
    if plot_inv:
        for time in times_inv:
            plt.plot([time,time],[0,15],'b',lw=2)
    if glacier == 'Helheim':
        plt.ylim([3,11])
    elif glacier == 'Kanger':
        plt.ylim([1.5,13])
    #plt.text(2008.25,3.8,'b',fontsize=16,fontname='Arial',fontweight='bold')

    # Plot surface elevations
    plt.subplot(gs[4:, :])
    ax = plt.gca()
    #for i in range(2000,2017):
    #    plt.plot([i,i],[-30,400],color='0.5',lw=2)
    # Set up points for legend
    if vbars == 'runoff':
        for i in range(0,len(year_runoff)):
            path = matplotlib.path.Path([[day1_runoff[i],-50],[day1_runoff[i],240],[day2_runoff[i],240],[day2_runoff[i],-50]])
            patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
            ax.add_patch(patch)
    elif vbars == 'seasonal':
        for i in range(0,len(years)):
            path = matplotlib.path.Path([[years[i]+0.25,-50],[years[i]+0.25,240],[years[i]+0.75,240],[years[i]+0.75,-50]])
            patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
            ax.add_patch(patch)  
    x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax.xaxis.set_major_formatter(x_formatter)
    plt.xticks(np.arange(2000,2017,2))
    labels=[]
    for i in np.arange(2000,2017,2):
        labels.append('Jan \n'+str(i))
    x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax.xaxis.set_major_formatter(x_formatter)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params('both', length=6, width=1.25, which='major')
    ax.tick_params('x', length=6, width=1.25, which='minor')
    ax.set_xticklabels(labels,fontsize=16,fontname='Arial')
    plt.xlim([time1,time2])
    #plt.legend(handles,labels,loc=3,borderpad=0.3,handleheight=0.1,fontsize=10,numpoints=1,handlelength=0.4,labelspacing=0.05,ncol=2,columnspacing=0.7,handletextpad=0.5)
  #plt.text(2008.25,-2,'c',fontsize=16,fontname='Arial',fontweight='bold')

    if glacier == 'Helheim':
        nonnan = np.where(~(np.isnan(zpt_dem[:,1])))[0]
        ax.plot(time_dem[nonnan],zpt_dem[nonnan,1],'o',color=coloptions[1],markersize=5,mec='k')
    elif glacier == 'Kanger':
        nonnan = np.where(~(np.isnan(zpt_dem[:,1])))[0]
        ax.plot(time_dem[nonnan],zpt_dem[nonnan,1],'o',color=coloptions[1],markersize=5,mec='k')
    if plot_inv:
        for time in times_inv:
            plt.plot([time,time],[0,500],'b',lw=1.5)
    ax.set_ylabel('Surface elevation \n (m)',fontsize=16,fontname='Arial')
    if glacier == 'Helheim':
        #ax.set_yticks(np.arange(60,180,20))
        ax.set_ylim([100,300])
        plt.yticks(np.arange(100,350,100),fontsize=16,fontname='arial')
    elif glacier == 'Kanger':
        ax.set_ylim([50,390])
        plt.yticks(np.arange(100,400,100),fontsize=16,fontname='arial')

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.02,wspace=0.02,top=0.96,right=0.98,left=0.13,bottom=0.09) 
    plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_all.pdf"),FORMAT='PDF',dpi=600)
    plt.close()

    # Image for plotting
    if glacier == "Helheim":
        imagetime = datelib.date_to_fracyear(2014,7,4)
        ximage,yimage,image = geotifflib.readrgb(os.path.join(os.getenv("DATA_HOME"),"Imagery/Landsat/Helheim/TIF/20140704140535_LC82330132014185LGN00.tif"))
    elif glacier == "Kanger":
        imagetime = datelib.date_to_fracyear(2014,7,6)
        ximage,yimage,image = geotifflib.readrgb(os.path.join(os.getenv("DATA_HOME"),"Imagery/Landsat/Kanger/TIF/20140706135251_LC82310122014187LGN00.tif"))

    extent = meshlib.shp_to_xy(os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/Glaciers/3D/"+glacier+"/glacier_extent_inversion_front.shp"))

    # Load velocity record
    xvel = np.arange(np.min(ximage),np.max(ximage),100)
    yvel = np.arange(np.min(yimage),np.max(yimage),100)
    vx,vy = vellib.inversion_3D(glacier,xvel,yvel,imagetime,dir_velocity_out='none',blur=False)
    vel = np.sqrt(vx**2+vy**2)
    del vx,vy

    # Load mask
    xmask,ymask,mask = masklib.load_grid(glacier,np.min(xvel),np.max(xvel),np.min(yvel),np.max(yvel),100,icefront_time=datelib.date_to_fracyear(2014,7,4))
    vel_masked = np.ma.masked_array(vel,mask)
  
    fig = plt.figure(figsize=(2.5,2.5))

    cx = cubehelix.cmap(start=1.2,rot=-1.1,reverse=True,minLight=0.1,sat=2)
    p=plt.imshow(vel_masked/1e3,extent=[np.min(xvel),np.max(xvel),np.min(yvel),np.max(yvel)],origin='lower',clim=[0,10],cmap=cx)
    ax = plt.gca()
    ax.imshow(image[:,:,0],extent=[np.min(ximage),np.max(ximage),np.min(yimage),np.max(yimage)],cmap='Greys_r',origin='lower',clim=[0,0.6])
    ax.imshow(vel_masked/1e3,extent=[np.min(xvel),np.max(xvel),np.min(yvel),np.max(yvel)],origin='lower',alpha=0.5,clim=[0,10],cmap=cx)
    ax.axis('equal')
    plt.plot(np.r_[extent[:,0],extent[0,0]],np.r_[extent[:,1],extent[0,1]],'k',lw=1.5)

    if glacier == 'Kanger':
        xmin = 468000.
        xmax = 498000.
        ymin = -2299000.
        ymax = -2264000.
    elif glacier == 'Helheim':
        xmin = 283000.
        xmax = 313000.
        ymin = -2587000.
        ymax = -2552000.

    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_ylim([ymin,ymax])
    ax.set_xlim([xmin,xmax])

    plt.tight_layout()
  
    xmin,xmax = plt.xlim()
    ymin,ymax = plt.ylim()
    matplotlib.rcParams['hatch.linewidth'] = 0.5
    path = matplotlib.path.Path([[0.6*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
  			[0.99*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
  			[0.99*(xmax-xmin)+xmin,0.72*(ymax-ymin)+ymin],
  			[0.6*(xmax-xmin)+xmin,0.72*(ymax-ymin)+ymin],
  			[0.6*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin]])
    patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1)
    ax.add_patch(patch)
    cbaxes = fig.add_axes([0.64, 0.92, 0.28, 0.02]) 
    cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=[0,5,10]) 
    ax.text(xmin+0.63*(xmax-xmin),ymin+0.81*(ymax-ymin),'Speed (km yr$^{-1}$)',fontsize=8)
    cb.ax.tick_params(labelsize=8)
    ax.plot([xmin+0.63*(xmax-xmin),xmin+0.63*(xmax-xmin)+5e3],[ymin+0.78*(ymax-ymin),ymin+0.78*(ymax-ymin)],'k',linewidth=1.5)
    ax.plot([xmin+0.63*(xmax-xmin),xmin+0.63*(xmax-xmin)],[ymin+0.78*(ymax-ymin),ymin+0.76*(ymax-ymin)],'k',linewidth=1.5)
    ax.plot([xmin+0.63*(xmax-xmin)+5e3,xmin+0.63*(xmax-xmin)+5e3],[ymin+0.78*(ymax-ymin),ymin+0.76*(ymax-ymin)],'k',linewidth=1.5)
    ax.text(xmin+0.67*(xmax-xmin)+5e3,ymin+0.75*(ymax-ymin),'5 km',fontsize=8)

    ax2 = fig.add_axes([0.04, 0.04, 0.2, 0.36])
    ax2.axis('off')
    m = Basemap(width=1600000,height=3000000,
            resolution='l',projection='stere',\
            lat_ts=60,lat_0=72,lon_0=-41.)
    m.drawcoastlines(linewidth=0.75)
    m.drawmapboundary(fill_color='lightblue')
    m.fillcontinents(color='w',lake_color='w')
    if glacier == 'Helheim':
        xg,yg = m([-38.3],[66.1])
    elif glacier == 'Kanger':
        xg,yg = m ([-33.0],[68.63333])
    m.plot(xg,yg,'ro',mec='k',mew=0.5,markersize=3.5)

    plt.subplots_adjust(wspace=0,hspace=0,top=0.98,right=0.98,left=0.02,bottom=0.02)

    plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_velocity_map.pdf"),transparent=True,FORMAT='PDF',dpi=600)
    plt.close()

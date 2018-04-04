import os, elmerreadlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.ndimage

glacier = 'Kanger'
output = 'zs'
temperature = 'modelT'
beta_date = '20120522'
beta_file = '1e13_SSA_DEM'
beta_suffix = beta_file+beta_date+'_'+temperature
slidinglaw = 'linear'
maindir = '/Users/kehrl/Models/'+glacier+'/3D/'
dirs = os.listdir(maindir)

rho_i = 917.
rho_w = 1025.
g = 9.8

if glacier == 'Kanger':
    fig = plt.figure(figsize=(7.5,11))
elif glacier == 'Helheim':
    fig = plt.figure(figsize=(7.5,15))
matplotlib.rc('font',family='Arial')
gs = matplotlib.gridspec.GridSpec(6,4)

if output == 'zs':
    x_orig,y_orig,grid_orig = elmerreadlib.input_file(maindir+'DEM'+beta_date+'_'+temperature+'/inputs/zsdem.xy')
elif output == 'velocity':
    x_orig,y_orig,u_orig = elmerreadlib.input_file(maindir+'DEM'+beta_date+'_'+temperature+'/inputs/udem.xy')
    x_orig,y_orig,v_orig = elmerreadlib.input_file(maindir+'DEM'+beta_date+'_'+temperature+'/inputs/vdem.xy')
    grid_orig = np.sqrt(u_orig**2+v_orig**2)
elif output == 'N':
    x_orig,y_orig,zs_orig = elmerreadlib.input_file(maindir+'DEM'+beta_date+'_'+temperature+'/inputs/zsdem.xy')
    x_orig,y_orig,zb_orig = elmerreadlib.input_file(maindir+'DEM'+beta_date+'_'+temperature+'/inputs/zbdem.xy')
    grid_orig = rho_i*g*(zs_orig-zb_orig)
    grid_orig[zb_orig < 0] = grid_orig[zb_orig < 0] + rho_w*g*zb_orig[zb_orig < 0]
    grid_orig[(grid_orig < 0) | (zs_orig < zb_orig)] = np.float('NaN')
    beta_orig = np.loadtxt(maindir+'DEM'+beta_date+'_'+temperature+'/inputs/beta_'+slidinglaw+'_'+beta_suffix+'.dat',skiprows=1)

if output == 'zs':
    vmin = -500.0
    vmax = 500.0
elif output == 'velocity':
    vmin = -20.0
    vmax = 20.0
cmap = 'RdYlBu_r'

n = 0
for dir in dirs:
    if dir.startswith('DEM') and dir.endswith(temperature):
        date = dir[3:11]

        if output == 'zs':
            x,y,grid = elmerreadlib.input_file(maindir+dir+'/inputs/zsdem.xy')
        elif output == 'velocity':
            x,y,u = elmerreadlib.input_file(maindir+dir+'/inputs/udem.xy')
            x,y,v = elmerreadlib.input_file(maindir+dir+'/inputs/vdem.xy')
            grid = np.sqrt(u**2+v**2)
        elif output == 'N':
            x,y,zs = elmerreadlib.input_file(maindir+dir+'/inputs/zsdem.xy')
            x,y,zb = elmerreadlib.input_file(maindir+dir+'/inputs/zbdem.xy')
            grid = rho_i*g*(zs-zb)
            grid[zb < 0] = grid[zb < 0] + rho_w*g*zb[zb < 0]
            grid[(zs < zb) | (grid < 0)] = np.float('NaN')
            
        if n == 0:
            extent = np.loadtxt(maindir+dir+'/inputs/mesh_extent.dat')
            xmin = np.floor(np.min(extent[:,0])/100)*100
            xmax = np.ceil(np.max(extent[:,0])/100)*100
            ymin = np.floor(np.min(extent[:,1])/100)*100
            ymax = np.ceil(np.max(extent[:,1])/100)*100
            try:
                hole1 = np.loadtxt(maindir+dir+"/inputs/mesh_hole1.dat")
                hole2 = np.loadtxt(maindir+dir+"/inputs/mesh_hole2.dat")
                holes=[hole1,hole2]
            except:
                holes = []
            path = matplotlib.path.Path(extent)
            xgrid,ygrid = np.meshgrid(x_orig,y_orig)
            inside = path.contains_points(np.column_stack([xgrid.flatten(),ygrid.flatten()]))
            mask = ~(np.reshape(inside,(len(y_orig),len(x_orig))))
            if len(holes) > 0:
                for i in range(0,len(holes)):
                    path = matplotlib.path.Path(holes[i])
                    inside = path.contains_points(np.column_stack([xgrid.flatten(),ygrid.flatten()]))
                    maskhole = (np.reshape(inside,(len(y_orig),len(x_orig))))
                    mask[maskhole == 1] = 1
        if (len(x_orig) != len(x)) or (len(y_orig) != len(y)):
            f = scipy.interpolate.RegularGridInterpolator((y,x),grid,bounds_error=False)
            grid = f((ygrid.flatten(),xgrid.flatten())).reshape((len(y_orig),len(x_orig)))
        
        if output == 'zs':
            grid_diff = grid-grid_orig
        elif output == 'velocity':
            grid_diff = (grid-grid_orig)/grid_orig*100.0
        elif output == 'N':
            grid_diff = grid/grid_orig

        grid_diff[mask == 1] = np.float('NaN')

        if output == 'N':
            if slidinglaw == 'linear':
                m = 1
            elif slidinglaw == 'weertman':
                m = 3
            nonnan = np.where(~(np.isnan(grid_diff)))
            xnonnan = xgrid[nonnan].flatten()
            ynonnan = ygrid[nonnan].flatten()
            gnonnan = grid_diff[nonnan].flatten()
            beta_new = (scipy.interpolate.griddata((xnonnan,ynonnan),gnonnan,(beta_orig[:,0],beta_orig[:,1]),\
                    method = 'nearest')**m)*beta_orig[:,2]

            fid = open(maindir+dir+'/inputs/beta_'+slidinglaw+'_'+beta_suffix+'_N.dat','w')
            fid.write('{0}\n'.format(len(beta_orig)))
            for i in range(0,len(beta_orig)):
                fid.write('{0} {1} {2}\n'.format(beta_orig[i,0],beta_orig[i,1],beta_new[i]))
            fid.close()

        
        plt.subplot(gs[n])
        if output == 'zs':
            im=plt.imshow(grid_diff,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',\
                    vmin=vmin,vmax=vmax,cmap=cmap,norm=matplotlib.colors.SymLogNorm(2))
        elif output == 'velocity':
            im=plt.imshow(grid_diff,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',\
                    vmin=vmin,vmax=vmax,cmap=cmap)
        elif output == 'N':
            im=plt.imshow(grid_diff,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',\
                vmin=0.75,vmax=1.25,cmap=cmap)
        plt.plot(np.r_[extent[:,0],extent[0,0]],np.r_[extent[:,1],extent[0,1]],'k')
        plt.axis('equal')
        plt.xlim([xmin,xmax])
        plt.ylim([ymin,ymax])
        plt.yticks([])
        plt.xticks([])
        if len(holes) > 0:
            for i in range(0,len(holes)):
                plt.plot(np.r_[holes[i][:,0],holes[i][0,0]],np.r_[holes[i][:,1],holes[i][0,1]],'k')
        if glacier == 'Kanger':
            plt.text(xmin+0.05*(xmax-xmin),ymin+0.05*(ymax-ymin),date,fontsize=10)
            if output == 'velocity':
                plt.text(xmax-0.25*(xmax-xmin),ymax-0.1*(ymax-ymin),\
                    'Change',fontsize=10)
                plt.text(xmax-0.25*(xmax-xmin),ymax-0.2*(ymax-ymin),\
                    '{0:>4.1f}%'.format(np.nanmean((grid_diff[grid>1000]))),fontsize=10)
        elif glacier == 'Helheim':
            plt.text(xmin+0.6*(xmax-xmin),ymin-0.03*(ymin-ymax),date,fontsize=10)

        n = n+1

if glacier == 'Kanger':
    cbar_ax1 = fig.add_axes([0.25,0.05,0.5,0.02])
elif glacier == 'Helheim':
    cbar_ax1 = fig.add_axes([0.25,0.03,0.5,0.015])
cb1 = fig.colorbar(im, cax=cbar_ax1,orientation='horizontal',extend='both')
if output == 'zs':
    cb1.set_label(r'Elevation change (m)',fontsize=10,fontname='Arial')
elif output == 'velocity':
    cb1.set_label(r'Velocity change (%)',fontsize=10,fontname='Arial')
elif output == 'N':
    cb1.set_label(r'$N/N_o$',fontsize=10,fontname='Arial')

plt.tight_layout()
if glacier == 'Kanger':
    plt.subplots_adjust(hspace=0.03,wspace=0.03,top=0.99,right=0.99,left=0.01,bottom=0.08)
elif glacier == 'Helheim':
    plt.subplots_adjust(hspace=0.02,wspace=0.03,top=0.995,right=0.99,left=0.01,bottom=0.05)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+'_misfit_'+output+'_'+beta_suffix+'.pdf'),FORMAT='PDF',dpi=400)
plt.close()


import os, elmerreadlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import inverselib

glacier = 'Kanger'
temperature = 'model'
maindir = '/Users/kehrl/Models/'+glacier+'/3D/'
dirs = os.listdir(maindir)
beta_date = 'average'
slidinglaw = 'linear'
cutoff = 1000

if temperature == 'model':
    temperature_text = 'modelT'
else:
    temperature_text = 'constantT'

if glacier == 'Kanger':
    fig = plt.figure(figsize=(7.5,11))
elif glacier == 'Helheim':
    fig = plt.figure(figsize=(7.5,15))
matplotlib.rc('font',family='Arial')
gs = matplotlib.gridspec.GridSpec(6,4)

if beta_date == '':
    no_beta_file = True
else:
    no_beta_file = False

vmin = -20
vmax = 20
cmap = 'RdYlBu_r'

# Get indices where velocity is always greater than the cutoff value
x_cutoff,y_cutoff,vsurfini_cutoff,ind_cutoff_grid, ind_cutoff  = inverselib.get_velocity_cutoff(glacier,velocity_cutoff=cutoff)

# Now create plots and get misfit values for velocities > cutoff        
n = 0
for dir in dirs:
    if dir.startswith('DEM') and dir.endswith(temperature_text):
        if no_beta_file:
            beta_date = dir[3:11]
            beta_suffix = '1e13_SSA_DEM'+beta_date+'_'+temperature_text+'_'+slidinglaw
        else:
            beta_suffix = '1e13_SSA_average_2000_2016_'+temperature_text+'_'+slidinglaw
       
        data = elmerreadlib.pvtu_file(maindir+dir+'/mesh2d/steady_'+beta_suffix+\
                '0001.pvtu',['vsurfini','ssavelocity'])
        surf = elmerreadlib.values_in_layer(data,'surf')
        if n == 0:
            # Mesh boundaries
            extent = np.loadtxt(maindir+dir+"/inputs/mesh_extent.dat")
            try:
                hole1 = np.loadtxt(maindir+dir+"/inputs/mesh_hole1.dat")
                hole2 = np.loadtxt(maindir+dir+"/inputs/mesh_hole2.dat")
                holes=[hole1,hole2]
            except:
                holes = []

        surf_misfit = np.zeros(surf.shape,dtype=np.dtype(surf.dtype.descr+[('misfit',np.float64)]))
        for name in surf.dtype.descr:
            surf_misfit[name[0]] = surf[name[0]]
        surf_misfit['misfit'] = (surf['ssavelocity']-surf['vsurfini'])/surf['vsurfini']*100

        x,y,grid = elmerreadlib.grid3d(surf_misfit,'misfit',holes,extent,dx=50)
        grid[ind_cutoff_grid] = np.float('NaN')
        grid[np.isnan(vsurfini_cutoff)] = np.float('NaN') 

        plt.subplot(gs[n])
        im=plt.imshow(grid,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',\
                    vmin=vmin,vmax=vmax,cmap=cmap)
        plt.plot(np.r_[extent[:,0],extent[0,0]],np.r_[extent[:,1],extent[0,1]],'k')
        if len(holes) > 0:
            for hole in holes:
                plt.plot(np.r_[hole[:,0],hole[0,0]],np.r_[hole[:,1],hole[0,1]],'k')
        plt.contour(x_cutoff,y_cutoff,vsurfini_cutoff,[cutoff],colors='k',linestyles='dashed')
        plt.axis('equal')
        xmin,xmax = plt.xlim()
        ymin,ymax = plt.ylim()
        plt.yticks([])
        plt.xticks([])
        #ind_cutoff = np.where(surf_misfit['vsurfini'] > cutoff)
        if glacier == 'Kanger':
            plt.text(xmax-0.25*(xmax-xmin),ymax-0.1*(ymax-ymin),\
                'MAR',fontsize=10)
            plt.text(xmax-0.25*(xmax-xmin),ymax-0.2*(ymax-ymin),\
                '{0:>4.1f}%'.format(np.nanmean(abs(surf_misfit['misfit'][ind_cutoff]))),fontsize=10)
            plt.text(xmin+0.05*(xmax-xmin),ymin+0.05*(ymax-ymin),dir[3:11],fontsize=10)
        elif glacier == 'Helheim':
            plt.text(xmax-0.14*(xmax-xmin),ymax-0.08*(ymax-ymin),\
                'MAR',fontsize=10)
            plt.text(xmax-0.16*(xmax-xmin),ymax-0.16*(ymax-ymin),\
                '{0:>4.1f}%'.format(np.nanmean(abs(surf_misfit['misfit'][ind_cutoff]))),fontsize=10)
            plt.text(xmin+0.6*(xmax-xmin),ymin-0.03*(ymin-ymax),dir[3:11],fontsize=10)

        n = n+1

if glacier == 'Kanger':
    cbar_ax1 = fig.add_axes([0.25,0.05,0.5,0.02])
elif glacier == 'Helheim':
    cbar_ax1 = fig.add_axes([0.25,0.03,0.5,0.015])
cb1 = fig.colorbar(im, cax=cbar_ax1,orientation='horizontal',extend='both')
cb1.set_label(r'Residual (%)',fontsize=10,fontname='Arial')

plt.tight_layout()
if glacier == 'Kanger':
    plt.subplots_adjust(hspace=0.03,wspace=0.03,top=0.99,right=0.99,left=0.01,bottom=0.08)
elif glacier == 'Helheim':
    plt.subplots_adjust(hspace=0.02,wspace=0.03,top=0.995,right=0.99,left=0.01,bottom=0.05)

if no_beta_file:
    beta_suffix = '1e13_SSA_'+temperature_text+'_'+slidinglaw
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+'_misfit_beta_'+beta_suffix+'.pdf'),FORMAT='PDF',dpi=400)
print "Saving as "+os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+'_misfit_beta_'+beta_suffix+'.pdf')
plt.close()


import os, sys, elmerreadlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.ndimage
import argparse

args = sys.argv

# Get inputs to file
parser = argparse.ArgumentParser()
parser.add_argument("-glacier",dest="glacier",required = True,
        help = "Name of glacier.")
parser.add_argument("-temperature", dest="temperature",required = True)

args, _ = parser.parse_known_args(sys.argv)
glacier = args.glacier
T = args.temperature

if T == 'model':
    temperature = 'modelT'
else:
    temperature = 'constantT'

maindir = '/Users/kehrl/Models/'+glacier+'/3D/'
dirs = os.listdir(maindir)

if glacier == 'Kanger':
    fig = plt.figure(figsize=(7.5,11))
elif glacier == 'Helheim':
    fig = plt.figure(figsize=(7.5,15))
matplotlib.rc('font',family='Arial')
gs = matplotlib.gridspec.GridSpec(6,4)

vmin = -500.0
vmax = 500.0
cmap = 'RdYlBu_r'

rho_i = 917.
g = 9.8

n = 0
for dir in dirs:
    if dir.startswith('DEM') and dir.endswith(temperature):
        
        x,y,u = elmerreadlib.input_file(maindir+dir+'/inputs/udem.xy')
        x,y,v = elmerreadlib.input_file(maindir+dir+'/inputs/vdem.xy')
        x,y,zs = elmerreadlib.input_file(maindir+dir+'/inputs/zsdem.xy')
        x,y,zb = elmerreadlib.input_file(maindir+dir+'/inputs/zbdem.xy')
        
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
        xgrid,ygrid = np.meshgrid(x,y)
        inside = path.contains_points(np.column_stack([xgrid.flatten(),ygrid.flatten()]))
        mask = ~(np.reshape(inside,(len(y),len(x))))
        if len(holes) > 0:
            for i in range(0,len(holes)):
                path = matplotlib.path.Path(holes[i])
                inside = path.contains_points(np.column_stack([xgrid.flatten(),ygrid.flatten()]))
                maskhole = (np.reshape(inside,(len(y),len(x))))
                mask[maskhole == 1] = 1

        # Calculate grad(h) & H for driving stress
        nx = len(x)
        ny = len(y)
        dhdx = np.zeros([ny,nx])
        dhdy = np.zeros([ny,nx])
        dhdx[:,:] = float('nan')
        dhdy[:,:] = float('nan')

        H = zs-zb

        dhdx[1:-1,1:-1] = (zs[1:-1,2:]-zs[1:-1,0:-2])/(x[2:]-x[0:-2])
        dhdy[1:-1,1:-1] = ((zs[2:,1:-1]-zs[0:-2,1:-1]).T/(y[2:]-y[0:-2])).T
        taud = -rho_i*g*H*(dhdx*u+dhdy*v)/np.sqrt(u**2+v**2)
        taud_blur = scipy.ndimage.filters.gaussian_filter(taud,sigma=2,truncate=4)
        taud_blur[mask == 1] = np.float('NaN')
        taud[mask == 1] = np.float('NaN')

        plt.subplot(gs[n])
        im=plt.imshow(taud/1e3,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',\
                    vmin=vmin,vmax=vmax,cmap=cmap)
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
            plt.text(xmax-0.2*(xmax-xmin),ymax-0.1*(ymax-ymin),\
                r'$\overline{\tau_d}$',fontsize=10)
            plt.text(xmax-0.3*(xmax-xmin),ymax-0.2*(ymax-ymin),\
                '{0} kPa'.format(int(np.nanmean(taud)/1e3)),fontsize=10)
            plt.text(xmin+0.05*(xmax-xmin),ymin+0.05*(ymax-ymin),dir[3:11],fontsize=10)
        elif glacier == 'Helheim':
            plt.text(xmax-0.1*(xmax-xmin),ymax-0.08*(ymax-ymin),\
                r'$\overline{\tau_d}$',fontsize=10)
            plt.text(xmax-0.12*(xmax-xmin),ymax-0.16*(ymax-ymin),\
                '{0}'.format(int(np.nanmean(taud)/1e3)),fontsize=10)
            plt.text(xmax-0.12*(xmax-xmin),ymax-0.24*(ymax-ymin),\
                'kPa',fontsize=10)
            plt.text(xmin+0.6*(xmax-xmin),ymin-0.03*(ymin-ymax),dir[3:11],fontsize=10)

        n = n+1

if glacier == 'Kanger':
    cbar_ax1 = fig.add_axes([0.25,0.05,0.5,0.02])
elif glacier == 'Helheim':
    cbar_ax1 = fig.add_axes([0.25,0.03,0.5,0.015])
cb1 = fig.colorbar(im, cax=cbar_ax1,orientation='horizontal',extend='both')
cb1.set_label(r'$\tau_d$ (kPa)',fontsize=10,fontname='Arial')

plt.tight_layout()
if glacier == 'Kanger':
    plt.subplots_adjust(hspace=0.03,wspace=0.03,top=0.99,right=0.99,left=0.01,bottom=0.08)
elif glacier == 'Helheim':
    plt.subplots_adjust(hspace=0.02,wspace=0.03,top=0.995,right=0.99,left=0.01,bottom=0.05)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+'_taud.pdf'),FORMAT='PDF',dpi=400)
plt.close()


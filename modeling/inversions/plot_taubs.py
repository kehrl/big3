import os, sys, elmerreadlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse

args = sys.argv

# Get inputs to file
parser = argparse.ArgumentParser()
parser.add_argument("-glacier",dest="glacier",required = True,
        help = "Name of glacier.")
parser.add_argument("-regpar", dest="regpar",required = False,
        default = '1e13', help = "Regularization parameter.")
parser.add_argument("-output", dest="output",required = False,
        default = 'taub', help = "Output (taub or betasquared).")
parser.add_argument("-temperature", dest="temperature",required = True)

args, _ = parser.parse_known_args(sys.argv)
glacier = args.glacier
regpar = args.regpar
output = args.output
T = args.temperature

if T == 'model':
    temperature = 'modelT'
else:
    temperature = 'constantT'

maindir = '/Users/kehrl/Models/'+glacier+'/3D/'
dirs = os.listdir(maindir)

fig = plt.figure(figsize=(7.5,15))
matplotlib.rc('font',family='Arial')
gs = matplotlib.gridspec.GridSpec(6,4)

if output == 'taub':
    vmin = 0.0
    vmax = 0.4e3
    cmap = 'viridis'
else:
    vmin = 1.0e-3
    vmax = 2.0
    cmap = 'viridis'

n = 0
for dir in dirs:
    if dir.startswith('DEM') and dir.endswith(temperature):
        try:
            subdirs = os.listdir(maindir+dir+'/mesh2d/inversion_adjoint')
            for subdir in subdirs:
                if (regpar in subdir) and not(subdir.endswith('.pdf')):
                    bed = elmerreadlib.pvtu_file(maindir+dir+'/mesh2d/inversion_adjoint/'+subdir+
                        '/adjoint_beta_ssa0001.pvtu',['ssavelocity','beta','vsurfini'])            
            if n == 0:
                # Mesh boundaries
                extent = np.loadtxt(maindir+dir+"/inputs/mesh_extent.dat")
                try:
                    hole1 = np.loadtxt(maindir+dir+"/inputs/mesh_hole1.dat")
                    hole2 = np.loadtxt(maindir+dir+"/inputs/mesh_hole2.dat")
                    holes=[hole1,hole2]
                except:
                    holes = []
            
            x,y,grid = elmerreadlib.grid3d(bed,output,holes,extent,dx=50)

            plt.subplot(gs[n])
            if output == 'taub':
                im=plt.imshow(grid*1e3,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',\
                    vmin=vmin,vmax=vmax,cmap=cmap)
            else:
                im=plt.imshow(grid*1e3,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',\
                    norm=matplotlib.colors.LogNorm(vmin=vmin,vmax=vmax),cmap=cmap)
            plt.axis('equal')
            xmin,xmax = plt.xlim()
            ymin,ymax = plt.ylim()
            plt.yticks([])
            plt.xticks([])
            if output == 'taub':
                if glacier == 'Kanger':
                    plt.text(xmax-0.2*(xmax-xmin),ymax-0.1*(ymax-ymin),\
                        r'$\overline{\tau_b}$',fontsize=10)
                    plt.text(xmax-0.3*(xmax-xmin),ymax-0.2*(ymax-ymin),\
                        '{0} kPa'.format(int(np.nanmean(grid)*1e3)),fontsize=10)
                elif glacier == 'Helheim':
                    plt.text(xmax-0.1*(xmax-xmin),ymax-0.08*(ymax-ymin),\
                        r'$\overline{\tau_b}$',fontsize=10)
                    plt.text(xmax-0.12*(xmax-xmin),ymax-0.16*(ymax-ymin),\
                        '{0}'.format(int(np.nanmean(grid)*1e3)),fontsize=10)
                    plt.text(xmax-0.12*(xmax-xmin),ymax-0.24*(ymax-ymin),\
                        'kPa',fontsize=10)
            if glacier == 'Kanger':
                plt.text(xmin+0.05*(xmax-xmin),ymin-0.05*(ymin-ymax),dir[3:11],fontsize=10)
            elif glacier == 'Helheim':
                plt.text(xmin+0.6*(xmax-xmin),ymin-0.03*(ymin-ymax),dir[3:11],fontsize=10)

            n = n+1
        except:
            pass

if glacier == 'Kanger':
    cbar_ax1 = fig.add_axes([0.25,0.05,0.5,0.02])
elif glacier == 'Helheim':
    cbar_ax1 = fig.add_axes([0.25,0.03,0.5,0.015])
if output == 'taub':
    cb1 = fig.colorbar(im, cax=cbar_ax1,orientation='horizontal',extend='max')
    cb1.set_label(r'$\tau_b$ (kPa)',fontsize=10,fontname='Arial')
else:
    cb1 = fig.colorbar(im, cax=cbar_ax1,orientation='horizontal',extend='both')
    cb1.set_label(r'$\beta$ (kPa yr / m)',fontsize=10,fontname='Arial')

plt.tight_layout()
if glacier == 'Kanger':
    plt.subplots_adjust(hspace=0.03,wspace=0.03,top=0.99,right=0.99,left=0.01,bottom=0.08)
elif glacier == 'Helheim':
    plt.subplots_adjust(hspace=0.02,wspace=0.03,top=0.995,right=0.99,left=0.01,bottom=0.05)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+'_'+output+'_'+temperature+'_'+regpar+'.pdf'),FORMAT='PDF',dpi=400)
plt.close()


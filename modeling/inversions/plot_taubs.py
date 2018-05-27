import os, sys, elmerreadlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse
import datelib

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
parser.add_argument("-model", dest="model", required = True)

args, _ = parser.parse_known_args(sys.argv)
glacier = args.glacier
regpar = args.regpar
output = args.output
T = args.temperature
model = args.model

if T == 'model':
    temperature = 'modelT'
    temperature_text = 'ModelT'
else:
    temperature = 'constantT'
    temperature_text = 'ConstantT'

maindir = os.path.join(os.getenv("MODEL_HOME"),glacier+'/3D/INV_'+model+'_'+temperature_text+'/')
dirs = os.listdir(maindir)

if glacier == 'Helheim':
    fig = plt.figure(figsize=(7.5,7.5))
elif glacier == 'Kanger':
    fig = plt.figure(figsize=(7.5,6.2))
matplotlib.rc('font',family='Arial')
gs = matplotlib.gridspec.GridSpec(4,6)

if output == 'taub':
    vmin = 0.0
    vmax = 0.4e3
    cmap = 'viridis'
else:
    vmin = 1.0e-3
    vmax = 2.0
    cmap = 'viridis'

n = 0
for dir in np.sort(dirs):
    if dir.startswith('DEM') and dir.endswith(temperature):
        subdirs = os.listdir(maindir+dir+'/mesh2d/inversion_adjoint')
        for subdir in subdirs:
            if (regpar in subdir) and not(subdir.endswith('.pdf')):
		if os.path.isfile(maindir+dir+'/mesh2d/inversion_adjoint/'+subdir+\
		        '/adjoint_beta_ssa0001.pvtu'):
		    data = elmerreadlib.pvtu_file(maindir+dir+'/mesh2d/inversion_adjoint/'+subdir+
                        '/adjoint_beta_ssa0001.pvtu',['ssavelocity','beta','vsurfini'])
	        else:
                    data = elmerreadlib.pvtu_file(maindir+dir+'/mesh2d/inversion_adjoint/'+subdir+
		        '/adjoint_beta0001.pvtu',['velocity','beta','vsurfini'])
	        bed = elmerreadlib.values_in_layer(data,'bed')
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
        if glacier == 'Kanger':
            ymin = ymin+1.2e3
            ymax = ymax+1.2e3
            plt.ylim(ymin,ymax)
        plt.yticks([])
        plt.xticks([])
        #if output == 'taub':
        #    if glacier == 'Kanger':
        #        plt.text(xmax-0.2*(xmax-xmin),ymax-0.1*(ymax-ymin),\
        #            r'$\overline{\tau_b}$',fontsize=10)
        #        plt.text(xmax-0.3*(xmax-xmin),ymax-0.2*(ymax-ymin),\
        #            '{0} kPa'.format(int(np.nanmean(grid)*1e3)),fontsize=10)
        #    elif glacier == 'Helheim':
        #        plt.text(xmax-0.1*(xmax-xmin),ymax-0.08*(ymax-ymin),\
        #            r'$\overline{\tau_b}$',fontsize=10)
        #        plt.text(xmax-0.12*(xmax-xmin),ymax-0.16*(ymax-ymin),\
        #            '{0}'.format(int(np.nanmean(grid)*1e3)),fontsize=10)
        #        plt.text(xmax-0.12*(xmax-xmin),ymax-0.24*(ymax-ymin),\
        #            'kPa',fontsize=10)
        plt.text(xmin+0.30*(xmax-xmin),ymax-0.07*(ymax-ymin),dir[9:11]+' '+datelib.month_to_string(int(dir[7:9]))+' '+dir[3:7],\
                    fontsize=10,bbox=dict(facecolor='w', edgecolor='none', pad=1))

        n = n+1

if glacier == 'Helheim':
    cbar_ax1 = fig.add_axes([0.337,0.055,0.34,0.015])
elif glacier == 'Kanger':
    cbar_ax1 = fig.add_axes([0.337,0.07,0.34,0.02])
if output == 'taub':
    cb1 = fig.colorbar(im, cax=cbar_ax1,orientation='horizontal',extend='max')
    cb1.set_label(r'$\tau_b$ (kPa)',fontsize=10,fontname='Arial')
else:
    cb1 = fig.colorbar(im, cax=cbar_ax1,orientation='horizontal',extend='both')
    cb1.set_label(r'$\beta$ (kPa yr / m)',fontsize=10,fontname='Arial')

if glacier == 'Helheim':
    plt.subplots_adjust(hspace=0.00,wspace=0.0,top=0.995,right=0.99,left=0.01,bottom=0.075)
elif glacier == 'Kanger':
    plt.subplots_adjust(hspace=0.0,wspace=0.0,top=0.995,right=0.99,left=0.01,bottom=0.095)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+'_'+output+'_'+model+'_'+temperature+'_'+regpar+'.pdf'),FORMAT='PDF',dpi=300)
plt.close()


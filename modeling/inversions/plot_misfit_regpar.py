# This code plots basal shear stress and velocity misfit for three different regularization parameters. 
#  
#
# LMK, UW, 5/30/2015

import numpy as np
import sys
import os
import matplotlib
import matplotlib.pyplot as plt
import elmerreadlib, geotifflib
import shapely.geometry
from pylab import *
import argparse

args = sys.argv

# Get inputs to file
parser = argparse.ArgumentParser()
parser.add_argument("-glacier",dest="glacier",required = True,
        help = "Name of glacier.")
parser.add_argument("-mesh", dest="meshname", required = True,
        help = "Name of mesh.")
parser.add_argument("-dim", dest="dimension",required = False,
	default='3D',help = "3D or Flowline. Default is 3D.")
parser.add_argument("-method", dest="method",required = False,
        default='adjoint', help = "Adjoint or robin. Default is adjoint.")
parser.add_argument("-regpar1", dest="regpar1",required = True,
        help = "First regularization parameter for plotting.")
parser.add_argument("-regpar2", dest="regpar2",required = True,
        help = "Second regularization parameter for plotting.")
parser.add_argument("-regpar3", dest="regpar3",required = True,
        help = "Third regularization parameter for plotting.")

args, _ = parser.parse_known_args(sys.argv)
RES = args.meshname
ver = args.dimension
method = args.method
glacier = args.glacier
regpar1 = args.regpar1
regpar2 = args.regpar2
regpar3 = args.regpar3

regpars = [regpar1,regpar2,regpar3]

# Input directory
DIR = os.path.join(os.getenv("MODEL_HOME"),glacier+"/"+ver+"/"+RES+"/")
DIRR = DIR+"mesh2d/inversion_"+method+"/"

if not(os.path.isdir(DIRR)):
    print DIRR
    sys.exit("That is not a model results directory.")

# Get mesh boundary
extent = np.loadtxt(DIR+"inputs/mesh_extent.dat")
try:
    hole1 = np.loadtxt(DIR+"inputs/mesh_hole1.dat")
    hole2 = np.loadtxt(DIR+"inputs/mesh_hole2.dat")
    holes=[hole1,hole2]
except:
    holes = []

if glacier == "Helheim":
    date = '20120316'
    ximage,yimage,image = geotifflib.readrgb(os.path.join(os.getenv("DATA_HOME"),"Imagery/Landsat/Helheim/"+\
        "TIF/20140704140535_LC82330132014185LGN00.tif"))
elif glacier == "Kanger":
    date = '20120213'
    ximage,yimage,image = geotifflib.readrgb(os.path.join(os.getenv("DATA_HOME"),"Imagery/Landsat/Kanger/"+\
        "TIF/20140706135251_LC82310122014187LGN00.tif"))

# Figure
if glacier == 'Helheim':
    fig=plt.figure(figsize=(5,4.4))
elif glacier == 'Kanger':
    fig=plt.figure(figsize=(5,3.2))
matplotlib.rc('font',family='sans-serif',size=10)
gs = matplotlib.gridspec.GridSpec(2,3)

# Get vtu files
dirs = os.listdir(DIRR)
for i in range(0,len(regpars)):
    for dir in dirs:
        if (dir.startswith('lambda_'+regpars[i])) and not(dir.endswith('.pdf')):
            try:
                vtudata = elmerreadlib.pvtu_file(DIRR+dir+'/adjoint_beta_ssa0001.pvtu',['vsurfini','ssavelocity','beta'])
                vname = 'ssavelocity'
            except:
                vtudata = elmerreadlib.pvtu_file(DIRR+dir+'/adjoint_beta0001.pvtu',['vsurfini','velocity','beta'])
                vname = 'velocity'
            surf = elmerreadlib.values_in_layer(vtudata,'surface')
            bed = elmerreadlib.values_in_layer(vtudata,'bed')
            surf_misfit = np.zeros(surf.shape,dtype=np.dtype(surf.dtype.descr+[('misfit',np.float64)]))
            for name in surf.dtype.descr:
                surf_misfit[name[0]] = surf[name[0]]
            surf_misfit['misfit'] = (surf[vname]-surf['vsurfini'])/surf['vsurfini']*100
    
    if i == 0:
        misfits_all = np.zeros([len(surf_misfit['misfit']),3])
        misfits_all[:,:] = np.float('NaN')
    misfits_all[:,i] = surf_misfit['misfit']

    plt.subplot(gs[0,i])
    plt.title(r'$\lambda$ = '+regpars[i],fontsize=10,fontname='Arial',color='k')
    x,y,taub = elmerreadlib.grid3d(bed,'taub',holes,extent)
    plt.imshow(image[:,:,0],extent=[ximage[0],ximage[-1],yimage[0],yimage[-1]],\
                cmap='Greys_r',origin='lower',clim=[0,0.6])
    im1 = plt.imshow(taub*1e3,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',vmin=0,vmax=400,cmap='viridis')
    xmin,xmax = plt.xlim()
    ymin,ymax = plt.ylim()
    plt.plot(np.r_[extent[:,0],extent[0,0]],np.r_[extent[:,1],extent[0,1]],'k',lw=0.5)
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    plt.xticks([])
    plt.yticks([])

    plt.subplot(gs[1,i])
    ax = plt.gca()
    x,y,misfit = elmerreadlib.grid3d(surf_misfit,'misfit',holes,extent)
    plt.imshow(image[:,:,0],extent=[ximage[0],ximage[-1],yimage[0],yimage[-1]],\
                cmap='Greys_r',origin='lower',clim=[0,0.6])
    im2 = plt.imshow(misfit,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',vmin=-10,vmax=10,cmap='RdBu_r')
    plt.xticks([])
    plt.yticks([])
    plt.plot(np.r_[extent[:,0],extent[0,0]],np.r_[extent[:,1],extent[0,1]],'k',lw=0.5)
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    if glacier == 'Helheim':
        bot = 0.83
    elif glacier == 'Kanger':
        bot = 0.73
    path = matplotlib.path.Path([[0.63*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.97*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.97*(xmax-xmin)+xmin,bot*(ymax-ymin)+ymin],
                        [0.63*(xmax-xmin)+xmin,bot*(ymax-ymin)+ymin],
                        [0.63*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin]])
    patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1,zorder=3)
    ax.add_patch(patch)
    MAPD = np.mean(abs(surf_misfit['misfit']))
    if glacier == 'Helheim':
        plt.text(xmax-0.35*(xmax-xmin),y[-1]-3.5e3,'MAR',fontsize=10)
        plt.text(xmax-0.35*(xmax-xmin),y[-1]-6.0e3,'{0:.1f}%'.format(MAPD),fontsize=10)
    elif glacier == 'Kanger':
        plt.text(xmax-0.35*(xmax-xmin),y[-1]-3.2e3,'MAR',fontsize=10)
        plt.text(xmax-0.35*(xmax-xmin),y[-1]-5.7e3,'{0:.1f}%'.format(MAPD),fontsize=10)

    #plt.subplot(gs[1,i])
    #ax = plt.gca()
    #values=plt.hist((surf_misfit['misfit']),range=[-15,15],bins=30)
    #plt.ylim([0,1300])
    #plt.xlim([-15,15])
    #plt.xticks([-10,0,10])
    #if i != 0:
    #    ax.set_yticklabels([])
    #else:
    #    plt.ylabel('Count')
    #if i == 1:
    #    plt.xlabel('Residual (%)',fontsize=10)

if glacier == 'Helheim':
    cbar_ax1 = fig.add_axes([0.86,0.49,0.02,0.45])
    cbar_ax2 = fig.add_axes([0.86,0.01,0.02,0.45])
elif glacier == 'Kanger':
    cbar_ax1 = fig.add_axes([0.86,0.49,0.02,0.43])
    cbar_ax2 = fig.add_axes([0.86,0.01,0.02,0.43])
cb1 = fig.colorbar(im1, cax=cbar_ax1,ticks=np.arange(0,450,100))
cb1.set_label(r'$\tau_b$ (kPa)',fontsize=10,fontname='Arial',color='k')
cb1.ax.yaxis.set_tick_params(color='k')
plt.setp(plt.getp(cb1.ax.axes, 'yticklabels'), color='k')

cb2 = fig.colorbar(im2, cax=cbar_ax2,ticks=np.arange(-8,12,4))
cb2.set_label(r'Residual (%)',fontsize=10,fontname='Arial',color='k')
cb2.ax.yaxis.set_tick_params(color='k')
plt.setp(plt.getp(cb2.ax.axes, 'yticklabels'), color='k')

plt.tight_layout()
if glacier == 'Helheim':
    plt.subplots_adjust(left=0.01, bottom=0.01, right=0.85, top=0.94, wspace=0.0, hspace=0.0)
elif glacier == 'Kanger':
    plt.subplots_adjust(left=0.01, bottom=0.01, right=0.85, top=0.92, wspace=0.0, hspace=0.0)
plt.savefig(DIR+"/figures/Regpars_taub_misfit_"+method+".pdf",DPI=600,transparent=False)
plt.close()

print "Saving as "+DIR+"figures/Regpars_taub_misfit_"+method+".pdf"

bins = np.arange(0,np.ceil(np.max(surf_misfit['vsurfini'])/100)*100,100)
e_bins = np.zeros([len(bins),3])
s_bins = np.zeros([len(bins),3])
e_bins[:] = np.float('NaN')
s_bins[:] = np.float('NaN')
for i in range(1,len(bins)):
    ind = np.where((surf_misfit['vsurfini'] < bins[i]) & (surf_misfit['vsurfini'] >= bins[i-1]))[0]
    e_bins[i-1,:] = np.mean(abs(misfits_all[ind,:]),axis=0)
    s_bins[i-1,:] = np.std(abs(misfits_all[ind,:]),axis=0)

fig = plt.figure(figsize=(3.75,2.5))
plt.ylim([0,10])
plt.plot([0,(bins[-1]+50)/1e3],[3,3],'k')
plt.xticks(np.arange(0,bins[-1]/1e3,1))
plt.xlim([0,(bins[-1]+50)/1e3])
colors = ['g','b','r']
for i in [0,1,2]:
    plt.plot((bins+50)/1e3,e_bins[:,i],c=colors[i],label=r'$\lambda$='+regpars[i])
    plt.fill_between((bins+50)/1e3, e_bins[:,i]-s_bins[:,i], e_bins[:,i]+s_bins[:,i],alpha=0.2,edgecolor=colors[i])
plt.legend(fontsize=10)
plt.ylabel('MAR (%)',fontsize=10)
plt.xlabel(r'$u^{obs}$ bin (km / yr)',fontsize=10)
plt.tight_layout()
plt.subplots_adjust(left=0.13, bottom=0.19, right=0.98, top=0.97, wspace=0.0, hspace=0.0)
plt.savefig(DIR+"/figures/Regpars_MAPE_misfit_"+method+".pdf",DPI=600,transparent=False)
plt.close()

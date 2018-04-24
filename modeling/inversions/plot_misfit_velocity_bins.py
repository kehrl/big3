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
parser.add_argument("-mesh",dest = "mesh",required = True,
        help = "Mesh name.")
parser.add_argument("-subpanel",dest="subpanel",required = False,
        help = "Subpanel label.", default="none")

args, _ = parser.parse_known_args(sys.argv)
glacier = args.glacier
subpanel = args.subpanel
mesh = args.mesh

# Input directory
DIRG = os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/")
DIR_SSA_MT = DIRG+"INV_SSA_ModelT/"
DIR_SSA_CT = DIRG+"INV_SSA_ConstantT/"
DIR_FS_MT = DIRG+"INV_FS_ModelT/"
DIR_FS_CT = DIRG+"INV_FS_ConstantT/"
DIRs = [DIR_FS_MT,DIR_FS_CT,DIR_SSA_MT,DIR_SSA_CT]
temperature_text = ['modelT','constantT','modelT','constantT']
modelnames = ['FS-MT','FS-CT','SSA-MT','SSA-CT']
regpars = ['1e12','1e12','1e13','1e13']

bins = np.arange(0,10000,100)
e_bins = np.zeros([len(bins),4])
s_bins_min = np.zeros([len(bins),4])
s_bins_max = np.zeros([len(bins),4])
e_bins[:,:] = np.float('NaN')
s_bins_min[:,:] = np.float('NaN')
s_bins_max[:,:] = np.float('NaN')

for j in range(0,len(DIRs)):
    DIRM = DIRs[j]+mesh+'_'+temperature_text[j]+'_Lcurve'

    # Get vtu files
    dirs = os.listdir(DIRM+'/mesh2d/inversion_adjoint/')
    for dir in dirs:
        if (dir.startswith('lambda_'+regpars[j])) and not(dir.endswith('.pdf')):
            try:
                vtudata = elmerreadlib.pvtu_file(DIRM+'/mesh2d/inversion_adjoint/'+\
                        dir+'/adjoint_beta_ssa0001.pvtu',['vsurfini','ssavelocity','beta'])
                vname = 'ssavelocity'
            except:
                vtudata = elmerreadlib.pvtu_file(DIRM+'/mesh2d/inversion_adjoint/'+\
                        dir+'/adjoint_beta0001.pvtu',['vsurfini','velocity','beta'])
                vname = 'velocity'
            surf = elmerreadlib.values_in_layer(vtudata,'surface')
            surf_misfit = np.zeros(surf.shape,dtype=np.dtype(surf.dtype.descr+[('misfit',np.float64)]))
            for name in surf.dtype.descr:
                surf_misfit[name[0]] = surf[name[0]]
            surf_misfit['misfit'] = (surf[vname]-surf['vsurfini'])/surf['vsurfini']*100
    
    for i in range(1,len(bins)):
        ind = np.where((surf_misfit['vsurfini'] < bins[i]) & (surf_misfit['vsurfini'] >= bins[i-1]))[0]
        if len(ind) > 0:
            e_bins[i-1,j] = np.mean(abs(surf_misfit['misfit'][ind]),axis=0)
            s_bins_min[i-1,j] = np.percentile(abs(surf_misfit['misfit'][ind]),5,axis=0)
            s_bins_max[i-1,j] = np.percentile(abs(surf_misfit['misfit'][ind]),95,axis=0)

fig = plt.figure(figsize=(2.9,2.0))
ax = plt.gca()
ax.tick_params(labelsize=8)
matplotlib.rc('font',family='Arial')
plt.ylim([0,16])
plt.yticks(np.arange(0,16,5))
plt.xticks(np.arange(0,bins[-1]/1e3,1),fontname='Arial')
plt.xlim([0,(bins[-1]+50)/1e3])
colors = ['b','c','red','orange']
linestyles = ['--','-','--','-']
for i in [3,2,1,0]:
    plt.plot((bins+50)/1e3,e_bins[:,i],linestyles[i],c=colors[i],label=modelnames[i])
    plt.fill_between((bins+50)/1e3, s_bins_min[:,i], s_bins_max[:,i],color=colors[i],alpha=0.2,edgecolor=colors[i])
plt.plot([0,(bins[-1]+50)/1e3],[3,3],'k',label = 'TSX velocity error')
plt.legend(labelspacing=0.25,handlelength=1.5,handletextpad=0.25,columnspacing=0.5,fontsize=8)
if subpanel != 'none':
    if glacier == 'Kanger':
        plt.text(0.6,13.8,subpanel,fontsize=8,fontweight='bold')
    elif glacier == 'Helheim':
        plt.text(0.3,13.8,subpanel,fontsize=8,fontweight='bold')
nonnan = np.where(~(np.isnan(e_bins)))[0]
plt.xlim([0,bins[nonnan[-1]]/1e3])
if glacier == 'Kanger':
    plt.text(0.7,14.7,'(c)',fontsize=8,fontname='Arial',fontweight='bold')
elif glacier == 'Helheim':
    plt.text(0.3,14.7,'(c)',fontsize=8,fontname='Arial',fontweight='bold')
plt.ylabel('Absolute Residual (%)',fontsize=8)
plt.xlabel(r'$u^{obs}$ bin (km yr$^{-1}$)',fontsize=8,fontname='Arial')
plt.tight_layout()
plt.subplots_adjust(left=0.15, bottom=0.2, right=0.98, top=0.98, wspace=0.0, hspace=0.0)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_residual_velocity_bins.pdf"),DPI=400,transparent=False)
plt.close()

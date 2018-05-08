# This code plots L-curves for SSA-MT and FS-MT for Kanger or Helheim, highlighting 
# the chosen regularization paramter.
#
# LMK, UW, 24 April 2018

import numpy as np
import sys
import os
import matplotlib
import matplotlib.pyplot as plt
import elmerreadlib
import shapely.geometry
from pylab import *
import argparse

args = sys.argv

# Get inputs to file
parser = argparse.ArgumentParser()
parser.add_argument("-glacier",dest="glacier",required = True,
        help = "Name of glacier.")

args, _ = parser.parse_known_args(sys.argv)
glacier = args.glacier

if glacier == 'Kanger':
    date = '20120522'
elif glacier == 'Helheim':
    date = '20120316'

# Input directory
DIR = os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/")
DIR_FS = DIR+"INV_FS_ModelT/DEM"+date+"_modelT_Lcurve/mesh2d/inversion_adjoint/"
DIR_SSA = DIR+"INV_SSA_ModelT/DEM"+date+"_modelT_Lcurve/mesh2d/inversion_adjoint/"

regpar_FS, regpar_SSA = [],[]
regpar_str_FS, regpar_str_SSA = [],[]
cost_tot_FS, cost_tot_SSA = [],[]
cost_sur_FS, cost_sur_SSA = [],[]
cost_bed_FS, cost_bed_SSA = [],[]
fid = open(DIR_FS+"summary.dat","r")
lines = fid.readlines()
for line in lines:
    p=line.split()
    if p[0] == 'Lambda':
        pass
    else:
	if p[0].startswith('5') or p[0].startswith('1'):
            regpar_str_FS.append(p[0])
	    regpar_FS.append(float(p[0]))
            cost_tot_FS.append(float(p[2]))
            cost_sur_FS.append(float(p[3]))
            cost_bed_FS.append(float(p[4]))
fid.close()
fid = open(DIR_SSA+"summary.dat","r")
lines = fid.readlines()
for line in lines:
    p=line.split()
    if p[0] == 'Lambda':
        pass
    else:
        if p[0].startswith('5') or p[0].startswith('1'):
            regpar_str_SSA.append(p[0])
            regpar_SSA.append(float(p[0]))
            cost_tot_SSA.append(float(p[2]))
            cost_sur_SSA.append(float(p[3]))
            cost_bed_SSA.append(float(p[4]))
fid.close()

regpar_FS = np.array(regpar_FS)
regpar_SSA = np.array(regpar_SSA)
cost_tot_FS = np.array(cost_tot_FS)
cost_tot_SSA = np.array(cost_tot_SSA)
cost_bed_FS = np.array(cost_bed_FS)
cost_bed_SSA = np.array(cost_bed_SSA)
cost_sur_FS = np.array(cost_sur_FS)
cost_sur_SSA = np.array(cost_sur_SSA)

ind = np.argsort(regpar_FS)
regpar_FS = regpar_FS[ind]
regpar_str_FS = [regpar_str_FS[x] for x in ind]
cost_tot_FS = cost_tot_FS[ind]
cost_bed_FS = cost_bed_FS[ind]
cost_sur_FS = cost_sur_FS[ind]

ind = np.argsort(regpar_SSA)
regpar_SSA = regpar_SSA[ind]
regpar_str_SSA = [regpar_str_SSA[x] for x in ind]
cost_tot_SSA = cost_tot_SSA[ind]
cost_bed_SSA = cost_bed_SSA[ind]
cost_sur_SSA = cost_sur_SSA[ind]

# Get area to get average "misfit"
area = (shapely.geometry.Polygon(np.loadtxt(DIR+'INV_FS_ModelT/DEM'+date+'_modelT_Lcurve/inputs/mesh_extent.dat'))).area
#np.sqrt(cost_sur*2/area)

ymin = np.min([cost_sur_SSA,cost_sur_FS])
ymax = np.max([cost_sur_SSA,cost_sur_FS])
ymin, ymax = ymin-0.04*(ymax-ymin), ymax+0.25*(ymax-ymin)

fig =  plt.figure(figsize=(3.5,2))
gs = matplotlib.gridspec.GridSpec(1,2) 
ax1 = plt.subplot(gs[0])
matplotlib.rc('font',family='Arial')
ax1.tick_params(axis='both',labelsize=8)
strings=["{:.0e}".format(i) for i in regpar_SSA]   
plt.plot(cost_bed_SSA,cost_sur_SSA,'ko--',lw=1.5,markersize=4,markerfacecolor='r')
ax1.set_xlabel(r'                                 $J_{reg}$',fontsize=8,fontname='Arial')
ax1.set_ylabel(r'$J_o$',fontsize=8,fontname='Arial')
plt.xticks(fontsize=8)
ind = np.where(regpar_SSA == 1e13)[0]
ax1.plot(cost_bed_SSA[ind],cost_sur_SSA[ind],'k*',markersize=9,markerfacecolor='y')
ax1.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.1E"))
ax1.set_ylim(ymin,ymax)
if glacier == 'Kanger':
    ax1.set_xlim([-0.05,1.35])
elif glacier == 'Helheim':
    ax1.set_xlim([-0.1,2.4])
xmin,xmax = plt.xlim()
for i in flipud(range(0,len(strings))):
    if regpar_SSA[i] > 1e9 and (strings[i].startswith('1') or strings[i].startswith('5e+14')):
        ax1.text(cost_bed_SSA[i]+0.01*(xmax-xmin),cost_sur_SSA[i]+0.13*(ymax-ymin),strings[i],fontsize=8,fontname='Arial',rotation=45)
if glacier == 'Kanger':
    ax1.text(xmin+0.05*(xmax-xmin),ymax-0.08*(ymax-ymin),'(a)',fontsize=8,fontweight='bold',fontname='Arial')
    ax1.text(xmin+0.20*(xmax-xmin),ymax-0.08*(ymax-ymin),'SSA-MT',fontsize=8)
elif glacier == 'Helheim':
    ax1.text(xmin+0.35*(xmax-xmin),ymax-0.08*(ymax-ymin),'(a)',fontsize=8,fontweight='bold',fontname='Arial')
    ax1.text(xmin+0.5*(xmax-xmin),ymax-0.08*(ymax-ymin),'SSA-MT',fontsize=8)

ax1 = plt.subplot(gs[1])
ax1.tick_params(labelsize=8)
plt.plot(cost_bed_FS,cost_sur_FS,'ko--',lw=1.5,markersize=4,markerfacecolor='b')
strings=["{:.0e}".format(i) for i in regpar_FS]
ind = np.where(regpar_SSA == 1e12)[0]
ax1.plot(cost_bed_FS[ind],cost_sur_FS[ind],'k*',markersize=9,markerfacecolor='y')
ax1.set_yticks([])
ax1.set_ylim(ymin,ymax)
#ax2 = ax1.twinx()
#ax2.set_ylim(ymin,ymax)
yticks_RMSE = np.arange(np.ceil(np.sqrt(ymin*2/area)/10)*10,np.floor(np.sqrt(ymax*2/area)/10)*10+1,10,dtype=int)
yticks_J = (yticks_RMSE**2.0)*area/2.0
#ax2.set_yticks(yticks_J)
#ax2.set_yticklabels(yticks_RMSE,fontsize=8,fontname='Arial')
#ax2.set_ylabel(r'RMSE (m yr$^{-1}$)',fontsize=8,fontname='Arial')
if glacier == 'Kanger':
    ax1.set_xlim([-1.25,36])
elif glacier == 'Helheim':
    ax1.set_xlim([-2.4,69])
xmin,xmax = plt.xlim()
for i in flipud(range(0,len(strings))):
    if (strings[i].startswith('1') or strings[i].startswith('5e+14')):
        ax1.text(cost_bed_FS[i]+0.01*(xmax-xmin),cost_sur_FS[i]+0.13*(ymax-ymin),strings[i],fontsize=8,fontname='Arial',rotation=45)
if glacier == 'Kanger':
    ax1.text(xmin+0.33*(xmax-xmin),ymax-0.08*(ymax-ymin),'(b)',fontsize=8,fontweight='bold',fontname='Arial')
    ax1.text(xmin+0.48*(xmax-xmin),ymax-0.08*(ymax-ymin),'FS-MT',fontsize=8)
elif glacier == 'Helheim':
    ax1.text(xmin+0.05*(xmax-xmin),ymax-0.08*(ymax-ymin),'(b)',fontsize=8,fontweight='bold',fontname='Arial')
    ax1.text(xmin+0.20*(xmax-xmin),ymax-0.08*(ymax-ymin),'FS-MT',fontsize=8)

plt.tight_layout()
plt.subplots_adjust(left=0.2, bottom=0.2, right=0.97, top=0.98, wspace=0.0, hspace=0.0)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_Lcurves.pdf"),DPI=300)
plt.close()

print "Saving as "+os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_Lcurves.pdf")

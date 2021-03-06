# This code plots the data-model misfit J_o through time for the snapshot inversions.
#  
#
# Laura Kehrl, UW, 04/10/2018
import numpy as np
import sys
import os
import matplotlib
import matplotlib.pyplot as plt
import elmerreadlib, geotifflib, datelib
import shapely.geometry
from pylab import *
import argparse

args = sys.argv

# Get inputs to file
parser = argparse.ArgumentParser()
parser.add_argument("-glacier",dest="glacier",required = True,
        help = "Name of glacier.")
parser.add_argument("-dim", dest="dimension",required = False,
	default='3D',help = "3D or Flowline. Default is 3D.")

args, _ = parser.parse_known_args(sys.argv)
dim = args.dimension
glacier = args.glacier

# Input directory
DIRG = os.path.join(os.getenv("MODEL_HOME"),glacier+'/'+dim+'/')
DIR_SSA_MT = DIRG+"INV_SSA_ModelT/"
DIR_SSA_CT = DIRG+"INV_SSA_ConstantT/"
DIR_FS_MT = DIRG+"INV_FS_ModelT/"
DIR_FS_CT = DIRG+"INV_FS_ConstantT/"

SSA_MT_dates, SSA_CT_dates, FS_MT_dates, FS_CT_dates = [], [], [], []
SSA_MT_time, SSA_CT_time, FS_MT_time, FS_CT_time = [], [], [], []
SSA_MT_cost, SSA_CT_cost, FS_MT_cost, FS_CT_cost = [], [], [], []
SSA_MT_rmse, SSA_CT_rmse, FS_MT_rmse, FS_CT_rmse = [], [], [], []

DIRs = [DIR_FS_CT, DIR_FS_MT, DIR_SSA_CT, DIR_SSA_MT]
dates = [FS_CT_dates, FS_MT_dates, SSA_CT_dates, SSA_MT_dates]
times = [FS_CT_time, FS_MT_time, SSA_CT_time, SSA_MT_time]
costs = [FS_CT_cost, FS_MT_cost, SSA_CT_cost, SSA_MT_cost]
rmses = [FS_CT_rmse, FS_MT_rmse, SSA_CT_rmse, SSA_MT_rmse]

# Correct regpar for FS and SSA
regpars = ['1e12','1e12','1e13','1e13']

for i in range(0,len(DIRs)):
    dirs = os.listdir(DIRs[i])
    for dir in dirs:
        area = (shapely.geometry.Polygon(np.loadtxt(DIRs[i]+dir+'/inputs/mesh_extent.dat'))).area    
        if dir.endswith('T'):
            fid = open(DIRs[i]+dir+'/mesh2d/inversion_adjoint/summary.dat','r')
            lines = fid.readlines()
            for line in lines:
                p = line.split()
                if p[0] == regpars[i]:
                    dates[i].append(dir[3:11])
                    times[i].append(datelib.date_to_fracyear(int(dir[3:7]),int(dir[7:9]),int(dir[9:11])))
                    costs[i].append(float(p[3]))
                    rmses[i].append(np.sqrt(2*float(p[3])/area))
	    fid.close()

fig = plt.figure(figsize=(6.5,1.5))
matplotlib.rc('font',family='Arial')
ax1 = plt.gca()
plt.plot(SSA_CT_time,SSA_CT_cost,'ks',markersize=4,markerfacecolor='orange',label='SSA-CT')
plt.plot(SSA_MT_time,SSA_MT_cost,'ko',markersize=4,markerfacecolor='red',label='SSA-MT')
plt.plot(FS_CT_time,FS_CT_cost,'ks',markersize=4,markerfacecolor='cyan',label='FS-CT')
plt.plot(FS_MT_time,FS_MT_cost,'ko',markersize=4,markerfacecolor='blue',label='FS-MT')
ax1.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1E'))
ax1.tick_params(axis='both',labelsize=8)
plt.ylabel(r'$J_o$',fontsize=8,fontname='Arial')
ax1.set_xlim([2001,2016])
ax2 = ax1.twinx()
mn, mx = ax1.get_ylim()
if mn < 0:
  mn = 0
ax2.set_ylim(mn,mx)
ax1.set_ylim(mn,mx)
if glacier == 'Kanger':
    rmse_tick = 25
elif glacier == 'Helheim':
    rmse_tick = 50
yticks_RMSE = np.arange(np.ceil(np.sqrt(mn*2/area)/rmse_tick)*rmse_tick,np.floor(np.sqrt(mx*2/area)/rmse_tick)*rmse_tick+1,rmse_tick,dtype=int)
yticks_J = (yticks_RMSE**2.0)*area/2.0
ax2.set_yticks(yticks_J)
ax2.set_yticklabels(yticks_RMSE,fontsize=8,fontname='Arial')
ax2.set_ylabel(r'RMSE (m yr$^{-1}$)',fontsize=8,fontname='Arial')
ax1.text(2001.2,mx-(mx-mn)*0.11,'(d)',fontsize=8,fontweight='bold',fontname='Arial')
ax1.legend(labelspacing=0.25,handlelength=1.5,handletextpad=0.25,columnspacing=0.5,fontsize=8)
plt.tight_layout()
plt.subplots_adjust(left=0.11, bottom=0.13, right=0.925, top=0.98, wspace=0.0, hspace=0.0)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_misfit_timeseries.pdf"),format="PDF",dpi=300)
plt.close()

# Check to make sure we only have one RMSE value for each snapshot inversion for the following stats
for i in range(0,len(times[0])):
    if (FS_CT_time[i] != FS_MT_time[i]) and (FS_CT_time[i] != SSA_CT_time[i]) and (FS_CT_time[i] != SSA_MT_time[i]):
        print "WARNING: Something is wrong with the RMSE value for "+str(FS_CT_time[i])

print "Average RMSE for "+glacier
print "FS-CT:  {0:.0f} m/yr".format(np.mean(FS_CT_rmse))
print "FS-MT:  {0:.0f} m/yr".format(np.mean(FS_MT_rmse))
print "SSA-CT: {0:.0f} m/yr".format(np.mean(SSA_CT_rmse))
print "SSA-MT: {0:.0f} m/yr".format(np.mean(SSA_MT_rmse))

# Find how many of the snapshot inversions have FS-CT as the closest match to the observations
n = 0
not_dates = []
for i in range(0,len(times[0])):
    cost = FS_CT_cost[i]
    if (cost < FS_MT_cost[i]) and (cost < SSA_CT_cost[i]) and (cost < SSA_MT_cost[i]):
        n = n + 1
    else:
        not_dates.append(FS_CT_dates[i])
print str(n)+"/"+str(len(times[0]))+" snapshot inversions have FS-CT as the closest match to the observations"
print "The exceptions are "+str(not_dates)


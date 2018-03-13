# This code plots an L-curve for the three dimensional model or flowline model. 
#  
#
# LMK, UW, 5/30/2015

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
parser.add_argument("-mesh", dest="meshname", required = True,
        help = "Name of mesh.")
parser.add_argument("-dim", dest="dimension",required = False,
	default='3D',help = "3D or Flowline.")
parser.add_argument("-method", dest="method",required = True,
        help = "adjoint or robin.")
args, _ = parser.parse_known_args(sys.argv)
RES = args.meshname
ver = args.dimension
method = args.method
glacier = args.glacier

# Input directory
DIR = os.path.join(os.getenv("MODEL_HOME"),glacier+"/"+ver+"/"+RES+"/")
DIRR = DIR+"mesh2d/inversion_"+method+"/"

if not(os.path.isdir(DIRR)):
    print DIRR
    sys.exit("That is not a model results directory.")

regpar = []
regpar_str = []
nsim = []
cost_tot = []
cost_sur = []
cost_bed = []
fid = open(DIRR+"summary.dat","r")
lines = fid.readlines()
for line in lines:
    p=line.split()
    if p[0] == 'Lambda':
        pass
    else:
        regpar_str.append(p[0])
        regpar.append(float(p[0]))
        nsim.append(float(p[1]))
        cost_tot.append(float(p[2]))
        cost_sur.append(float(p[3]))
        cost_bed.append(float(p[4]))

regpar = np.array(regpar)
nsim = np.array(nsim)
cost_tot = np.array(cost_tot)
cost_bed = np.array(cost_bed)
cost_sur = np.array(cost_sur)

ind1 = np.argsort(regpar)
ind = ind1#[np.where((regpar >= 1e11) & (regpar <= 1e14))[0]]
regpar = regpar[ind]
regpar_str = [regpar_str[x] for x in ind]
nsim = nsim[ind]
cost_tot = cost_tot[ind]
cost_bed = cost_bed[ind]
cost_sur = cost_sur[ind]

# Get vtu files
dirs = os.listdir(DIRR)
misfit = np.zeros([len(regpar),])
for i in range(0,len(regpar)):
    for dir in dirs:
        if (dir.startswith('lambda_'+regpar_str[i])) and not(dir.endswith('.pdf')):
            try:
                vtudata = elmerreadlib.pvtu_file(DIRR+dir+'/adjoint_beta_ssa0001.pvtu',['vsurfini','ssavelocity'])
                surf = elmerreadlib.values_in_layer(vtudata,'surface')
                misfit[i] = np.sqrt(np.mean((surf['vsurfini 1']-surf['ssavelocity 1'])**2+\
                        (surf['vsurfini 2']-surf['ssavelocity 2'])**2))
            except:
                vtudata = elmerreadlib.pvtu_file(DIRR+dir+'/adjoint_beta0001.pvtu',['vsurfini','velocity'])
                surf = elmerreadlib.values_in_layer(vtudata,'surface')
                misfit[i] = np.sqrt(np.mean((surf['vsurfini 1']-surf['velocity 1'])**2+\
                        (surf['vsurfini 2']-surf['velocity 2'])**2))

# Get area to get average "misfit"
area = (shapely.geometry.Polygon(np.loadtxt(DIR+'inputs/mesh_extent.dat'))).area
#np.sqrt(cost_sur*2/area)


fig =  plt.figure(figsize=(3.5,3))
ax1 = plt.gca()
matplotlib.rc('font',family='sans-serif',size=10)
strings=["{:.0e}".format(i) for i in regpar]   
plt.plot(cost_bed,cost_sur,'ko--',linewidth=1.5)
ax1.set_xlabel('$J_{reg}$',fontsize=10)
ax1.set_ylabel(r'$J_o$',fontsize=10)
plt.xticks(fontsize=10)
ax1.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.1E"))
ax1.set_ylim([np.min(cost_sur)-(np.max(cost_sur)-np.min(cost_sur))/10,np.max(cost_sur)+(np.max(cost_sur)-np.min(cost_sur))/4])
ymin,ymax = ax1.get_ylim()
ax2 = ax1.twinx()
mn, mx = ax1.get_ylim()
ax2.set_ylim(mn,mx)
#ax2.set_ylim(np.sqrt(mn*2/area), np.sqrt(mx*2/area))
if (mx-mn) < 25:
    yticks_RMSE = np.arange(np.ceil(np.sqrt(mn*2/area)/10)*10,np.floor(np.sqrt(mx*2/area)/10)*10+1,10,dtype=int)
else:
    yticks_RMSE = np.arange(np.ceil(np.sqrt(mn*2/area)/5)*5,np.floor(np.sqrt(mx*2/area)/5)*5+1,5,dtype=int)
yticks_J = (yticks_RMSE**2.0)*area/2.0
ax2.set_yticks(yticks_J)
ax2.set_yticklabels(yticks_RMSE)
ax2.set_ylabel('RMSE (m/yr)')
for i in range(0,len(strings)):
    if strings[i].startswith('1') or strings[i].startswith('5'):
        ax1.text(cost_bed[i]+0.0,cost_sur[i]+0.13*(ymax-ymin),strings[i],fontsize=10,rotation=45)
if np.max(cost_bed) > 40:
    plt.xlim([-2,np.max(cost_bed)+(np.max(cost_bed)-np.min(cost_bed))/3])
elif np.max(cost_bed) > 12:
    plt.xlim([-1,np.max(cost_bed)+(np.max(cost_bed)-np.min(cost_bed))/3])
else:
    plt.xlim([-0.02,np.max(cost_bed)+(np.max(cost_bed)-np.min(cost_bed))/3])

plt.tight_layout()
plt.subplots_adjust(left=0.26, bottom=0.15, right=0.84, top=0.98, wspace=0.0, hspace=0.0)
#plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Lcurve_"+glacier+"_"+RES+".pdf"),format='PDF')
plt.savefig(DIR+"/figures/Lcurve_"+method+".pdf")
plt.close()

print "Saving as "+DIR+"figures/Lcurve_"+method+".pdf"

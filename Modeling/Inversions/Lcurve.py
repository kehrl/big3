# This code plots an L-curve for the three dimensional model or flowline model. 
# Inputs is : lcurve.py Flowline DEM20120624 
#
# LMK, UW, 5/30/2015

import sys
import os
import matplotlib.pyplot as plt
from pylab import *
import argparse

args = sys.argv

# Get inputs to file
parser = argparse.ArgumentParser()
parser.add_argument("-mesh", dest="meshname", required = True,
        help = "Name of mesh.")
parser.add_argument("-dim", dest="dimension",required = True,
		help = "3D or Flowline.")
args, _ = parser.parse_known_args(sys.argv)
RES = args.meshname
ver = args.dimension

# Input directory
DIR = os.path.join(os.getenv("MODEL_HOME"),"Helheim/Results/"+ver+"/")
DIRR = DIR+RES+"/Inversion/"

if not(os.path.isdir(DIRR)):
  sys.exit("That is not a model results directory.")

regpar =[]
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

ind = np.argsort(cost_bed)
regpar = regpar[ind]
nsim = nsim[ind]
cost_tot = cost_tot[ind]
cost_bed = cost_bed[ind]
cost_sur = cost_sur[ind]


#plt.xlim([0,np.max(cost_sur)])  
plt.figure(figsize=(5,5))
strings=["{:.0e}".format(i) for i in regpar]   
plt.plot(cost_bed,cost_sur,'ko--',linewidth=1.5)
plt.xlabel('$J_{reg}$',fontsize=14)
plt.ylabel('$J_o$',fontsize=14)
i=0; plt.text(cost_bed[i]+cost_bed[i]*0.5,cost_sur[i],strings[i],fontsize=12) 
for i in np.arange(1,len(nsim)-4,1):
  if (i==0) or (strings[i][0]=='1'):
    plt.text(cost_bed[i],cost_sur[i]+cost_sur[i]*0.05,strings[i],fontsize=12) 
plt.text(cost_bed[-1]-cost_bed[-1]*0.5,cost_sur[-1]-cost_sur[-1]*0.15,strings[-1],fontsize=12) 
#plt.yticks([1e10,5e10,1e11],fontsize=10)
plt.xticks(fontsize=10)
plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
plt.gca().tick_params('both', length=10, width=1.5, which='major',labelsize=10)
plt.gca().tick_params('both', length=5, width=1.5, which='minor',labelsize=10)
plt.gca().yaxis.set_minor_formatter(FormatStrFormatter("%.0E"))
plt.gca().xaxis.set_major_formatter(FormatStrFormatter("%.0E"))
plt.gca().yaxis.set_major_formatter(FormatStrFormatter("%.0E"))
plt.subplots_adjust(left=0.18, bottom=None, right=None, top=None, wspace=0.0, hspace=0.0)
plt.gca()
#plt.ylim([1e10,1e11+1e10])
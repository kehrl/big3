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
parser.add_argument("-method", dest="method",required = True,
                help = "adjoint or robin.")
args, _ = parser.parse_known_args(sys.argv)
RES = args.meshname
ver = args.dimension
method = args.method

# Input directory
DIR = os.path.join(os.getenv("MODEL_HOME"),"Helheim/"+ver+"/"+RES+"/")
DIRR = DIR+"mesh2d/inversion_"+method+"/"

if not(os.path.isdir(DIRR)):
  print DIRR
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
plt.figure(figsize=(4.5,4))
matplotlib.rc('font',family='sans-serif',size=10)
strings=["{:.0e}".format(i) for i in regpar]   
plt.plot(cost_bed,cost_sur,'ko--',linewidth=1.5)
plt.xlabel('Model norm',fontsize=10)
plt.ylabel(r'Misfit',fontsize=10)
for i in range(0,len(strings)):
  plt.text(cost_bed[i]+0.15,cost_sur[i]+0.25e12,strings[i],fontsize=10,rotation=45) 
plt.xticks(fontsize=10)
#plt.gca().set_xscale('log')
#plt.gca().set_yscale('log')
#plt.gca().tick_params('both', length=10, width=1.5, which='major',labelsize=10)
#plt.gca().tick_params('both', length=5, width=1.5, which='minor',labelsize=10)
#plt.gca().yaxis.set_minor_formatter(FormatStrFormatter("%.0E"))
#plt.gca().xaxis.set_major_formatter(FormatStrFormatter("%.0E"))
plt.gca().yaxis.set_major_formatter(FormatStrFormatter("%.1E"))
plt.ylim([np.min(cost_sur)-(np.max(cost_sur)-np.min(cost_sur))/10,np.max(cost_sur)+(np.max(cost_sur)-np.min(cost_sur))/4])
plt.xlim([-0.4,np.max(cost_bed)+(np.max(cost_bed)-np.min(cost_bed))/4])
plt.tight_layout()
plt.subplots_adjust(left=0.2, bottom=0.11, right=0.98, top=0.97, wspace=0.0, hspace=0.0)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Lcurve_"+RES+".pdf"),format='PDF')
#plt.close()

# This code plots an L-curve for the three dimensional model. 
#
# LMK, UW, 5/30/2015

import sys
import os
import matplotlib.pyplot as plt
from pylab import *

args = sys.argv

DIR = os.path.join(os.getenv("HOME"),"Models/Helheim/Results/3D/")
RES = args[1]
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
      regpar.append(p[0])
      nsim.append(float(p[1]))
      cost_tot.append(float(p[2]))
      cost_sur.append(float(p[3]))
      cost_bed.append(float(p[4]))
     
  plt.plot(cost_bed,cost_sur,'ko--')
  plt.xlabel('$J_{reg}$',fontsize=24)
  plt.ylabel('$J_o$',fontsize=24)
  for i in range(0,len(nsim)):
    plt.text(cost_bed[i]+0.05,cost_sur[i]+0.25e9,regpar[i],fontsize=20) 
  plt.xticks([0,1,2,3,4],fontsize=16)
  plt.yticks(fontsize=16)
  #subplots_adjust(left=None, bottom=0.15, right=0.98, top=0.95, wspace=-0.1, hspace=None)

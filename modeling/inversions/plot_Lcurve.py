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
parser.add_argument("-method", dest="method",required = False,
        help = "adjoint or robin.",default='adjoint')
parser.add_argument("-highlight",dest="regpar_highlight", required = False,
        help = "Choose a regularization parameter to highlight with a red circle.",default="none")
parser.add_argument("-subpanel",dest="subpanel",required = False,
	help = "Subpanel label.",default="none")
parser.add_argument("-modelname",dest="modelname",required = False,
        help = "Model name.",default="none")

args, _ = parser.parse_known_args(sys.argv)
RES = args.meshname
ver = args.dimension
method = args.method
glacier = args.glacier
regpar_highlight = args.regpar_highlight
subpanel = args.subpanel
modelname = args.modelname

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
	if p[0].startswith('5') or p[0].startswith('1'):
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

ind = np.argsort(regpar)
regpar = regpar[ind]
regpar_str = [regpar_str[x] for x in ind]
nsim = nsim[ind]
cost_tot = cost_tot[ind]
cost_bed = cost_bed[ind]
cost_sur = cost_sur[ind]

# Get vtu files
dirs = os.listdir(DIRR)
misfit = np.zeros([len(regpar),])
cutoffs = np.arange(0,1001,250)
MAPE = np.zeros([len(regpar),len(cutoffs)])
for i in range(0,len(regpar)):
    for dir in dirs:
        if (dir.startswith('lambda_'+regpar_str[i])) and not(dir.endswith('.pdf')):
            try:
                vtudata = elmerreadlib.pvtu_file(DIRR+dir+'/adjoint_beta_ssa0001.pvtu',['vsurfini','ssavelocity'])
                surf = elmerreadlib.values_in_layer(vtudata,'surface')
                misfit[i] = np.sqrt(np.mean((surf['vsurfini 1']-surf['ssavelocity 1'])**2+\
                        (surf['vsurfini 2']-surf['ssavelocity 2'])**2))
		APE = abs((surf['vsurfini']-surf['ssavelocity'])/surf['vsurfini'])*100
	    except:
                vtudata = elmerreadlib.pvtu_file(DIRR+dir+'/adjoint_beta0001.pvtu',['vsurfini','velocity'])
                surf = elmerreadlib.values_in_layer(vtudata,'surface')
                misfit[i] = np.sqrt(np.mean((surf['vsurfini 1']-surf['velocity 1'])**2+\
                        (surf['vsurfini 2']-surf['velocity 2'])**2))
                APE = abs((surf['vsurfini']-surf['velocity'])/surf['vsurfini'])*100
            for j in range(0,len(cutoffs)):
	        MAPE[i,j] = np.mean(APE[[surf['vsurfini'] > cutoffs[j]]])

# Get area to get average "misfit"
area = (shapely.geometry.Polygon(np.loadtxt(DIR+'inputs/mesh_extent.dat'))).area
#np.sqrt(cost_sur*2/area)


fig =  plt.figure(figsize=(3.5,3))
ax1 = plt.gca()
matplotlib.rc('font',family='sans-serif',size=10)
strings=["{:.0e}".format(i) for i in regpar]   
plt.plot(cost_bed,cost_sur,'ko--',linewidth=1.5)
ax1.set_xlabel(r'$J_{reg}$',fontsize=10)
ax1.set_ylabel(r'$J_o$',fontsize=10)
plt.xticks(fontsize=10)
if not(regpar_highlight == 'none'):
    ind = np.where(regpar == float(regpar_highlight))[0]
    ax1.plot(cost_bed[ind],cost_sur[ind],'ko',markerfacecolor='r')
ax1.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.1E"))
ax1.set_ylim([np.min(cost_sur)-(np.max(cost_sur)-np.min(cost_sur))/25,np.max(cost_sur)+(np.max(cost_sur)-np.min(cost_sur))/5])
ymin,ymax = ax1.get_ylim()
ax2 = ax1.twinx()
mn, mx = ax1.get_ylim()
ax2.set_ylim(mn,mx)
#ax2.set_ylim(np.sqrt(mn*2/area), np.sqrt(mx*2/area))
yticks_RMSE = np.arange(np.ceil(np.sqrt(mn*2/area)/10)*10,np.floor(np.sqrt(mx*2/area)/10)*10+1,10,dtype=int)
yticks_J = (yticks_RMSE**2.0)*area/2.0
ax2.set_yticks(yticks_J)
ax2.set_yticklabels(yticks_RMSE)
ax2.set_ylabel('RMSE (m/yr)')
if np.max(cost_bed) > 40:
    plt.xlim([-2,np.max(cost_bed)+(np.max(cost_bed)-np.min(cost_bed))/4])
elif np.max(cost_bed) > 12:
    plt.xlim([-1,np.max(cost_bed)+(np.max(cost_bed)-np.min(cost_bed))/4])
elif np.max(cost_bed) > 1:
    plt.xlim([-0.15,np.max(cost_bed)+(np.max(cost_bed)-np.min(cost_bed))/4])
else:
    plt.xlim([-0.02,np.max(cost_bed)+(np.max(cost_bed)-np.min(cost_bed))/4])
xmin,xmax = plt.xlim()
for i in flipud(range(0,len(strings))):
    if strings[i].startswith('1') or strings[i].startswith('5'):
        if i != len(strings)-1 and (abs(cost_bed[i+1]-cost_bed[i]) < 0.05*(xmax-xmin) and \
                abs(cost_sur[i+1]-cost_sur[i]) < 0.05*(ymax-ymin)):
            if regpar[i] > 1e10:
                ax1.text(cost_bed[i]+0.04*(xmax-xmin),cost_sur[i]+0.10*(ymax-ymin),strings[i],fontsize=8,rotation=45)
        else:
            ax1.text(cost_bed[i]+0.01*(xmax-xmin),cost_sur[i]+0.105*(ymax-ymin),strings[i],fontsize=8,rotation=45)
if subpanel != 'none':
    ax1.text(xmin+0.33*(xmax-xmin),ymax-0.07*(ymax-ymin),subpanel,fontsize=10,fontweight='bold')
if modelname != 'none':
    ax1.text(xmin+0.4*(xmax-xmin),ymax-0.07*(ymax-ymin),modelname,fontsize=10)

plt.tight_layout()
plt.subplots_adjust(left=0.26, bottom=0.15, right=0.84, top=0.98, wspace=0.0, hspace=0.0)
plt.savefig(DIR+"/figures/Lcurve_"+method+".pdf")
plt.close()

fig =  plt.figure(figsize=(3.5,3))
ax1 = plt.gca()
matplotlib.rc('font',family='sans-serif',size=10)
plt.plot([-0.5,len(regpar)],[3,3],'k-')
for j in range(0,len(cutoffs)):
    plt.plot(MAPE[:,j],'o--',linewidth=1.5,label=r'$>$'+str(cutoffs[j]))
ax1.set_xticks(range(0,len(regpar)))
ax1.set_xlim([-0.2,len(regpar)-0.8])
ax1.set_xticklabels(strings,rotation=45,fontsize=10,ha='right')
ax1.set_xlabel('$\lambda$',fontsize=10)
ax1.set_ylabel(r'MAR (%)',fontsize=10)
ax1.set_ylim([1,15])
plt.legend(ncol=2,labelspacing=0.25,handlelength=1.5,handletextpad=0.25,columnspacing=0.5)

plt.tight_layout()
plt.subplots_adjust(left=0.15, bottom=0.24, right=0.98, top=0.98, wspace=0.0, hspace=0.0)
plt.savefig(DIR+"/figures/Lcurve_"+method+"_MAPE.pdf")
plt.close()

print "Saving as "+DIR+"figures/Lcurve_"+method+".pdf"

# This code plots the cost function through for individual inversions time.
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
parser.add_argument("-method", dest="method",required = False,
        default='adjoint', help = "Adjoint or robin. Default is adjoint.")
parser.add_argument("-regpar", dest="regpar",required = False,
        default = '1e13', help = "Regularization parameter.")

args, _ = parser.parse_known_args(sys.argv)
dim = args.dimension
method = args.method
glacier = args.glacier
regpar = args.regpar

# Input directory
DIRG = os.path.join(os.getenv("MODEL_HOME"),glacier+'/'+dim+'/')

dirs = os.listdir(DIRG)

model_dates = []
model_times = []
model_J_surf = []
constant_dates = []
constant_times = []
constant_J_surf = []
firsttime = True
for dir in dirs:
    if dir.endswith('modelT'):
        # Get area to get average "misfit"
        area = (shapely.geometry.Polygon(np.loadtxt(DIRG+dir+'/inputs/mesh_extent.dat'))).area
        try:
            fid = open(DIRG+dir+'/mesh2d/inversion_'+method+'/summary.dat','r')
            lines = fid.readlines()
            for line in lines[1:]:
                p = line.split()
                if p[0] == regpar:
                    model_dates.append(dir[3:11])
                    model_times.append(datelib.date_to_fracyear(int(dir[3:7]),int(dir[7:9]),int(dir[9:11])))
                    model_J_surf.append(float(p[3]))
	except:
	    pass
    if dir.endswith('constantT'):
        # Get area to get average "misfit"
        area = (shapely.geometry.Polygon(np.loadtxt(DIRG+dir+'/inputs/mesh_extent.dat'))).area
        try:
            fid = open(DIRG+dir+'/mesh2d/inversion_'+method+'/summary.dat','r')
            lines = fid.readlines()
            for line in lines[1:]:
                p = line.split()
                if p[0] == regpar:
                    constant_dates.append(dir[3:11])
                    constant_times.append(datelib.date_to_fracyear(int(dir[3:7]),int(dir[7:9]),int(dir[9:11])))
                    constant_J_surf.append(float(p[3]))
        except:
            pass

fig = plt.figure(figsize=(6,3))
ax1 = plt.gca()
matplotlib.rc('font',family='sans-serif',size=10)
plt.plot(model_times,model_J_surf,'ko',label='SSA-MT')
plt.plot(constant_times,constant_J_surf,'ro',label='SSA-CT')
ax1.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1E'))
plt.ylabel(r'$J_o$',fontsize=10)
ax2 = ax1.twinx()
mn, mx = ax1.get_ylim()
if mn < 0:
  mn = 0
ax2.set_ylim(mn,mx)
ax1.set_ylim(mn,mx)
yticks_RMSE = np.arange(np.ceil(np.sqrt(mn*2/area)/25)*25,np.floor(np.sqrt(mx*2/area)/25)*25+1,25,dtype=int)
yticks_J = (yticks_RMSE**2.0)*area/2.0
ax2.set_yticks(yticks_J)
ax2.set_yticklabels(yticks_RMSE)
ax2.set_ylabel('RMSE (m/yr)')
ax1.legend()
plt.tight_layout()
plt.subplots_adjust(left=0.16, bottom=0.1, right=0.9, top=0.98, wspace=0.0, hspace=0.0)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_misfit_timeseries_"+regpar+".pdf"),format="PDF",dpi=600)
plt.close()

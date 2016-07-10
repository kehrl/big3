import sys
import os
import matplotlib.pyplot as plt
from pylab import *
import argparse
import elmerreadlib,glaclib,colorlib

args = sys.argv

# Get inputs to file
parser = argparse.ArgumentParser()
parser.add_argument("-dir", dest="dir", required = True,
        help = "Directory where all the model results are located.")
parser.add_argument("-glacier", dest="glacier",required = True,
                help = "Glacier name.")

args, _ = parser.parse_known_args(sys.argv)
DIR = args.dir
glacier = args.glacier

# Input directory
method = 'robin'
regpar = '1e10'
DIR = os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/"+DIR+"/")

# Load flowline
x,y,zb,dists = glaclib.load_flowline(glacier,dx=100)

DEMDIRs = os.listdir(DIR)

# Get distinct colors for plotting
NUM_COLORS = len(DEMDIRs)
cm = plt.get_cmap('gist_rainbow')
fig1 = plt.figure(1)
ax1 = plt.gca()
ax1.set_color_cycle([cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
fig2 = plt.figure(2)
ax2 = plt.gca()
ax2.set_color_cycle([cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
fig3 = plt.figure(3)
ax3 = plt.gca()
ax3.set_color_cycle([cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
fig4 = plt.figure(4)
ax4 = plt.gca()
ax4.set_color_cycle([cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
fig5 = plt.figure(5)
ax5 = plt.gca()
ax5.set_color_cycle([cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])

n = 0
bed = {}
for DEMDIR in DEMDIRs:
  invdir = DIR+DEMDIR+"/inversion_"+method+"/"
  regpardirs = os.listdir(invdir)
  for regpardir in regpardirs:
    if (regpardir.startswith('lambda_'+regpar)) and not(regpardir.endswith('pdf')):
      finaldir = invdir+regpardir
  files = os.listdir(finaldir)
  finaliter = 0
  for file in files:
    if file.endswith('.pvtu'):
      if int(file[-9:-5]) > finaliter:
         finaliter = int(file[-9:-5])
  data = elmerreadlib.pvtu_file(finaldir+'/robin_beta'+'{:04d}'.format(finaliter)+'.pvtu',['velocity','beta'])
  surf = elmerreadlib.values_in_layer(data,'surf')
  bed = elmerreadlib.values_in_layer(data,'bed')
  flowsurf = elmerreadlib.grid_to_flowline_surface(surf,x,y)
  flowbed = elmerreadlib.grid_to_flowline_surface(bed,x,y)

  plt.figure(1)
  ax1.plot(dists/1e3,flowsurf['velocity'],label=DEMDIR[3:])
  
  plt.figure(2)
  ax2.plot(dists/1e3,flowbed['beta']**2,label=DEMDIR[3:])
  
  plt.figure(3)
  ax3.plot(dists/1e3,flowbed['taub']*1e3,label=DEMDIR[3:])

  plt.figure(4)
  alpha = np.diff(flowsurf['z'])/np.diff(dists)
  taud = -1*917.0*9.8*(flowsurf['z'][1:]-flowbed['z'][1:])*alpha
  ax4.plot(dists[1:]/1e3,taud/1e3,label=DEMDIR[3:])

  plt.figure(5)
  ax5.plot(dists[1:]/1e3,flowbed['taub'][1:]*1e3/(taud/1e3),label=DEMDIR[3:])

  plt.figure(6)
  plt.scatter(bed['x'],bed['y'],c=bed['taub']*1e3,vmin=0,vmax=500,lw=0)
  plt.colorbar()
  plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/inversion_"+DEMDIR[3:]+".pdf"),format='PDF')
  plt.close()

plt.figure(1)
plt.legend(loc=2,ncol=2)
plt.ylabel('Surface velocity (m/yr)')
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/inversion_surface_velocity.pdf"),format='PDF')

plt.figure(2)
plt.legend(ncol=2)
plt.ylabel('Beta^2')
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/inversion_betas.pdf"),format='PDF')

plt.figure(3)
plt.legend(ncol=2)
plt.ylabel('tau_b (kPa)')
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/inversion_taub.pdf"),format='PDF')

plt.figure(4)
plt.legend(ncol=2)
plt.ylabel('tau_d (kPa)')
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/inversion_taud.pdf"),format='PDF')

plt.figure(5)
plt.legend(ncol=2)
plt.ylabel('tau_b/tau_d')
plt.ylim([0,2])
nonnan = np.where(~(np.isnan(flowbed['taub'])))[0]
plt.plot([dists[nonnan[0]]/1e3,dists[nonnan[-1]]/1e3],[1,1],'k',lw=2)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/inversion_taudb.pdf"),format='PDF')

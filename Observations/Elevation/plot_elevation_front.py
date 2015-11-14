# This file pulls all the available TDX, WV DEMs and then plots them in a big figure.
#
#
# LMK, UW, 9/1/2015

import os
import sys
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"Util/Modules/"))
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"BigThreeGlaciers/Tools/"))
import glacier_flowline, elevation, icefronts, fracyear, flotation, dem_shading
import numpy as np
import matplotlib.pyplot as plt
import matplotlib


# Get arguments
args = sys.argv
glacier = args[1][:] # Options: Kanger, Helheim

# Plot extent
if glacier == 'Helheim':
  xmin = 304000.0
  xmax = 314000.0
  ymin = -2582500.0
  ymax = -2572500.0
elif glacier == 'Kanger':
  xmin = 485000.0
  xmax = 498000.0
  ymin = -2298000.0
  ymax = -2285000.0

# Load worldview DEMs
xwv,ywv,zwv,timewv = elevation.worldview_grid(glacier,xmin,xmax,ymin,ymax,years='all',resolution=32,verticaldatum='geoid')

# Calculate height above flotation at CreSIS radar picks
xf,yf,zabovefloat = flotation.extent(xwv,ywv,zwv,timewv,glacier,rho_i=917.0,rho_sw=1020.0,bedsource='cresis',verticaldatum='geoid')


# Length of plot
N = len(timewv)
years = np.arange(np.floor(np.min(timewv)),np.floor(np.max(timewv))+1)

# Find dimensions for plot
ncol = 0
for year in years:
  ndem = len(np.where(np.floor(timewv)==year)[0])
  if ndem > ncol:
    ncol = ndem
nrow = len(years)

# Make plot
plt.figure(figsize=(ncol*2.8,nrow*3))
cmap = matplotlib.cm.get_cmap('gray_r',4)
gs = matplotlib.gridspec.GridSpec(nrow,ncol)
gs.update(right=0.94,left=0.02,wspace=0.04,hspace=0.04)
for j in range(0,nrow):
  year = years[j]
  for i in range(0,len(np.where(np.floor(timewv)==year)[0])):
    plt.subplot(gs[j,i])
    wvind = np.where(np.floor(timewv)==year)[0][i]
    date = fracyear.fracyear_to_date(timewv[wvind])
    shadeddem = dem_shading.set_shade(zwv[:,:,wvind],0,220)
    plt.imshow(shadeddem,extent=[np.min(xwv),np.max(xwv),np.min(ywv),np.max(ywv)],clim=[0,200],origin='lower')
    plt.scatter(xf,yf,c=zabovefloat[:,wvind],s=3.0**2.0,cmap='gist_gray_r',edgecolors='none',vmin=-10,vmax=10)
    plt.xlim([np.min(xwv),np.max(xwv)])
    plt.ylim([np.min(ywv),np.max(ywv)])
    plt.xticks([])
    plt.yticks([])
    plt.text(xmin+500,ymax-1.25e3,str(date[0])+'-'+str(date[1])+'-'+str(int(np.floor(date[2]))),backgroundcolor='w',fontsize=8)
    
plt.savefig("/Users/kehrl/Bigtmp/"+glacier+"_dems.png",format='PNG',dpi=150)
plt.close()
# The Morlighem bed DEM does not extend to the terminus, so we need to come up 
# with a way to extend it. I'm first going to try to use the 2001 CreSIS line 
# and some kind of interpolating function in the transverse direction to come 
# up with a bed profile.
#
# LMK, UW, 5/18/2015

import os
import sys
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
sys.path.append(os.path.join(os.getenv("HOME"),"Code/BigThreeGlaciers/Tools"))
import geotiff
import numpy as np
import bed, dist, elevation
import matplotlib
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate
import elmer_mesh as mesh


##########
# Inputs #
##########

glacier = 'Helheim'

if glacier == 'Helheim':
  xmin = 290000.0
  xmax = 320000.0
  ymin = -2584000.0
  ymax = -2561000.0
elif glacier == 'Kanger':
  xmin = 480000.0
  xmax = 506000.0
  ymin = -2302000.0
  ymax = -2280000.0

# Directories
DIRX=os.path.join(os.getenv("HOME"),"Data/ShapeFiles/Glaciers/3D/"+glacier+"/")

# Mesh extent
#exterior = mesh.shp_to_xy(DIRX+"glacier_extent_normal")
#exterior_nofront = mesh.shp_to_xy(DIRX+"glacier_extent_nofront")

# Load CreSIS radar picks
cresis_2001 = bed.cresis('2001',glacier,verticaldatum='geoid')
cresis_all = bed.cresis('all',glacier,verticaldatum='geoid')
ind = np.where(cresis_all[:,2] < 0)[0]

# Load CreSIS grid
xCre,yCre,zCre = bed.cresis_grid(glacier,verticaldatum='geoid')

# Morlighem bed
xMor,yMor,zMor = bed.morlighem_grid(xmin,xmax,ymin,ymax,'geoid')
cresis_all_zMor = bed.morlighem_pts(cresis_all[ind,0],cresis_all[ind,1],glacier,'geoid')
cresis_2001_zMor = bed.morlighem_pts(cresis_2001[:,0],cresis_2001[:,1],glacier,'geoid')

# Make some plots
plt.figure(1)
plt.title('Morlighem grid')
plt.imshow(zMor,extent=[np.min(xMor),np.max(xMor),np.min(yMor),np.max(yMor)],origin='lower',cmap='RdBu_r',vmin=-1000,vmax=1000)
plt.scatter(cresis_all[ind,0],cresis_all[ind,1],c=cresis_all[ind,2],cmap='RdBu_r',vmin=-1000,vmax=1000)
plt.scatter(cresis_2001[:,0],cresis_2001[:,1],c=cresis_2001[:,2],cmap='RdBu_r',vmin=-1000,vmax=1000)
plt.colorbar()
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])

#plt.figure()
#plt.title('Cresis grid')
#plt.imshow(zCre,extent=[np.min(xCre),np.max(xCre),np.min(yCre),np.max(yCre)],origin='lower',cmap='RdBu_r',vmin=-1000,vmax=1000)
#plt.scatter(cresis_all[ind,0],cresis_all[ind,1],c=cresis_all[ind,2],cmap='RdBu_r',vmin=-1000,vmax=1000)
#plt.scatter(cresis_2001[:,0],cresis_2001[:,1],c=cresis_2001[:,2],cmap='RdBu_r',vmin=-1000,vmax=1000)
#plt.colorbar()
#plt.xlim([xmin,xmax])
#plt.ylim([ymin,ymax])

#plt.figure()
#plt.title('Cresis - Morlighem')
#plt.imshow(zMor,extent=[np.min(xMor),np.max(xMor),np.min(yMor),np.max(yMor)],origin='lower',cmap='RdBu_r',vmin=-1000,vmax=1000)
#cmap = plt.cm.jet
#bounds = np.linspace(-300,300,11)
#norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
#plt.scatter(cresis_all[ind,0],cresis_all[ind,1],c=cresis_all[ind,2]-cresis_all_zMor,cmap=cmap,norm=norm,edgecolors='none')
#plt.scatter(cresis_2001[:,0],cresis_2001[:,1],c=cresis_2001[:,2]-cresis_2001_zMor,cmap=cmap,norm=norm,edgecolors='none')
#plt.xlim([xmin,xmax])
#plt.ylim([ymin,ymax])
#plt.colorbar()

##############################
# Look at tranverse profile #
##############################

coloroptions = ['b','r','g','y']
labels = [['A','B'],['C','D'],['E','F'],['G','H']]
if glacier == 'Helheim':
  transects = [[309222,-2574650,308537,-2580150],[307964,-2575731,307695,-2579225],
               [306915,-2575288,307152,-2580378],[305387,-2574680,305798,-2580490]]
  addons =    [[309260,-2574010,308517,-2580680],[308023,-2574268,307665,-2580538],
               [306915,-2574059,307152,-2580378],[305299,-2573920,305828,-2580810]]
if glacier == 'Kanger':
  transects = [[493249,-2290191,490480,-2295650],[492543,-2290179,490218,-2294330],
               [491340,-2289179,489240,-2293783]]         
  addons =    [[493502,-2289560,490324,-2295990],[492731,-2288736,489727,-2295370],
               [491293,-2289272,488834,-2294751]]

for j in range(0,len(transects)):
  # Find radar pick indices
  ind1 = np.argmin((cresis_all[:,0]-transects[j][0])**2+(cresis_all[:,1]-transects[j][1])**2)
  ind2 = np.argmin((cresis_all[:,0]-transects[j][2])**2+(cresis_all[:,1]-transects[j][3])**2)
  if ind1 < ind2:
    cut_cresis = cresis_all[ind1:ind2,:]
    sortind = (np.argsort(cut_cresis[:,1]))
    cut_cresis = cut_cresis[sortind,:]
    other = np.row_stack([cresis_all[0:ind1,:],cresis_all[ind2:-1,:]])
  else: 
    cut_cresis = cresis_all[ind2:ind1,:]
    sortind = (np.argsort(cut_cresis[:,1]))
    cut_cresis = cut_cresis[sortind,:]
    other = np.row_stack([cresis_all[0:ind2,:],cresis_all[ind1:-1,:]])

  # Distances along transect
  xtemp = np.append(cut_cresis[:,0],[addons[j][0],addons[j][2]])
  ytemp = np.append(cut_cresis[:,1],[addons[j][1],addons[j][3]])
  sortind = ((np.argsort(ytemp)))
  xtemp = xtemp[sortind]
  ytemp = ytemp[sortind]
  dtemp = dist.transect(xtemp,ytemp)
  dists = dtemp[1:-1]

  # Get morlighem bed
  cut_dMor = np.linspace(0,dtemp[-1],dtemp[-1]/150.0)
  cut_xMor= np.interp(cut_dMor,dtemp,xtemp)
  cut_yMor = np.interp(cut_dMor,dtemp,ytemp)
  cut_zMor = bed.morlighem_pts(cut_xMor,cut_yMor,glacier,'geoid')
  cut_cresis_grid = bed.cresis_grid_pts(cut_xMor,cut_yMor,glacier,'geoid')
  del xtemp,ytemp,dtemp

  # Get 2001 Cresis point
  mindist = 100.0
  for i in range(0,len(cut_cresis)):
    pdist = np.min(np.sqrt((cut_cresis[i,0]-cresis_2001[:,0])**2+(cut_cresis[i,1]-cresis_2001[:,1])**2))
    if pdist < mindist:
      mindist = pdist
      pind = np.argmin(np.sqrt((cut_cresis[i,0]-cresis_2001[:,0])**2+(cut_cresis[i,1]-cresis_2001[:,1])**2))
      cut_2001 = [dists[i],cresis_2001[pind,2]]
  #del pind,mindist,pdist

  # Get other Cresis points
  cut_cresis_other = []
  for i in range(0,len(other)):
    mindist = 50.0
    pdist = np.min(np.sqrt((cut_cresis[:,0]-other[i,0])**2+(cut_cresis[:,1]-other[i,1])**2))
    if pdist < mindist:
      pind = np.argmin(np.sqrt((cut_cresis[:,0]-other[i,0])**2+(cut_cresis[:,1]-other[i,1])**2))
      cut_cresis_other.append([dists[pind],other[i,2]])
  cut_cresis_other = np.array(cut_cresis_other)  
  
  # Get worldview surface elevations
  zs,time = elevation.worldview_at_pts(cut_xMor,cut_yMor,glacier,years='all',verticaldatum='geoid')
  
  plt.figure(1)
  plt.plot(cut_xMor,cut_yMor,color=coloroptions[j],linewidth=1.5)
  plt.text(cut_xMor[0],cut_yMor[0],labels[j][0],backgroundcolor='w')
  plt.text(cut_xMor[-1],cut_yMor[-1],labels[j][1],backgroundcolor='w')
       
  plt.figure(j+2)
  plt.title(labels[j])
  rho_i = 910.0; rho_sw = 1020.0;
  for i in range(0,len(time)):
    iceind = np.where(zs[i,:] < 40.0)
    zs[i,iceind] = 'NaN'
    zb = rho_i*zs[i,:]/(rho_i-rho_sw)
    if i == 0 and glacier == 'Helheim':
      plt.plot(cut_dMor,zb,color='0.7',label='flotation bed',linewidth=3)
    else:
      plt.plot(cut_dMor,zb,color='0.7')
  plt.plot(cut_dMor,cut_zMor,'k',label='Morlighem',linewidth=1.5)
  plt.plot(cut_dMor,cut_cresis_grid,'k--',label='Cresis grid',linewidth=1.5)
  plt.plot(dists,cut_cresis[:,2],'o',color=coloroptions[j],label='Cresis transect')
  plt.plot(cut_cresis_other[:,0],cut_cresis_other[:,1],'cs',label='Cresis other')
  try:
    plt.plot(cut_2001[0],cut_2001[1],'ms',label='Cresis 2001')
  except:
    print "No 2001 points"
  plt.legend(loc=2)
  plt.ylim([np.min(cut_cresis[:,2])-100,100])


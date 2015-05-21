# The Morlighem bed DEM does not extend to the terminus, so we need to come up 
# with a way to extend it. I'm first going to try to use the 2001 CreSIS line 
# and some kind of interpolating function in the transverse direction to come 
# up with a bed profile.
#
# LMK, UW, 5/18/2015

import os
import sys
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Helheim/Tools"))
import geotiff
import numpy as np
import helheim_bed
import scipy.interpolate as interpolate
import elmer_mesh as mesh


##########
# Inputs #
##########

# Directories
DIRX=os.path.join(os.getenv("HOME"),"Data/ShapeFiles/Glaciers/3D/Helheim/")

# Mesh extent
exterior = np.array(mesh.shp_to_xy(DIRX+"glacier_extent_normal"))
exterior_nofront = np.array(mesh.shp_to_xy(DIRX+"glacier_extent_nofront"))

# Load CreSIS radar picks
bed2001 = helheim_bed.cresis('2001','geoid')

#####################
# Valley side walls #
#####################

# Find valley side walls so that we can fit a line through them. We can then easily look
# at the transverse bed profile at different locations along the fjord.
ind = np.where(np.array(exterior_nofront[0,:]) > 306000.0)
exterior_nofront = exterior_nofront[:,ind[0]]

# Find points along north wall and fit a line through them.
ind = np.where(exterior_nofront[1,:] > -2577000.0)
north_wall = exterior_nofront[:,ind[0]]
PN = np.polyfit(north_wall[0],north_wall[1],1)

# Find points along south wall and fit a line through them.
ind = np.where(exterior_nofront[1,:] < -2577000.0)
south_wall = exterior_nofront[:,ind[0]]
PS = np.polyfit(south_wall[0],south_wall[1],1)

PM = np.mean([PS[0],PN[0]])

# Points along north wall
xN = np.linspace(np.min(north_wall[0]),np.max(north_wall[0]),(np.max(north_wall[0])-np.min(north_wall[0]))/100)
yN = PN[0]*xN+PN[1]

# Corresponding points along south wall
xS = (yN+(1/PM)*xN-PS[1])/(PS[0]+1/PM)
yS = PS[0]*xS+PS[1]

###########################
# Transverse bed profiles #
###########################

xt = np.zeros([50,len(xN)])
yt = np.zeros([50,len(xN)])
zt = np.zeros([50,len(xN)])
for i in range(0,len(xN)):
  xt[:,i] = np.linspace(xS[i],xN[i],50)
  yt[:,i] = np.linspace(yS[i],yN[i],50)
  zt[:,i] = helheim_bed.morlighem_pts(xt[:,i],yt[:,i],'geoid')
zt[zt==-9999]='NaN'


plt.subplot(121)
plt.imshow(bedMor,extent=[np.min(xMor),np.max(xMor),np.min(yMor),np.max(yMor)],origin='lower',cmap='RdBu_r')
plt.clim([-1.5e3,1.5e3])
plt.ylim([-2.59e6,-2.565e6])
plt.xlim([295000,315000])
for i in range(10,30):
  plt.subplot(121)
  plt.plot(xt[:,i],yt[:,i],'k')
  plt.plot(xtC[:,i],ytC[:,i],'r')
  plt.plot(exterior[0],exterior[1],'k',linewidth=2)
  plt.subplot(122)
  plt.cla()
  #if i > 0:
  #  plt.plot(yt[:,i-1],zt[:,i-1],'k')
  plt.plot(yt[:,i],zt[:,i],'k',label="Morlighem")
  plt.plot(ytC[:,i],ztC[:,i],'r',label="CreSIS grid")

  mindist = 100
  for j in range(0,len(xt[:,i])):
    dist = np.min(np.sqrt((bed2001[:,0]-xt[j,i])**2+(bed2001[:,1]-yt[j,i])**2))
    if dist < mindist:
      mindist = dist
      ind = j
  plt.plot(bed2001[ind,1],bed2001[ind,2],'bo',label="Cresis 2001")
  plt.legend()
  plt.waitforbuttonpress()


##########################################
# Interpolate transverse profile forward #
##########################################

# Choose a transverse profile from previous section
ind = 21
xint = xt[:,ind]
yint = yt[:,ind]
zint = zt[:,ind]
 
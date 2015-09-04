import os
import sys
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules/"))
sys.path.append(os.path.join(os.getenv("HOME"),"Code/BigThreeGlaciers/Tools/"))
import elevation, icefronts, elmer_mesh, coords, glacier_extent
import numpy as np
import scipy.interpolate
from matplotlib import path

def height(zb,rho_i=917.0,rho_sw=1024.0):

  # Calculate flotation height for bed elevations zb

  # Set up output
  float = np.zeros(len(zb))
  float[:] = 'NaN' # no flotation height for bed elevations above sea level
  
  # Find places where bed is below sea level
  ind = np.where(zb < 0)[0]
  
  # Find flotation height
  float[ind] = (1-rho_sw/rho_i)*zb[ind]

  return float

def extent(xs,ys,zs,ztime,glacier,rho_i=917.0,rho_sw=1020.0,bedsource='cresis',verticaldatum='geoid'):

  # This function takes an array of surface elevations and finds the height above flotation,
  # according to the Morlighem bed DEM or CreSIS radar picks. Right now it's only set up to 
  # work for cresis radar picks.
  
  # Inputs:
  # xs,ys,zs,time from elevation.worldview_grid
  # glacier: name of glacier
  # rho_i: ice density
  # rho_sw: seawater density
  # bedsource: cresis or morlighem
  
  # Outputs:
  # xf,yf,zabovefloat: height above (or below flotation) at points xf,yf
  
  # Find bed elevations
  if bedsource=='cresis':
    cresis = bed.cresis('all',glacier,verticaldatum)
    cresis2001 = bed.cresis('2001',glacier,verticaldatum)
    cresis = np.row_stack([cresis[cresis[:,2]<-200.0,:],cresis2001[cresis2001[:,2]<-200.0,:]])
  else:
    print "need to work on morlighem"
  
  # Find flotation height
  zfloat = height(cresis[:,2],rho_i,rho_sw)
  xf = cresis[:,0]
  yf = cresis[:,1]
  
  # Calculate height above flotation through time
  N = len(ztime)
  zabovefloat=np.zeros([len(zfloat),N])
  zabovefloat[:,:] = 'NaN'
  for i in range(0,N): 
  
    # Get glacier extent so we're only checking if the glacier is floating where there is glacier ice
    xglacier,yglacier = glacier_extent.date(glacier,ztime[i])
    
    # Check what points are located on glacier ice
    box = path.Path(np.column_stack([xglacier,yglacier]))
  
    # Find indices in the grid that fall within the fluxbox
    inside = box.contains_points(np.column_stack([cresis[:,0:2]]))
    
    # Current surface elevation
    dem = scipy.interpolate.RegularGridInterpolator([xs,ys],zs[:,:,i],bounds_error = False,method='linear',fill_value=float('nan')) 
    zs_on_zb = dem(cresis[:,0:2])
    
    # Find height above flotation
    zabovefloat[inside,i] = zs_on_zb[inside]-zfloat[inside]
    
  return xf,yf,zabovefloat
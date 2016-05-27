# This module calculates flotation height and height above flotation according to different bed products.
#
# height(zb,rho,rho_sw): flotation height for bed elevations zb 
# shelfbase(zs,rho_i,rho_sw): ice shelf base height for freeboard zs
# extent(xs,ys,zs,ztime,glacier,rho_i,rho_sw,bedsource,verticaldatum): calculates height above flotation 
#	through time according to bed elevations from cresis and surface elevations from worldview grids
#
# LMK, UW, 9/2/2015

import os
import glaclib, bedlib
import numpy as np
import scipy.interpolate
from matplotlib import path

def shelfbase(zs,rho_i=917.0,rho_sw=1020.0):
  
  '''
  D = shelfbase(zs,rho_i=917.0,rho_sw=1020.0)
  
  Calculate the flotation depth given a surface profile
  
  Inputs:
  zs: array of surface elevations relative to sea level
  rho_i,rho_sw: ice and seawater densities
  
  Outputs:
  D: array of flotation depths using the input surface elevations as freeboard
  
  '''
  
  D = -rho_i/(rho_sw-rho_i)*zs

  return D

def icebottom(zb,zs,rho_i=917.0,rho_sw=1020.0):

  '''
  bottom = icebottom(zb,zs,rho_i=917.0,rho_sw=1020.0)
  
  Take a surface and bed topography (relative to sea level) and calculate the bottom of the
  glacier. In locations where the ice is floating, the ice bottom elevation is the flotation 
  depth. In locations where the glacier is grounded, the ice bottom elevation is the bed
  elevation.
  
  Inputs:
  zb,zs: array of bed and surface elevations
  rho_i,rho_sw: ice and seawater densities
  
  Outputs:
  bottom: array of ice bottom elevations
  '''

  # Find flotation depth
  floatdepth = shelfbase(zs,rho_i,rho_sw)
  
  # Find locations where the flotation depth is less deep (i.e., greater) than the bed elevation (these are locations
  # where the ice is floating)
  ind = np.where(floatdepth > zb)[0]
  
  bottom = np.array(zb)
  bottom[ind] = floatdepth[ind]
  
  # Also set nan bed depths to flotation thickness
  ind = np.where(np.isnan(zb))[0]
  bottom[ind] = floatdepth[ind]

  return bottom
  
#########################################################################################
  
def height(zb,rho_i=917.0,rho_sw=1020.0):

  '''
  float = height(zb,rho_i=917.0,rho_sw=1020.0)
  
  Calculate flotation height for bed elevations zb
  
  Inputs:
  zb: array of bed elevations relative to sea level
  rho_i, rho_sw: ice and seawater densities
  
  Outputs:
  float: flotation height (above this elevation the ice is grounded; below this elevation
  the ice is floating)
  '''

  # Set up output
  try:
    floatheight = np.zeros(len(zb))
    floatheight[:] = float('NaN') # no flotation height for bed elevations above sea level
  
    # Find places where bed is below sea level
    nonnan = np.where(~(np.isnan(zb)))[0]
    ind = np.where(zb[nonnan] < 0)[0]
  
    # Find flotation height
    floatheight[nonnan[ind]] = (1-rho_sw/rho_i)*zb[nonnan[ind]]
    
  except:
    if zb < 0:
      floatheight = (1-rho_sw/rho_i)*zb
    else:
      floatheight = float('nan')

  return floatheight

def extent(xs,ys,zs,ztime,glacier,rho_i=917.0,rho_sw=1020.0,bedsource='cresis',verticaldatum='geoid'):

  '''
  This function finds the height above (or below) flotation for an array of surface elevations,
  according to the Morlighem bed DEM or CreSIS radar picks. Right now it's only set up to 
  work for cresis radar picks.
  
  Inputs:
  xs,ys,zs,time from zslib.dem_grid
  glacier: name of glacier
  rho_i: ice density
  rho_sw: seawater density
  bedsource: cresis or morlighem
  
  Outputs:
  xf,yf,zabovefloat: height above (or below flotation) at points xf,yf
  '''
  
  # Find bed elevations
  if bedsource=='cresis':
    cresis = bedlib.cresis('all',glacier,verticaldatum,cleanup=True)
    cresis2001 = bedlib.cresis('2001',glacier,verticaldatum)
    cresis = np.row_stack([cresis[cresis[:,2]<-50.0,:],cresis2001[cresis2001[:,2]<-50.0,:]])
  else:
    print "need to work on morlighem"
  
  # Find flotation height
  zfloat = height(cresis[:,2],rho_i=rho_i,rho_sw=rho_sw)
  xf = cresis[:,0]
  yf = cresis[:,1]
  
  # Calculate height above flotation through time
  try:
    N = len(ztime)
  except: 
    N = 1
  zabovefloat=np.zeros([len(zfloat),N])
  zabovefloat[:,:] = float('NaN')
  for i in range(0,N): 
  
    # Get glacier extent so we're only checking if the glacier is floating where there is glacier ice
    if N == 1: # no iteration if only one DEM
      extent = glaclib.load_extent(glacier,ztime)
      dem = scipy.interpolate.RegularGridInterpolator([ys,xs],zs[:,:],bounds_error = False,method='linear',fill_value=float('nan')) 
    else:
      extent = glaclib.load_extent(glacier,ztime[i])
      dem = scipy.interpolate.RegularGridInterpolator([ys,xs],zs[:,:,i],bounds_error = False,method='linear',fill_value=float('nan')) 
    
    # Check what points are located on glacier ice
    box = path.Path(extent[:,0:2])
  
    # Find indices in the grid that fall within the fluxbox
    inside = box.contains_points(np.column_stack([cresis[:,0:2]]))
    
    # Current surface elevation
    zs_on_zb = dem((cresis[:,1],cresis[:,0]))
    
    # Find height above flotation
    zabovefloat[inside,i] = zs_on_zb[inside]-zfloat[inside]
    
  return xf,yf,zabovefloat
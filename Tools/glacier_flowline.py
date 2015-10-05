# This is simply a script to keep the chosen flowline for each glacier 
# consistent between scripts.
#
# LMK, UW, 7/2/2015

import sys
import os
import numpy as np
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"BigThreeGlaciers/Tools"))
import elmer_mesh as mesh
import icefronts, bed
import dist
import scipy.signal as signal

def load(glacier,shapefilename='center_flowline',filt_len='none',verticaldatum='geoid',bedmodel='aniso',bedsmoothing=4):

  '''
  x,y,zb_filt,dists = load(glacier,shapefilename='center_flowline')
  
  Load glacier flowline. This script is mostly to keep everything consistent (distance along flowline,
  chosen bed profile,etc.
  
  Inputs:
  glacier: glacier name
  shapefilename: shapefile to use for the flowline
  filt_len: filter length (in meters) for the bed profile
  verticaldatum: geoid or ellipsoid
  
  Outputs:
  x,y: x,y coordinates of flowline
  zb_filt: bed profile for flowline
  '''

  if glacier == "Helheim":
    file1 = 'helheim_'
    file_flowline_in = os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/Glaciers/Flowlines/Helheim/"+file1+shapefilename)
  elif glacier == 'Kanger':
    file1 = 'kanger_'
    file_flowline_in = os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/Glaciers/Flowlines/Kanger/"+file1+shapefilename)
  else:
    sys.exit("Unknown glacier.") 
  
  flowline = mesh.shp_to_xy(file_flowline_in)
  d = dist.transect(flowline[:,0],flowline[:,1])

  # Set uniform spacing between nodes along flowline
  dists_old = np.linspace(0,np.max(d),np.max(d)/50)
  x = np.interp(dists_old,d,flowline[:,0])
  y = np.interp(dists_old,d,flowline[:,1])
  
  # Get new distances
  dists = dist.transect(x,y)
  
  # Get terminus positions so that we can set distance along flowline 
  # relative to the average terminus position
  terminus_val, terminus_time = icefronts.distance_along_flowline(x,y,dists,glacier,type='icefront')

  # Get average terminus position
  time1 = 2008.0
  time2 = 2016.0
  indt = np.where((terminus_time > time1) & (terminus_time < time2))
  terminus_time = terminus_time[indt[0]]
  terminus_val = terminus_val[indt[0]]

  # Average terminus position
  terminus = np.mean(terminus_val)

  # Set dists relative to average terminus position
  dists = dists-terminus

  # Find bed elevation
  if glacier == 'Helheim':
    zb = bed.smith_at_pts(x,y,glacier,model=bedmodel,smoothing=bedsmoothing,verticaldatum=verticaldatum)
  elif glacier == 'Kanger':
    zb = bed.morlighem_pts(x,y,glacier,verticaldatum)
    zb_cresis1 = bed.cresis('2009a',glacier,verticaldatum)
    zb_cresis2 = bed.cresis('2008',glacier,verticaldatum)
    ind1 = np.where((x > 490970) & (x < 491300))[0]
    ind2 = np.where(x > 491300)[0]
    zb[ind1] = np.interp(x[ind1],zb_cresis1[:,0],zb_cresis1[:,2])
    zb[ind2] = np.interp(x[ind2],zb_cresis2[:,0],zb_cresis2[:,2])

  if filt_len != 'none':
    cutoff=(1/filt_len)/(1/(np.diff(dists[1:3])*2))
    b,a=signal.butter(4,cutoff,btype='low')
    zb_filt = signal.filtfilt(b,a,zb)
  else:
    zb_filt = zb

  return x,y,zb_filt,dists
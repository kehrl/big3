# This is simply a script to keep the chosen flowline for each glacier 
# consistent between scripts.
#
# LMK, UW, 7/2/2015

import sys
import os
import numpy as np
sys.path.append(os.path.join(os.getenv("HOME"),"Code/BigThreeGlaciers/Tools"))
import elmer_mesh as mesh
import icefronts, bed
import dist
import scipy.signal as signal

def load(glacier):

  if glacier == "Helheim":
    file_flowline_in = os.path.join(os.getenv("HOME"),"Data/ShapeFiles/Glaciers/Flowlines/Helheim/helheim_center_flowline")
  elif glacier == 'Kanger':
    file_flowline_in = os.path.join(os.getenv("HOME"),"Data/ShapeFiles/Glaciers/Flowlines/Kanger/kanger_center_flowline")
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
  time1 = 2009.0
  time2 = 2015.0
  indt = np.where((terminus_time > time1) & (terminus_time < time2))
  terminus_time = terminus_time[indt[0]]
  terminus_val = terminus_val[indt[0]]

  # Average terminus position
  terminus = np.mean(terminus_val)

  # Set dists relative to average terminus position
  dists = dists-terminus

  # Find bed elevation
  
  zb = bed.morlighem_pts(x,y,glacier,'geoid')
  if glacier == 'Helheim':
    zb_cresis = bed.cresis('2001',glacier,'geoid')
    ind = np.where(x > 309000)[0]
    zb[ind] = np.interp(x[ind],zb_cresis[:,0],zb_cresis[:,2])
  elif glacier == 'Kanger':
    zb_cresis1 = bed.cresis('2009a',glacier,'geoid')
    zb_cresis2 = bed.cresis('2008',glacier,'geoid')
    ind1 = np.where((x > 490970) & (x < 491300))[0]
    ind2 = np.where(x > 491300)[0]
    zb[ind1] = np.interp(x[ind1],zb_cresis1[:,0],zb_cresis1[:,2])
    zb[ind2] = np.interp(x[ind2],zb_cresis2[:,0],zb_cresis2[:,2])

  filt_len = 1000.0
  cutoff=(1/filt_len)/(1/(np.diff(dists[1:3])*2))
  b,a=signal.butter(4,cutoff,btype='low')
  zb_filt = signal.filtfilt(b,a,zb)

  return x,y,zb_filt,dists

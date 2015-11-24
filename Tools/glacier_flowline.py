# This is simply a script to keep the chosen flowline for each glacier 
# consistent between scripts.
#
# LMK, UW, 7/2/2015

import sys
import os
import numpy as np
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"BigThreeGlaciers/Tools"))
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"Util/Modules"))
import elmer_mesh as mesh
import icefronts, bed
import dist
import scipy.signal as signal

def load(glacier,shapefilename='center_flowline',filt_len=2.0e3,verticaldatum='geoid',bedmodel='aniso',bedsmoothing=4):

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
  dists_old = np.linspace(0,np.max(d),np.max(d)/20)
  x = np.interp(dists_old,d,flowline[:,0])
  y = np.interp(dists_old,d,flowline[:,1])
  
  # Get new distances
  dists = dist.transect(x,y)

  # Get average terminus position
  time1 = 2008.0
  time2 = 2016.0
  
  # Get terminus positions so that we can set distance along flowline 
  # relative to the average terminus position
  terminus_val, terminus_time = icefronts.distance_along_flowline(x,y,dists,glacier,type='icefront',time1=time1,time2=time2)

  # Average terminus position
  terminus = np.mean(terminus_val)

  # Set dists relative to average terminus position
  dists = dists-terminus

  # Find bed elevation
  if glacier == 'Helheim':
    zb = bed.smith_at_pts(x,y,glacier,model=bedmodel,smoothing=bedsmoothing,verticaldatum=verticaldatum)
  elif glacier == 'Kanger':
    cresis = bed.cresis('all',glacier,verticaldatum=verticaldatum)
    cutdist = 100.
    dcresis = []
    zcresis = []
    tcresis = []
    for i in range(0,len(cresis[:,0])):
      mindist = np.min(np.sqrt((cresis[i,0]-x)**2+(cresis[i,1]-y)**2))
      if mindist < cutdist:
        minind = np.argmin(np.sqrt((cresis[i,0]-x)**2+(cresis[i,1]-y)**2))
        dcresis.append(dists[minind])
        zcresis.append(cresis[i,2])
        tcresis.append(cresis[i,4])

    ind = np.argsort(dcresis)
    dcresis = np.array(dcresis)[ind]
    zcresis = np.array(zcresis)[ind]
    zb = np.interp(dists,dcresis,zcresis)

  if filt_len != 'none':
    ind = np.where(~(np.isnan(zb)))[0]
    zb_filt = np.zeros_like(zb)
    zb_filt[:] = float('nan')
    cutoff=(1/filt_len)/(1/(np.diff(dists[1:3])*2))
    b,a=signal.butter(4,cutoff,btype='low')
    zb_filt[ind] = signal.filtfilt(b,a,zb[ind])
  else:
    zb_filt = zb
    
  if (glacier == 'Kanger') and (shapefilename=='flowline_flightline'):
    ind = np.where(dists > 3000.)[0]
    zb_filt[ind] = float('nan')
  elif (glacier == 'Helheim') and (shapefilename=='flowline_flightline'):
    ind = np.where(dists > 3900.)[0]
    zb_filt[ind] = float('nan')
    

  return x,y,zb_filt,dists
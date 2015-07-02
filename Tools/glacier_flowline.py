# This is simply a script to keep the chosen flowline for each glacier 
# consistent between scripts.
#
# LMK, UW, 7/2/2015

import sys
import os
import numpy as np
sys.path.append(os.path.join(os.getenv("HOME"),"Code/BigThreeGlaciers/Tools"))
import elmer_mesh as mesh
import dist


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
  dists = np.linspace(0,np.max(d),np.max(d)/100)
  x = np.interp(dists,d[:,0],flowline[:,0])
  y = np.interp(dists,d[:,0],flowline[:,1])
  del d,flowline

  return x,y,dists
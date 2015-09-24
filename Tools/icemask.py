# This module produces an ice mask where ice/water is zero and land is one. 
#
# LMK, UW, 9/15/2015

import os
import sys
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
import elmer_mesh
import matplotlib.path as path
import numpy as np

def load(glacier,xmin,xmax,ymin,ymax,dx,ice=0):

  # Directory for holes for mask
  DIR = os.path.join(os.getenv("HOME"),'Data/Shapefiles/IceMasks/'+glacier+'/')
  files = os.listdir(DIR)
  
  # Set up grid
  nx = np.ceil((xmax-xmin)/dx)+1
  x = np.linspace(xmin,xmin+dx*nx,nx)
  ny = np.ceil((ymax-ymin)/dx)+1
  y = np.linspace(ymin,ymin+dx*ny,ny)
  
  if ice == 1:
    mask = np.ones([ny,nx])
  else: 
    mask = np.zeros([ny,nx])
  xgrid,ygrid = np.meshgrid(x,y)
  
  xflatten = xgrid.flatten()
  yflatten = ygrid.flatten()
  points = np.column_stack([xflatten,yflatten])
  maskflatten = mask.flatten()
  
  for file in files:
    if file.endswith('.shp'):
      hole = elmer_mesh.shp_to_xy(DIR+file)
      holepath = path.Path(hole[:,0:2])
      cond = holepath.contains_points(points)
      ind = np.where(cond==True)[0]
      if ice == 1:
        maskflatten[ind]=0
      else: 
        maskflatten[ind]=1
          
  mask = np.reshape(maskflatten,[ny,nx])    

  return x,y,mask
'''
This module deals with our ice masking products. You can produce a ice mask grid or 
load all the points that define where land is located. 

x,y,mask = load_grid(glacier,xmin,xmax,ymin,ymax,dx,ice=0)
pts = load_points(glacier)

LMK, UW, 9/15/2015
'''

import os
import sys
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"Util/Modules"))
import elmer_mesh
import matplotlib.path as path
import numpy as np

def load_grid(glacier,xmin,xmax,ymin,ymax,dx,ice=0):

  '''
  Load ice mask grid.
  
  x,y,mask = load_grid(glacier,xmin,xmax,ymin,ymax,dx,ice=0)
  
  Inputs:
  glacier: glacier name
  xmin,xmax,ymin,ymax: output dimensions for grid
  dx: grid spacing
  ice: set ice grid points to 0 or 1
  
  Outputs:
  x,y: grid points
  mask: mask grid of 0,1's 
  '''

  # Directory for holes for mask
  DIR = os.path.join(os.getenv("DATA_HOME"),'Shapefiles/IceMasks/'+glacier+'/')

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

def load_points(glacier):

  '''
  Load all points that define where land is located.
  
  pts = load_points(glacier)
  
  Inputs:
  glacier: glacier name
  
  Outputs:
  pts: 2d array of points that define land extent

  '''

  DIR = os.path.join(os.getenv("DATA_HOME"),'Shapefiles/IceMasks/'+glacier+'/')

  files = os.listdir(DIR)
  
  n=0
  for file in files:
    if file.endswith('.shp'):
      hole = elmer_mesh.shp_to_xy(DIR+file)
      if n==0:
        pts = hole[:,0:2]  
      else:
        pts = np.row_stack([pts,hole[:,0:2]])
      n=1
      
  return pts
'''
This module deals with our ice masking products. You can produce a ice mask grid or 
load all the points that define where land is located. 

x,y,mask = load_grid(glacier,xmin,xmax,ymin,ymax,dx,ice=0)
pts = load_points(glacier)

LMK, UW, 9/15/2015
'''

import os
import meshlib, icefrontlib
import matplotlib.path as path
import numpy as np

def load_grid(glacier,xmin,xmax,ymin,ymax,dx,ice=0,icefront_time='none',icefront_type='all'):

  '''
  Load ice mask grid.
  
  x,y,mask = load_grid(glacier,xmin,xmax,ymin,ymax,dx,ice=0,icefront_time='none')
  
  Inputs:
  glacier: glacier name
  xmin,xmax,ymin,ymax: output dimensions for grid
  dx: grid spacing
  ice: set ice grid points to 0 or 1
  icefront_time: sometimes we will want the fjord to be part of the mask, so 
  	this option also masks the fjord according to the ice front for the given time
  
  
  Outputs:
  x,y: grid points
  mask: mask grid of 0,1's 
  '''

  # Directory for holes for mask
  DIR = os.path.join(os.getenv("DATA_HOME"),'ShapeFiles/IceMasks/'+glacier+'/')
  files = os.listdir(DIR)
  
  # Set up grid
  nx = np.ceil((xmax-xmin)/dx)+1
  x = np.linspace(xmin,xmin+dx*(nx-1),nx)
  ny = np.ceil((ymax-ymin)/dx)+1
  y = np.linspace(ymin,ymin+dx*(ny-1),ny)
  
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
    if file.endswith('.shp') and 'hole' in file:
      hole = meshlib.shp_to_xy(DIR+file)
      holepath = path.Path(hole[:,0:2])
      cond = holepath.contains_points(points)
      ind = np.where(cond==True)[0]
      if ice == 1:
        maskflatten[ind]=0
      else: 
        maskflatten[ind]=1
 
  # If we want to mask out icefree fjord
  if icefront_time !='none':
    xfjord = []
    yfjord = []
    xfront,yfront,time = icefrontlib.near_time(icefront_time,glacier,type=icefront_type)
    south = meshlib.shp_to_xy(DIR+'icemask_southfjordwall.shp')
    north = meshlib.shp_to_xy(DIR+'icemask_northfjordwall.shp')
    if north[0,1] > north[-1,0]:
      north = np.flipud(north)
    if south[0,1] < south[-1,0]:
      south = np.flipud(south) 
    if yfront[0] > yfront[-1]:
      xfront = np.flipud(xfront)
      yfront = np.flipud(yfront)  
    xfjord = xfront
    yfjord = yfront
    nind = np.argmin((xfront[-1]-north[:,0])**2+(yfront[-1]-north[:,1])**2)
    xfjord = np.r_[xfjord,north[nind:,0]]
    yfjord = np.r_[yfjord,north[nind:,1]]
    sind = np.argmin((xfront[0]-south[:,0])**2+(yfront[0]-south[:,1])**2)
    xfjord = np.r_[xfjord,south[0:sind+1,0]]
    yfjord = np.r_[yfjord,south[0:sind+1,1]]
    fjordpath = path.Path(np.column_stack([xfjord,yfjord]))
    cond = fjordpath.contains_points(points)
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

  DIR = os.path.join(os.getenv("DATA_HOME"),'ShapeFiles/IceMasks/'+glacier+'/')

  files = os.listdir(DIR)
  
  n=0
  for file in files:
    if (file.endswith('.shp')) and ('hole' in file):
      hole = meshlib.shp_to_xy(DIR+file)
      if n==0:
        pts = hole[:,0:2]  
      else:
        pts = np.row_stack([pts,hole[:,0:2]])
      n=1
      
  return pts

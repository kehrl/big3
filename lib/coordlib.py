import numpy as np
from osgeo import osr


def convert(xin,yin,epsg_in,epsg_out):

  srs_in = osr.SpatialReference() #create an input Spatial Reference
  srs_in.ImportFromEPSG(epsg_in) #import from EPSG:4326

  srs_out = osr.SpatialReference() #create an output Spatial Reference
  srs_out.ImportFromEPSG(epsg_out) #import from EPSG:3338

  ct = osr.CoordinateTransformation(srs_in, srs_out)

  try:
    nx,ny=np.shape(xin)
    xout=np.zeros_like(xin)
    yout=np.zeros_like(yin)
    zout=np.zeros_like(xin)
    for j in range(0,ny):
      for i in range(0,nx):
        (xout[i,j],yout[i,j],zout[i,j]) = ct.TransformPoint(float(xin[i,j]),float(yin[i,j])) #Calculate our transformed points
  except:
    try:
      xout=np.zeros_like(xin)
      yout=np.zeros_like(yin)
      zout=np.zeros_like(xin)
      for i in range(0,len(xin)):
        (xout[i],yout[i],zout[i]) = ct.TransformPoint(float(xin[i]),float(yin[i])) #Calculate our transformed points
    except:
      (xout,yout,zout) = ct.TransformPoint(xin,yin)
  
  return (xout,yout)
  
def geoidheight(xin,yin,zin,epsg=3413):
  
  import os
  import sys
  import numpy as np
  from subprocess import Popen, STDOUT, PIPE
  
  # Length of coordinates
  N = len(xin)
  
  srs_in = osr.SpatialReference() #create an input Spatial Reference
  srs_in.ImportFromEPSG(epsg)

  srs_out = osr.SpatialReference() #create an output Spatial Reference
  srs_out.ImportFromProj4("+proj=longlat +ellps=WGS84 +no_defs +geoidgrids=egm08_25.gtx") 
  
  ct = osr.CoordinateTransformation(srs_in, srs_out)
  
  try:
    nx,ny=np.shape(xin)
    xout=np.zeros_like(xin)
    yout=np.zeros_like(yin)
    zout=np.zeros_like(xin)
    for j in range(0,ny):
      for i in range(0,nx):
        (xout[i,j],yout[i,j],zout[i,j]) = ct.TransformPoint(float(xin[i,j]),float(yin[i,j]),float(zin[i,j])) #Calculate our transformed points
  except:
    try:
      xout=np.zeros_like(xin)
      yout=np.zeros_like(yin)
      zout=np.zeros_like(xin)
      for i in range(0,len(xin)):
        (xout[i],yout[i],zout[i]) = ct.TransformPoint(float(xin[i]),float(yin[i]),float(zin[i])) #Calculate our transformed points
    except:
      (xout,yout,zout) = ct.TransformPoint(xin,yin,zin)
  
  return zout
  
def sort_xy(xin,yin):

  import numpy as np

  # Sort points
  ind = np.argmin(xin)
  x2=[] 
  x2.append(xin[ind])
  y2=[]
  y2.append(yin[ind])
  indices=range(0,len(xin))
  indices.remove(ind)
  n=0
  while len(indices) > 1:
    min_d=1000000000
    d=[]
    for i in indices:
      d=(np.sqrt((x2[n]-xin[i])**2+(y2[n]-yin[i])**2))
      if (d < min_d):
        min_d=d
        ind=i
    x2.append(xin[ind])
    y2.append(yin[ind])
    indices.remove(ind)
    n=n+1

  return x2,y2
  
def sortind_xy(xin,yin):

  import numpy as np

  # Sort points
  ind = np.argmin(xin)
  sortindices=[]
  sortindices.append(ind)
  indices=range(0,len(xin))
  indices.remove(ind)
  x2=[] 
  x2.append(xin[ind])
  y2=[]
  y2.append(yin[ind])
  n=0
  while len(indices) > 0:
    min_d=1000000000
    d=[]
    for i in indices:
      d=(np.sqrt((x2[n]-xin[i])**2+(y2[n]-yin[i])**2))
      if (d < min_d):
        min_d=d
        ind=i
    sortindices.append(ind)
    x2.append(xin[ind])
    y2.append(yin[ind])
    indices.remove(ind)
    n=n+1

  return sortindices
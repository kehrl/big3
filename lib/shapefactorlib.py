# This function creates the inputs necessary for the ShapeFactor function in Elmersolver.
#
# Functions:
# glacierwidth(flowline,rightside,leftside) - calculate glacier width
# trapezoid - calculate shapefactor given a trapezoid shape for the glacier
# lateralconvergence- calculate nodal mass input for lateral convergence


def glacierwidth(flowline,file_rightside_in,file_leftside_in,filt_len):
  print "\n## Calculating glacier width ##"
  
  import os
  import sys
  import math
  import numpy as np
  import meshlib
  import scipy.interpolate as interpolate
  import scipy.signal as signal
  import matplotlib.pyplot as plt
  
  R = len(flowline[:,0])
  width=[0]*R
  width_interped=[0]*R
  
  # Load channel boundaries for calculating width
  rightside = meshlib.shp_to_xy(file_rightside_in)
  leftside = meshlib.shp_to_xy(file_leftside_in)
  del file_rightside_in, file_leftside_in
  
  # Set up channel boundaries as functions
  f_left = interpolate.interp1d(leftside[:,0],leftside[:,1],'linear')
  f_right = interpolate.interp1d(rightside[:,0],rightside[:,1],'linear')
  
  # Go through all points in the flowline
  no_interp=0 
  for i in range(1,R-1):
    # Surrounding points
    x1 = flowline[i-1,1]
    y1 = flowline[i-1,2]
    x2 = flowline[i,1]
    y2 = flowline[i,2]
    x3 = flowline[i+1,1]
    y3 = flowline[i+1,2]
  
    # Find intersection between perpendicular line and the "leftside" and "rightside" of the channel
    if y1 == y3: 
      # In case the perpendicular line is vertical, which will make the slope calculation fail
      yR = f_right(x2)
      yL = f_left(x2)
      width[i]=yR-yL
      del yR, yL
    else: 
      # In all other cases, we first find the slope of the perpendicular line, find the line 
      # with that slope, subtract the two functions, and find the minimum distance (the root)
      m_para = (y3-y1)/(x3-x1)
      m_perp = -1/m_para
      b = y2-m_perp*x2
  
      xL = np.linspace(min(leftside[:,0]),max(leftside[:,0]),math.ceil((max(leftside[:,0])-min(leftside[:,0]))/5))
      f_perp_L=interpolate.interp1d(xL,m_perp*xL+b,'linear') 
      xR = np.linspace(min(rightside[:,0]),max(rightside[:,0]),math.ceil((max(rightside[:,0])-min(rightside[:,0]))/5))
      f_perp_R=interpolate.interp1d(xR,m_perp*xR+b,'linear')
      
      ind_left=[]
      mindist=200
      dists=f_perp_L(xL)-f_left(xL)
      for j in range(0,len(dists)):
        if abs(dists[j]) < mindist:
          mindist=abs(dists[j])
          ind_left = j
      del dists, mindist
      ind_right=[]
      mindist=200
      dists=f_perp_R(xR)-f_right(xR)
      for j in range(0,len(dists)):
        if abs(dists[j]) < mindist:
          mindist=abs(dists[j])
          ind_right = j
      if not ind_right or not ind_left:
        no_interp=no_interp+1
      try:
        width[i]=(math.sqrt((xL[ind_left]-xR[ind_right])**2+(f_left(xL[ind_left])-f_right(xR[ind_right]))**2))     
      except:
        pass
      del xL,xR, ind_left, ind_right, mindist, dists, m_para, m_perp
    del x1, x2, x3, y1, y2, y3
  print "Could not find an intersection point for width estimation at",no_interp,"points"
    
  m=[]
  n=[]
  indices=[]
  for i in range(0,R):
    if width[i] is not 0:
      m.append(flowline[i,0])
      n.append(width[i])
      indices.append(i)
  
  fwidth=interpolate.interp1d(m,n)
  width_interped[indices[0]:indices[len(indices)-1]]=fwidth(flowline[indices[0]:indices[len(indices)-1],0])
  for i in range(0,indices[0]):
    width_interped[i] = width_interped[indices[0]]
  for i in range(indices[len(indices)-1],R):
    width_interped[i] = width_interped[indices[len(indices)-2]]
  
  print "Filtering with a butterworth filter with a cutoff of", filt_len, "m"
  cutoff=(1/filt_len)/(1/(np.diff(flowline[1:3,0])*2))
  b,a=signal.butter(4,cutoff,btype='low')
  width_filtered=signal.filtfilt(b,a,width_interped)
  
  return width_filtered
  
def trapezoid(width,height):
  import numpy as np
  
  W=width
  H=height
  nx=len(W)
  
  # Cross-sectional area
  A=(2.0/3.0)*W*H+(1.0/6.0)*W*H # Rectangle + side triangles
  
  # Perimeter in contact with the ice
  P=(2.0/3.0)*W+2.0*np.sqrt(H**2.0+(1.0/36.0)*W**2.0)
  
  # Shape factor, f
  f = A/(H*P)
  
  return f

def dwdx(coord1,width):
  import numpy as np
  
  dw = np.zeros_like(coord1)
  dw[0:-1] = np.diff(width)/np.diff(coord1)
  dw[-1]=dw[-2] 

  return dw
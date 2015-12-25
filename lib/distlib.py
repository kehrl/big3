# Compute distances between points or along transects
#
# LMK, UW, 11/2/2014

def transect(x,y):
  # Calculates distances along a transect
  
  import numpy as np
  
  n=len(x)
  
  dists=np.zeros(n)
  for i in range(1,n):
    dists[i]=dists[i-1]+np.sqrt((x[i-1]-x[i])**2+(y[i-1]-y[i])**2)

  return dists
  
def between_pts(x1,y1,x2,y2):
  # Calculates distances between two points x1,y1 and x2,y2
  
  import numpy as np
  
  dist = np.sqrt((x2-x1)**2+(y2-y1)**2)
    
  return dist
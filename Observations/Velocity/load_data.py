#Export the TSX velocities as text files for use in GMT and Matlab. I eventually decided
#that it makes more sense to just do this in Matlab. LMK, UW, 04/01/2014

import os
import math
import sys
import numpy as np
sys.path.append(os.path.join(os.getenv("HOME"),"Dropbox/Code/Modules"))
import geodat


input=os.path.join(os.getenv("HOME"),"Data/Velocity/RADARSAT/Helheim/winter05-06/mosaicOffsets")
output=os.path.join(os.getenv("HOME"),"Data/Velocity/RADARSAT/Helheim/Outputs/vel_winter05-06.txt")

(x,y,vx,vy,ex,ey)=geodat.readvelocity(input)

(I,J) = np.where(vx!=-2.0e+9)
(imin,imax) = (min(I),max(I))
(jmin,jmax) = (min(J),max(J))

vx = vx[imin:imax+1,jmin:jmax+1]
vy = vy[imin:imax+1,jmin:jmax+1]
ex = ex[imin:imax+1,jmin:jmax+1]
ey = ey[imin:imax+1,jmin:jmax+1]
x = x[jmin:jmax+1]
y = y[imin:imax+1]
nx = len(x)
ny = len(y)

v=np.zeros([ny,nx])
for i in range(0,ny):
  for j in range(0,nx):
    if vx[i,j] == -2.0e+9:
      v[i,j]='NaN'
    else:
      v[i,j]=math.sqrt(vx[i,j]**2+vy[i,j]**2)

fid = open(output,'w')
for i in range(0,ny):
  for j in range(0,nx):
    fid.write('{} {} {} {} {} {} {}\n'.format(x[j],y[i],v[i,j],vx[i,j],vy[i,j],ex[i,j],ey[i,j]))
fid.close()
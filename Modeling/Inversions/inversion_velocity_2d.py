def inversion_velocity_2d(dist,coords,file_velocity_in1,file_velocity_in2,file_velocity_out)
  import os
  import sys
  import numpy as np
  sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
  import scipy.io
  from scipy import interpolate

  # Large velocity map, in case the small map isn't big enough
  xmin=np.min(coords[:,0])
  xmax=np.max(coords[:,0])
  ymin=np.min(coords[:,1])
  ymax=np.max(coords[:,1])
  
  vel2 = scipy.io.loadmat(file_velocity_in2)
  xind1=(abs(xmin-vel2['x'])).argmin()-2
  xind2=(abs(xmax-vel2['x'])).argmin()+2
  yind1=(abs(ymin-vel2['y'])).argmin()-2
  yind2=(abs(ymax-vel2['y'])).argmin()+2
  nx=xind2-xind1
  ny=yind2-yind1
  
  x2,y2=np.meshgrid(vel2['x'][0][xind1:xind2],vel2['y'][0][yind1:yind2])
  x2=np.reshape(x2,[nx*ny,1])
  y2=np.reshape(y2,[nx*ny,1])
  v2=np.reshape(vel2['v'][yind1:yind2,xind1:xind2],[nx*ny,1])
  
  nans=np.isnan(v2)
  x2=x2[~np.isnan(v2)]
  y2=y2[~np.isnan(v2)]
  v2=v2[~np.isnan(v2)]
  coords2=np.column_stack([x2,y2])

  # Individual velocity map
  vel1 = scipy.io.loadmat(file_velocity_in1)
  x1,y1=np.meshgrid(vel1['x'],vel1['y'])
  x1=np.reshape(x1,[len(x1[:,0])*len(x1[0,:]),1])
  y1=np.reshape(y1,[len(x1),1])
  v1=np.reshape(vel1['v'],[len(y1),1])
  coords1=np.column_stack([x1,y1])

  # Interpolate velocities along flowline
  vint2=interpolate.griddata(coords2,v2,coords,method='linear')
  vint1=interpolate.griddata(coords1,v1,coords,method='linear')
  
  # Combine velocity records
  vfinal=vint1
  for i in range(0,len(vfinal)):
    if np.isnan(vfinal[i]):
      vfinal[i]=vint2[i]
  
  # Write out the velocity data
  fid = open(file_velocity_out,'w')


  R=len(vfinal)
  fid.write('{0}\n'.format(R))

  for j in range(0,R):
     fid.write('{} {}\n'.format(dist[j],vfinal[j]))
  fid.close()

  return 0
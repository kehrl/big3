import os
import sys
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules/"))
sys.path.append(os.path.join(os.getenv("HOME"),"Code/BigThreeGlaciers/Tools/"))
import elevation, icefronts, elmer_mesh, coords
import numpy as np
from matplotlib.path import Path
import dist

def date(glacier,time):
  
  '''
  xextent,yextent = date(glacier,time)
  
  Find the glacier extent that is closest to "time".
  
  Inputs:
  glacier: glacier name
  time: time when we want the glacier extent:
  
  xextent,yextent: 1-D arrays of x and y coordinates that define glacier extent for that date
  '''
  
  # Glacier extent with no ice front
  extent = elmer_mesh.shp_to_xy(os.path.join(os.getenv("HOME"),"Data/ShapeFiles/Glaciers/3D/"+glacier+"/glacier_extent_nofront"))
  xextent = extent[:,0]
  yextent = extent[:,1]
  
  dold = dist.transect(xextent,yextent)
  dnew = np.arange(0,dold[-1],20.0)
  xextent = np.interp(dnew,dold,xextent)
  yextent = np.interp(dnew,dold,yextent)
  
  # Terminus coordinates
  xterminus,yterminus,time_terminus = icefronts.near_time(time,glacier)
  if yterminus[-1] > yterminus[0]:
    xterminus = np.flipud(xterminus)
    yterminus = np.flipud(yterminus)
  
  glacierperimeter = Path(extent[:,0:2])
  ind = glacierperimeter.contains_points(np.column_stack([xterminus,yterminus]))
  xterminus = xterminus[ind]
  yterminus = yterminus[ind]

  ind1 = np.argmin((xterminus[0]-xextent)**2+(yterminus[0]-yextent)**2)
  ind2 = np.argmin((xterminus[-1]-xextent)**2+(yterminus[-1]-yextent)**2)
  if ind2 > ind1:
    xextent = np.r_[np.delete(xextent[0:ind2],range(ind1+1,ind2)),xterminus,xextent[ind2+1:]]
    yextent = np.r_[np.delete(yextent[0:ind2],range(ind1+1,ind2)),yterminus,yextent[ind2+1:]]
  else:
    xextent = np.r_[np.delete(xextent[0:ind1],range(ind2+1,ind1)),xterminus,xextent[ind1+1:]]
    yextent = np.r_[np.delete(yextent[0:ind1],range(ind2+1,ind1)),yterminus,yextent[ind1+1:]]
    
  
  return xextent,yextent
  
  
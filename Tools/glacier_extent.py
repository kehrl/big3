import os
import sys
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules/"))
sys.path.append(os.path.join(os.getenv("HOME"),"Code/BigThreeGlaciers/Tools/"))
import elevation, icefronts, elmer_mesh, coords
import numpy as np

def date(glacier,time):
  
  # Glacier extent with no ice front
  extent = elmer_mesh.shp_to_xy(os.path.join(os.getenv("HOME"),"Data/ShapeFiles/Glaciers/3D/"+glacier+"/glacier_extent_nofront"))
  xextent = extent[:,0]
  yextent = extent[:,1]
  
  # Terminus coordinates
  xterminus,yterminus,time_terminus = icefronts.near_time(time,glacier)
  if yterminus[0] > yterminus[-1]:
    xterminus = np.flipud(xterminus)
    yterminus = np.flipud(yterminus)
  
  ind1 = np.argmin((xterminus[0]-xextent)**2+(yterminus[0]-yextent)**2)
  ind2 = np.argmin((xterminus[-1]-xextent)**2+(yterminus[-1]-yextent)**2)
  if ind1 > ind2:
    xextent = np.r_[np.delete(xextent,range(ind2+1,ind1)),xterminus]
    yextent = np.r_[np.delete(yextent,range(ind2+1,ind1)),yterminus]
  else:
    xextent = np.r_[np.delete(xextent,range(ind1+1,ind2)),xterminus]
    yextent = np.r_[np.delete(yextent,range(ind1+1,ind2)),yterminus]
    
  xextent,yextent = coords.sort_xy(xextent,yextent)
  
  return xextent,yextent
  
  
import os
import shapefile
import jdcal
from shapely.geometry import LineString
import numpy as np

def distance_along_flowline(x,y,dists):

  DIRI=os.path.join(os.getenv("HOME"),"Data/ShapeFiles/IceFronts/Helheim/")

  files = os.listdir(DIRI)

  terminus_val = []
  terminus_time = []

  lineflow = LineString(np.row_stack([np.column_stack([x,y]),[315715,-2577820]]))
  n = 0
  for file in files:
    if file.endswith('.shp') and (not "moon" in file):
      # Load shapefile
      sf = shapefile.Reader(DIRI+file)
      shapes = sf.shapes()
      termpts = np.array(shapes[0].points[:])
      lineterm = LineString(termpts)
    
      # Find intersection
      intersect = (lineflow.intersection(lineterm))
      if not(intersect.is_empty):
        flowind = (abs(x-intersect.x)).argmin()
        if x[flowind] > intersect.x: #make sure that we have the smaller value
          flowind=flowind-1;
    
      # Terminus position
      terminus_val.append( dists[flowind]+((intersect.x-x[flowind])**2+(intersect.y-y[flowind])**2)**(0.5) )
    
      # Time of that terminus position
      if ("TSX" in file) or ("moon" in file):
        year = float(file[0:4])
        day = float(file[5:8])
        day1 = jdcal.gcal2jd(year,1,1) # Get length of year
        day2 = jdcal.gcal2jd(year+1,1,1)
        terminus_time.append( year + day/(day2[1]+day2[0]-day1[0]-day1[1]))
      elif ("ASTER" in file) or ("Landsat" in file):
        year = float(file[0:4])
        day = jdcal.gcal2jd(year,float(file[5:7]),float(file[8:10]))
        day2 = jdcal.gcal2jd(year,12,31)
        day1 = jdcal.gcal2jd(year-1,12,31)
        doy = day[1]+day[0]-day1[1]-day1[0]
        terminus_time.append( year + doy/(day2[1]+day2[0]-day1[0]-day1[1]))

  terminus_time = np.array(terminus_time)
  terminus_val = np.array(terminus_val)

  sortind=np.argsort(terminus_time,0)
  terminus_time = terminus_time[sortind]
  terminus_val = terminus_val[sortind]
  
  return terminus_val, terminus_time
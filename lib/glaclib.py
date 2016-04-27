import os
import sys
import icefrontlib, meshlib, distlib, bedlib, geotifflib, datelib
import numpy as np
from matplotlib.path import Path
import scipy.interpolate
import scipy.signal

def load_flowline(glacier,shapefilename='center_flowline',filt_len=2.0e3,verticaldatum='geoid',bedsource='cresis',bedmodel='aniso',bedsmoothing=4):

  '''
  x,y,zb_filt,dists = load(glacier,shapefilename='center_flowline')
  
  Load glacier flowline. This script is mostly to keep everything consistent (distance 
  along flowline, chosen bed profile,etc.
  
  Inputs:
  glacier: glacier name
  shapefilename: shapefile to use for the flowline
  filt_len: filter length (in meters) for the bed profile
  verticaldatum: geoid or ellipsoid
  
  Outputs:
  x,y: x,y coordinates of flowline
  zb_filt: bed profile for flowline
  '''

  if glacier == "Helheim":
    file1 = 'helheim_'
    file_flowline_in = os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/Glaciers/Flowlines/Helheim/"+file1+shapefilename)
  elif glacier == 'Kanger':
    file1 = 'kanger_'
    file_flowline_in = os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/Glaciers/Flowlines/Kanger/"+file1+shapefilename)
  else:
    sys.exit("Unknown glacier.") 
  
  flowline = meshlib.shp_to_xy(file_flowline_in)
  d = distlib.transect(flowline[:,0],flowline[:,1])

  # Set uniform spacing between nodes along flowline
  dists_old = np.linspace(0,np.max(d),np.max(d)/20)
  x = np.interp(dists_old,d,flowline[:,0])
  y = np.interp(dists_old,d,flowline[:,1])
  
  # Get new distances
  dists = distlib.transect(x,y)

  # Get average terminus position
  time1 = 2008.0
  time2 = 2016.0
  
  # Get terminus positions so that we can set distance along flowline 
  # relative to the average terminus position
  terminus_val, terminus_time = icefrontlib.distance_along_flowline(x,y,dists,glacier,type='icefront',time1=time1,time2=time2)

  # Average terminus position
  terminus = np.mean(terminus_val)

  # Set dists relative to average terminus position
  dists = dists-terminus

  # Find bed elevation
  if (glacier == 'Helheim') and (bedsource == 'smith'):
    zb = bedlib.smith_at_pts(x,y,glacier,model=bedmodel,smoothing=bedsmoothing,verticaldatum=verticaldatum)
  elif bedsource == 'morlighem':
    zb = bedlib.morlighem_pts(x,y,verticaldatum=verticaldatum)
  elif (glacier == 'Kanger' and (bedsource != 'morlighem')) or (bedsource == 'cresis'):
    cresis = bedlib.cresis('all',glacier,verticaldatum=verticaldatum)
    if glacier == 'Helheim':
      cresis2001 = bedlib.cresis('2001',glacier,verticaldatum=verticaldatum)
      cresis = np.row_stack([cresis,cresis2001])
    cutdist = 200.
    dcresis = []
    zcresis = []
    tcresis = []
    for i in range(0,len(cresis[:,0])):
      mindist = np.min(np.sqrt((cresis[i,0]-x)**2+(cresis[i,1]-y)**2))
      if mindist < cutdist:
        minind = np.argmin(np.sqrt((cresis[i,0]-x)**2+(cresis[i,1]-y)**2))
        dcresis.append(dists[minind])
        zcresis.append(cresis[i,2])
        tcresis.append(cresis[i,4])

    ind = np.argsort(dcresis)
    dcresis = np.array(dcresis)[ind]
    zcresis = np.array(zcresis)[ind]
    zb = np.interp(dists,dcresis,zcresis)


  if filt_len != 'none':
    ind = np.where(~(np.isnan(zb)))[0]
    zb_filt = np.zeros_like(zb)
    zb_filt[:] = float('nan')
    cutoff=(1/filt_len)/(1/(np.diff(dists[1:3])*2))
    b,a=scipy.signal.butter(4,cutoff,btype='low')
    zb_filt[ind] = scipy.signal.filtfilt(b,a,zb[ind])
  else:
    zb_filt = zb
    
  if (glacier == 'Kanger') and (shapefilename=='flowline_flightline'):
    ind = np.where(dists > 3000.)[0]
    zb_filt[ind] = float('nan')
  elif (glacier == 'Helheim') and (shapefilename=='flowline_flightline'):
    ind = np.where(dists > 3900.)[0]
    zb_filt[ind] = float('nan')
    

  return x,y,zb_filt,dists

def load_extent(glacier,time,nofront_shapefile='glacier_extent_nofront'):
  
  '''
  xextent,yextent = date(glacier,time)
  
  Find the glacier extent that is closest to "time".
  
  Inputs:
  glacier: glacier name
  time: time (fractional year) when we want the glacier extent
  
  Outputs:
  xextent,yextent: 1-D arrays of x and y coordinates that define glacier extent for that date
  '''
  
  # Glacier extent with no ice front
  extent = meshlib.shp_to_xy(os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/Glaciers/3D/"+glacier+"/"+nofront_shapefile))
  
  if extent[1,1] > extent[0,1]:
    extent = np.flipud(extent)
  
  xextent = extent[:,0]
  yextent = extent[:,1]
  bound = extent[:,2]
  
  dold = distlib.transect(xextent,yextent)
  dnew = np.arange(0,dold[-1],20.0)
  xextent = np.interp(dnew,dold,xextent)
  yextent = np.interp(dnew,dold,yextent)

  f = scipy.interpolate.interp1d(dold,bound,kind='nearest')
  bound = f(dnew)
  
  # Terminus coordinates
  xterminus,yterminus,time_terminus = icefrontlib.near_time(time,glacier)
  if yterminus[-1] > yterminus[0]:
    xterminus = np.flipud(xterminus)
    yterminus = np.flipud(yterminus)
    bound = np.flipud(bound)
  
  glacierperimeter = Path(extent[:,0:2])
  ind = glacierperimeter.contains_points(np.column_stack([xterminus,yterminus]))
  xterminus = xterminus[ind]
  yterminus = yterminus[ind]
  boundterminus = np.ones(len(xterminus))*2.0

  ind1 = np.argmin((xterminus[0]-xextent)**2+(yterminus[0]-yextent)**2)
  ind2 = np.argmin((xterminus[-1]-xextent)**2+(yterminus[-1]-yextent)**2)
  if ind2 > ind1:
    if yextent[ind2] < yextent[ind1]:
      xextent = np.r_[np.delete(xextent[0:ind2],range(ind1+1,ind2)),np.flipud(xterminus),xextent[ind2+1:]]
      yextent = np.r_[np.delete(yextent[0:ind2],range(ind1+1,ind2)),np.flipud(yterminus),yextent[ind2+1:]]
      bound = np.r_[np.delete(bound[0:ind2],range(ind1+1,ind2)),boundterminus,bound[ind2+1:]]
    else:
      xextent = np.r_[np.delete(xextent[0:ind2],range(ind1+1,ind2)),xterminus,xextent[ind2+1:]]
      yextent = np.r_[np.delete(yextent[0:ind2],range(ind1+1,ind2)),yterminus,yextent[ind2+1:]]
      bound = np.r_[np.delete(bound[0:ind2],range(ind1+1,ind2)),boundterminus,bound[ind2+1:]]
  else:
    if yextent[ind2] < yextent[ind1]:
      xextent = np.r_[np.delete(xextent[0:ind1],range(ind2+1,ind1)),np.flipud(xterminus),xextent[ind1+1:]]
      yextent = np.r_[np.delete(yextent[0:ind1],range(ind2+1,ind1)),np.flipud(yterminus),yextent[ind1+1:]]
      bound = np.r_[np.delete(bound[0:ind1],range(ind2+1,ind1)),np.flipud(boundterminus),bound[ind1+1:]]
    else:
      xextent = np.r_[np.delete(xextent[0:ind1],range(ind2+1,ind1)),xterminus,xextent[ind1+1:]]
      yextent = np.r_[np.delete(yextent[0:ind1],range(ind2+1,ind1)),yterminus,yextent[ind1+1:]]
      bound = np.r_[np.delete(bound[0:ind1],range(ind2+1,ind1)),boundterminus,bound[ind1+1:]]
  
  return np.column_stack([xextent,yextent,bound])

def load_satimages(glacier,xmin,xmax,ymin,ymax,time1=-np.inf,time2=np.inf,data='all'):

  '''
  Load satellite images for a particular glacier for the grid defined by xmin,xmax,ymin,
  ymax, over the time interval time1 to time2.
  '''
  
  DIRLANDSAT = os.path.join(os.getenv("DATA_HOME"),"Imagery/Landsat/"+glacier+"/TIF/")
  DIRWV = os.path.join(os.getenv("DATA_HOME"),"Imagery/Worldview/"+glacier+"/")
  DIRTSX = os.path.join(os.getenv("DATA_HOME"),"Mosaics/"+glacier+"/")
  
  # Find files to load
  dirs = []
  times = []
  types = []
  images = []
  
  # Check landsat files
  if ('LANDSAT' in data) or (data=='all'):
    files = os.listdir(DIRLANDSAT)
    for file in files:
      if file.endswith('.tif'):
        filetime = datelib.date_to_fracyear(float(file[0:4]),float(file[4:6]),float(file[6:8]))
        if (filetime >= time1) and (filetime <= time2):
          types.append('Landsat')
          dirs.append(DIRLANDSAT+file)
          times.append(filetime)
          images.append(geotifflib.read(DIRLANDSAT+file))
  
  # Check TSX files
  if ('TSX' in data) or (data == 'all'):  
    files = os.listdir(DIRTSX)
    for file in files:
      if file.endswith('.tif'):
        if glacier == 'Helheim':
          filetime = datelib.doy_to_fracyear(float(file[14:18]),float(file[19:22]))
        elif glacier == 'Kanger':
          filetime = datelib.doy_to_fracyear(float(file[11:15]),float(file[16:19]))        
        if (filetime >= time1) and (filetime <= time2):
          types.append('TSX')
          dirs.append(DIRLANDSAT+file)
          times.append(filetime)
          images.append(geotifflib.read(DIRTSX+file))
  
  sortind = np.argsort(times)
  images_sorted = []
  types_sorted = []
  times_sorted = []
  for ind in sortind:
    images_sorted.append(images[ind])
    types_sorted.append(types[ind])
    times_sorted.append(times[ind])
    
          
  return images_sorted,times_sorted,types_sorted
     
          
       
    
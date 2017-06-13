import os
import sys
import icefrontlib, meshlib, distlib, bedlib, geotifflib, datelib, vellib
import numpy as np
from shapely.geometry import LineString
from matplotlib.path import Path
import scipy.interpolate
import scipy.signal

def load_flowline(glacier,shapefilename='center_flowline',filt_len=2.0e3,verticaldatum='geoid',bedsource='cresis',bedmodel='aniso',bedsmoothing=4,dx=20):

  '''
  x,y,zb_filt,dists = load(glacier,shapefilename='center_flowline')
  
  Load glacier flowline. This script is mostly to keep everything consistent (distance 
  along flowline, chosen bed profile,etc.
  
  Inputs:
  glacier: glacier name
  shapefilename: shapefile to use for the flowline
  filt_len: filter length (in meters) for the bed profile
  verticaldatum: geoid or ellipsoid
  dx: distance between points
  
  Outputs:
  x,y: x,y coordinates of flowline
  zb_filt: bed profile for flowline
  '''

  if glacier == 'Helheim':
    file1 = 'helheim_'
  elif glacier == 'Kanger':
    file1 = 'kanger_'
  elif glacier == 'Midgaard':
    file1 = 'midgaard_'
  else:
    sys.exit("Unknown glacier.") 
  file_flowline_in = os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/Glaciers/Flowlines/"+glacier+"/"+file1+shapefilename)

  
  flowline = meshlib.shp_to_xy(file_flowline_in)
  d = distlib.transect(flowline[:,0],flowline[:,1])

  # Set uniform spacing between nodes along flowline
  dists_old = np.linspace(0,np.max(d),np.max(d)/dx)
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
  extent = load_extent(glacier,time)
  
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

def load_extent_timeseries(glacier,time1,time2,dt,nofront_shapefile='glacier_extent_nofront',datatypes=['Landsat','TSX','WV']):

  '''
  time, xextents, yxextents, bounds = load_extent_timeseries(glacier,time1,
     time2,dt,nofront_shapefile='glacier_extent_nofront',
     datatypes=['Landsat','TSX','WV'])
     
  Interpolates ice-front positions from picked ice-front positions to create 
  a timeseries of meshes to be used in the terminus-driven model.
  
  Inputs:
  glacier           : glacier name (Helheim, Kanger)
  time1             : fractional start time for timeseries
  time2             : fractional end time for timeseries
  dt                : timestep
  nofront_shapefile : nofront shapefile name for mesh extent
  datatypes         : satellite image types for picked ice fronts
  
  Outputs:
  time     : interpolated time between time1,time2 with timestep dt
  xextents : 2-d array of x-coordinates of extents
  yextents : 2-d array of y-coordinates of extents
  bounds   : boundary numbers for extents
  '''

  # Glacier extent with no ice front
  extent = meshlib.shp_to_xy(os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/Glaciers/3D/"+glacier+"/"+nofront_shapefile))
  if extent[1,1] > extent[0,1]:
    extent = np.flipud(extent)
  
  xextent = extent[:,0]
  yextent = extent[:,1]
  bound_extent = extent[:,2]

  # Interpolate glacier extent to a finer grid to make ice-front interpolation easier
  dextent = distlib.transect(xextent,yextent)
  #dextent = np.arange(0,dold[-1],20.0)
  #xextent = np.interp(dextent,dold,xextent)
  #yextent = np.interp(dextent,dold,yextent)
  #f = scipy.interpolate.interp1d(dold,bound,kind='nearest')
  #bound = f(dextent)
  extent = LineString(np.column_stack([extent[:,0:2]]))

  # Load all ice front positions for that time period
  termx,termy,termt = icefrontlib.load_all(time1-0.5,time2+0.5,glacier,type='icefront',datatypes=datatypes)

  # In case we have multiple ice front picks for the same day, we want to
  # use only one of those for continuity
  [junk,ind] = np.unique(termt,return_index=True)
  termt = np.array(termt)[ind]
  termx = termx[:,ind]
  termy = termy[:,ind]

  # Load a velocity profile to use as interpolation direction to figure out ice-front position
  # on timesteps that fall between picked ice fronts
  x,y,u = geotifflib.read(os.path.join(os.getenv("DATA_HOME"),"Velocity/TSX/"+glacier+"/TIF/all-2008-2016_vx.tif"))
  x,y,v = geotifflib.read(os.path.join(os.getenv("DATA_HOME"),"Velocity/TSX/"+glacier+"/TIF/all-2008-2016_vy.tif"))  

  fu = scipy.interpolate.RegularGridInterpolator((y,x),u)
  fv = scipy.interpolate.RegularGridInterpolator((y,x),v)

  # Get ice-front position for each timestep
  time = np.arange(time1,time2+dt,dt)
  timeseries_x = np.zeros([len(termx[:,0]),len(time)])
  timeseries_y = np.zeros([len(termx[:,0]),len(time)])
  xextents = np.zeros([len(termx[:,0])+len(xextent),len(time)])
  yextents = np.zeros([len(termx[:,0])+len(xextent),len(time)])
  bounds = np.zeros([len(termx[:,0])+len(xextent),len(time)])
  for i in range(0,len(time)):
    # Find picked ice-front positions for before and after timestep for the interpolation
    ind = np.argmin(abs(time[i]-termt))
    if termt[ind] < time[i]:
      ind1 = ind
      ind2 = ind+1
    else:
      ind1 = ind-1
      ind2 = ind
    
    # Fractional time between ind1,ind2 to use for interpolation
    frac = (time[i]-termt[ind1])/(termt[ind2]-termt[ind1])
    
    # Get picked ice-front positions that we will use for the interpolation
    nonnan = np.where(~(np.isnan(termx[:,ind1])))[0]
    termx1 = termx[nonnan,ind1]
    termy1 = termy[nonnan,ind1]
    nonnan = np.where(~(np.isnan(termx[:,ind2])))[0]
    termx2 = termx[nonnan,ind2]
    termy2 = termy[nonnan,ind2]   
    
    if termy1[-1] > termy1[0]:
      termx1 = np.flipud(termx1)
      termy1 = np.flipud(termy1)
    if termy2[-1] > termy2[0]:
      termx2 = np.flipud(termx2)
      termy2 = np.flipud(termy2)
  
    # Get locations where interpolate ice front intersects the glacier extent
    # First, get intersection pts for two closest ice-front positions in time
    term1 = LineString(np.column_stack([termx1,termy1]))
    intersect = extent.intersection(term1)
    try:
      if len(intersect) == 2:
        if intersect[0].y > intersect[1].y:
          top1 = [intersect[0].x,intersect[0].y]
          bot1 = [intersect[1].x,intersect[1].y]
        else:
          top1 = [intersect[1].x,intersect[1].y]
          bot1 = [intersect[0].x,intersect[0].y]
      else:
        print "Need to look at date ", datelib.fracyear_to_date(termt[ind1])
    except:
      print "Need to look at date ", datelib.fracyear_to_date(termt[ind1])
    
    term2 = LineString(np.column_stack([termx2,termy2]))
    intersect = extent.intersection(term2)
    try:
      if len(intersect) == 2:
        if intersect[0].y > intersect[1].y:
          top2 = [intersect[0].x,intersect[0].y]
          bot2 = [intersect[1].x,intersect[1].y]
        else:
          top2 = [intersect[1].x,intersect[1].y]
          bot2 = [intersect[0].x,intersect[0].y]
      else:
        print "Need to look at date ", datelib.fracyear_to_date(termt[ind2])
    except:
      print "Need to look at date ", datelib.fracyear_to_date(termt[ind2])
    # Now find new intersection points
    if top1[0] < top2[0]: # advancing on this side
      ind_top = np.where((xextent > top1[0]) & (xextent < top2[0]) & (abs(top1[1]-yextent) < 500.))[0]
      sortind = np.argsort(xextent[ind_top])
      xtops = np.r_[top1[0],xextent[ind_top[sortind]],top2[0]]
      ytops = np.r_[top1[1],yextent[ind_top[sortind]],top2[1]]
      dtops = distlib.transect(xtops,ytops)
      dtop = dtops[-1]*frac
    elif top1[0] > top2[0]: # retreating on this side
      ind_top = np.where((xextent < top1[0]) & (xextent > top2[0]) & (abs(top1[1]-yextent) < 500.))[0]
      sortind = np.argsort(xextent[ind_top])
      xtops = np.r_[top2[0],xextent[ind_top[sortind]],top1[0]]
      ytops = np.r_[top2[1],yextent[ind_top[sortind]],top1[1]]
      dtops = distlib.transect(xtops,ytops)
      dtop = dtops[-1]*(1-frac)
    else:
      print "not advancing or retreating on top"
    xtop = np.interp(dtop,dtops,xtops)
    ytop = np.interp(dtop,dtops,ytops) 
      
    if bot1[0] < bot2[0]: # advancing on this side
      ind_bot = np.where((xextent > bot1[0]) & (xextent < bot2[0]) & (abs(bot1[1]-yextent) < 500.))[0]
      sortind = np.argsort(xextent[ind_bot])
      xbots = np.r_[bot1[0],xextent[ind_bot[sortind]],bot2[0]]
      ybots = np.r_[bot1[1],yextent[ind_bot[sortind]],bot2[1]]
      dbots= distlib.transect(xbots,ybots)
      dbot = (dbots[-1])*frac
    elif bot1[0] > bot2[0]: # retreating on this side
      ind_bot = np.where((xextent < bot1[0]) & (xextent > bot2[0]) & (abs(bot1[1]-yextent) < 500.))[0]
      sortind = np.argsort(xextent[ind_bot])
      xbots = np.r_[bot2[0],xextent[ind_bot[sortind]],bot1[0]]
      ybots = np.r_[bot2[1],yextent[ind_bot[sortind]],bot1[1]]
      dbots= distlib.transect(xbots,ybots)
      dbot = (dbots[-1])*(1-frac)
    else:
      print "not advancing or retreating on bot"

    xbot = np.interp(dbot,dbots,xbots)
    ybot = np.interp(dbot,dbots,ybots)

    # Now that we know the bottom and top points (extent of the ice front), we can start
    # calculating the shape, again based on linear interpolation
    # May need to change next expression to find indices between the sidewalls for Kanger, 
    # but this should work for Helheim
    ind_term1 = np.where((termy1 > bot1[1]) & (termy1 < top1[1]))[0]
    icefront_x = []
    icefront_y = []
    icefront_x.append(xtop)
    icefront_y.append(ytop)
    for j in ind_term1:
      # Get velocities to create a line, to interpolate between ice fronts
      uj = fu((termy1[j],termx1[j]))
      vj = fv((termy1[j],termx1[j]))
    
      # Create flowline that intersects that point of the ice front
      xunit = uj/(np.sqrt(uj**2+vj**2))
      yunit = vj/(np.sqrt(uj**2+vj**2))
      b = termy1[j] - (yunit/xunit)*termx1[j]
      flowlinex = np.arange(-3000.,3005.,10) + termx1[j]
      flowliney = (yunit/xunit)*flowlinex + b
      flowline = LineString(np.column_stack([flowlinex,flowliney]))
    
      # Find where flowline intersects the next ice-front position
      intersect = flowline.intersection(term2)
      add = False
      try:   
        if len(intersect) > 0: 
          ind = np.argmin(abs([intersect[k].x for k in range(0,len(intersect))]-termx1[j]))
          term2_flowline = [intersect[ind].x,intersect[ind].y]
          add = True
      except:
        try:
          term2_flowline = [intersect.x,intersect.y]
          add = True
        except:
          pass
      if add == True:
        dflow = distlib.between_pts(termx1[j],termy1[j],term2_flowline[0],term2_flowline[1])
        dmid = frac*dflow
        xmid = np.interp(dmid,[0,dflow],[termx1[j],term2_flowline[0]])
        ymid = np.interp(dmid,[0,dflow],[termy1[j],term2_flowline[1]])
        icefront_x.append(xmid)
        icefront_y.append(ymid)
    
    icefront_x.append(xbot)
    icefront_y.append(ybot)
    
    # Try sorting icefront to get rid of potential tangles
    icefront_x_old = np.asarray(icefront_x)
    icefront_y_old = np.asarray(icefront_y)
    icefront_x = np.zeros_like(icefront_x_old)
    icefront_y = np.zeros_like(icefront_y_old)
    icefront_x[0] = icefront_x_old[0]
    icefront_y[0] = icefront_y_old[0]
    ind = range(1,len(icefront_x_old))
    for k in range(1,len(icefront_x_old)):
      mindist = 10000.
      for j in range(0,len(ind)):
        dist = distlib.between_pts(icefront_x[k-1],icefront_y[k-1],icefront_x_old[ind[j]],icefront_y_old[ind[j]])
        if dist < mindist:
          mindist = dist
          minind = ind[j]
      icefront_x[k] = icefront_x_old[minind]
      icefront_y[k] = icefront_y_old[minind]
      ind.remove(minind)
    
    # Save icefront in timeseries 
    timeseries_x[0:len(icefront_x),i] = icefront_x
    timeseries_y[0:len(icefront_y),i] = icefront_y
    
    # Now create mesh extent and BC numbers using the interpolated ice front
    boundterminus = np.ones(len(icefront_x))*2.0
    ind1 = np.where((xextent > xbot) & (abs(yextent - ybot) < 1.0e3))[0][0] 
    ind2 = np.where((xextent > xtop) & (abs(yextent - ytop) < 1.0e3))[0][-1]

    extent_x = np.r_[xextent[0:ind1],np.flipud(icefront_x),xextent[ind2+1:]]
    extent_y = np.r_[yextent[0:ind1],np.flipud(icefront_y),yextent[ind2+1:]]
    bound = np.r_[bound_extent[0:ind1],boundterminus,bound_extent[ind2+1:]]
    
    xextents[0:len(extent_x),i] = extent_x
    yextents[0:len(extent_x),i] = extent_y
    bounds[0:len(extent_x),i] = bound    
    
  return  time, xextents, yextents, bounds

def load_satimages(glacier,xmin,xmax,ymin,ymax,time1=-np.inf,time2=np.inf,data='all'):

  '''
  images,times,types = load_satimages(glacier,xmin,xmax,ymin,ymax,time1=-np.inf,time2=np.inf,data='all')

  Load satellite images for a particular glacier for the grid defined by xmin,xmax,ymin,
  ymax, over the time interval time1 to time2. If time1 == time2, the code will find the 
  sat image closest to the chosen time.
  
  '''
  
  DIRLANDSAT = os.path.join(os.getenv("DATA_HOME"),"Imagery/Landsat/"+glacier+"/TIF/")
  DIRWV = os.path.join(os.getenv("DATA_HOME"),"Imagery/Worldview/"+glacier+"/")
  DIRTSX = os.path.join(os.getenv("DATA_HOME"),"Mosaics/"+glacier+"/")
  
  # Find files to load
  dirs = []
  if time1 == time2:
    times = 0.0
  else:
    times = []
  types = []
  images = []
  
  # Check landsat files
  if ('LANDSAT' in data) or (data=='all'):
    files = os.listdir(DIRLANDSAT)
    for file in files:
      if file.endswith('.tif'):
        filetime = datelib.date_to_fracyear(float(file[0:4]),float(file[4:6]),float(file[6:8]))
        if (time1 == time2): 
          if (abs(filetime - times) < abs(time1 - times)):
            # Last constrain is to prevent unnecessary loading of files, 
            types = 'Landsat'
            dirs = DIRLANDSAT+file
            times = filetime
        elif (filetime >= time1) and (filetime <= time2):
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
        if (time1 == time2):
          if (abs(filetime - times) <= abs(time1 - times)) and file.endswith('_1-20mgeo.tif'):
            types = 'TSX'
            dirs = DIRTSX+file
            times = filetime
        elif (filetime >= time1) and (filetime <= time2):
          types.append('TSX')
          dirs.append(DIRLANDSAT+file)
          times.append(filetime)
          images.append(geotifflib.read(DIRTSX+file))
  
  if time1 != time2:
    sortind = np.argsort(times)
    images_sorted = []
    types_sorted = []
    times_sorted = []
    for ind in sortind:
      images_sorted.append(images[ind])
      types_sorted.append(types[ind])
      times_sorted.append(times[ind])
  else:
    images_sorted = geotifflib.read(dirs)
    times_sorted = times
    types_sorted = types
    
          
  return images_sorted,times_sorted,types_sorted        
       
    

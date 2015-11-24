# This script loads different profiles for the Helheim bed.

# LMK, UW, 10/06/2014

import os
import sys
import numpy as np
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"Util/Modules"))
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"BigThreeGlaciers/Tools"))
import coords, elevation, flotation
import scipy.interpolate
import geotiff

def cresis(year,glacier,verticaldatum='geoid',cleanup=True):

  '''
  cresis = cresis(year,glacier,verticaldatum='geoid',cleanup=True)
  
  Loads CreSIS data. You can either select a year to load (look at options below) or 
  load all years (year='all'). If you load all years, there is an option to cleanup the data
  using the option cleanup=True.
  
  Inputs:
  year: load selected years or 'all' years (see options below)
  glacier: glacier name
  verticaldatum: ellipsoid or geoid
  cleanup: True or False
  
  Outputs:
  cresis: four column array of x,y,zb, and date for radar picks
  
  '''

  if (year == '2001'):
    if glacier == 'Helheim':
      file = os.path.join(os.getenv("DATA_HOME"),"Bed/Cresis/Helheim/helheim_20010521_good_wgs84.csv")
      ind = range(3180,3297)
    elif glacier == 'Kanger':
      file = os.path.join(os.getenv("DATA_HOME"),"Bed/Cresis/Kanger/Data_20010520_01.csv")
      ind=range(4387,4475)
      
    data = np.loadtxt(file,skiprows=1,delimiter=',')
    y=data[:,0]
    x=data[:,1]
    H=data[:,3]
    x2,y2 = coords.convert(x,y,4326,3413)
    
    surf = elevation.atm('2001','ellipsoid')
    if glacier == 'Helheim':
      zs = scipy.interpolate.griddata(surf['20010521'][:,0:2],surf['20010521'][:,2],np.column_stack([x2,y2]),method='nearest')
      dates = np.ones(len(x2))*20010521
    elif glacier == 'Kanger':
      zs = scipy.interpolate.griddata(surf['20010520'][:,0:2],surf['20010520'][:,2],np.column_stack([x2,y2]),method='nearest')
      dates = np.ones(len(x2))*20010520
    
    zb = zs-H
  
  else:
    if glacier == 'Helheim':
      print "Using data set Helheim_2006_2014_Composite"
      file=os.path.join(os.getenv("DATA_HOME"),"Bed/Cresis/Helheim/Helheim_2006_2014_Composite/flightlines/Helheim_2006_2014_Composite_Flightlines.txt")
    elif glacier == 'Kanger':
      print "Using data set Kanger_2006_2014_Composite"
      file=os.path.join(os.getenv("DATA_HOME"),"Bed/Cresis/Kanger/Kangerdlugssuaq_2006_2014_Composite/flightlines/Kangerdlugssuaq_2006_2014_Composite_Flightlines.txt")
    else: 
      print "unknown glacier"
    
    fid = open(file,"r")
    x=[]
    y=[]
    zs=[]
    zb=[]
    years=[]
    type=[]
    dates = []
    lines = fid.readlines()
    for line in lines[1:-1]:
      fields = line.split()
      try:
        years.append(float(fields[8][0:4]))
      except:
        years.append(0)
      x.append(float(fields[5]))
      y.append(float(fields[4]))
      zb.append(float(fields[1]))
      zs.append(float(fields[0]))
      type.append(fields[3])
      dates.append(float(fields[8]))
    years=np.array(years)
    x=np.array(x)
    y=np.array(y)
    zb=np.array(zb)
    zs=np.array(zs)
    dates=np.array(dates)
    
    x2,y2 = coords.convert(x,y,4326,3413)
    
    if glacier == 'Helheim':
      # Radar picks to toss, if cleanup is set to true. The indices are for the Helheim 
      # 006-2014 composite product.
      if cleanup==True:
        badind = np.r_[range(22406,22513),range(26969,27020),range(29300,29336),range(30253,30408),range(33224,33440),range(33552,33620),range(37452,37531),\
      		range(37640,37675),range(37819,37836),range(44127,44207),range(46942,47030),range(53595,53663),\
      		range(53713,53793),range(53974,53987),range(56646,56726),range(64006,64013),range(61237,61529),range(61745,62000),range(68541,68810),\
      		range(69202,69475),range(75645,75904),range(77285,77538),range(77728,77970)] 
        # Find locations where ice is floating so we can toss those out
        floatind = np.where(coords.geoidheight(x2,y2,zs) < flotation.height(coords.geoidheight(x2,y2,zb),rho_i=917.,rho_sw=1025.))[0]
        badind = np.union1d(badind,floatind)
      else: 
        badind = []
        floatind = []
      if year == '2006a':
        ind = range(22285,22406)
      elif year == '2006b':
        ind = range(26856,26969)
      elif year == '2008a':
        ind = range(29200,29335)
      elif year == '2008b':
        ind = range(30121,30253)
      elif year == '2009':
        ind = range(53500,53595)
      elif year == '2011':
        ind = range(61240,61265)
      elif year == '2012':
        ind = range(63900,64013)
      elif year == '2013':
        ind = range(68400,68542)
      elif year == '2014':
        ind = range(77100,77285)
      elif year == 'all':
        ind = []
        for i in range(0,len(type)):
          if ("ICESat" not in type[i]) and (i not in badind):
            ind.append(i)
    elif glacier == 'Kanger':
      # Radar picks to toss, if cleanup is set to true. The indices are for the Kanger
      # 2006-2014 composite product.
      if cleanup==True:
        badind = np.r_[range(28502,28581),range(28672,28756),range(31032,31216),
        		range(29927,30075),range(33768,33847),range(28851,29063),
        		range(39150,39210),range(40471,40550),range(23853,23860),
        		range(36055,36095)] 
        # Find locations where ice is floating so we can toss those out
        floatind = np.where(coords.geoidheight(x2,y2,zs) < flotation.height(coords.geoidheight(x2,y2,zb),rho_i=910.,rho_sw=1025.))[0]
        badind = np.union1d(badind,floatind)
      else: 
        badind = []
        floatind = []
      if year == '2008':
        ind = range(23600,23853)
      elif year == '2009a':
        ind = range(30800,31032)
      elif year == '2009b':
        ind = range(29800,29927)
      elif year == '2012':
        ind = np.arange(36130,36065,-1)
      elif year == '2013':
        ind = np.arange(39370,39210,-1)
      elif year == '2014':
        ind = range(40445,40517)
      elif year == 'all':
        ind = []
        for i in range(0,len(type)):
          if 'ICESat' not in type[i] and (i not in badind):
            ind.append(i)
    else:
      print "Unrecognized CreSIS profile"
  
  
  # Select what reference we want for the elevation  
  if verticaldatum == "geoid":
    zb = coords.geoidheight(x2,y2,zb)
    zs = coords.geoidheight(x2,y2,zs)
  elif verticaldatum == "ellipsoid":
    zb = zb
    zs = zs
  else:
    sys.exit("Unknown vertical datum, exiting...")
    
  return np.column_stack([x2[ind],y2[ind],zb[ind],zs[ind],dates[ind]])

#########################################################################################

def cresis_grid(glacier,verticaldatum='geoid'):

  '''
  x,y,grid = cresis_grid(glacier,verticaldatum='geoid')

  Load gridded CreSIS product for "glacier". The gridded products aren't great, so you 
  probably don't actually want to use this function.
  
  Inputs:
  glacier: glacier name
  verticaldatum: ellipsoid or geoid
  
  Outputs:
  x,y: coordinates for grid
  grid: bed elevation grid
  '''

  # Read the ASCII grid. Why they use ASCII grids, no one knows. Ugh.
  if glacier == 'Helheim':
    print "Using data set Helheim_2006_2014_Composite"
    file = os.path.join(os.getenv("DATA_HOME"),"Bed/Cresis/Helheim/Helheim_2006_2014_Composite/grids/helheim_2006_2014_composite_bottom.txt")
  elif glacier == 'Kanger':
    print "Using data set Kanger_2006_2014_Composite"
    file = os.path.join(os.getenv("DATA_HOME"),"Bed/Cresis/Kanger/Kangerdlugssuaq_2006_2014_Composite/grids/kangerdlugssuaq_2006_2014_composite_bottom.txt")

  # Get file info
  fid = open(file)
  ny = int(fid.readline().split()[-1])
  nx = int(fid.readline().split()[-1])
  xllcorner = float(fid.readline().split()[-1])
  yllcorner = float(fid.readline().split()[-1])
  cellsize = float(fid.readline().split()[-1])
  nodata = float(fid.readline().split()[-1])
  fid.close()
  del fid
  
  # Set up x,y
  x = np.linspace(xllcorner+cellsize/2,xllcorner+cellsize*(nx-1),nx)
  y = np.linspace(yllcorner+cellsize/2,yllcorner+cellsize*(ny-1),ny)
  
  # Get data
  grid = np.flipud(np.loadtxt(file,skiprows=6))
  grid[grid==nodata] = 'NaN'
  
  if verticaldatum == 'geoid':
    x_grid,y_grid = np.meshgrid(x,y)
    geoidheight = coords.geoidheight(x_grid.flatten(),y_grid.flatten(),grid.flatten())
    grid = np.reshape(geoidheight,(ny,nx))
  elif verticaldatum == 'ellipsoid':
    grid = grid
  else:
    print "Unknown datum, defaulting to geoid"

  return x,y,grid

#########################################################################################

def cresis_grid_pts(xpts,ypts,glacier,verticaldatum='geoid',method='linear'):

  '''
  bed = cresis_grid_pts(xpts,ypts,glacier,verticaldatum='geoid')
  
  Get the CreSIS gridded product at specific points.
  
  Inputs: 
  xpts,ypts: points where you want bed elevation
  glacier: glacier name
  verticaldatum: geoid or ellipsoid
  method: interpolation method (same options as RegularGridInterpolator)
  
  Outputs:
  bed: bed elevations at xpts,ypts
  '''

  # Get Cresis grid
  x,y,z = cresis_grid(glacier,verticaldatum)

  # Create interpolation function
  f = scipy.interpolate.RegularGridInterpolator((y,x),z,method=method,bounds_error=False)
  
  # Get points
  bed = f(np.column_stack([ypts,xpts]))

  return bed

#########################################################################################

def morlighem_pts(xpts,ypts,verticaldatum='geoid',method='linear'):
  
  '''
  bed = morlighem_pts(xpts,ypts,glacier,verticaldatum='geoid')
  
  Interpolate morlighem gridded bed to select x,y points.
  
  Inputs:
  xpts,ypts: locations where you want bed elevation
  glacier: glacier name
  verticaldatum: ellipsoid or geoid
  method: interpolation method (same options as RegularGridInterpolator)
  
  Outputs:
  bed: interpolated bed elevations at xpts,ypts
  '''
  
  # Set dimensions of grid to load morlighem bed
  xmin = np.min(xpts)-2.0e3
  xmax = np.max(xpts)+2.0e3
  ymin = np.min(ypts)-2.0e3
  ymax = np.max(ypts)+2.0e3
  
  # Load bed DEM
  xbed,ybed,zbed = morlighem_grid(xmin,xmax,ymin,ymax,verticaldatum=verticaldatum)
  
  f = scipy.interpolate.RegularGridInterpolator((ybed,xbed),zbed,method=method,bounds_error=False)
  zb = f(np.column_stack([ypts,xpts]))
      
  return zb

#########################################################################################
  
def morlighem_grid(xmin=-np.inf,xmax=np.inf,ymin=-np.inf,ymax=np.Inf,verticaldatum='geoid'):

  '''
  xb,yb,bed = morlighem_grid(xmin,xmax,ymin,ymax,verticaldatum='geoid')
  
  Export morlighem gridded bed.
  
  Inputs:
  xmin,xmax,ymin,ymax: extent of desired grid
  verticaldatum: geoid or ellipsoid
  
  Outputs:
  xb,yb: coordinates for grid
  bed: gridded bed
  '''

  # Load Bed DEM
  file = os.path.join(os.getenv("DATA_HOME"),"Bed/Morlighem_2014/MCdataset-2015-04-27.tif")
  [xb,yb,zb]=geotiff.read(file,xmin,xmax,ymin,ymax)
  zb[zb==-9999] = 'NaN'
  
  # Morlighem bed DEM is given as elevation above mean sea level (at geoid). So we need
  # to correct only if we want the ellipsoid height.
  if verticaldatum == "ellipsoid":
    # Load Geoid 
    file = os.path.join(os.getenv("DATA_HOME"),"Bed/Morlighem_2014/geoid.tif")
    [xg,yg,zg]=geotiff.read(file,xmin,xmax,ymin,ymax)
    bed = zb+zg
  elif verticaldatum == "geoid":
    bed = zb
  else:
    sys.exit("Unknown vertical datum, defaulting to geoid height")
  
  return xb,yb,bed

#########################################################################################
  
def smith_grid(glacier,xmin=-np.Inf,xmax=np.Inf,ymin=-np.Inf,ymax=np.Inf,grid='unstructured',model='aniso',smoothing=3,verticaldatum='geoid'):

  '''
               
  xpts,ypts,zpts = smith_grid(glacier, xmin=-np.Inf ,xmax=np.Inf, 
  			ymin=-np.Inf, ymax=np.Inf, grid='unstructured',
  			model='aniso', smoothing=1, verticaldatum='geoid')

  Output Ben Smith's bed, either as an unstructured mesh (grid='unstructured') 
  or as a structured grid (grid='structured').
  
  Inputs:
  glacier: glacier name
  xmin,xmax,ymin,ymax: coordinates for grid
  grid: unstructured or structured grid
  model: anisotropic or isotropic bed interpolation
  smoothing: smoothing level (1-8)
  verticaldatum: geoid or ellipsoid
  
  Outputs:
  xpts,ypts,zpts: either unstructured points or a grid, depending on your choice
  '''

  # Load Ben Smith's interpolated bed. The input data set is in EPSG 3413, relative to WGS84
  # ellipsoid.
  if model == 'iso':
    file = os.path.join(os.getenv("DATA_HOME"),"Bed/Smith_2015/"+glacier+"/iso_models.txt")
  else:
    file = os.path.join(os.getenv("DATA_HOME"),"Bed/Smith_2015/"+glacier+"/aniso_models.txt")

  # Load data
  data = np.loadtxt(file)
  x = data[:,0]
  y = data[:,1]
  z = data[:,1+smoothing]

  if verticaldatum == 'geoid':
    z = coords.geoidheight(x,y,z)
  elif verticaldatum == 'ellipsoid':
    z = z
  else:
    print "Unknown datum, defaulting to geoid"

  # Create structured grid, if desired  
  if grid == 'structured':
    xpts = np.arange(xmin,xmax+100,100.)
    ypts = np.arange(ymin,ymax+100,100.)
    xgrid,ygrid = np.meshgrid(xpts,ypts)
    zpts_flat = scipy.interpolate.griddata((y,x),z,(ygrid.flatten(),xgrid.flatten()),method='linear',fill_value=float('nan'))
    zpts = np.reshape(zpts_flat,(len(ypts),len(xpts)))
  else:
    # Find points that fall within the desired spatial extent
    ind = np.where((x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax))[0]
    xpts = x[ind]
    ypts = y[ind]
    zpts = z[ind]

  return xpts,ypts,zpts

#########################################################################################

def smith_at_pts(xpts,ypts,glacier,model='aniso',smoothing=1,verticaldatum='geoid',method='linear'):

  '''
  zbed_interp = smith_at_pts(xpts,ypts,glacier,model='aniso'
  				,smoothing=1,verticaldatum='geoid',method='linear')
  
  Interpolate Ben Smith's interpolated bed elevations to point locations.
  
  Inputs:
  xpts,ypts: points where you want bed elevation
  glacier: glacier name
  model: anisotropic ('aniso') or isotropic ('iso') interpolation
  smoothing: 1 to 8
  verticaldatum: geoid or ellipsoid
  method: linear or nearest
  
  Outputs:
  zbed_interp: interpolated beds at xpts,ypts
  '''

  # Load smith bed, but don't correct the vertical datum yet
  xbed,ybed,zbed = smith_grid(glacier,grid='unstructured',model=model,smoothing=smoothing,verticaldatum='ellipsoid')

  # Interpolate smith bed to points
  zbed_interp = scipy.interpolate.griddata((xbed,ybed),zbed,(xpts,ypts),method=method)

  # Now fix the vertical datum. We didn't fix it in the call to "smith_grid" because we would 
  # have had to calculate the geoid height at more points. This saves a bit of time.
  if verticaldatum == 'geoid':
    zbed_interp = coords.geoidheight(xpts,ypts,zbed_interp)
  elif verticaldatum == 'ellipsoid':
    zbed_interp = zbed_interp
  else:
    sys.exit("Unknown vertical datum.")


  return zbed_interp
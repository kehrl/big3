import numpy as np
import matplotlib.mlab as mlab
import math
from scipy.interpolate import griddata
from matplotlib.path import Path
from multiprocessing import Pool
import os, shutil

def saveline(DIR,runname,variables):
  '''
  data = saveline(DIR,runname)
  
  Read results from a saveline call in Elmer. This command will return model results on 
  all boundaries. You can pull a single boundary from "saveline_boundary."
  
  Inputs:
  DIR: directory name
  runname: model run name 
  
  '''
  
  # Find all partitions of the SaveLine results
  names=os.listdir(DIR)
  files=[]
  datanames = []
  for name in names:
    if name.startswith(runname) and not name.endswith('names') and ('.dat' in name):
      files.append(name)
    if name.startswith(runname) and name.endswith('names'):
      datanames.append(name)
  
  # Find the variable names in the ".name" file and initialize them   
  fid = open(DIR+datanames[0],"r")
  lines=fid.readlines()
  for i in range(0,4):
    lines.remove(lines[0])
  varnames=[]
  types = []
  columns = []
  n = 0
  for line in lines:
    if line[5:-1] == 'coordinate 1':
      varnames.append('x')
      types.append(np.float64) 
      columns.append(n)  
    elif line[5:-1] == 'coordinate 2':
      varnames.append('y')
      types.append(np.float64)  
      columns.append(n)        
    elif line[5:-1] == 'coordinate 3':
      varnames.append('z')
      types.append(np.float64)  
      columns.append(n)  
    elif line[5:-1] == 'Node index':
      varnames.append('Node Number')  
      types.append(np.int64)    
      columns.append(n) 
    elif line[5:-1] == 'Boundary condition':
      varnames.append('BC')  
      types.append(np.int64)    
      columns.append(n) 
    for variable in variables:
      if line[5:-1].startswith(variable):
        varnames.append(line[5:-1])
        types.append(np.float64)
        columns.append(n)     
    n=n+1

  # Adding some variables for ease of use   
  if ('vsurfini 1') in varnames:
    varnames.append('vsurfini')
    types.append(np.float64)    
  if ('velocity 1') in varnames:
    varnames.append('velocity')
    types.append(np.float64) 
  if ('beta' in varnames) and ('velocity 1' in varnames):
    varnames.append('taub')
    types.append(np.float64) 
  
  # Get number of points 
  n=0
  for file in files:
    fid = open(DIR+file,"r")
    lines = fid.readlines()
    n = n+len(lines)
  data = np.empty(n,dtype = zip(varnames, types))  
      
  # Read in data
  n = 0
  for file in files:
    datafile = np.loadtxt(DIR+file)
    for i in range(0,len(columns)):
      data[varnames[i]][n:len(datafile[:,0])+n] = datafile[:,columns[i]]
    n = n+len(datafile[:,0])
  
  # Calculate velocity magnitude for 3D simulations
  if ('velocity 3' in varnames):
    data['velocity'] = np.sqrt(data['velocity 1']**2+data['velocity 2']**2)
  elif ('velocity 2' in varnames):
    data['velocity'] = data['velocity 1']
  if ('vsurfini 2' in varnames):
    data['vsurfini'] = np.sqrt(data['vsurfini 1']**2+data['vsurfini 2']**2)
  elif ('vsurfini 1' in varnames):
    data['vsurfini'] = data['vsurfini 1']
   
  # Calculate basal shear stress, if it is an inversion and we have a variable beta
  if 'taub' in varnames:
    data['taub']= (data['beta']**2*data['velocity'])
         
  return data
  
def saveline_boundary(DIR,runname,bound,variables):
  '''
  subset = saveline_boundary(DIR,runname,bound)
  
  Get model results for a single boundary saved with the Elmer saveline solver. 
  
  Inputs:
  DIR: directory name for model results
  runname: name of model run name for the .dat files
  bound: boundary number
  '''
  
  # Load data
  data = saveline(DIR,runname,variables)
  
  # Get desired indices  
  ind = np.where(data['BC'] == bound)[0]
  
  points = []
  for i in ind:
    points.append(data[i])
  subset = np.asarray(points)
       
  return subset

def grid3d(data,variable,holes,extent,dx=50):

  '''
  x,y,masked = grid3d(data,variable,holes,extent,dx=50)
  
  Get grid from saveline_boundary results for a given variable. The result is masked using
  the mesh extent and holes.
  
  Inputs:
  data: data from saveline_boundary
  variable: name of variable that you want to grid
  holes: list of mesh holes or holes=[]
  extent: mesh extent
  dx: grid spacing for output grid
  
  Outputs:
  x: list of x coords for grid
  y: list of y coords for grid
  masked: masked grid
  '''
  
  xmin=math.floor(min(data['x'])/100)*100
  xmax=math.ceil(max(data['x'])/100)*100
  ymin=math.floor(min(data['y'])/100)*100
  ymax=math.ceil(max(data['y'])/100)*100

  x=np.linspace(xmin,xmax,(xmax-xmin)/dx+1)
  y=np.linspace(ymin,ymax,(ymax-ymin)/dx+1)
  
  xx,yy=np.meshgrid(x,y)
  
  xy_i=np.column_stack([data['x'],data['y']])
  
  zz=griddata((data['x'],data['y']),data[variable],(xx,yy),method='linear')
  
  # Create mask
  nx = len(x)
  ny = len(y)
  xR = xx.flatten()
  yR = yy.flatten()
  mask = np.ones((ny,nx),dtype="bool")
  
  exterior = Path(np.column_stack((extent[:,0],extent[:,1])))
  pts = exterior.contains_points(np.column_stack((xR,yR)))
  bool = ~(np.reshape(pts,(ny,nx)))
  mask[bool==0] = 0
  
  # Mask out points inside holes
  for i in range(0,len(holes)):
    hole = Path(np.column_stack((holes[i][:,0],holes[i][:,1])))
    pts = hole.contains_points(np.column_stack((xR,yR)))
    bool = (np.reshape(pts,(ny,nx)))
    mask[bool==1] = 1
    
  masked = np.ma.masked_array(zz,mask)

  return x,y,masked
  
def result_file(variables, mesh_dir, result_fn, maps=[], debug=False, **kwargs):
  '''
  Pull variables from a .result file.
  
  modelresult = result_file(variables, result_fn, maps=[], debug=False, **kwargs)
  
  Inputs:
  variables: name of variables to export
  result_fn: name of result file
  '''
  
  parts = None
  if 'partitions' in kwargs:
    parts = kwargs['partitions']
    if os.path.isdir(mesh_dir + '/partitioning.' + str(parts)):
      parts_dir = mesh_dir + '/partitioning.' + str(parts)
      print str(parts) + ' partition mesh found in ' + parts_dir
    elif os.path.isdir(mesh_dir + '/../partitioning.' + str(parts)):
      parts_dir = mesh_dir + '/../partitioning.' + str(parts)
      print str(parts) + ' partition mesh found in ' + parts_dir
    else:
      print 'There is no mesh with the given partitions, trying to find with other partitioning'
      parts = None
  # try to figure out the number of partitions. Capped at 256.
  if parts is None:
    k = 2
    while k < 257:
      # Try result directory
      if os.path.isdir(mesh_dir + '/partitioning.' + str(k)):
        parts = k
        parts_dir = mesh_dir + '/partitioning.' + str(k)
        print str(parts) + ' partition mesh found in ' + parts_dir
        break
      else:
        k = k + 1
    if not k < 257:
      print 'No mesh partitions found'
      return None
  
  varnames = ['Node Number', 'x', 'y', 'z']
  types = [np.int64, np.float64, np.float64, np.float64]
  for i in range(len(variables)):  # Prepare to receive data about variables
    types.append(np.float64)  # data must be a float
    # node and coordinates fixed, name the desired variables
    varnames.append(variables[i])
  if maps is not None:
    for mp in maps:
      types.append(np.float64)
      varnames.append(variables[mp[1]] + mp[0].split('.')[-1])
  
  if debug:
	  print 'Varnames are ', varnames
  
  # Make a pool of workers to handle each partition separately
  pool = Pool()
  results = [pool.apply_async(result_file_partition, (i, parts_dir, varnames, types, variables, maps, mesh_dir, result_fn, debug)) for i in range(parts)]
  
  # Get, stack, and sort the results
  return np.sort(np.hstack(tuple([result.get() for result in results])), order=['x', 'y', 'z'])

def result_file_partition(i, parts_dir, varnames, types, variables, maps, mesh_dir, result_fn, debug=False):

  if debug and not i==0:
    debug=False
  
  file_len = bufcount(parts_dir + '/part.' + str(i + 1) + '.nodes')
  data = np.empty(file_len, dtype = zip(varnames, types))

  if debug:
    print 'For partition ', i, ' file length is ', file_len
	
  with open(parts_dir + '/part.' + str(i + 1) + '.nodes') as fin:
    points = fin.readlines()
	
  for k, pt in enumerate(points):
    points[k] = tuple([map(float, pt.split())[j] for j in [
                0, 2, 3, 4]] + [0 for j in range(len(variables) + len(maps))])
  data[:] = points
    
  if debug:
    print 'After going through partition files, data is:'
    print data

  if os.path.exists(result_fn + '.' + str(i)):
    rfn = result_fn + '.' + str(i)
  else:
    rfn = mesh_dir+result_fn + '.' + str(i)
  if not os.path.exists(rfn):
    # allows for use of dummy file name to quickly just get coordinates
    if os.path.exists(mesh_dir+result_fn + str(i)):
      found = True
      rfn = mesh_dir + result_fn + str(i)
      if debug:
        print 'Found result file ', rfn
    else:
      print(rfn + ' not found')
      print('Assuming dummy filename')
      found = False
  else:
    if debug:
      print 'Found result file ', rfn
    found = True
	
  if found:
    with open(rfn) as f:
      dat = f.readlines()
	  
    # get all the indices for each of the variables
    var_indices = [None for L in variables]
    for L, var in enumerate(variables):
      if debug:
        print 'For variable ', var, ' indices are ', [k for k, x in enumerate(dat) if var in x and var==x.rstrip('\n')]
      try:
        var_indices[L] = [k for k, x in enumerate(dat) if var in x and var==x.rstrip('\n')][-1]
      except IndexError:
        raise IndexError('No lines found that exactly match variable '+var)
      if debug:
        print 'For variable ', var, ' index is ', var_indices[L]
      try:
        data[var][:] = map(
           float, dat[(var_indices[L] + 2):(var_indices[L] + 2 + file_len)])
        if debug:
          print 'After variable ', var, ' data is:'
          print data
        if maps is not None:
          for mp in maps:
            if mp[1] == L:
              data[var + mp[0].split('.')[-1]][:] = map(eval(mp[0]), map(float, dat[var_indices[L] + 2:var_indices[L] + 2 + file_len]))
      except ValueError:
        if i == 0:
          print 'Variable \"' + var + '\" appears to be wonky, attempting to fetch regardless'
        var_indices[L] = var_indices[L] + file_len
        try:
          data[var][:] = map(
                        float, dat[var_indices[L] + 2:var_indices[L] + 2 + file_len])
          if maps is not None:
            for mp in maps:
              if mp[1] == L:
                data[var + mp[0].split('.')[-1]][:] = map(
                        eval(mp[0]), map(float, dat[var_indices[L] + 2:var_indices[L] + 2 + file_len]))
        except:
          print 'Failed to fetch'
          return None

  return data 
  
def bufcount(filename):
    '''
    Gets the length of a file
    '''
    
    f = open(filename)
    lines = 0
    buf_size = 1024 * 1024
    read_f = f.read  # loop optimization

    buf = read_f(buf_size)
    while buf:
        lines += buf.count('\n')
        buf = read_f(buf_size)

    return lines
  
def pvtu_file(file,variables):

  from paraview import numpy_support, simple  
  
  # Load vtu file
  reader = simple.XMLPartitionedUnstructuredGridReader(FileName=file)
  vtudata = simple.servermanager.Fetch(reader)
  
  # Get number of nodes
  n = vtudata.GetNumberOfPoints()
  
  # Set up output array
  varnames = ['Node Number','x','y','z']
  types = [np.int64, np.float64, np.float64, np.float64]
  for var in variables:
    if var == 'velocity':
      opts = ['velocity 1','velocity 2','velocity 3','velocity']
      for opt in opts:
        varnames.append(opt)
        types.append(np.float64)
    elif var == 'vsurfini':
      opts = ['vsurfini 1','vsurfini 2','vsurfini']
      for opt in opts:
        varnames.append(opt)
        types.append(np.float64)
    else:
      types.append(np.float64)
      varnames.append(var)
  if ('beta' in varnames) and ('velocity 1' in varnames):
    varnames.append('taub')
    types.append(np.float64) 
  data = np.empty(n,dtype = zip(varnames, types))  

  # Get coordinates
  x = np.zeros(n)
  y = np.zeros(n)
  z = np.zeros(n)
  for i in range(0,n):
    x[i],y[i],z[i] = vtudata.GetPoint(i)
  
  data['x'][:] = x
  data['y'][:] = y
  data['z'][:] = z
  
  # Get variables
  for var in variables:
    if var == 'velocity':
      velocity = numpy_support.vtk_to_numpy(vtudata.GetPointData().GetArray(var))
      data['velocity 1'][:] = velocity[:,0]
      data['velocity 2'][:] = velocity[:,1]
      data['velocity 3'][:] = velocity[:,2]
      data['velocity'][:] = np.sqrt(data['velocity 1']**2+data['velocity 2']**2)
    elif var == 'vsurfini': 
      data['vsurfini 1'][:] = numpy_support.vtk_to_numpy(vtudata.GetPointData().GetArray('vsurfini 1'))
      data['vsurfini 2'][:] = numpy_support.vtk_to_numpy(vtudata.GetPointData().GetArray('vsurfini 2'))
      data['vsurfini'][:] = np.sqrt(data['vsurfini 1']**2+data['vsurfini 2']**2)
    else:
      data[var][:] = numpy_support.vtk_to_numpy(vtudata.GetPointData().GetArray(var))
  if 'taub' in varnames:
    data['taub'][:] = data['beta']**2*np.sqrt(data['velocity 1']**2+data['velocity 2']**2) 
 
  # Some of the nodes are shared/repeated between partitions, so we need to remove those
  var1,ind1 = np.unique(x,return_index = True)
  var2,ind2 = np.unique(y,return_index = True)
  var3,ind3 = np.unique(z,return_index = True)
  
  ind = np.union1d(np.union1d(ind1,ind2),ind3)
  
  data = data[ind]

  return data

def grid_to_cross_section(data,xf,yf,df,variable,dz=20.,method='linear'):

  '''
  Make cross section from model
    
  Inputs:
  data: structured array with required fields x,y,z
  xf,yf,df: points for cross-section
  method: presently this module only takes the closest point from the model
  maxdist: how far away it is acceptable to take points in the nearest neighbor interpolation
  '''

  import numpy as np
  from scipy.interpolate import griddata

  # Do a quick check that we have the necessary variables
  if not set(('Node Number', 'x', 'y')).issubset(set(data.dtype.names)):
    print 'Not a valid dataset, does not contain Node Number, x, and y variables, returning None'
    return None

  # This approach is a bit involved, but it's probably the best way to do this. 
  # First we find the bed and surface for the cross-section, then we interpolate onto 
  # a grid that extends across the full range of z values. Finally we mask out values 
  # outside of the model domain.
  bed = values_in_layer(data,'bed')
  sur = values_in_layer(data,'surf') 
  
  zbed = griddata((bed['x'],bed['y']),bed['z'],(xf,yf),method='linear')
  zsur = griddata((sur['x'],sur['y']),sur['z'],(xf,yf),method='linear')
    
  zmin = np.round(np.nanmin(zbed)-100,decimals=-2)
  zmax = np.round(np.nanmax(zsur)+100,decimals=-2)
  z = np.arange(zmin,zmax,dz)
    
  # Make grid for interpolation
  dgrid,zgrid = np.meshgrid(df,z)
  dgrid_flat = dgrid.flatten()
  xgrid_flat = np.interp(dgrid_flat,df,xf)
  ygrid_flat = np.interp(dgrid_flat,df,yf)
  zgrid_flat = zgrid.flatten()
    
  outflat = griddata((data['x'],data['y'],data['z']),data[variable],(xgrid_flat,ygrid_flat,zgrid_flat),method='linear')
  output = np.reshape(outflat,(len(z),len(df)))

  # Delete points that fall outside of the model domain
  zsur_int = np.column_stack([df,zsur])
  zbed_int = np.column_stack([df,zbed])
  domain = np.row_stack([zsur_int,np.flipud(zbed_int)])
  ind = np.where(np.isnan(domain[:,1]))[0]
  domain = np.delete(domain,ind,axis=0)
  domainpath = Path(domain)
  pts = domainpath.contains_points(np.column_stack((dgrid_flat,zgrid_flat)))
  bool = ~(np.reshape(pts,(len(z),len(df))))

  output[bool==1] = float('nan')

  
  return df,z,output


# def cross_section_grid(cross_section,variable,dz=50.,dd=100.,method=method):
# 
#   import scipy.interpolate
# 
#   # Set up grid
#   zmin = np.round(np.min(cross_section['z'])-100,decimals=-2)
#   zmax = np.round(np.max(cross_section['z'])+100,decimals=-2)
#   z = np.arange(zmin,zmax,dz)
# 
#   dmin = np.round(np.min(cross_section['dist'])-100,decimals=-2)
#   dmax = np.round(np.max(cross_section['dist'])+100,decimals=-2)
#   d = np.arange(dmin,dmax,dd)
# 
#   dgrid,zgrid = np.meshgrid(d,z)
# 
#   # Now interpolate points
#   outflat = scipy.interpolate.griddata((cross_section['dist'],cross_section['z']),cross_section[variable],(dgrid.flatten(),zgrid.flatten()),method=method)
#   output = np.reshape(outflat,(len(z),len(d)))
#   
#   # Get glacier geometry so we can mask out values outside of the glacier
#   dgeo = []
#   zsur = []
#   zbed = []
#   for val in np.unique(cross_section['dist']):
#     points = cross_section[cross_section['dist'] == val]
#     zsur.append(np.max(points['z']))
#     zbed.append(np.min(points['z']))
#     dgeo.append(val)
#     
#   zsur_int = np.column_stack([d,np.interp(d,dgeo,zsur)])
#   zbed_int = np.column_stack([d,np.interp(d,dgeo,zbed)])
# 
#   exterior = Path(np.row_stack([zsur_int,np.flipud(zbed_int)]))
#   pts = exterior.contains_points(np.column_stack((dgrid.flatten(),zgrid.flatten())))
#   bool = ~(np.reshape(pts,(len(d),len(z))))
#   mask = np.ones((len(d),len(z)),dtype="bool")
#   mask[bool==0] = 0
#   
#   masked = np.ma.masked_array(output,mask)
#   
#   return masked

def depth_averaged_variables(data):
  '''
  Pull depth averaged values from pointwise data list
    
  Inputs:
  data: structured array with required fields x,y,z
  '''

  # Do a quick check that we have the necessary variables
  if not set(('Node Number', 'x', 'y')).issubset(set(data.dtype.names)):
    print 'Not a valid dataset, does not contain Node Number, x, and y variables, returning None'
    return None
  # Try to get things right with what variables we are returning
  varnames = list(data.dtype.names)
  varnames.remove('Node Number')
  
  # Do the actual sorting
  points = []
  for x_val in np.unique(data['x']):
    x_points = data[data['x'] == x_val]
    for y_val in np.unique(x_points['y']):
      y_points = x_points[x_points['y'] == y_val]
      pt = y_points[np.argmin(y_points['z'])]
      sorted_list = np.sort(y_points, order='z')
      weights = (np.pad(sorted_list['z'], (1,0), 'edge')[0:-1]-np.pad(sorted_list['z'], (0,1), 'edge')[1:] ) / (sorted_list['z'][0] - sorted_list['z'][-1]) / 2.0
      for var in varnames:
        # The following will assume that the points are evenly spaced in the vertical direction
        # pt[var] = np.sum(y_points[var]) / len(y_points[var])
        # Instead, do things correctly to weight by percentage of the column
        pt[var] = np.sum(sorted_list[var] * weights)
      points.append(pt)
  points = np.asarray(points)
  
  return points[varnames]

def values_in_layer(data,layer='surface'):
  '''
  Pull values at either the surface or bed
    
  Inputs:
  data: structured array with required fields x,y,z
  '''

  # Do a quick check that we have the necessary variables
  if not set(('Node Number', 'x', 'y')).issubset(set(data.dtype.names)):
    print 'Not a valid dataset, does not contain Node Number, x, and y variables, returning None'
    return None
  # Try to get things right with what variables we are returning
  varnames = list(data.dtype.names)
  
  if layer.startswith('surf'):
    ind = -1
  elif layer == 'bed':
    ind = 0
  
  # Do the actual sorting
  points = []
  for x_val in np.unique(data['x']):
    x_points = data[data['x'] == x_val]
    for y_val in np.unique(x_points['y']):
      y_points = x_points[x_points['y'] == y_val]
      pt = y_points[np.argmin(y_points['z'])]
      sorted_list = np.sort(y_points, order='z')
      for var in varnames:
        # The following will assume that the points are evenly spaced in the vertical direction
        # pt[var] = np.sum(y_points[var]) / len(y_points[var])
        # Instead, do things correctly to weight by percentage of the column
        pt[var] = sorted_list[var][ind]
      points.append(pt)
  points = np.asarray(points)
  
  return points[varnames]

def values_in_column(data,x,y):
  '''
  Pull data at a particular point
    
  Inputs:
  x,y: list of desired coordinates
  data: structured array with required fields x,y,z
  '''

  # Do a quick check that we have the necessary variables
  if not set(('Node Number', 'x', 'y')).issubset(set(data.dtype.names)):
    print 'Not a valid dataset, does not contain Node Number, x, and y variables, returning None'
    return None

  varnames = list(data.dtype.names)
  
  # Do the actual sorting
  points = []
  try:
    n = len(x)
  except:
    n = 1
    x2 = np.zeros(1)
    y2 = np.zeros(1)
    x2[0] = x
    y2[0] = y
    x = np.array(x2)
    y = np.array(y2)
    del x2,y2
  print "Getting profiles for",n," points"
  
  # Get unique coordinates for extruded mesh
  junk,ind1 = np.unique(data['x'],return_index = True)
  junk,ind2 = np.unique(data['y'],return_index = True)
  ind = np.union1d(ind1,ind2)
  xnode = data['x'][ind]
  ynode = data['y'][ind]
  
  # Grab closest coordinate to desired point  
  for i in range(0,n):
    minind = np.argmin(np.sqrt((x[i]-xnode)**2+(y[i]-ynode)**2))
    print "distance is",np.min(np.sqrt((x[i]-xnode)**2+(y[i]-ynode)**2))
    ind = np.where((data['x'] == xnode[minind]) & (data['y'] == ynode[minind]))[0]
    sortind = ind[np.argsort(data['z'][ind])]
    for j in sortind:
      pt = data[j]
      pt['Node Number'] = i
      points.append(pt)
  points = np.asarray(points)
  
  return points

def grid_to_flowline_surface(surf,xf,yf):
  '''
  flow = grid_to_flowline(data,x,y)
  
  Takes results and interpolates them to a flowline.
  
  Inputs:
  data: model results from saveline_boundary
  x,y: coordinates for flowline
  
  Outputs:
  flow: same as "data", except now only for values interpolated to the flowline
  '''
  
  import numpy as np
  from scipy.interpolate import griddata


  # Do a quick check that we have the necessary variables
  if not set(('Node Number', 'x', 'y')).issubset(set(surf.dtype.names)):
    print 'Not a valid dataset, does not contain Node Number, x, and y variables, returning None'
    return None
  # Try to get things right with what variables we are returning
  varnames = list(surf.dtype.names)
  varnames.remove('Node Number') 
  types = []
  for var in varnames:
    types.append(np.float64)
 
  # Do the actual sorting
  points = np.empty(len(xf), dtype = zip(varnames,types))
  nonnan = np.intersect1d(np.where(~(np.isnan(surf['x'])))[0],np.where(~(np.isnan(surf['y'])))[0])
  for var in varnames:
    points[var] = griddata((surf['y'][nonnan],surf['x'][nonnan]),surf[var][nonnan],(yf,xf),method='linear')
  points = np.asarray(points)
  
  return points

# def grid_to_flowline_old(data,x,y):
#   '''
#   flow = grid_to_flowline(data,x,y)
#   
#   Takes results from saveline_boundary and interpolates them to a flowline.
#   
#   Inputs:
#   data: model results from saveline_boundary
#   x,y: coordinates for flowline
#   
#   Outputs:
#   flow: same as "data", except now only for values interpolated to the flowline
#   '''
#   
#   import numpy as np
#   from scipy.interpolate import griddata
#   
#   gridx=data['coord1']
#   gridy=data['coord2']
#   
#   variables=data.keys()
#   
#   flow={}
#   for variable in variables:
#     flow[variable]=griddata(np.column_stack([gridx,gridy]),data[variable],np.column_stack([x,y]),method='linear')
#   
#   return flow
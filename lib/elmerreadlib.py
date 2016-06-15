import numpy as np
import matplotlib.mlab as mlab
import math
from scipy.interpolate import griddata
from matplotlib.path import Path
from multiprocessing import Pool
import os

def saveline(DIR,runname):
  '''
  data = saveline(DIR,runname)
  
  Read results from a saveline call in Elmer. This command will return model results on 
  all boundaries. You can pull a single boundary from "saveline_boundary."
  
  Inputs:
  DIR: directory name
  runname: model run name 
  
  '''

  import os
  import math
  import numpy as np
  
  os.chdir(DIR)
  # Find all partitions of the SaveLine results
  names=os.listdir(DIR)
  files=[]
  varnames=[]
  for name in names:
    if name.startswith(runname) and not name.endswith('names') and ('.dat' in name):
      files.append(name)
    if name.startswith(runname) and name.endswith('names'):
      varnames.append(name)

  # Find the variable names in the ".name" file and initialize them 
  fid = open(varnames[0],"r")
  os.rename(varnames[0],runname+".dat.names")

  lines=fid.readlines()
  for i in range(0,4):
    lines.remove(lines[0])
  variables=[]
  for line in lines:
    p=line.split()
    variable=''
    for i in range(1,len(p)):
      if p[i]=='coordinate':
        p[i]='coord'
      elif p[i]=='step':
        p[i]=''
      elif p[i]=='Time':
        p[i]='timestep'
      elif p[i]=='Iteration':
        p[i]='iterstep'
      elif p[i]=='index':
        p[i]=''
      elif p[i]=='condition':
        p[i]=''
      elif p[i]=='Boundary':
        p[i]='bound'
      elif p[i]=='Node':
        p[i]='node'
      elif p[i]=='velocity':
        p[i]='vel'
      else:
        pass
      variable=variable+p[i]
    variables.append(variable)
      
  R=len(variables)
  
  data={}
  for i in range(0,R):
    data[variables[i]]=[]
  
  # Read in data
  for file in files:
    fid = open(file,"r")
    lines = fid.readlines()
    for line in lines:
      p = line.split()
      if len(p) == R:
        for i in range(0,R):
            data[variables[i]].append(float(p[i]))
      else: 
        print "The number of variables and length of the file do not agree. \n Compare the initialized variables to reported values in the .names file" 
    fid.close()
  
  # Calculate velocity magnitude for 3D simulations
  if ('velod4' in variables) or ('velocityb4' in variables):
    data['velmes']=[]
    data['vel']=[]
    for i in range(0,len(data['vel1'])):
      data['vel'].append(math.sqrt(data['vel1'][i]**2+data['vel2'][i]**2))
      data['velmes'].append(math.sqrt(data['vsurfini1'][i]**2+data['vsurfini2'][i]**2))
  
  # Calculate basal shear stress
  if 'beta' in variables:
    data['taub']=[]
    if ('velod4' in variables) or ('velocityb4' in variables):
      for i in range(0,len(data['vel1'])):
        data['taub'].append(data['beta'][i]**2*data['vel'][i])
    else:
      for i in range(0,len(data['vel1'])):
        data['taub'].append(data['beta'][i]**2*data['vel1'][i])
    
  # Save model results for future use
  if len(files) > 1:
    fid = open(DIR+runname+".dat","w")
    for i in range(0,len(data[variables[0]])):
      for variable in variables:
        fid.write('{} '.format(data[variable][i]))
      fid.write('\n')
    fid.close()
    #for i in range(0,len(varnames)):
      #if varnames[i] is not runname+".dat.names":
        #os.remove(varnames[i])
      #if files[i] is not runname+".dat": 
        #os.remove(files[i])
  
  for variable in data.keys():
    data[variable] = np.array(data[variable])
         
  return data
  
def saveline_boundary(DIR,runname,bound):
  '''
  subset = saveline_boundary(DIR,runname,bound)
  
  Get model results for a single boundary saved with the Elmer saveline solver. 
  
  Inputs:
  DIR: directory name for model results
  runname: name of model run name for the .dat files
  bound: boundary number
  '''
  
  data=saveline(DIR,runname)
    
  subset={}
  variables=data.keys()
  R1=len(variables)
  R2=len(data['node'])
  for i in range(0,R1):
    subset[variables[i]]=[]
  
  ind = np.where(data['bound']==float(bound))
  
  subset={}
  for variable in variables:
    subset[variable]=data[variable][ind]
  
  #for i in range(0,R2):
  #  if data['bound'][i]==float(bound):
  #    for j in range(0,R1):
  #      try:
  #        subset[variables[j]].append(data[variables[j]][i])
  #      except:
  #        subset[variables[j]].append([])
       

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
  
  xmin=math.floor(min(data['coord1'])/100)*100
  xmax=math.ceil(max(data['coord1'])/100)*100
  ymin=math.floor(min(data['coord2'])/100)*100
  ymax=math.ceil(max(data['coord2'])/100)*100

  x=np.linspace(xmin,xmax,(xmax-xmin)/dx+1)
  y=np.linspace(ymin,ymax,(ymax-ymin)/dx+1)
  
  xx,yy=np.meshgrid(x,y)
  
  xy_i=np.column_stack([data['coord1'],data['coord2']])
  
  zz=griddata((data['coord1'],data['coord2']),data[variable],(xx,yy),method='linear')
  
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

def grid_to_flowline(data,x,y):
  '''
  flow = grid_to_flowline(data,x,y)
  
  Takes results from saveline_boundary and interpolates them to a flowline.
  
  Inputs:
  data: model results from saveline_boundary
  x,y: coordinates for flowline
  
  Outputs:
  flow: same as "data", except now only for values interpolated to the flowline
  '''
  
  import numpy as np
  from scipy.interpolate import griddata
  
  gridx=data['coord1']
  gridy=data['coord2']
  
  variables=data.keys()
  
  flow={}
  for variable in variables:
    flow[variable]=griddata(np.column_stack([gridx,gridy]),data[variable],np.column_stack([x,y]),method='linear')
  
  return flow
  
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
      opts = ['velocity 1','velocity 2','velocity 3']
      for opt in opts:
        varnames.append(opt)
        types.append(np.float64)
    else:
      types.append(np.float64)
      varnames.append(var)
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
    else:
      data[var][:] = numpy_support.vtk_to_numpy(vtudata.GetPointData().GetArray(var))
  
  # Some of the nodes are shared/repeated between partitions, so we need to remove those
  var1,ind1 = np.unique(x,return_index = True)
  var2,ind2 = np.unique(y,return_index = True)
  var3,ind3 = np.unique(z,return_index = True)
  
  ind = np.union1d(np.union1d(ind1,ind2),ind3)
  
  data = data[ind]

  return data

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

def values_in_column(x,y,data):
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
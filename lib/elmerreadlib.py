import numpy as np
import matplotlib.mlab as mlab
import math
from scipy.interpolate import griddata
from matplotlib.path import Path
import shapely
from multiprocessing import Pool
import os, shutil, sys

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
    if name.startswith(runname) and not('names' in name) and ('.dat' in name):
      files.append(name)
    if name.startswith(runname) and name.endswith('names'):
      datanames.append(name)
  
  # Find the variable names in the ".name" file and initialize them   
  fid = open(DIR+datanames[0],"r")
  lines=fid.readlines()
  varnames=[]
  types = []
  columns = []
  for line in lines:
    p = line.split()
    if ' '.join(p[1:]) == 'coordinate 1':
      varnames.append('x')
      types.append(np.float64) 
      columns.append(int(p[0][0:-1])-1)  
    elif ' '.join(p[1:]) == 'coordinate 2':
      varnames.append('y')
      types.append(np.float64)  
      columns.append(int(p[0][0:-1])-1)        
    elif ' '.join(p[1:]) == 'coordinate 3':
      varnames.append('z')
      types.append(np.float64)  
      columns.append(int(p[0][0:-1])-1)  
    elif ' '.join(p[1:]) == 'Node index':
      varnames.append('Node Number')  
      types.append(np.int64)    
      columns.append(int(p[0][0:-1])-1) 
    elif ' '.join(p[1:]) == 'Boundary condition':
      varnames.append('BC')  
      types.append(np.int64)    
      columns.append(int(p[0][0:-1])-1)
    for variable in variables:
      if ' '.join(p[1:]).startswith(variable):
        varnames.append(' '.join(p[1:]))
        types.append(np.float64)
        columns.append(int(p[0][0:-1])-1)     

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
    varnames.append('betasquared')
    types.append(np.float64)
  if ('velocity' in varnames) and ('vsurfini' in varnames): 
    varnames.append('misfit')
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
  #elif ('velocity 2' in varnames):
  #  data['velocity'] = data['velocity 1']
  if ('vsurfini 2' in varnames):
    data['vsurfini'] = np.sqrt(data['vsurfini 1']**2+data['vsurfini 2']**2)
  #elif ('vsurfini 1' in varnames):
  #  data['vsurfini'] = data['vsurfini 1']
   
  # Calculate basal shear stress, if it is an inversion and we have a variable beta
  if 'taub' in varnames:
    data['taub'] = (data['beta']**2.0*data['velocity'])

  # Calculate betasquared, if it is an inversion and we have a variable beta
  if 'betasquared' in varnames:
    data['betasquared'] = data['beta']**2

  # Calculate misfit
  if 'misfit' in varnames:
    data['misfit'] = data['velocity']-data['vsurfini']
         
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

def grid3d(data,variable,holes=[],extent=[],dx=50,clip_boundary=0):

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
  clip_boundary: mask grid elements that are within a certain distance of the mesh boundaries
  
  Outputs:
  x: list of x coords for grid
  y: list of y coords for grid
  masked: masked grid
  '''
  
  xmin = math.floor(min(data['x'])/100)*100
  xmax = math.ceil(max(data['x'])/100)*100
  ymin = math.floor(min(data['y'])/100)*100
  ymax = math.ceil(max(data['y'])/100)*100

  x = np.linspace(xmin,xmax,(xmax-xmin)/dx+1)
  y = np.linspace(ymin,ymax,(ymax-ymin)/dx+1)
  
  xx,yy = np.meshgrid(x,y)
  
  xy_i = np.column_stack([data['x'],data['y']])
  
  zz = griddata((data['x'],data['y']),data[variable],(xx,yy),method='linear')
  
  # Create mask
  nx = len(x)
  ny = len(y)
  xR = xx.flatten()
  yR = yy.flatten()
  mask = np.ones((ny,nx),dtype="bool")
  
  if len(extent) > 0:
    exterior = Path(np.column_stack((extent[:,0],extent[:,1])))
    pts = exterior.contains_points(np.column_stack((xR,yR)))
    bool = ~(np.reshape(pts,(ny,nx)))
    mask[bool==0] = 0
    
    if clip_boundary > 0.0:
      line = shapely.geometry.LineString(np.row_stack([extent,extent[0,:]]))
      for i in range(0,nx):
        for j in range(0,ny):
          if line.distance(shapely.geometry.Point(x[i],y[j])) <= clip_boundary:    
            mask[j,i] = 1
  
  if len(holes) > 0:
    # Mask out points inside holes
    for i in range(0,len(holes)):
      hole = Path(np.column_stack((holes[i][:,0],holes[i][:,1])))
      pts = hole.contains_points(np.column_stack((xR,yR)))
      bool = (np.reshape(pts,(ny,nx)))
      mask[bool==1] = 1

      if clip_boundary > 0.0:
        for k in range(0,nx):
          for j in range(0,ny):
            if np.min(np.sqrt((x[k]-holes[i][:,0])**2 + (y[j]-holes[i][:,1])**2)) <= clip_boundary:
              mask[j,k] = 1
  
  if (len(holes) > 0) or (len(extent) > 0):  
    zz[mask] = float('nan')
    masked = np.ma.masked_array(zz,mask)
  else:
    masked = zz

  return x,y,masked
  
def result_file(mesh_dir, result_fn, variables, maps=[], debug=False, **kwargs):
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
  data = np.sort(np.hstack(tuple([result.get() for result in results])), order=['x', 'y', 'z'])
  
  # Add some additional variables for ease of use
  if ('ssavelocity 1' in varnames) and ('ssavelocity 2' in varnames):
    varnames.append('ssavelocity')
    types.append(np.float64)
  if ('beta' in varnames) and ('ssavelocity' in varnames): 
    varnames.append('taub')
    types.append(np.float64)
    varnames.append('betasquared')
    types.append(np.float64)
  if ('vsurfini 1' in varnames) and ('vsurfini 2' in varnames):
    varnames.append('vsurfini')
    types.append(np.float64) 
  
  data2 = np.empty(len(data),dtype = zip(varnames,types))
  for varname in varnames:
    if varname == 'ssavelocity':
      data2[varname] = np.sqrt(data['ssavelocity 1']**2+data['ssavelocity 2']**2)
    elif varname == 'taub':
      data2[varname] = (data['beta']*data['beta'])*np.sqrt(data['ssavelocity 1']**2+data['ssavelocity 2']**2)
    elif varname == 'betasquared':
      data2[varname] = data['beta']**2
    elif varname == 'vsurfini':
      data2[varname] = np.sqrt(data['vsurfini 1']**2+data['vsurfini 2']**2)
    else:
      data2[varname] = data[varname]
      
  return data2

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
  
def pvtu_file(file,variables,reader='none',returnreader=False):

  if os.path.isfile(file):
    pass
    tarfile = False
  elif os.path.isfile(file+'.tar.gz'):
    tarfile = True
    cwd = os.getcwd()
    os.chdir(os.path.dirname(file))
    os.system('tar -xzf '+file+'.tar.gz')
  else:
    sys.exit("File "+file+" does not exist.")  

  try:
    from paraview import numpy_support, simple
    # Load vtu file
    if reader == 'none':
      reader = simple.XMLPartitionedUnstructuredGridReader(FileName=file)
    reader.FileName = file
    vtudata = simple.servermanager.Fetch(reader) 
  except:
    try:
      import vtk
      from vtk.util import numpy_support
      if reader == 'none':
        reader = vtk.vtkXMLPUnstructuredGridReader()
      reader.SetFileName(file)
      reader.Update()
      vtudata = reader.GetOutput()
    except:
      sys.exit("You do not have the necessary modules (vtk or paraview) to import vtu files.")

  
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
      del opt, opts
    elif var == 'ssavelocity':
      opts = ['ssavelocity 1','ssavelocity 2','ssavelocity 3','ssavelocity']
      for opt in opts:
        varnames.append(opt)
        types.append(np.float64)
      del opt, opts
    elif var.endswith('update'):
      opts = [var+' 1',var+' 2',var+' 3']
      for opt in opts:
        varnames.append(opt)
        types.append(np.float64)
    elif var == 'vsurfini':
      opts = ['vsurfini 1','vsurfini 2','vsurfini']
      for opt in opts:
        varnames.append(opt)
        types.append(np.float64)
      del opt, opts
    else:
      types.append(np.float64)
      varnames.append(var)
  if ('beta' in varnames) and (('velocity 1' in varnames) or ('ssavelocity 1' in varnames)):
    varnames.append('taub')
    types.append(np.float64) 
    varnames.append('betasquared')
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
      data['velocity 1'] = velocity[:,0]
      data['velocity 2'] = velocity[:,1]
      data['velocity 3'] = velocity[:,2]
      data['velocity'] = np.sqrt(data['velocity 1']**2+data['velocity 2']**2)
    elif var == 'ssavelocity':
      ssavelocity = numpy_support.vtk_to_numpy(vtudata.GetPointData().GetArray(var))
      data['ssavelocity 1'] = ssavelocity[:,0]
      data['ssavelocity 2'] = ssavelocity[:,1]
      data['ssavelocity 3'] = ssavelocity[:,2]
      data['ssavelocity'] = np.sqrt(data['ssavelocity 1']**2+data['ssavelocity 2']**2)
    elif var == 'vsurfini': 
      try:
        vsurfini = numpy_support.vtk_to_numpy(vtudata.GetPointData().GetArray(var))
        data['vsurfini 1'] = vsurfini[:,0]
        data['vsurfini 2'] = vsurfini[:,1]
        data['vsurfini'] = np.sqrt(data['vsurfini 1']**2+data['vsurfini 2']**2)
      except: 
        try:
          data['vsurfini 1'] = numpy_support.vtk_to_numpy(vtudata.GetPointData().GetArray('vsurfini 1'))
          data['vsurfini 2'] = numpy_support.vtk_to_numpy(vtudata.GetPointData().GetArray('vsurfini 2'))
          data['vsurfini'] = np.sqrt(data['vsurfini 1']**2+data['vsurfini 2']**2)
        except:
          # To account for bug in steady.sif, which won't save "vsurfini 1" and "vsurfini 2" 
          # at the same time, so I've saved "vsurfini 2" as vsurfini2
          data['vsurfini 1'] = numpy_support.vtk_to_numpy(vtudata.GetPointData().GetArray('vsurfini 1'))
          data['vsurfini 2'] = numpy_support.vtk_to_numpy(vtudata.GetPointData().GetArray('vsurfini2'))
          data['vsurfini'] = np.sqrt(data['vsurfini 1']**2+data['vsurfini 2']**2)
    elif var.endswith('update'):
      update = numpy_support.vtk_to_numpy(vtudata.GetPointData().GetArray(var))
      data[var+' 1'] = update[:,0]
      data[var+' 2'] = update[:,1]
      data[var+' 3'] = update[:,2]
    else:
      data[var] = numpy_support.vtk_to_numpy(vtudata.GetPointData().GetArray(var))
  if ('taub' in varnames) and ('ssavelocity' in varnames):
    data['taub'] = data['beta']**2.0*np.sqrt(data['ssavelocity 1']**2.0+data['ssavelocity 2']**2.0)
  elif ('taub' in varnames) and ('velocity' in varnames):
    data['taub'] = data['beta']**2.0*np.sqrt(data['velocity 1']**2.0+data['velocity 2']**2.0) 
  if ('betasquared' in varnames):
    data['betasquared'] = data['beta']**2

  # Some of the nodes are shared/repeated between partitions, so we need to remove those
  var1,ind1 = np.unique(x,return_index = True)
  var2,ind2 = np.unique(y,return_index = True)
  var3,ind3 = np.unique(z,return_index = True)
  
  ind = np.union1d(np.union1d(ind1,ind2),ind3)
  
  data_final = data[ind]

  if tarfile:
    i = int(file[-9:-5])
    os.system('rm '+file[0:-10]+'*{0:04d}.'.format(i)+'*vtu')
    os.chdir(cwd)
  
  vtudata.ReleaseData()
  del ind,var1,var2,var3,ind1,ind2,ind3,varnames,vtudata,data,x,y,z,var,types,variables
  del n

  if returnreader:
    return data_final,reader
  else:
    return data_final

def pvtu_timeseries_grid(x,y,DIR,fileprefix,variables,inputsdir,layer='surface',debug=False,crop_mesh=True,t1=1,t2=np.Inf):

  from scipy.interpolate import griddata
  import numpy as np
  import matplotlib.path
  import gc

  # First get number of timesteps
  files = os.listdir(DIR)

  xgrid,ygrid = np.meshgrid(x,y)

  totsteps = 0
  for file in files:
    if file.startswith(fileprefix) and file.endswith('.pvtu'):
      timestep = int(file[len(fileprefix):-5])
      numfilelen = len(file)-len('.pvtu')-len(fileprefix)
      if timestep > totsteps:
        totsteps = timestep
    elif file.startswith(fileprefix) and file.endswith('.pvtu.tar.gz'):
      timestep = int(file[-16:-12])
      numfilelen = len(file)-len('.pvtu.tar.gz')-len(fileprefix)
      if timestep > totsteps:
        totsteps = timestep     
  if totsteps == 0:
    sys.exit("Check that file "+DIR+fileprefix+" actually exists.")
 
  if t2 > totsteps:
    t2 = totsteps
  
  print "Loading "+str(t2-t1+1)+" out of "+str(totsteps)+" timesteps"

  if layer == 'surface':
    freesurfacevar = 'zs top'
  elif layer == 'bed' and 'zs bottom' not in variables:
    freesurfacevar = 'zs bottom'
  if freesurfacevar not in variables:
    variables.append(freesurfacevar)
  
  if crop_mesh:
    mesh_extent_x = np.loadtxt(inputsdir+'/mesh_timeseries_x.dat')
    mesh_extent_y = np.loadtxt(inputsdir+'/mesh_timeseries_y.dat')
  try:
    mesh_extent_x = np.loadtxt(inputsdir+'/mesh_timeseries_x.dat')
    mesh_extent_y = np.loadtxt(inputsdir+'/mesh_timeseries_y.dat')

    path = matplotlib.path.Path(np.column_stack([mesh_extent_x[:,0],mesh_extent_y[:,0]]))
    inmesh = path.contains_points(np.column_stack([xgrid.flatten(),ygrid.flatten()]))
    inmesh = inmesh.reshape(len(y),len(x))
  except:
    mesh_extent = np.loadtxt(inputsdir+'/mesh_extent.dat')

    path = matplotlib.path.Path(np.column_stack([mesh_extent_x[:,0],mesh_extent_y[:,1]]))
    inmesh = path.contains_points(np.column_stack([xgrid.flatten(),ygrid.flatten()]))
    inmesh = inmesh.reshape(len(y),len(x))    
  try: 
    mesh_hole1 = np.loadtxt(inputsdir+'/mesh_hole1.dat')
    mesh_hole2 = np.loadtxt(inputsdir+'/mesh_hole2.dat')
    holes = True

    path = matplotlib.path.Path(np.column_stack([mesh_hole1[:,0],mesh_hole1[:,1]]))
    inhole1 = path.contains_points(np.column_stack([xgrid.flatten(),ygrid.flatten()]))
    inhole1 = inhole1.reshape(len(y),len(x))
    path = matplotlib.path.Path(np.column_stack([mesh_hole2[:,0],mesh_hole2[:,1]]))
    inhole2 = path.contains_points(np.column_stack([xgrid.flatten(),ygrid.flatten()]))
    inhole2 = inhole2.reshape(len(y),len(x))
    del path
  except:
    holes = False

  for i in range(0,t2-t1+1):
    t = i+t1
    # Get filename
    pvtufile = '{0}{2:0{1}d}{3}'.format(fileprefix,numfilelen,t,'.pvtu')
    if debug:
      print "Loading file "+pvtufile
      
    # Get data
    if i==0:
      reader='none'
    data,reader = pvtu_file(DIR+pvtufile,variables,reader=reader,returnreader=True)
    surf = data[data[freesurfacevar] != 0]
    del data
    # If first timestep, set up output variable name
    if i==0:
      varnames = list(surf.dtype.names)
      varnames.remove('Node Number')
      varnames.append('dh')
      types = []
      for var in varnames:
        types.append(np.float64)
      datagrid = np.zeros([len(y),len(x),t2-t1+1], dtype=zip(varnames,types)) 
      del types

    if crop_mesh:
      ind = np.where(mesh_extent_x[:,t-1] != 0)
      path = matplotlib.path.Path(np.column_stack([mesh_extent_x[:,t-1],mesh_extent_y[:,t-1]]))
      inmesh = path.contains_points(np.column_stack([xgrid.flatten(),ygrid.flatten()]))
      inmesh = inmesh.reshape(len(y),len(x))
      del path

    for var in varnames:
      if var == 'dh':
        if i > 0:
          z_old = griddata((surf_last['x'],surf_last['y']),surf_last['z'],(surf['x'],surf['y']))
          dh = surf['z']-z_old
          datagrid[var][:,:,i] = griddata((surf['x'],surf['y']),dh,(xgrid,ygrid))
          del surf_last
        surf_last = np.array(surf)
      else:
        datagrid[var][:,:,i] = griddata((surf['x'],surf['y']),surf[var],(xgrid,ygrid))
      
      #if crop_mesh: 
      datagrid[var][~inmesh,i] = float('nan')
      
      if holes:
        datagrid[var][inhole1,i] = float('nan')
        datagrid[var][inhole2,i] = float('nan')   

    del surf,pvtufile
    if crop_mesh:
      del inmesh,ind
    gc.collect()
    
  return datagrid

def pvtu_flowline(x,y,pvtufile,variables,layer='surface'):

  from scipy.interpolate import griddata
  import numpy as np
  import matplotlib.path


  data = pvtu_file(pvtufile,variables)
  surf = values_in_layer(data,layer=layer)
  del data
  
  varnames = list(surf.dtype.names)
  varnames.remove('Node Number')
  types = []
  for var in varnames:
    types.append(np.float64)
  dataflow = np.zeros([len(x),], dtype=zip(varnames,types)) 

  for var in varnames:
    dataflow[var][:] = np.float('nan')
    dataflow[var] = griddata((surf['x'],surf['y']),surf[var],(x,y))
    
  del surf
     
  return dataflow


def pvtu_timeseries_flowline(x,y,DIR,fileprefix,variables,inputsdir='none',layer='surface',debug=False,t1=1,t2=np.Inf):

  from scipy.interpolate import griddata
  import numpy as np
  import matplotlib.path

  # First get number of timesteps
  files = os.listdir(DIR)

  totsteps = 0
  for file in files:
    if file.startswith(fileprefix) and file.endswith('.pvtu'):
      timestep = int(file[len(fileprefix):-5])
      numfilelen = len(file)-len('.pvtu')-len(fileprefix)
      if timestep > totsteps:
        totsteps = timestep
    elif file.startswith(fileprefix) and file.endswith('.pvtu.tar.gz'):
      timestep = int(file[-16:-12])
      numfilelen = len(file)-len('.pvtu.tar.gz')-len(fileprefix)
      if timestep > totsteps:
        totsteps = timestep     
  if totsteps == 0:
    sys.exit("Check that file "+DIR+fileprefix+" actually exists.")
  
  if t2 > totsteps:
    t2 = totsteps
  
  print "Loading "+str(t2-t1+1)+" out of "+str(totsteps)+" timesteps"

  if layer == 'surface':
    freesurfacevar = 'zs top'
  elif layer == 'bed' and 'zs bottom' not in variables:
    freesurfacevar = 'zs bottom'
  if freesurfacevar not in variables:
    variables.append(freesurfacevar)

  if not(inputsdir == 'none'):
    mesh_extent_x = np.loadtxt(inputsdir+'/mesh_timeseries_x.dat')
    mesh_extent_y = np.loadtxt(inputsdir+'/mesh_timeseries_y.dat')

  for i in range(0,t2-t1+1):
    t = i+t1
    # Get filename
    pvtufile = '{0}{2:0{1}d}{3}'.format(fileprefix,numfilelen,t,'.pvtu')
    if debug:
      print "Loading file "+pvtufile
    # Get data

    if i==0:
      reader = 'none'
    data,reader = pvtu_file(DIR+pvtufile,variables,reader=reader,returnreader=True)
    surf = data[data[freesurfacevar] != 0]
    del data
    # If first timestep, set up output variable name
    if i==0:
      varnames = list(data.dtype.names)
      varnames.remove('Node Number')
      types = []
      for var in varnames:
        types.append(np.float64)
      dataflow = np.zeros([len(x),t2-t1+1], dtype=zip(varnames,types)) 
      for var in varnames:
        dataflow[var][:,i] = float('nan')
    
    if not(inputdirs == 'none'):
      ind = np.where(mesh_extent_x[:,t-1] != 0)
      path = matplotlib.path.Path(np.column_stack([mesh_extent_x[:,t-1],mesh_extent_y[:,t-1]]))
      inmesh = path.contains_points(np.column_stack([x,y]))
    else:
      inmesh = np.arange(0,len(x))

    for var in varnames:
      dataflow[var][inmesh,i] = griddata((surf['x'],surf['y']),surf[var],(x[inmesh],y[inmesh]))
    
    del surf
     
  return dataflow

def pvtu_timeseries_grounding_line(DIR,fileprefix,debug=False,t1=1,t2=np.Inf):

  from scipy.interpolate import griddata
  import numpy as np
  import matplotlib.path

  # First get number of timesteps
  files = os.listdir(DIR)

  totsteps = 0
  for file in files:
    if file.startswith(fileprefix) and file.endswith('.pvtu'):
      timestep = int(file[len(fileprefix):-5])
      numfilelen = len(file)-len('.pvtu')-len(fileprefix)
      if timestep > totsteps:
        totsteps = timestep
    elif file.startswith(fileprefix) and file.endswith('.pvtu.tar.gz'):
      timestep = int(file[-16:-12])
      numfilelen = len(file)-len('.pvtu.tar.gz')-len(fileprefix)
      if timestep > totsteps:
        totsteps = timestep     
  if totsteps == 0:
    sys.exit("Check that file "+DIR+fileprefix+" actually exists.")
  
  if t2 > totsteps:
    t2 = totsteps
  
  print "Loading "+str(t2-t1+1)+" out of "+str(totsteps)+" timesteps"


  variables = ['zs bottom','groundedmask']
  freesurfacevar = 'zs bottom'

  for i in range(0,t2-t1+1):
    t = i+t1
    # Get filename
    pvtufile = '{0}{2:0{1}d}{3}'.format(fileprefix,numfilelen,t,'.pvtu')
    if debug:
      print "Loading file "+pvtufile
    # Get data

    if i==0:
      reader = 'none'
    data,reader = pvtu_file(DIR+pvtufile,variables,reader=reader,returnreader=True)
    surf = data[data[freesurfacevar] != 0]
    # If first timestep, set up output variable name
    if i==0:
      varnames = ['x','y']
      types = []
      for var in varnames:
        types.append(np.float64)
      GLs = np.zeros([500,t2-t1+1], dtype=zip(varnames,types)) 
      for var in varnames:
        GLs[var][:,i] = float('nan')

    ind = np.where(surf['groundedmask'] == 0)
    GLs['x'][0:len(surf['x'][ind].flatten()),i] = surf['x'][ind].flatten()
    GLs['y'][0:len(surf['y'][ind].flatten()),i] = surf['y'][ind].flatten()

    del ind,data,surf 
  return GLs


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
  for x_val in np.unique((data['x'])):
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

def input_file(file,dim=2):
  '''
  x,y,grid = input_file(file,dim=2)
  
  Load input grid file for elmersolver.
  
  Inputs:
  file: filename
  dim: dimension of file, currently only works for 2
  
  Outputs:
  x,y,grid values of file
  '''

  fid = open(file,"r")
  lines = fid.readlines()
  nx = int(lines[0].split()[0])
  ny = int(lines[1].split()[0])
  if dim==3:
    nz = int(lines[2].split()[0])

  x = np.zeros(nx)
  y = np.zeros(ny) 
  if dim==2:
    grid = np.zeros([ny,nx])

    n = 0
    for i in range(0,nx):
      for j in range(0,ny):
        x[i],y[j],grid[j,i] = lines[2+n].split()
        n = n+1
  elif dim==3:
    grid = np.zeros([ny,nx,nz])

    n = 0
    for i in range(0,nx):
      for j in range(0,ny):
        p = (lines[3+n]).split()
        x[i] = p[0]
        y[j] = p[1]
        for k in range(0,nz):
          grid[j,i,k] = p[2+k]
        n = n+1

  fid.close()

  return x,y,grid


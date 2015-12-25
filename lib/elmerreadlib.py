# This module reads the results from Elmer.

def saveline(DIR,runname):
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
  if 'velod4' in variables:
    data['velmes']=[]
    data['vel']=[]
    for i in range(0,len(data['vel1'])):
      data['vel'].append(math.sqrt(data['vel1'][i]**2+data['vel2'][i]**2))
      data['velmes'].append(math.sqrt(data['vsurfini1'][i]**2+data['vsurfini2'][i]**2))
  
  # Calculate basal shear stress
  if 'beta' in variables:
    data['taub']=[]
    if 'velod4' in variables:
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
  import numpy as np
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

def grid3d(data,variable,holes,extent):

  import matplotlib.mlab as mlab
  import math
  import numpy as np
  from scipy.interpolate import griddata
  from matplotlib.path import Path
  
  xmin=math.floor(min(data['coord1'])/100)*100
  xmax=math.ceil(max(data['coord1'])/100)*100
  ymin=math.floor(min(data['coord2'])/100)*100
  ymax=math.ceil(max(data['coord2'])/100)*100

  x=np.linspace(xmin,xmax,(xmax-xmin)/50+1)
  y=np.linspace(ymin,ymax,(ymax-ymin)/50+1)
  
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
  # Takes a 3D model output and calculates values along a flowline given 
  # by "x" and "y"
  
  import numpy as np
  from scipy.interpolate import griddata
  
  gridx=data['coord1']
  gridy=data['coord2']
  
  variables=data.keys()
  
  flow={}
  for variable in variables:
    flow[variable]=griddata(np.column_stack([gridx,gridy]),data[variable],np.column_stack([x,y]),method='linear')
  
  return flow
import os
import sys
import numpy as np
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
import coords

def atm(year):
  
  if year == '2001':
    DIR = os.path.join(os.getenv("HOME"),'/Users/kehrl/Data/Elevation/ATM/Helheim/20010521')
    files = os.listdir(DIR)
  
    os.chdir(DIR)
    y=[]
    x=[]
    z=[]
    for file in files:
      if file.endswith('nadir5seg'):
        data=np.loadtxt(file)
        xfile = data[:,2]
        yfile = data[:,1]
        zfile = data[:,3]
        x=np.hstack([x,xfile])
        y=np.hstack([y,yfile])
        z=np.hstack([z,zfile])
      
        x2,y2 = coords.convert(x-360,y,4326,3413)
        surf = np.column_stack([x2,y2,z])
    
  return surf
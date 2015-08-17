# This file takes the Microsoft Excel spreadsheets with potential TanDEM DEMs and selects 
# the best image pairs for DEM production.
#
#
# LMK, UW, 8/16/2015

import openpyxl
import os,sys
import numpy as np
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
import fracyear

# Choose glacier (Kanger, Helheim, Jak)
glacier = 'jak'

# Ranges of acceptable values:
hoa_lims = [20.0,200.0]
ebase_lims = [30.0,200.0]

# Choose relative orbits
if glacier == 'kanger':
  orbits = [163,65]
elif glacier == 'helheim':
  orbits = [148]
elif glacier == 'jak':
  orbits = [42,118]
  
# Load excel spreadsheet
excelfile=os.path.join(os.getenv("HOME"),"Data/Elevation/TanDEM/OrderLists/"+glacier+"_TDM.xlsx")
data = openpyxl.load_workbook(excelfile)[glacier]

# Values that we want from the file
N = len(data.columns[0])-1
hoa = np.zeros(N) # height of ambiguity
ebase = np.zeros(N) # effective baseline
scenenum = np.zeros(N) # scenenumber
date = np.zeros([N,3]) # Date
relorbit = np.zeros(N) # Relative orbit 
absorbit = np.zeros(N) # Absolute orbit
location = []

# Load those values from the file
for column in data.columns:
  if column[0].value == 'Start Date':
    n = 0
    for cell in column[1:]:
      date[n,:] = ([int(cell.value[0:4]),int(cell.value[5:7]),int(cell.value[8:10])])  
      n = n+1
  elif column[0].value == 'Relative Orbit':
    n = 0
    for cell in column[1:]:
      relorbit[n] = (int(cell.value))
      n = n+1
  elif column[0].value == 'Effective Baseline':
    n = 0
    for cell in column[1:]:
      ebase[n] = (float(cell.value))
      n = n+1
  elif column[0].value == 'Scene Number':
    n = 0
    for cell in column[1:]:
      scenenum[n] = (int(cell.value))
      n = n+1
  elif column[0].value == 'Absolute Orbit':
    n = 0
    for cell in column[1:]:
      absorbit[n] = (int(cell.value))
      n = n+1
  elif column[0].value == 'Height of Ambiguity':
    n = 0
    for cell in column[1:]:
      hoa[n] = (float(cell.value))
      n = n+1
  elif column[0].value == 'FOOTPRINT':
    for cell in column[1:]:
      location.append (str(cell.value))
del column, cell, n

# Select candidates for DEM production
cand_orbit = []
cand_ebase = []
cand_hoa = []
for i in range(0,len(date)):
  if relorbit[i] in orbits:
    cand_orbit.append(i)
  if (ebase[i] > ebase_lims[0]) and (ebase[i] < ebase_lims[1]):
    cand_ebase.append(i)
  if (abs(hoa[i]) > hoa_lims[0]) and (abs(hoa[i]) < hoa_lims[1]):    
    cand_hoa.append(i)

# Find intersection of candidate indices
cand_ind = list(set(cand_orbit).intersection(cand_ebase,cand_hoa))

# Get potential image pairs
hoa_cand = hoa[cand_ind]
ebase_cand = ebase[cand_ind]
relorbit_cand = relorbit[cand_ind]
absorbit_cand = absorbit[cand_ind]
date_cand = date[cand_ind,:]
scenenum_cand = scenenum[cand_ind]
fracyear_cand = fracyear.date_to_fracyear(date[cand_ind,0],date[cand_ind,1],date[cand_ind,2])
location_cand = np.array(location)[cand_ind]
print len(cand_ind)

fid = open(glacier+'_candidate_TDM.dat','w')
fid.write('Date\t\t AbsoluteOrbit\t RelativeOrbit\t SceneNumber\t HeightOfAmbiguity\t EffectiveBaseline\t Footprint\n') 
for i in range(0,len(cand_ind)):
  fid.write('{0}-{1:02d}-{2:02d}\t {3:3d}\t\t\t {4:3d}\t\t\t {5}\t\t\t\t {6:.4f}\t\t\t {7:.4f}\t\t\t {8}\n'.format(int(date_cand[i,0]),(int(date_cand[i,1])),(int(date_cand[i,2])),int(absorbit_cand[i]),int(relorbit_cand[i]),int(scenenum_cand[i]),hoa_cand[i],ebase_cand[i],location_cand[i]))
fid.close() 

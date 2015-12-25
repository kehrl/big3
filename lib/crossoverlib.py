'''
Find and process crossover points.

LMK, UW, 10/01/2015
'''

import os
import sys
import distlib, coordlib
import numpy as np
import shapefile
from shapely.geometry import LineString

def find(x,y,z,maxdist=200.0):

  '''
  cross,ycross,zcross1,zcross2,zcrossdiff = find(x,y,z,maxdist=200.0)
  
  Find crossovers and their differences for an array of x,y,z points.
  
  Inputs:
  x,y,z : input coordinates
  maxdist : maximum acceptable distance between points in crossovers
  
  Outputs:
  xcross,ycross : locations of crossovers
  zcross1,zcross2 : heights at crossovers
  zcrossdiff: difference between zcross1 and zcross2
  '''

  # Number of points
  N = len(x)
  
  # Minimum number of indices between two points to be considered a crossover.
  minind = 10
  
  # Get a list of potential crossover indices
  crossind = []
  for i in range(0,N):
    ind = (np.where(np.sqrt((x[i]-x)**2+(y[i]-y)**2) < maxdist))[0]
    if np.max(ind)-i > minind or i-np.min(ind) > minind: 
      ind = np.delete(ind,np.where(abs(ind-i) < minind)[0])
      for j in range(0,len(ind)):
        if (j == 0) or abs(ind[j]-ind[j-1]) > 1:
          if [ind[j],i] not in crossind:
            crossind.append([i,ind[j]]) 
  
  # OK. Now we go through and create line segments around each of these potential 
  # crossover points. If the lines actually intersect, then we actually calculate and
  # save the crossover point.
  
  # Number of potential crossover points
  M = len(crossind)
  
  xcross=[]
  ycross=[]
  zcross1=[]
  zcross2=[]
  n=0
  for i in range(0,M):
    # Indices for the line segments
    if (crossind[i][0] - 5 < 0): 
      start1 = 0   
    else:
      start1 = crossind[i][0]-5
    if  (crossind[i][1] - 5 < 0):
      start2 = 0
    else:
      start2 = crossind[i][1]-5
    if (crossind[i][0] + 5 > len(x) - 1):
      fin1 = len(x)-1
    else:
      fin1 = crossind[i][0]+5
    if (crossind[i][1] + 5 > len(x) - 1):
      fin2 = len(x)-1
    else:
      fin2 = crossind[i][1]+5
    ind1 = np.arange(start1,fin1,1)
    ind2 = np.arange(start2,fin2,1)
    sortind1 = coordlib.sortind_xy(x[ind1],y[ind1])
    sortind2 = coordlib.sortind_xy(x[ind2],y[ind2])
    
    # line segments
    line1_x = x[ind1[sortind1]]
    line1_y = y[ind1[sortind1]]
    line1_z = z[ind1[sortind1]]
    line2_x = x[ind2[sortind2]]
    line2_y = y[ind2[sortind2]]
    line2_z = z[ind2[sortind2]]
    
    line1 = LineString(np.column_stack([line1_x,line1_y]))
    line2 = LineString(np.column_stack([line2_x,line2_y]))

    intersections = (line1.intersection(line2))
    if not(intersections.is_empty):
      xpt = []
      ypt = []
      if intersections.geometryType() is 'Point':
        xpt.append(intersections.x)
        ypt.append(intersections.y)
      else:
        for intersection in intersections:
          xpt.append(intersection.x)
          ypt.append(intersection.y)
      for j in range(0,len(xpt)):
        
        ind1_1 = np.argmin((xpt[j]-line1_x)**2+(ypt[j]-line1_y)**2)
        ind1_2 = np.argmin((xpt[j]-line2_x)**2+(ypt[j]-line2_y)**2)
        d1_1 = np.sqrt((xpt[j]-line1_x[ind1_1])**2+(ypt[j]-line1_y[ind1_1])**2)
        d1_2 = np.sqrt((xpt[j]-line2_x[ind1_2])**2+(ypt[j]-line2_y[ind1_2])**2)
        
        if (d1_1 < maxdist) and (d1_2 < maxdist) and (ind1_1 < len(line1_x)-1) and (ind1_1 > 0) and (ind1_2 < len(line2_x)-1) and (ind1_2 > 0) and xpt[j] not in xcross:
          # Save crossover location
          xcross.append(xpt[j])
          ycross.append(ypt[j])

          # Find height for first line at intersection

          if (distlib.between_pts(line1_x[ind1_1-1],line1_y[ind1_1-1],xpt[j],ypt[j]) < \
          		distlib.between_pts(line1_x[ind1_1+1],line1_y[ind1_1+1],xpt[j],ypt[j])) and \
          		(distlib.between_pts(line1_x[ind1_1-1],line1_y[ind1_1-1],xpt[j],ypt[j]) < \
          		distlib.between_pts(line1_x[ind1_1-1],line1_y[ind1_1-1],line1_x[ind1_1],line1_y[ind1_1])):
            ind2_1 = ind1_1-1
          else:
            ind2_1 = ind1_1+1
          d2_1 = np.sqrt((xpt[j]-line1_x[ind2_1])**2+(ypt[j]-line1_y[ind2_1])**2)
          zcross1.append((line1_z[ind2_1]-line1_z[ind1_1])/(d1_1+d2_1)*d1_1+line1_z[ind1_1])
          
          # Find height for second line at intersection
          if (distlib.between_pts(line2_x[ind1_2-1],line2_y[ind1_2-1],xpt[j],ypt[j]) < \
          		distlib.between_pts(line2_x[ind1_2+1],line2_y[ind1_2+1],xpt[j],ypt[j])) and \
          		(distlib.between_pts(line2_x[ind1_2-1],line2_y[ind1_2-1],xpt[j],ypt[j]) < \
          		distlib.between_pts(line2_x[ind1_2-1],line2_y[ind1_2-1],line2_x[ind1_2],line2_y[ind1_2])):
            ind2_2 = ind1_2-1
          else:
            ind2_2 = ind1_2+1
          d2_2 = np.sqrt((xpt[j]-line2_x[ind2_2])**2+(ypt[j]-line2_y[ind2_2])**2)
          zcross2.append((line2_z[ind2_2]-line2_z[ind1_2])/(d1_2+d2_2)*d1_2+line2_z[ind1_2]) 
            

  xcross = np.array(xcross)
  ycross = np.array(ycross)
  zcross1 = np.array(zcross1)
  zcross2 = np.array(zcross2)
  zcrossdiff = abs(zcross1-zcross2)

  return xcross,ycross,zcross1,zcross2,zcrossdiff
import numpy as np
from scipy.spatial import cKDTree
import sys
import os
import shutil

def arrhenius(T):
  '''
  This function takes a temperature and outputs the flow parameter 
  following Cuffey and Paterson 2010, pg. 72.
  '''
  
  A=np.zeros_like(T)
  for i in range(0,len(T)):
    if (T[i] > 273.15):
      A[i]=2.4E-24
    elif (T[i] > 263.15):
      A[i]=3.5E-25* np.exp(-115.0E03/8.314 * (1/T[i] - 1/263.15))
    else:
      A[i]=3.5E-25 * np.exp(-60.0E03/8.314 * (1/T[i] - 1/263.15))
    	
  return A 

def load_kristin(glacier,x,y,type='A',dir='none'):

  # Load file if it already exists
  fileA = os.path.join(os.getenv("DATA_HOME"),"Climate/IceTemperature/"+glacier+"/flowA.xyz")
  fileT = os.path.join(os.getenv("DATA_HOME"),"Climate/IceTemperature/"+glacier+"/flowT.xyz")
  if (type=='A') and (os.path.isfile(fileA)):
    shutil.copy(fileA,dir+"flowA.xyz")
    output = 'success'
  elif (type=='T') and (os.path.isfile(fileT)):
    shutil.copy(fileT,dir+"flowT.xyz")  
    output = 'success' 
  else:
    # Get modeled temperatures from Kristin's work
    if glacier == 'Helheim':
      kristin_file=os.path.join(os.getenv("DATA_HOME"),"Climate/IceTemperature/Helheim/helheim_TA.xyz")
    else:
      sys.exit('Unknown glacier')  
  
    # Get grid values for gridding Kristin's values
    nx = len(x)
    ny = len(y)
  
    dx = x[1]-x[0]
    dy = y[1]-y[0]
  
    xmin = x[0]
    xmax = x[-1]
    ymin = y[0]
    ymax = y[-1]

    # Variables for loading Kristin's values
    xT = []
    yT = []
    zT = []
    tT = []
    aT = []

    # Read in the data from a file
    fid = open(kristin_file, "r")
    line = fid.readline().split(',')
    line = fid.readline().split(',')

    while(line[0] != ""):
      xT.append(float(line[0]))
      yT.append(float(line[1]))
      zT.append(float(line[2]))
      tT.append(float(line[3]))
      aT.append(float(line[4][:-1]))
      line = fid.readline().split(',')
    fid.close()

    # Condense the list so that all the temperature / viscosity from the same
    # (x, y)-point is contained in one list, and the data in a list of lists
    X = []
    Y = []
    Z = []
    T = []
    A = []
    
    X.append(xT[0])
    Y.append(yT[0])
    Z.append([zT[0]])
    T.append([tT[0]])
    A.append([aT[0]])
    
    k = 0
    
    for n in range(1, len(xT)):
      if xT[n] == xT[n - 1] and yT[n] == yT[n - 1]:
        Z[k].append(zT[n])
        T[k].append(tT[n])
        A[k].append(aT[n])
      else:
        X.append(xT[n])
        Y.append(yT[n])
        Z.append([zT[n]])
        T.append([tT[n]])
        A.append([aT[n]])
        k += 1
        
    nn = len(X)
    
    X = np.asarray(X)
    Y = np.asarray(Y)
    
    for k in range(nn):
      Z[k] = np.asarray(Z[k])
      T[k] = np.asarray(T[k])
      A[k] = np.asarray(A[k])

    # Set the number of vertical layers
    nz = 20
    
    temp = np.zeros((ny,nx,nz)) # Kelvin
    rate = np.zeros((ny,nx,nz))

		# Make a gridded data set from the model output

    # For each point in the grid,
    for i in range(ny):
		  for j in range(nx):

				# Find five closest points
				dists = np.sqrt((x[j]-X)**2+(y[i]-Y)**2)
				ind = dists.argsort()[0:5]
			
				# Initialize the weights to 1; if there aren't any
				# other nearby model points to interpolate from, we'll
				# keep the temperature value at -10C
				weights = 0.0

				# For all the nearby model points,
				for l in ind:
					xp = X[l]
					yp = Y[l]

					# find the distance to the current point and the
					# appropriate weight
					r = np.sqrt( (x[j] - xp)**2 + (y[i] - yp)**2 )
					w = (1000/(r+dx))**3
					weights += w

					# For each point within the current vertical column,
					for k in range(nz):
						# find which point within the nearby vertical
						# column to interpolate from
						m = (k * (len(Z[l]) - 1)) / (nz - 1)

						# Add up the value to the running average
						rate[i, j, k] += w * A[l][m]
						temp[i, j, k] += w * T[l][m]

				# Normalize the running average by the weight sum
				rate[i,j,:] /= weights  
				temp[i,j,:] /= weights  

    if type == 'A':
      output = rate
    elif type == 'T':
      output = temp
  
    if dir != 'none':
      fid = open(dir+"flowA.xyz", "w")
      fid.write("{0}\n{1}\n{2}\n".format(len(x), len(y), len(output[0,0,:])))

      for j in range(len(x)):
        for i in range(len(y)):
          fid.write("{0} {1} ".format(x[j], y[i]))
          for k in range(len(output[0,0,:])):
            fid.write("{0} ".format(output[i, j, k]))
          fid.write("\n")
      fid.close()
      
  return output
import numpy as np
import sys
import os
import shutil
import distlib, elmerreadlib
import scipy
from scipy.spatial import cKDTree
import scipy.interpolate
import math

def arrhenius(T,paterson_version=2010):
  '''
  A = arrhenius(T)
  
  Inputs:
  T: temperature in Kelvin
  
  Outputs: 
  A: Arrhenius parameter
  
  This function takes a temperature in Kelvin and outputs the flow parameter 
  following Cuffey and Paterson 2010, pg. 72.
  
  '''
  
  R = 8.314         # gas constant
  Tstar = 263.15    # temperature where activation energy flips
  Qbelow = 60.0E3
  Qabove = 115.0E3
  Ao = 3.5E-25
  
  try:
    n = len(T)
  except:
    n = 1
    T = np.ones(n)*T
  
  A = np.zeros(n)
  for i in range(0,len(T)):
    if (T[i] > 273.15):
      A[i] = 2.4E-24
    elif (T[i] > Tstar):
      A[i] = Ao* np.exp(-1*Qabove/R * (1/T[i] - 1/Tstar))
    else:
      A[i] = Ao * np.exp(-Qbelow/R * (1/T[i] - 1/Tstar))
    	
  return A 

def arrhenius_reverse(A,paterson_version=2010):
  '''
  A = arrhenius(T)
  
  Inputs:
  T: temperature in Kelvin
  
  Outputs: 
  A: Arrhenius parameter
  
  This function takes a temperature in Kelvin and outputs the flow parameter 
  following Cuffey and Paterson 2010, pg. 72.
  
  '''
  
  R = 8.314         # gas constant
  Tstar = 263.15    # temperature where activation energy flips
  Qbelow = 60.0E3
  Qabove = 115.0E3
  Ao = 3.5E-25
  
  try:
    n = len(A)
  except:
    n = 1
    A = np.ones(n)*A
  
  T = np.zeros(n)
  for i in range(0,len(T)):
    if (A[i] >= 2.4E-24):
      T[i] = 273.15
    elif (A[i] > 3.5E-25):
      T[i] = 1/((-R/Qabove)*np.log(A[i]/Ao)+1/Tstar)
    else:
      T[i] = 1/((-R/Qbelow)*np.log(A[i]/Ao)+1/Tstar)
      	
  return T

def load_temperature_model(glacier,x,y,modelfile='none',outputdir='none',method='linear'):

  # Choose file
  if (modelfile == 'none') and (glacier == 'Helheim'):
    modelfile = os.path.join(os.getenv("MODEL_HOME"),"Helheim/3D/BASIN20120316/mesh2d/temperature/temperature_20170906/temperature0005.pvtu")
  elif (modelfile == 'none') and (glacier == 'Kanger'):
    modelfile = os.path.join(os.getenv("MODEL_HOME"),"Kanger/3D/BASIN20120213/mesh2d/temperature/temperature_20170901/temperature0006.pvtu")

  print "Pulling temperatures from "+modelfile

  # What variables to load 
  variables = ['temp homologous','velocity']
 
  # Get temperatures from model
  data = elmerreadlib.pvtu_file(modelfile,variables)

  # Get info about output grid
  nx = len(x)
  ny = len(y)
  dx = x[1]-x[0]
  dy = y[1]-y[0]

  # Get temperatures, velocities in a column
  X = []
  Y = []
  Z = []
  T = []
  U = []
  V = []
    
  for x_val in np.unique(data['x']):
    x_points = data[data['x'] == x_val]
    for y_val in np.unique(x_points['y']):
      X.append(x_val)
      Y.append(y_val)
      y_points = x_points[x_points['y'] == y_val]
      sorted_list = np.sort(y_points, order='z')
      Z.append(sorted_list['z'])
      T.append(sorted_list['temp homologous'])
      U.append(sorted_list['velocity 1'])
      V.append(sorted_list['velocity 2']) 
        
  nn = len(X)
  X = np.asarray(X)
  Y = np.asarray(Y)
    
  for k in range(nn):
    Z[k] = np.asarray(Z[k])
    T[k] = np.asarray(T[k])

  # Set the number of vertical layers
  nz = len(Z[0])
    
  temp = np.zeros((ny,nx,nz))
  temp[:,:,:] = 0.
  u = np.zeros((ny,nx,nz))
  u[:,:,:] = 0.
  v = np.zeros((ny,nx,nz))
  v[:,:,:] = 0.

  if method=='linear':
    Tnew = np.zeros([len(X),len(T[0])])
    Unew = np.zeros([len(X),len(U[0])])
    Vnew = np.zeros([len(X),len(V[0])])
    for i in range(0,len(T)):
      Tnew[i,:] = T[i]    
      Unew[i,:] = U[i]
      Vnew[i,:] = V[i]

    # Get flattened grid for interpolation
    xgrid,ygrid = np.meshgrid(x,y)
    xflat = xgrid.flatten()
    yflat = ygrid.flatten()
    ind = np.where((X >= (x[0]-5e3)) & (X <= (x[-1]+5e3)) & (Y >= (y[0]-5e3)) & (Y <= (y[-1]+5e3)))
    
    for k in range(0,nz):
      temp[:,:,k] = np.reshape(scipy.interpolate.griddata((X[ind],Y[ind]),\
                  Tnew[ind,k][0],(xflat,yflat),method='nearest'),(ny,nx))
      u[:,:,k] = np.reshape(scipy.interpolate.griddata((X[ind],Y[ind]),\
                  Unew[ind,k][0],(xflat,yflat),method='nearest'),(ny,nx))
      v[:,:,k] = np.reshape(scipy.interpolate.griddata((X[ind],Y[ind]),\
                  Vnew[ind,k][0],(xflat,yflat),method='nearest'),(ny,nx))      

      temp_linear = np.reshape(scipy.interpolate.griddata((X[ind],Y[ind]),\
                  Tnew[ind,k][0],(xflat,yflat),method='linear') ,(ny,nx))
      u_linear = np.reshape(scipy.interpolate.griddata((X[ind],Y[ind]),\
                  Unew[ind,k][0],(xflat,yflat),method='linear') ,(ny,nx))
      v_linear = np.reshape(scipy.interpolate.griddata((X[ind],Y[ind]),\
                  Vnew[ind,k][0],(xflat,yflat),method='linear') ,(ny,nx))

      
      nonnan = np.where(~(np.isnan(temp_linear)))
      temp[nonnan[0],nonnan[1],k] = temp_linear[nonnan[0],nonnan[1]]
      v[nonnan[0],nonnan[1],k] = u_linear[nonnan[0],nonnan[1]]
      u[nonnan[0],nonnan[1],k] = v_linear[nonnan[0],nonnan[1]] 

  else:
    # Make a KD-tree so we can do range searches fast
    tree = cKDTree(np.column_stack([X,Y]))

    # Make a gridded data set from the model output
    # For each point in the grid,
    for i in range(ny):
      for j in range(nx):
		    
        L = tree.query_ball_point( (x[j], y[i]), 1500. )
		    
        # Initialize the weights to 0
        weights = 0.0
		    
        # For all the nearby model points,
        if len(L) > 0:
          for l in L:
            xp = X[l]
            yp = Y[l]
		      
            # find the distance to the current point and the
	    # appropriate weight
	    r = np.sqrt( (x[j] - xp)**2 + (y[i] - yp)**2 )
            w = (1500./(r+dx)**2.0)
            weights += w
		      
	    # For each point within the current vertical column,
            for k in range(nz):
              # find which point within the nearby vertical
              # column to interpolate from
	      m = (k * (len(Z[l]) - 1)) / (nz - 1)
		        
	      # Add up the value to the running average
	      temp[i, j, k] += w * T[l][m]
              u[i, j, k] += w * U[l][m]
              v[i, j, k] += w * V[l][m]
		    
          # Normalize the running average by the weight sum
          temp[i,j,:] /= weights
          u[i,j,:] /= weights
          v[i,j,:] /= weights  
    
        else:
          # If no points within that distance, set to constant temperature value
          temp[i,j,:] = -10.0
          u[i,j,:] = -2.0e9
          v[i,j,:] = -2.0e9
          #L = np.argmin(np.sqrt((X-x[j])**2+(Y-y[i])**2))
          #temp[i,j,:] = T[L]
  
  Agrid = arrhenius(273.15+temp.flatten()).reshape(ny,nx,nz)
  output = temp 
 
  if outputdir != 'none':
    fidT = open(outputdir+"modelT.xyz", "w")
    fidT.write("{0}\n{1}\n{2}\n".format(len(x), len(y), len(output[0,0,:])))
    fidU = open(outputdir+"modelU.xyz", "w")
    fidU.write("{0}\n{1}\n{2}\n".format(len(x), len(y), len(u[0,0,:])))        
    fidV = open(outputdir+"modelV.xyz", "w")
    fidV.write("{0}\n{1}\n{2}\n".format(len(x), len(y), len(v[0,0,:])))

    for j in range(len(x)):
      for i in range(len(y)):
        fidT.write("{0} {1} ".format(x[j], y[i]))
        fidU.write("{0} {1} ".format(x[j], y[i]))        
        fidV.write("{0} {1} ".format(x[j], y[i]))
        for k in range(len(output[0,0,:])):
          fidT.write("{0} ".format(output[i, j, k]))
          fidU.write("{0} ".format(u[i, j, k]))
          fidV.write("{0} ".format(v[i, j, k]))
        fidT.write("\n") 
        fidU.write("\n")
        fidV.write("\n")
    fidT.close()
    fidU.close()
    fidV.close()
          
  return temp,Agrid,u,v

def load_kristin(glacier,x,y,type='A',dir='none'):

  # Load file if it already exists
  fileA = os.path.join(os.getenv("DATA_HOME"),"Climate/IceTemperature/"+glacier+"/modelA.xyz")
  fileT = os.path.join(os.getenv("DATA_HOME"),"Climate/IceTemperature/"+glacier+"/modelT.xyz")
  if (type=='A') and (os.path.isfile(fileA)):
    shutil.copy(fileA,dir+"modelA.xyz")
    output = 'success'
  elif (type=='T') and (os.path.isfile(fileT)):
    shutil.copy(fileT,dir+"modelT.xyz")  
    output = 'success' 
  else:
    # Get modeled temperatures from Kristin's work
    if glacier == 'Helheim':
      kristin_file=os.path.join(os.getenv("DATA_HOME"),"Climate/IceTemperature/Helheim/xyzTAhelheim_2016.txt")
    elif glacier == 'Kanger':
      kristin_file=os.path.join(os.getenv("DATA_HOME"),"Climate/IceTemperature/Kanger/xyzTAkanger_2016.txt")
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

    Xn = []
    Yn = []
    Zn = []
    Tn = []
    An = []
    fracs = []
    
    D = distlib.transect(X,Y)
    for i in range(1,len(D)):
      if D[i] - D[i-1] < 5000.:
        dtot = D[i] - D[i-1]
        dn = np.arange(D[i-1],D[i],500)
        for j in range(0,len(dn)):
          Xn.append(float(np.interp(dn[j],[D[i-1],D[i]],[X[i-1],X[i]])))
          Yn.append(float(np.interp(dn[j],[D[i-1],D[i]],[Y[i-1],Y[i]])))
          frac = (D[i]-dn[j])/dtot
          fracs.append(frac)
          Zn.append(frac*np.asarray(Z[i-1])+(1-frac)*np.asarray(Z[i]))
          Tn.append(frac*np.asarray(T[i-1])+(1-frac)*np.asarray(T[i]))
          An.append(frac*np.asarray(A[i-1])+(1-frac)*np.asarray(A[i]))
      else:
        Xn.append(X[i])
        Yn.append(Y[i])
        Zn.append(Z[i])
        Tn.append(T[i])
        An.append(A[i])
    Xn = X
    Yn = Y
    Zn = Z
    Tn = T
    An = A

    # Set the number of vertical layers
    nz = 21
    
    temp = np.zeros((ny,nx,nz)) # Kelvin
    temp[:,:,:] = 263.
    rate = np.zeros((ny,nx,nz))
    rate[:,:,:] = arrhenius([263.])

    # Make a KD-tree so we can do range searches fast
    #tree = cKDTree(np.column_stack([Xn,Yn]))

		# Make a gridded data set from the model output

    # For each point in the grid,
    for i in range(ny):
      for j in range(nx):
		    
        #L = tree.query_ball_point( (x[j], y[i]), 5000 )
		    
		    # Find five closest points
		    dists = np.sqrt((x[j]-Xn)**2+(y[i]-Yn)**2)
		    allind = np.where(dists < 5000.)[0]
		    sortind = dists[allind].argsort()
		    ind = allind[sortind[0:10]]
		    #ind = ind[dists.sort < 5000.]
		    
		    # Initialize the weights to 0
		    weights = 1.0
		    
		    # For all the nearby model points,
		    if len(ind) > 0:
		      for l in ind:
		        xp = Xn[l]
		        yp = Yn[l]
		      
		        # find the distance to the current point and the
		        # appropriate weight
		        r = np.sqrt( (x[j] - xp)**2 + (y[i] - yp)**2 )
		        w = (5000./(r+dx))**3
		        weights += w
		      
		        # For each point within the current vertical column,
		        for k in range(nz):
		          # find which point within the nearby vertical
		          # column to interpolate from
		          m = (k * (len(Zn[l]) - 1)) / (nz - 1)
		        
		          # Add up the value to the running average
		          rate[i, j, k] += w * An[l][m]
		          temp[i, j, k] += w * Tn[l][m]
		    
		      # Normalize the running average by the weight sum
		      rate[i,j,:] /= weights  
		      temp[i,j,:] /= weights
		      if np.max(temp[i,j,:])>300.:
		        print i,j
    
    smoothtemp = np.zeros_like(temp)
    smoothrate = np.zeros_like(rate)
    for j in range(0,nz):
      smoothtemp[:,:,j] = scipy.ndimage.filters.gaussian_filter(temp[:,:,j],sigma=6,truncate=6)
      smoothrate[:,:,j] = scipy.ndimage.filters.gaussian_filter(rate[:,:,j],sigma=6,truncate=6)
    
    # Calculate flow law parameter for temperatures, using parameterization from Paterson 2010
    newrate = np.zeros_like(temp)
    for i in range(0,nx):
      for j in range(0,ny):
        newrate[j,i,:] = arrhenius(np.array(smoothtemp[j,i,:]))
    
    if dir != 'none':
      fidA = open(dir+"modelA.xyz", "w")
      fidA.write("{0}\n{1}\n{2}\n".format(len(x), len(y), len(newrate[0,0,:])))
      
      fidT = open(dir+"modelT.xyz", "w")
      fidT.write("{0}\n{1}\n{2}\n".format(len(x), len(y), len(temp[0,0,:])))

      for j in range(len(x)):
        for i in range(len(y)):
          fidA.write("{0} {1} ".format(x[j], y[i]))
          fidT.write("{0} {1} ".format(x[j], y[i]))
          for k in range(len(newrate[0,0,:])):
            fidA.write("{0} ".format(newrate[i, j, k]))
            fidT.write("{0} ".format(temp[i, j, k]))
          fidA.write("\n")
          fidT.write("\n")
      fidT.close()
      fidA.close()
      
  return output
  
def steadystate_vprofile(H,Ts,bdot_myr,levels=20):

  '''
  To estimate temperature at the ice divide, I use the steady-state analytical 
  solution from Cuffey and Paterson, pg. 410. These values are used as our boundary condition
  for calculating temperatures.

  T = Ts + zstar*(sqrt(pi)/2)*dTdz*(erf(z/zstar)-erf(H/zstar))

  where zstar = sqrt(2*alpha*H/bdot)

  Inputs:
  H: grid of heights (m)
  Ts: grid of Ts (K)
  bdot_myr: grid of accumulation rates (m/yr) -- will need to be converted to m/s for these calculations


  Output:
  Tgrid: 3D grid of temperatures assuming steady state vertical advection/diffusion and no horizontal motion
  
  
  '''

  # Constants
  Tref = 263.
  c = 152.5+7.122*Tref # empirical formula for specific heat capacity, assuming T = 263.
  k = 9.828*np.exp(-5.7*1e-3*Tref)
  rho = 917. # ice density
  alpha = k/(rho*c) # Thermal diffusivity (m2 s-1)
  dTdz = -0.06/2.2 # Geothermal heat flux (60 mW) / thermal conductivity (K m-1)

  ny,nx = np.shape(H)
  Tgrid = np.zeros([ny,nx,levels])
  Tgrid[:,:,:] = float('nan')
  
  bdot = bdot_myr/(365.25*24*60*60) # change units from m/yr to m/s
  
  for j in range(0,ny):
    for i in range(0,nx):
      # If we are in the ablation zone the solution won't work, or if ice height is unknown
      # In that case we just set the temperature at all depths to the surface temperature
      if (bdot[j,i] < 0) | (np.isnan(H[j,i])) | (H[j,i] < 0):
        Tgrid[j,i,:] = Ts[j,i]
      else:
        zstar = np.sqrt(2*alpha*H[j,i]/bdot[j,i])
        z = np.linspace(0,H[j,i],levels)
        for k in range(0,levels):
          Tgrid[j,i,k] = Ts[j,i] + zstar*(np.sqrt(np.pi)/2)*dTdz*(math.erf(z[k]/zstar)-math.erf(H[j,i]/zstar))
  
  return Tgrid
  

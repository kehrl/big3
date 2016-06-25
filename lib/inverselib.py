def guess_beta(x,y,zs,zb,u,v,frac):   
   
  import numpy as np 
  import math
   
  # Inputs:
  # x    : list of horizontal coordinates of the grid
  # y    : list of vertical coordinates of the grid
  # zs   : array of ice surface elevations
  # zb   : array of ice bed elevations
  # u    : array of ice velocities in the x-direction
  # v    : array of ice velocities in the y-direction
  # frac : optional parameter; the fraction of the driving stress that the
  #        basal shear stress is assumed to support
    
  # Outputs:
  # beta : basal sliding coefficient, under the shallow ice approximation
  # ub   : basal sliding velocity in the x-direction
  # vb   : basal sliding velocity in the y-direction

  # Constants
  yearinsec = 365.25 * 24 * 60 * 60
  rho = 917 * 1.0e-6 / yearinsec**2
  g = 9.81 * yearinsec**2
  A = 3.5e-25*yearinsec*1.0e18 # For -10 deg ice
  E = 3.0 #Enhancement factor

  nx = len(x)
  ny = len(y)
  dx = x[1] - x[0]
  dy = y[1] - y[0]

  #############################
  # Compute the surface slope #
  #############################

  # Smooth the ice surface elevation
  for i in range(1, ny - 1):
    for j in range(1, nx - 1):
      zs[i, j] = (4 * zs[i, j] +
      	zs[i + 1, j] + zs[i - 1, j] +
        zs[i, j + 1] + zs[i, j - 1]) / 8.0
  
  dsdx = np.zeros((ny, nx))
  dsdy = np.zeros((ny, nx))

  for i in range(1, ny - 1):
    for j in range(1, nx - 1):
       dsdx[i, j] = 0.5 * (zs[i, j + 1] - zs[i, j - 1]) / dx
       dsdy[i, j] = 0.5 * (zs[i + 1, j] - zs[i - 1, j]) / dx

  # Compute the magnitude of the surface slope, the average slope and the 
  # standard deviation of the slope
  ds = np.sqrt(dsdx**2 + dsdy**2)
  avg = np.average(ds)
  stddev = np.std(ds)

  # The maximum slope we'll allow, for the purposes of coming up with a
  # bed sliding velocity, is the average + 1 standard deviation.
  mva = avg + 0.25 * stddev

  # Scale the surface slope at any point where it's too steep
  for i in range(ny):
    for j in range(nx):
      if ds[i, j] > mva:
        dsdx[i, j] = mva / ds[i, j] * dsdx[i, j]
        dsdy[i, j] = mva / ds[i, j] * dsdy[i, j]
        ds[i, j] = mva


  #---------------------------------------------------------
  # Compute the bed sliding velocities & friction parameter
  ub = np.zeros((ny, nx))
  vb = np.zeros((ny, nx))
  beta = np.zeros((ny, nx)) 

  q = 0.0
  speed = 0.0

  for i in range(1, ny - 1):
    for j in range(1, nx - 1):
      if not(np.isnan(u[i,j])) and not(np.isnan(zb[i,j])) and not(np.isnan(zs[i,j])):
        alpha = frac
        h = max(zs[i, j] - zb[i, j], 0.0)
        # internal deformation from SIA, Paterson pg. 310, eq 8.35
        q = E * A * (rho * g * h)**3 * ds[i, j]**3 / 2
        
        # Measured surface speed
        speed = np.sqrt(u[i, j]**2 + v[i, j]**2)
        
        # Inferred basal speed if bed supports ~1/2 of driving stress
        basal_speed = speed - alpha**3*h*q
        if basal_speed <= 0.0:
          basal_speed = min(10.0, 0.1 * speed)
          alpha = ((speed - basal_speed) / (h*q))**(1.0/3)

        # The basal sliding velocities are assumed to have the same
        # direction as the surface velocities, only with lower speed
        # according to a rough SIA-like approximation.
        ub[i, j] = basal_speed/speed * u[i, j]
        vb[i, j] = basal_speed/speed * v[i, j]

        # Since we've already guessed the sliding speed and the
        # x-z strain rate from the SIA, the boundary condition
        #     tau_xz = -beta**2 * u    
        # gives us the value of beta consistent with the guesses
        # we've already made.
        
        # Note that this is the square root of beta, so that we can use it as an initial 
        # guess for beta in inversions.
        beta[i, j] = (2*alpha**3*q / (A*basal_speed**3))**(1.0/6)
        if np.isnan(beta[i,j]):
          beta[i,j] = 1.0e-2
      else:
        ub[i, j] = -2.0e+9
        vb[i, j] = -2.0e+9
        beta[i, j] = -2.0e+9
      
  return ub,vb,beta
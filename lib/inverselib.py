def guess_beta(x,y,zs,zb,u,v,frac,A=3.5e-25*365.25*24*60*60*1.0e18):   
  ''' 
  ub,vb,beta = guess_beta(x,y,zs,zb,u,v,frac,A=3.5e-25*365.25*24*60*60*1.0e18)
  
  Guess a basal friction coefficient beta for a particular glacier geometry assuming
  the bed supports a fraction (frac) of the driving stress.

  Inputs:
  x    : list of horizontal coordinates of the grid
  y    : list of vertical coordinates of the grid
  zs   : array of ice surface elevations
  zb   : array of ice bed elevations
  u    : array of ice velocities in the x-direction
  v    : array of ice velocities in the y-direction
  frac : the fraction of the driving stress that the
          basal shear stress is assumed to support
  A    : flow law parameter
    
  Outputs:
  beta : basal sliding coefficient, under the shallow ice approximation
  ub   : basal sliding velocity in the x-direction
  vb   : basal sliding velocity in the y-direction
  '''

  import numpy as np
  import math

  # Constants
  yearinsec = 365.25 * 24 * 60 * 60
  rho = 917 * 1.0e-6 / yearinsec**2
  g = 9.81 * yearinsec**2
  #A = 3.5e-25*yearinsec*1.0e18 # For -10 deg ice
  E = 3.0 #Enhancement factor

  nx = len(x)
  ny = len(y)
  dx = x[1] - x[0]
  dy = y[1] - y[0]

  # Check if arrhenius parameter is a constant or matrix
  if len(A) < 2:
    A = np.ones([ny,nx])*A



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
        q = E * A[i,j] * (rho * g * h)**3 * ds[i, j]**3 / 2
        
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
        beta[i, j] = (2*alpha**3*q / (A[i,j]*basal_speed**3))**(1.0/6)
        if np.isnan(beta[i,j]):
          beta[i,j] = 1.0e-2
      else:
        ub[i, j] = -2.0e+9
        vb[i, j] = -2.0e+9
        beta[i, j] = -2.0e+9
      
  return ub,vb,beta

def get_velocity_cutoff(glacier,velocity_cutoff=1000,temperature='model',model_dir='',SSA=True,sign='over'):
    '''
    x_grid,y_grid,vsurfini_grid,ind_cutoff_grid,ind_cutoff = get_velocity_cutoff(glacier,
    velocity_cutoff=1000,temperature='model',model_dir='',SSA=True,sign='over')
    
    Find nodes and grid indices where the velocity remains above a particular cutoff value, since
    inversion results seem to be better for higher velocities.

    Inputs:
    glacier         : glacier name
    velocity_cutoff : cutoff value for velocity (m/yr)
    temperature     : model temperature
    model_dir       : add directory if the model results are NOT located in $MODEL_HOME/glacier/3D/
    SSA             : if the model is SSA
    sign            : find indices where values remain "over" or "under" the cutoff
    
    Outputs:
    x_grid          : x values for vsurfini_grid
    y_grid          : y values for vsurfini_grid
    vsurfini_grid   : grid of minimum velocities through time
    ind_cutoff_grid : grid indices that are below the velocity cutoff
    ind_cutoff      : node indices that remain above the velocity cutoff
    '''

    import os, scipy, elmerreadlib
    import numpy as np
    import matplotlib.pyplot as plt

    # Get model result directories
    maindir = os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/"+model_dir)
    dirs = os.listdir(maindir)

    # Set up various parameters for pulling the correct model results
    if temperature == 'model':
        temperature_text = 'modelT'
    else:
        temperature_text = 'constantT'
    if SSA:
        model_text = '1e13_SSA'
    else:
        model_text = '1e12_FS'

    # Check to make sure we have an acceptable value for "sign". If not, exit code.
    if not(sign == 'under') and not(sign == 'over'):
        sys.exit("Unacceptable value for sign of "+sign)

    # Get indices where velocity is always greater than the cutoff value
    n = 0
    for dir in dirs:
        if dir.startswith('DEM') and dir.endswith(temperature_text):
            beta_date = dir[3:11]
            beta_suffix = model_text+'_DEM'+beta_date+'_'+temperature_text+'_'
            data = elmerreadlib.pvtu_file(maindir+dir+'/mesh2d/steady_'+beta_suffix+'linear0001.pvtu',\
                ['vsurfini'])
            surf = elmerreadlib.values_in_layer(data,'surf')
            
            if n == 0:
                surf_min = np.zeros(surf.shape,dtype=np.dtype(surf.dtype.descr))
                for name in surf.dtype.descr:
                    surf_min[name[0]] = surf[name[0]]
                extent = np.loadtxt(maindir+dir+'/inputs/mesh_extent.dat')
                if glacier == 'Helheim':
                    hole1 = np.loadtxt(maindir+dir+"/inputs/mesh_hole1.dat")
                    hole2 = np.loadtxt(maindir+dir+"/inputs/mesh_hole2.dat")
                    holes=[hole1,hole2]
                else:
                    holes = []

            # Save minimum velocity through time
            if sign == 'over':
                ind = np.where(surf_min['vsurfini'] >= surf['vsurfini'])[0]
            elif sign == 'under':
                ind = np.where(surf_min['vsurfini'] <= surf['vsurfini'])[0]
            if len(ind) > 0:
                surf_min['vsurfini'][ind] = surf['vsurfini'][ind]

            n = n+1

    # Grid and filter minimum velocity. Filtering removes some of the spurious single grid cells
    # that remain above the cutoff value.
    x_grid,y_grid,vsurfini_grid = elmerreadlib.grid3d(surf_min,'vsurfini',extent=extent,holes=holes)
    vsurfini_grid = scipy.ndimage.filters.gaussian_filter(vsurfini_grid,sigma=2.5,truncate=4)

    # Find grid indices that are below the cutoff value
    if sign == 'over':
        ind_cutoff_grid = np.where(vsurfini_grid <= velocity_cutoff)
    elif sign == 'under':
        ind_cutoff_grid = np.where(vsurfini_grid >= velocity_cutoff)

    # Find nodes that remain above the cutoff value
    if sign == 'over':
        ind_cutoff = np.where(surf_min['vsurfini'] >= velocity_cutoff)[0]
    elif sign == 'under':
        ind_cutoff = np.where(surf_min['vsurfini'] <= velocity_cutoff)[0]

    return x_grid, y_grid, vsurfini_grid, ind_cutoff_grid, ind_cutoff

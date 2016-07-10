'''
This module determines flux through a flux gate.

fluxgate_thinning(glacier,fluxgate_filename,bedsource): uses two fluxgates to determine the volume
	change in a box and thereby calculated thinning rates according to the velocity record
dem_thinning(glacier,x,y,zs,time,fluxgate_filename): calculates the measured thinning rate in 
	a box according to the worldview grids (this should give values comparable to fluxgate_thinning,
	plus or minus surface accumulation/melt)
fluxbox_geometry(glacier,fluxgate_filename):
'''

import os
import sys
import shapefile
import distlib, bedlib, zslib, vellib, meshlib, coordlib, masklib, climlib
import shapely.geometry
import numpy as np
import scipy.interpolate, scipy
from matplotlib import path

def fluxgate(glacier,fluxgate_filename,bedsource='smith',dl=20.0,timing='velocity'):

  '''
  
  time, sumQ, Hbar, ubar, error = fluxgate(glacier,fluxgate_filename,bedsource='smith',dl=20.0,timing='velocity')
  
  Inputs:
  glacier: glacier name
  fluxgate_filename: fluxgate filename
  bedsource: should be use CreSIS radar transects or morlighem/smith bed DEM to define the bed elevation?
  dl: size of blocks for integrating ice flux across fluxgate
  timing: calculate flux either when we have velocities ('velocity') or surface elevaions ('elevation')
  '''
  
  # Get Shapefiles for flux gates and for glacier extent (so we can calculate surface area accurately)
  DIR = os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/FluxGates/"+glacier+"/")
  gate_sf = shapefile.Reader(DIR+fluxgate_filename)
   
  # Get end points for flux gates
  gate_pts = np.array(gate_sf.shapes()[0].points)
  sortind = np.argsort(gate_pts[:,1])
  gate_pts = gate_pts[sortind]

  # Calculate length of flux gates (i.e., the "width" of glacier)
  L = distlib.transect(gate_pts[:,0],gate_pts[:,1])
  
  # Get coordinates of points along flux gate
  l = np.linspace(0,L[-1],np.ceil(L[-1]/dl)+1)
  dl = l[1]-l[0]
  l = l[0:-1]+dl/2
  x = np.interp(l,L,gate_pts[:,0])
  y = np.interp(l,L,gate_pts[:,1])
  
  # Get bed elevation along flux gates
  if bedsource == 'morlighem':
    print "Using Morlighem bed DEM, should maybe be changed to CreSIS radar transects or Smith bed DEM"
    zb = bedlib.morlighem_pts(x,y,verticaldatum='ellipsoid')
  elif bedsource == 'smith':
    print "Using Smith bed DEM"
    zb = bedlib.smith_at_pts(x,y,glacier,verticaldatum='ellipsoid')
  else:
    print "Using CreSIS radar products"
    cresis_all = bedlib.cresis('all',glacier,verticaldatum='ellipsoid')

    # Find radar pick indices for flux gate     
    zb = np.zeros_like(x)
    # Find points that are within 200~m of CreSIS transect for interpolation
    ind = []
    for i in range(0,len(x)):
      d = np.min(np.sqrt((x[i]-cresis_all[:,0])**2+(y[i]-cresis_all[:,1])**2))
      if d < 200.0:
        ind.append(i)       
    zb[ind] = scipy.interpolate.griddata(cresis_all[:,0:2],cresis_all[:,2],(x[ind],y[ind]))
    f = scipy.interpolate.interp1d(l[ind],zb[ind],kind='cubic',bounds_error=False)
    fex = extrap1d(f)
    zb[zb == 0] = fex(l[zb==0])
    filt_len=500.0
    cutoff=(1/filt_len)/(1/((l[1]-l[0])*2))
    b,a=scipy.signal.butter(4,cutoff,btype='low')
    zb = scipy.signal.filtfilt(b,a,zb)
    
  # Get surface elevations
  zs_all,zs_error,ztime=zslib.dem_at_pts(x,y,glacier,verticaldatum='ellipsoid',method='linear')
  
  # Get velocities
  vpt,tpt,ept,vxpt,vypt = vellib.velocity_at_eulpoints(x,y,glacier,data='TSX',xy_velocities='True')
  
  # Get normal to fluxgate so we can calculate flux through it
  xperp = np.zeros(len(l))
  yperp = np.zeros(len(l))
  for i in range(0,len(l)-1):
    xperp[i+1] = -(-y[i+1]+y[i])/np.sqrt((y[i+1]-y[i])**2+(x[i+1]-x[i])**2)
    yperp[i+1] = -(x[i+1]-x[i])/np.sqrt((y[i+1]-y[i])**2+(x[i+1]-x[i])**2)
  xperp[0] = xperp[1]
  yperp[0] = yperp[1]

  # Find dates with complete velocity profiles. We will use these profiles to create a 
  # "shape factor" to fill in discontinuous velocity records.
  ind=[]
  for i in range(0,len(tpt)):
    nans = np.where(np.isnan(vxpt[i,:]))[0]
    if len(nans) < 1:
      ind.append(i)
  
  # Velocity shape factor
  sf = np.nanmean(abs(xperp*vxpt[ind,:]+yperp*vypt[ind,:]),0)/np.max(np.nanmean(abs(xperp*vxpt[ind,:]+yperp*vypt[ind,:]),0))
  vperp = (xperp*vxpt+yperp*vypt)
  vperp[:,0] = sf[0]*np.nanmax(vperp,1)
  vperp[:,-1] = sf[-1]*np.nanmax(vperp,1)
  ind = np.where(vperp < 0)
  vperp[ind] = 0.0

  midind = np.argmax(np.nanmean(vperp,axis=0))

  if timing == 'velocity':
    # Find DEMs where more than 90% of the data points are nonnan.
    ind_DEM = []
    for i in range(0,len(ztime)):
      nonnan = np.where(~(np.isnan(zs_all[i,:])))[0]
      if len(nonnan) > (9./10.)*len(l):
        ind_DEM.append(i)
  
    Q = np.zeros_like(vperp) # actual fluxes
    Q[:,:] = float('nan')
    allH = np.zeros_like(vperp) # actual fluxes
    allH[:,:] = float('nan')
    allU = np.zeros_like(vperp) # actual fluxes
    allU[:,:] = float('nan')
    time = tpt
    error = np.zeros_like(time)
    error[:] = float('nan')
    ubar = np.zeros_like(time)
    ubar[:] = float('NaN')
    width = np.zeros_like(time)
    width[:] = float('NaN')
    Hbar = np.zeros([len(time),2])
    Hbar[:,:] = float('NaN')
    Across = np.zeros_like(time)
    Across[:] = float('nan')
    for i in range(0,len(time)):
      # Find places where we have no surface velocities
      nans_vperp = np.where(np.isnan(vperp[i,:]))[0]
      if len(nans_vperp) < 1.0/3.0*len(l):
        # If the number of locations without surface velocities is small, let's use the known 
        # surface velocities to interpolate.
        nonnans_vperp = np.where(~(np.isnan(vperp[i,:])))[0]
        vperp_errors = np.zeros(len(l))
              
        # Get surface elevation for that timestep
        zs = np.zeros(len(l))
        zs_error_ind = np.zeros(len(l))

        # Linearly interpolate surface elevations from the available ones.
        zs = np.zeros_like(l)
        for j in range(0,len(zs)):
          nonnan = np.where(~(np.isnan(zs_all[ind_DEM,j])))[0]
          newind = [ind_DEM[k] for k in nonnan]
          zs[j] = np.interp(tpt[i],ztime[newind],zs_all[newind,j])
        zs_error_ind = zs_error
        
        # Interpolate surface velocities using a shape factor for the 
        # cross-sectional velocity profile
        f = np.arange(np.nanmin(np.nanmax(vperp[:,:],axis=1))-500,np.nanmax(np.nanmax(vperp[:,:],axis=1))+500,10)
        mindiff = 1e20
        for j in range(0,len(f)):
          diff = np.sqrt(np.mean(sf[nonnans_vperp]*f[j]-vperp[i,nonnans_vperp])**2)
          if diff < mindiff:
            mindiff = float(diff)
            fbest = int(j)
        vperp[i,nans_vperp] = sf[nans_vperp]*f[fbest]
        
        # Set velocity errors
        vperp_errors[nonnans_vperp] = vperp[i,nonnans_vperp]*0.03
        vperp_errors[nans_vperp] = vperp[i,nans_vperp]*0.1
      
        # Let's try filtering data before calculating ice flux
        filt_len = 300.
        
        cutoff=(1/filt_len)/(1/(np.diff(l[1:3])*2))
        b,a=scipy.signal.butter(4,cutoff,btype='low')
        zs_filt = scipy.signal.filtfilt(b,a,zs)
        vperp_filt = scipy.signal.filtfilt(b,a,vperp[i,:])
        
        # Calculate fluxes, only where zs > zb 
        Q[i,:] = 0.0
        ind = np.where((zs_filt-zb > 0) & (vperp_filt > 0))[0]
        Q[i,ind] = ((vperp_filt[ind])*(zs_filt[ind]-zb[ind])*dl)
        
        # Get average surface elevations and ice flow velocities for fluxgate
        sf_flux = Q[i,:]/np.nanmax(Q[i,:])
        #Hbar[i,0] = np.nansum(sf_flux*(zs_filt-zb))/np.sum(sf_flux)
        Hbar[i,0] = np.mean(zs_filt[midind-200:midind+201]-zb[midind-200:midind+201])
        Hbar[i,1] = np.mean(zs_error_ind[midind-200:midind+201])/np.sqrt(2000./300.)
        ubar[i] = np.mean(vperp_filt[midind-200:midind+201])
        #ubar[i] = np.nansum(sf_flux*(vperp_filt))/np.sum(sf_flux)
        Across[i] = np.sum((zs_filt[ind]-zb[ind])*dl)
        
        allH[i,ind] = zs_filt[ind]-zb[ind]
        allU[i,ind] = (vperp_filt[ind])
        
        # Calculate errors
        ubar_error = np.sqrt(1/(np.sum(np.mean(1/vperp_errors)**2*np.ones(int(L[1]/500.)))))
        hbar_error = np.sqrt(1/(np.sum(np.mean(1/zs_error_ind)**2*np.ones(int(L[1]/32.)))))
        error[i] = np.sum(Q[i,:])*np.sqrt((ubar_error/ubar[i])**2+(hbar_error/Hbar[i,0])**2)
      
  elif timing == 'elevation':
    time = ztime
    Q = np.zeros([len(time),len(l)]) # actual fluxes
    Q[:,:] = float('NaN')
    error = np.zeros_like(time)
    error[:] = float('nan')
    ubar = np.zeros_like(time)
    ubar[:] = float('NaN')
    Hbar = np.zeros([len(time),2])
    Hbar[:,:] = float('NaN')
    Across = np.zeros_like(time)
    Across[:] = float('nan')
    for i in range(0,len(time)):
      nonnans = np.where(~(np.isnan(zs_all[i,:])))[0]
      nans = np.where((np.isnan(zs_all[i,:])))[0]
      if len(nonnans) > (9./10)*len(l):
        zs_ind = zs_all[i,:]
        zs_ind[nans] = np.interp(l[nans],l[nonnans],zs_all[i,nonnans])      
        ind = np.argmin(abs(time[i]-tpt))
        if (time[i]-tpt[ind]) < 1/12.:
          nonnans = np.where(~(np.isnan(vperp[ind,:])))[0]
          nans = np.where(np.isnan(vperp[ind,:]))[0]
          
          vperp_ind = vperp[ind,:]
          vperp_ind[nans] = sf[nans]*np.nanmax(vperp[i,:]) #

          # Calculate fluxes  
          Q[i,:] = (vperp_ind*(zs_ind-zb)*dl)  
      
          # We don't want negative fluxes so let's toss them out.
          Q[i,zs_ind < zb] = 0
          
          # Get average surface elevations and ice flow velocities for fluxgate
          ind = np.where(zs_ind-zb > 0)[0]
          Hbar[i,0] = np.mean(zs_ind[midind-200:midind+201]-zb[midind-200:midind+201])
          Hbar[i,1] = zs_error[i]/np.sqrt(2000./200.)
          ubar[i] = np.mean(vperp_ind[midind-200:midind+201])
          Across[i] = np.sum((zs_ind[ind]-zb[ind])*dl)

  # Add up fluxes along glacier width
  sumQ = np.sum(Q,1)

  return time, sumQ, Hbar, ubar, error, L[1]


def fluxgate_thinning(glacier,fluxgate_filename,bedsource='smith',dl=20.0,timing='elevation',rho_i=900.):

  '''
  Inputs:
  glacier: glacier name
  fluxgate_filename: fluxgate filename, minus the "_in" or "_out" which defines the extent of the box
  bedsource: should be use CreSIS radar transects or morlighem/smith bed DEM to define the bed elevation?
  dl: size of blocks for integrating ice flux across fluxgate
  timing: calculate dH at time steps when we have ice flow velocity measurements ('velocity') 
  	or when we have surface elevation measurements ('elevation')
  rho_i: ice density for calculating dH/dt due to surface mass balance at the same time steps	
  
  '''

  # Get fluxbox geometry
  fluxbox_x,fluxbox_y = fluxbox_geometry(glacier,fluxgate_filename)

  # Get surface mass balance for flux box 
  xrac,yrac,zsrac,timerac = climlib.racmo_at_pts(np.mean(fluxbox_x),np.mean(fluxbox_y),'smb',filt_len='none')
  zsrac=zsrac*365.25/rho_i

  # Get fluxbox geometry
  fluxbox_x,fluxbox_y = fluxbox_geometry(glacier,fluxgate_filename)

  # Get flux through fluxgates
  time_in,Q_in,hbar_in,ubar_in,error_in,width_in = fluxgate(glacier,fluxgate_filename+'_in',bedsource=bedsource,dl=dl,timing=timing)
  time_out,Q_out,hbar_out,ubar_out,error_out,width_out = fluxgate(glacier,fluxgate_filename+'_out',bedsource=bedsource,dl=dl,timing=timing)

  # Get surface area
  shape = shapely.geometry.Polygon(np.column_stack([fluxbox_x,fluxbox_y]))
  A = shape.area
  
  time = []
  hbar = []
  ubar = []
  Q = []
  smb = []
  error_Q = []
  for i in range(0,len(time_in)):
    if time_in[i] in time_out:
      j = np.argmin(abs(time_in[i]-time_out))
      time.append(time_in[i])
      hbar.append([hbar_in[i,0],hbar_out[j,0],hbar_in[i,1],hbar_out[j,1]])
      ubar.append([ubar_in[i],ubar_out[j]])
      Q.append([Q_in[i],Q_out[j]])
      error_Q.append(error_in[i]+error_out[j])
      minind = np.where(abs(time_in[i]-timerac) < 5/365.)[0]
      #smb.append(np.mean(np.diff(zsrac[minind])/np.diff(timerac[minind])))
      smb.append(np.mean(zsrac[minind]))
  hbar = np.array(hbar)
  ubar = np.array(ubar)
  time = np.array(time)
  smb = np.array(smb)
  error_Q = np.array(error_Q)
  Q = np.array(Q)

  dH = np.column_stack([(Q[:,0]-Q[:,1])/A,(error_Q)/A])
  
  widths = np.column_stack([width_in,width_out])
  
  return time, dH, Q, hbar, ubar, smb, widths,A

def dem_thinning(glacier,x,y,zs,time,zserror,fluxgate_filename,type='rate'):
  
  '''
  Calculate the measured thinning rates in the fluxboxes so that we can compare the measured
  values to the inferred values from the change in flux (fluxgate_thinning)
  '''
  
  # Grid points
  xgrid,ygrid = np.meshgrid(x,y)
  
  # Flux box vertices
  fluxbox_x,fluxbox_y = fluxbox_geometry(glacier,fluxgate_filename)
  box = path.Path(np.column_stack([fluxbox_x,fluxbox_y]))
  
  # Find indices in the grid that fall within the fluxbox
  inside = box.contains_points(np.column_stack([xgrid.flatten(),ygrid.flatten()]))
  ninside = len(xgrid.flatten()[inside])
  
  # Calculate thinning rates from DEMs
  nt = len(time)
  dH_time_mid = []
  dH_time_extent = []
  dH_ave = []
  dH_error = []
  if type=='rate':
    for i in range(0,nt):
      for j in range(i+1,nt):
        if i != j:
          # Surface elevations inside the fluxbox at times 1 and 2
          zs_inside_t1 = zs[:,:,i].flatten()[inside]
          zs_inside_t2 = zs[:,:,j].flatten()[inside]
  
          # Find locations where we have thinning rates
          nonnan = np.where(~(np.isnan(zs_inside_t2-zs_inside_t1)))[0]
          if (len(nonnan) > 0.9*ninside) and ((time[j]-time[i]) > 1/12.0) and ((time[j]-time[i]) < 1/3.):
            # dH_time_extent.append(float(len(nonnan))/float(ninside))
            # Calculate time as halfway between the two times
            dH_time_mid.append((time[i]+time[j])/2) # time
            dH_time_extent.append(abs(time[j]-time[i])/2) # time interval
    
            # Only calculate the average thinning rates if the measured thinning rates cover a 
            # sufficient area of the fluxbox and surface elevation profiles are more than a month apart.
            dH_ave.append(np.nanmean(zs_inside_t2-zs_inside_t1)/(time[j]-time[i])) # average thinning rate
            dH_error.append(np.sqrt(zserror[i]**2+zserror[j]**2)/(time[j]-time[i])) # error according to co-registration 
    
    dH_time = np.column_stack([dH_time_mid,dH_time_extent])
    dH = np.column_stack([dH_ave,dH_error])
  elif type == 'absolute':
    dH = np.zeros([nt,2])
    dH_time = np.zeros([nt,1])
    for i in range(0,nt):
      # Calculate time as halfway between the two times
      dH_time[i,0] = time[i] # time
    
      # Surface elevations inside the fluxbox at times 1 and 2
      zs_inside = zs[:,:,i].flatten()[inside]
    
      # Find locations where we have thinning rates
      nonnan = np.where(~(np.isnan(zs_inside)))[0]
    
      # Only calculate the average thinning rates if the measured thinning rates cover a 
      # sufficient area of the fluxbox and surface elevation profiles are more than a month apart.
      if len(nonnan) > 0.9*ninside:
        dH[i,0] = np.nanmean(zs_inside)
        dH[i,1] = np.nanstd(zs_inside)
      else:
        dH[i,0] = float('NaN')
        dH[i,1] = float('NaN')
    
    
  return dH_time,dH
  
def fluxbox_geometry(glacier,fluxgate_filename):

  '''
  Find the geometry of the fluxbox by combining the shapefiles for the fluxgates and for the glacier 
  extent.
  '''
  
  # Get Shapefiles for flux gates and for glacier extent (so we can calculate surface area accurately)
  DIR = os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/FluxGates/"+glacier+"/")
  gate_in_sf = shapefile.Reader(DIR+fluxgate_filename+"_in")
  gate_out_sf = shapefile.Reader(DIR+fluxgate_filename+"_out")
  
  # Glacier extent shapefile
  glacier_extent = masklib.load_points(glacier)
  #glacier_extent = meshlib..shp_to_xy(os.path.join(os.getenv("DATA_HOME"),
  #		"ShapeFiles/Glaciers/3D/"+glacier+"/glacier_extent_normal"))
  #if glacier == 'Helheim':
  #  glacier_hole1 = meshlib..shp_to_xy(os.path.join(os.getenv("DATA_HOME"),
  #  	"ShapeFiles/Glaciers/3D/Helheim/glacier_hole1"))
  #  glacier_hole2 = meshlib..shp_to_xy(os.path.join(os.getenv("DATA_HOME"),
  #  	"ShapeFiles/Glaciers/3D/Helheim/glacier_hole2"))
  #  glacier_extent = np.row_stack([glacier_extent,glacier_hole1,glacier_hole2])
  
  # Get end points for flux gates
  gate_in_pts = np.array(gate_in_sf.shapes()[0].points)
  if gate_in_pts[-1,1] < gate_in_pts[0,1]:
    gate_in_pts = np.flipud(gate_in_pts)
  gate_out_pts = np.array(gate_out_sf.shapes()[0].points)
  if gate_out_pts[-1,1] < gate_out_pts[0,1]:
    gate_out_pts = np.flipud(gate_out_pts)
  
    # Calculate length of flux gates (i.e., the "width" of glacier)
  L_in = distlib.transect(gate_in_pts[:,0],gate_in_pts[:,1])
  L_out = distlib.transect(gate_out_pts[:,0],gate_out_pts[:,1])
  
  dl = 50.
  
  # Get coordinates of points along upstream flux gate
  l_in = np.linspace(0,L_in[-1],np.ceil(L_in[-1]/dl)+1)
  dl_in = l_in[1]-l_in[0]
  l_in = l_in[0:-1]+dl_in/2
  x_in = np.interp(l_in,L_in,gate_in_pts[:,0])
  y_in = np.interp(l_in,L_in,gate_in_pts[:,1])
  
  # Get coordinates of points along downstream flux gate
  l_out = np.linspace(0,L_out[-1],np.ceil(L_out[-1]/dl)+1)
  dl_out = l_out[1]-l_out[0]
  l_out = l_out[0:-1]+dl_out/2
  x_out = np.interp(l_out,L_out,gate_out_pts[:,0])
  y_out = np.interp(l_out,L_out,gate_out_pts[:,1])
  
  gate_in_pts = np.column_stack([x_in,y_in])
  gate_out_pts = np.column_stack([x_out,y_out])
  
  del l_out,l_in,dl_in,dl_out,x_in,x_out,y_in,y_out
  
  # Find glacier geometry (i.e., points in glacier extent) that lie between the flux gates, 
  # so we can calculate glacier area.
  
  fluxbox_x = [] # boundary points for area between flux gates 
  fluxbox_x.append(gate_in_pts[:,0].T.tolist())
  fluxbox_y = []
  fluxbox_y.append(gate_in_pts[:,1].T.tolist())
  for i in [0,-1]:
    ind_in = np.argmin((gate_in_pts[i,0]-glacier_extent[:,0])**2+(gate_in_pts[i,-1]-glacier_extent[:,1])**2) # get closest index
    ind_out = np.argmin((gate_out_pts[i,0]-glacier_extent[:,0])**2+(gate_out_pts[i,-1]-glacier_extent[:,1])**2) # get closest index
  
    # Check to make sure we didn't grab a point outside the flux box
    # The x-value of glacier extent near the upstream flux gate needs to be greater than the
    # x-value of the flux gate to be in the box. Check to see if it's smaller and if so correct it.
    if gate_in_pts[i,0] > glacier_extent[ind_in,0]: 
      if gate_in_pts[i,0] < glacier_extent[ind_in-1,0]:
        ind_in = ind_in-1
      else: 
        ind_in = ind_in+1 
    # The opposite is true for the downstream gate. 
    if gate_out_pts[i,0] < glacier_extent[ind_out,0]: 
      if gate_out_pts[i,0] > glacier_extent[ind_out-1,0]:
        ind_out = ind_out-1
      else: 
        ind_out = ind_out+1 
  
    # Add points from glacier extent to flux box
    if ind_in < ind_out:
      fluxbox_x.append(glacier_extent[ind_in:ind_out+1,0].tolist())
      fluxbox_y.append(glacier_extent[ind_in:ind_out+1,1].tolist())
    else:
      fluxbox_x.append(glacier_extent[ind_out:ind_in+1,0].tolist())  
      fluxbox_y.append(glacier_extent[ind_out:ind_in+1,1].tolist()) 
  fluxbox_x.append(gate_out_pts[:,0].T.tolist())
  fluxbox_y.append(gate_out_pts[:,1].T.tolist())
  fluxbox_x = sum(fluxbox_x,[])
  fluxbox_y = sum(fluxbox_y,[])
  fluxbox_x,fluxbox_y = coordlib.sort_xy(fluxbox_x,fluxbox_y)
  fluxbox_x.append(fluxbox_x[0])
  fluxbox_y.append(fluxbox_y[0])  
  
  return fluxbox_x,fluxbox_y

def compare_thinning_rates(demtime,demdH,fluxtime,fluxdH,smbdH):

  '''
  Compare the calculated thinning rates from the DEMs and from the fluxgate method.
  '''
  
  def nearest_interp(xi, x, y):
    idx = np.abs(x - xi[:,None])
    return y[idx.argmin(axis=1)]
  
  # Set up time
  dt = (1/365.25)

  starttime = demtime[0,0]-demtime[0,1]
  endtime = demtime[-1,0]+demtime[-1,1]
  Nt = int(np.ceil(endtime-starttime)/dt)+1
  
  time = np.arange(starttime,starttime+Nt*dt,dt)
  
  # Add surface mass balance to fluxgate thinning rates so that we can compare these rates
  # directly to the dem difference rates
  fluxsmbdH = fluxdH[:,0]+smbdH
  
  # Set up variables
  dH_time = np.array(demtime)
  dH_dem = np.array(demdH)
  dH_flux = np.zeros_like(demdH)
  dH_flux[:,:] = float('nan')
  
  # Iterate through DEM times for comparison
  for i in range(0,len(demtime[:,0])):
    # Check to make sure there is at least one estimated thinning rate from the 
    # velocities within the time period when we have a DEM thinning rate.
    inside = np.where((fluxtime > demtime[i,0]-demtime[i,1]) & (fluxtime < demtime[i,0]+demtime[i,1]))[0]
    if (len(inside) > 1):
      if len(np.where(~(np.isnan(fluxdH[inside])))[0]) > 1:
        ind = np.where((time > demtime[i,0]-demtime[i,1]) & (time < demtime[i,0]+demtime[i,1]))[0]
        nonnan = np.where(~(np.isnan(fluxdH[:,0])))[0]
        #values = nearest_interp(time[ind],fluxtime[nonnan],fluxdH[nonnan,0])
        values = np.interp(time[ind],fluxtime,fluxsmbdH)
        errors = np.interp(time[ind],fluxtime,fluxdH[:,1])
        dH_flux[i,0] = np.nanmean(values)  
        dH_flux[i,1] = np.nanmean(errors) 

  
  return dH_time,dH_flux,dH_dem
  
def extrap1d(interpolator):
    
  xs = interpolator.x
  ys = interpolator.y

  def pointwise(x):
    if x < xs[0]:
      return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
    elif x > xs[-1]:
      return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
    else:
      return interpolator(x)

  def ufunclike(xs):
    return scipy.array(map(pointwise, scipy.array(xs)))

  return ufunclike
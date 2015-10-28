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
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"Util/Modules"))
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"BigThreeGlaciers/Tools"))
import shapefile
import fracyear, dist, bed, elevation, velocity, elmer_mesh, coords, dist, icemask
import shapely.geometry
import numpy as np
import scipy.interpolate, scipy
from matplotlib import path

def fluxgate_thinning(glacier,fluxgate_filename,bedsource='smith',dl=20.0):

  '''
  Inputs:
  glacier: glacier name
  fluxgate_filename: fluxgate filename, minus the "_in" or "_out" which defines the extent of the box
  bedsource: should be use CreSIS radar transects or morlighem/smith bed DEM to define the bed elevation?
  dl: size of blocks for integrating ice flux across fluxgate
  '''
  
  # Get Shapefiles for flux gates and for glacier extent (so we can calculate surface area accurately)
  DIR = os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/FluxGates/"+glacier+"/")
  gate_in_sf = shapefile.Reader(DIR+fluxgate_filename+"_in")
  gate_out_sf = shapefile.Reader(DIR+fluxgate_filename+"_out")
   
  # Get end points for flux gates
  gate_in_pts = np.array(gate_in_sf.shapes()[0].points)
  sortind = np.argsort(gate_in_pts[:,1])
  gate_in_pts = gate_in_pts[sortind]
  gate_out_pts = np.array(gate_out_sf.shapes()[0].points)
  sortind = np.argsort(gate_out_pts[:,1])
  gate_out_pts = gate_out_pts[sortind]
  
  # Get fluxbox geometry
  fluxbox_x,fluxbox_y = fluxbox_geometry(glacier,fluxgate_filename)

  # Calculate length of flux gates (i.e., the "width" of glacier)
  L_in = dist.transect(gate_in_pts[:,0],gate_in_pts[:,1])
  L_out = dist.transect(gate_out_pts[:,0],gate_out_pts[:,1])
  
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
  
  # Get bed elevation along flux gates
  if bedsource == 'morlighem':
    print "Using Morlighem bed DEM, should maybe be changed to CreSIS radar transects or Smith bed DEM"
    zb_in = bed.morlighem_pts(x_in,y_in,verticaldatum='ellipsoid')
    zb_out = bed.morlighem_pts(x_out,y_out,verticaldatum='ellipsoid')
  elif bedsource == 'smith':
    print "Using Smith bed DEM"
    zb_in = bed.smith_at_pts(x_in,y_in,glacier,verticaldatum='ellipsoid')
    zb_out = bed.smith_at_pts(x_out,y_out,glacier,verticaldatum='ellipsoid')
  else:
    print "Using CreSIS radar products"
    cresis_all = bed.cresis('all',glacier,verticaldatum='ellipsoid')

    # Find radar pick indices for upstream flux gate
    zb_in = np.zeros_like(x_in)
    ind1 = np.argmin((cresis_all[:,0]-x_in[0])**2+(cresis_all[:,1]-y_in[0])**2)
    ind2 = np.argmin((cresis_all[:,0]-x_in[-1])**2+(cresis_all[:,1]-y_in[-1])**2)
    if ind1 < ind2:
      cut_cresis = cresis_all[ind1:ind2,:]
    else: 
      cut_cresis = cresis_all[ind2:ind1,:]
    # Find points that are within 500~m of CreSIS transect for interpolation
    ind = []
    for i in range(0,len(x_out)):
      d = np.min(np.sqrt((x_in[i]-cresis_all[:,0])**2+(y_in[i]-cresis_all[:,1])**2))
      if d < 50.0:
        ind.append(i)      
    zb_in[ind] = scipy.interpolate.griddata(cresis_all[:,0:2],cresis_all[:,2],(x_in[ind],y_in[ind]),method='cubic')
    f = scipy.interpolate.interp1d(l_in[ind],zb_in[ind],kind='cubic',bounds_error=False)
    fex = extrap1d(f)
    zb_in[zb_in == 0] = fex(l_in[zb_in==0])
    filt_len=1000.0
    cutoff=(1/filt_len)/(1/((l_in[1]-l_out[0])*2))
    b,a=scipy.signal.butter(4,cutoff,btype='low')
    zb_in = scipy.signal.filtfilt(b,a,zb_in)

    # Find radar pick indices for downstream flux gate
    zb_out = np.zeros_like(x_out)
    ind1 = np.argmin((cresis_all[:,0]-x_out[0])**2+(cresis_all[:,1]-y_out[0])**2)
    ind2 = np.argmin((cresis_all[:,0]-x_out[-1])**2+(cresis_all[:,1]-y_out[-1])**2)
    if ind1 < ind2:
      cut_cresis = cresis_all[ind1:ind2,:]
    else: 
      cut_cresis = cresis_all[ind2:ind1,:]
    # Find points that are within 500~m of CreSIS transect for interpolation
    ind = []
    for i in range(0,len(x_out)):
      d = np.min(np.sqrt((x_out[i]-cresis_all[:,0])**2+(y_out[i]-cresis_all[:,1])**2))
      if d < 50.0:
        ind.append(i)       
    zb_out[ind] = scipy.interpolate.griddata(cresis_all[:,0:2],cresis_all[:,2],(x_out[ind],y_out[ind]),method='cubic')
    f = scipy.interpolate.interp1d(l_out[ind],zb_out[ind],kind='cubic',bounds_error=False)
    fex = extrap1d(f)
    zb_out[zb_out == 0] = fex(l_out[zb_out==0])
    filt_len=1000.0
    cutoff=(1/filt_len)/(1/((l_out[1]-l_out[0])*2))
    b,a=scipy.signal.butter(4,cutoff,btype='low')
    zb_out_filt = scipy.signal.filtfilt(b,a,zb_out)

    
  # Get surface elevations, eventually we'll probably want to let the DEM vary so that 
  # we're getting the correct cross sectional area, but for now let's just use a single 
  # surface elevation profile
  zs_all,zs_std,ztime=elevation.worldview_at_pts(np.r_[x_in,x_out],np.r_[y_in,y_out],glacier,verticaldatum='ellipsoid',method='linear')
  zs_all_in = np.array(zs_all[:,0:len(x_in)])
  zs_std_in = np.array(zs_std[:,0:len(x_in)])
  zs_all_out = np.array(zs_all[:,len(x_in):])
  zs_std_out = np.array(zs_all[:,len(x_in):])
  for i in range(0,len(ztime)):
    nonnan1 = np.where(zs_all_in[i,:])[0]
    nonnan2 = np.where(zs_all_out[i,:])[0]
    if len(nonnan1) < 0.9*len(l_in) or len(nonnan2) < 0.9*len(l_out):
      zs_in_all[i,:] = float('nan')
      zs_out_all[i,:] = float('nan')
  #zs_in = np.nanmedian(zs_all_in,axis=0)
  zs_in = elevation.dem_continuous_flowline(x_in,y_in,glacier,'20110712',verticaldatum='ellipsoid',fillin='true')
  #zs_out = np.nanmedian(zs_all_out,axis=0)
  zs_out = elevation.dem_continuous_flowline(x_out,y_out,glacier,'20110712',verticaldatum='ellipsoid',fillin='true')
  
  # Get velocities
  vpt_in,tpt_in,ept_in,vxpt_in,vypt_in = velocity.velocity_at_eulpoints(x_in,y_in,glacier,data='all',xy_velocities='True')
  vpt_out,tpt_out,ept_out,vxpt_out,vypt_out = velocity.velocity_at_eulpoints(x_out,y_out,glacier,data='all',xy_velocities='True')
  
  # Time for dH/dt calculations
  time = tpt_in
  dH = np.zeros([len(time),2])
  
  # Anticipated error between surface elevations and/or bed DEM at fluxgate_in and 
  # fluxgate_out (total error is 2*error)
  error = 5.0/2.
  
  # Get flux through upstream gate
  xperp = np.zeros(len(l_in))
  yperp = np.zeros(len(l_in))
  for i in range(0,len(l_in)-1):
    xperp[i+1] = (-y_in[i+1]+y_in[i])/np.sqrt((y_in[i+1]-y_in[i])**2+(x_in[i+1]-x_in[i])**2)
    yperp[i+1] = (x_in[i+1]-x_in[i])/np.sqrt((y_in[i+1]-y_in[i])**2+(x_in[i+1]-x_in[i])**2)
  vperp_in = abs(xperp*vxpt_in+yperp*vypt_in)
  vperp_in[:,0] = 0.0
  vperp_in[:,-1] = 0.0
  Q_in = np.zeros_like(vperp_in) # actual fluxes
  Q_in[:,:] = 'NaN'
  Q_in_max = np.zeros_like(vperp_in) # large possible flux given the errors
  Q_in_max[:,:] = 'NaN'
  for i in range(0,len(vperp_in[:,0])):
    # Find places where we have no surface velocities
    nans = np.where(np.isnan(vperp_in[i,:]))[0]
    if len(nans) < 1.0/4.0*len(l_in):
      # If the number of locations without surface velocities is small, let's use the known 
      # surface velocities to interpolate.
      nonnans = np.where(~(np.isnan(vperp_in[i,:])))[0]
      
      #f = scipy.interpolate.interp1d(l_in[nonnans],vperp_in[i,nonnans],kind='cubic')
      #vperp_in[i,nans] = f(l_in[nans])
      #vperp_in[i,vperp_in[i,:] < 0] = np.interp(l_in[vperp_in[i,:] < 0],l_in[nonnans],vperp_in[i,nonnans])
      vperp_in[i,nans] = np.interp(l_in[nans],l_in[nonnans],vperp_in[i,nonnans])
      
      # Calculate fluxes
      Q_in[i,:] = (vperp_in[i,:]*(zs_in-zb_in)*dl_in)
      Q_in_max[i,:] = (vperp_in[i,:]*(zs_in+error-zb_in)*dl_in)
      
      # We don't want negative fluxes so let's toss them out.
      Q_in[i,zs_in < zb_in] = 0
      Q_in_max[i,zs_in < zb_in] = 0  

  # Get flux through downstream gate
  xperp = np.zeros(len(l_out))
  yperp = np.zeros(len(l_out))
  for i in range(0,len(l_out)-1):
    xperp[i+1] = (-y_out[i+1]+y_out[i])/np.sqrt((y_out[i+1]-y_out[i])**2+(x_out[i+1]-x_out[i])**2)
    yperp[i+1] =  (x_out[i+1]-x_out[i])/np.sqrt((y_out[i+1]-y_out[i])**2+(x_out[i+1]-x_out[i])**2)
  vperp_out = abs(xperp*vxpt_out+yperp*vypt_out)
  vperp_out[:,0] = 0.0
  vperp_out[:,-1] = 0.0
  Q_out = np.zeros_like(vperp_out)
  Q_out[:,:] = 'NaN'
  Q_out_max = np.zeros_like(vperp_out)
  Q_out_max[:,:] = 'NaN'
  for i in range(0,len(vperp_out[:,0])):
    # Find places where we have no surface velocities
    nans = np.where(np.isnan(vperp_out[i,:]))[0]
    if len(nans) < 1.0/4.0*len(l_out):
      # If the number of locations without surface velocities is small, let's use the known 
      # surface velocities to interpolate. Otherwise the fluxes are set to nan's.
      nonnans = np.where(~(np.isnan(vperp_out[i,:])))[0]
      
      #f = scipy.interpolate.interp1d(l_out[nonnans],vperp_out[i,nonnans],kind='cubic')
      #vperp_out[i,nans] = f(l_out[nans])
      vperp_out[i,nans] = np.interp(l_out[nans],l_out[nonnans],vperp_out[i,nonnans])
      #vperp_out[i,vperp_out < 0] = np.interp(l_out[vperp_out[i,:]<0],l_out[nonnans],vperp_out[i,nonnans])
      #vperp_out[i,nans] = np.nanmean(vperp_out[i,:])*scale_out[nans]
      
      # Calculate fluxes
      Q_out[i,:] = (vperp_out[i,:]*(zs_out-zb_out)*dl_out)
      Q_out_max[i,:] = (vperp_out[i,:]*(zs_out-error-zb_out)*dl_out)
    
      # We don't want negative fluxes
      Q_in[i,zs_in < zb_in] = 0
      Q_in_max[i,zs_in < zb_in] = 0
  
  # Get surface area
  shape = shapely.geometry.Polygon(np.column_stack([fluxbox_x,fluxbox_y]))
  A = shape.area
  
  # Calculate the change in flux between the two gates
  dQ = np.sum(Q_in,1)-np.sum(Q_out,1)
  dQ_max = np.sum(Q_in_max,1)-np.sum(Q_out_max,1)

  dH[:,0] = dQ/A # anticipated dH/dt
  dH[:,1] = (dQ_max-dQ)/A # error in dH/dt
  
  return time, dH

def dem_thinning(glacier,x,y,zs,time,fluxgate_filename,type='rate'):
  
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
  
  # Calculate average thinning rates
  nt = len(time)
  j = 0
  dH_time_mid = []
  dH_time_extent = []
  dH_ave = []
  dH_error = []
  if type=='rate':
    for i in range(1,nt):
      # Surface elevations inside the fluxbox at times 1 and 2
      zs_inside_t1 = zs[:,:,j].flatten()[inside]
      zs_inside_t2 = zs[:,:,i].flatten()[inside]
  
      # Find locations where we have thinning rates
      nonnan = np.where(~(np.isnan(zs_inside_t2-zs_inside_t1)))[0]    
      
      if (len(nonnan) > 0.8*ninside) and (time[i]-time[j] > 11/365.0):
        #dH_time_extent.append(float(len(nonnan))/float(ninside))
        # Calculate time as halfway between the two times
        dH_time_mid.append((time[i]+time[j])/2) # time
        dH_time_extent.append((time[i]-time[j])/2) # time interval
    
        # Only calculate the average thinning rates if the measured thinning rates cover a 
        # sufficient area of the fluxbox and surface elevation profiles are more than a month apart.
        dH_ave.append(np.nanmean(zs_inside_t2-zs_inside_t1)/(time[i]-time[j])) # average thinning rate
        dH_error.append(np.nanstd(zs_inside_t2-zs_inside_t1)) # standard deviation of thinning rates in box
        j = i
      else:
        if len(np.where(~(np.isnan(zs_inside_t1)))[0]) < 0.8*ninside:
          j = i
    
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
  glacier_extent = icemask.load_points(glacier)
  #glacier_extent = elmer_mesh.shp_to_xy(os.path.join(os.getenv("DATA_HOME"),
  #		"ShapeFiles/Glaciers/3D/"+glacier+"/glacier_extent_normal"))
  #if glacier == 'Helheim':
  #  glacier_hole1 = elmer_mesh.shp_to_xy(os.path.join(os.getenv("DATA_HOME"),
  #  	"ShapeFiles/Glaciers/3D/Helheim/glacier_hole1"))
  #  glacier_hole2 = elmer_mesh.shp_to_xy(os.path.join(os.getenv("DATA_HOME"),
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
  L_in = dist.transect(gate_in_pts[:,0],gate_in_pts[:,1])
  L_out = dist.transect(gate_out_pts[:,0],gate_out_pts[:,1])
  
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
  fluxbox_x,fluxbox_y = coords.sort_xy(fluxbox_x,fluxbox_y)
  fluxbox_x.append(fluxbox_x[0])
  fluxbox_y.append(fluxbox_y[0])  
  
  return fluxbox_x,fluxbox_y

def compare_thinning_rates(demtime,demdH,fluxtime,fluxdH,smbtime,smbdH,rho_i=900.0):

  '''
  Compare the calculated thinning rates from the DEMs and from the fluxgate method.
  '''
  
  def nearest_interp(xi, x, y):
    idx = np.abs(x - xi[:,None])
    return y[idx.argmin(axis=1)]
  
  # Set up time
  dt = (0.2/365.25)

  starttime = demtime[0,0]-demtime[0,1]
  endtime = demtime[-1,0]+demtime[-1,1]
  Nt = int(np.ceil(endtime-starttime)/dt)+1
  
  time = np.arange(starttime,starttime+Nt*dt,dt)
  
  dH_time = np.array(demtime)
  dH_dem = np.array(demdH)
  dH_flux = np.zeros_like(demdH)
  dH_smb = np.zeros(len(demdH))
  dH_flux[:,:] = float('nan')
  for i in range(0,len(demtime[:,0])):
    # Check to make sure there is at least one estimated thinning rate from the 
    # velocities within the time period when we have a DEM thinning rate.
    inside = np.where((fluxtime > demtime[i,0]-demtime[i,1]) & (fluxtime < demtime[i,0]+demtime[i,1]))[0]
    if (len(inside) > 0):
      if len(np.where(~(np.isnan(fluxdH[inside])))[0]) > 0:
        ind = np.where((time > demtime[i,0]-demtime[i,1]) & (time < demtime[i,0]+demtime[i,1]))[0]
        nonnan = np.where(~(np.isnan(fluxdH[:,0])))[0]
        #values = nearest_interp(time[ind],fluxtime[nonnan],fluxdH[nonnan,0])
        values = np.interp(time[ind],fluxtime,fluxdH[:,0])
        errors = np.interp(time[ind],fluxtime,fluxdH[:,1])
        dH_flux[i,0] = np.nanmean(values)  
        dH_flux[i,1] = np.nanmean(errors) 
    
    values = np.interp(time[ind],smbtime,smbdH)
    dH_smb[i] = np.mean(values)*365.25/rho_i
  
  return dH_time,dH_flux,dH_dem,dH_smb
  
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
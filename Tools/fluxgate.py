# This module determines flux through a flux gate.
#
# fluxgate_thinning(glacier,fluxgate_filename,bedsource): uses two fluxgates to determine the volume
#	change in a box and thereby calculated thinning rates according to the velocity record
# dem_thinning(glacier,x,y,zs,time,fluxgate_filename): calculates the measured thinning rate in 
#	a box according to the worldview grids (this should give values comparable to fluxgate_thinning,
#	plus or minus surface accumulation/melt)
# fluxbox_geometry(glacier,fluxgate_filename):

import os
import sys
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
import shapefile
import fracyear, dist, bed, elevation, velocity, elmer_mesh, coords
import shapely.geometry
import numpy as np
import scipy.interpolate
from matplotlib import path

def fluxgate_thinning(glacier,fluxgate_filename,bedsource='morlighem'):

  # Inputs:
  # glacier: glacier name
  # fluxgate_filename: fluxgate filename, minus the "_in" or "_out" which defines the extent of the box
  # bedsource: should be use CreSIS radar transects or morlighem bed DEM to define the bed elevation?
  
  # Size of blocks for integrating ice flux across fluxgate
  dl = 50.0 

  # Get Shapefiles for flux gates and for glacier extent (so we can calculate surface area accurately)
  DIR = os.path.join(os.getenv("HOME"),"Data/ShapeFiles/FluxGates/"+glacier+"/")
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
  L_in = dist.between_pts(gate_in_pts[0,0],gate_in_pts[0,1],gate_in_pts[1,0],gate_in_pts[1,1])
  L_out = dist.between_pts(gate_out_pts[0,0],gate_out_pts[0,1],gate_out_pts[1,0],gate_out_pts[1,1])
  
  # Get coordinates of points along upstream flux gate
  l_in = np.linspace(0,L_in,np.ceil(L_in/dl)+1)
  dl_in = l_in[1]-l_in[0]
  l_in = l_in[0:-1]+dl_in/2
  x_in = np.interp(l_in,[0,L_in],gate_in_pts[:,0])
  y_in = np.interp(l_in,[0,L_in],gate_in_pts[:,1])
  
  # Get coordinates of points along downstream flux gate
  l_out = np.linspace(0,L_out,np.ceil(L_out/dl)+1)
  dl_out = l_out[1]-l_out[0]
  l_out = l_out[0:-1]+dl_out/2
  x_out = np.interp(l_out,[0,L_out],gate_out_pts[:,0])
  y_out = np.interp(l_out,[0,L_out],gate_out_pts[:,1])
  
  # Get bed elevation along flux gates
  if bedsource == 'morlighem':
    print "Using Morlighem bed DEM, should maybe be changed to CreSIS radar transects"
    zb_in = bed.morlighem_pts(x_in,y_in,glacier,verticaldatum='geoid')
    zb_out = bed.morlighem_pts(x_out,y_out,glacier,verticaldatum='geoid')
  else:
    print "Using CreSIS radar products"
    cresis_all = bed.cresis('all',glacier,verticaldatum='geoid')

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
    zb_in[0] = cut_cresis[0,2]
    zb_in[-1] = cut_cresis[-1,2]     
    zb_in[ind] = scipy.interpolate.griddata(cresis_all[:,0:2],cresis_all[:,2],(x_in[ind],y_in[ind]),method='nearest')
    zb_in[zb_in==0] = np.interp(y_in[zb_in==0],y_in[zb_in!=0],zb_in[zb_in!=0])

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
    zb_out[0] = cut_cresis[0,2]
    zb_out[-1] = cut_cresis[-1,2]     
    zb_out[ind] = scipy.interpolate.griddata(cresis_all[:,0:2],cresis_all[:,2],(x_out[ind],y_out[ind]),method='nearest')
    zb_out[zb_out==0] = np.interp(y_out[zb_out==0],y_out[zb_out!=0],zb_out[zb_out!=0])
    
  # Get surface elevations, eventually we'll probably want to let the DEM vary so that 
  # we're getting the correct cross sectional area, but for now let's just use a single 
  # surface elevation profile
  if glacier == 'Helheim':
    zs_in = elevation.dem_continuous_flowline(x_in,y_in,glacier,'20120624',verticaldatum='geoid',fillin='true')
    zs_out = elevation.dem_continuous_flowline(x_out,y_out,glacier,'20120624',verticaldatum='geoid',fillin='true')
  elif glacier == 'Kanger':
    zs_in = elevation.dem_continuous_flowline(x_in,y_in,glacier,'20110712',verticaldatum='geoid',fillin='true')
    zs_out = elevation.dem_continuous_flowline(x_out,y_out,glacier,'20110712',verticaldatum='geoid',fillin='true')
  
  # Get velocities
  vpt_in,tpt_in,ept_in,vxpt_in,vypt_in = velocity.velocity_at_eulpoints(x_in,y_in,glacier,data='all',xy_velocities='True')
  vpt_out,tpt_out,ept_out,vxpt_out,vypt_out = velocity.velocity_at_eulpoints(x_out,y_out,glacier,data='all',xy_velocities='True')
  
  # Time for dH/dt calculations
  time = tpt_in
  dH = np.zeros([len(time),2])
  
  # Anticipated error between surface elevations and/or bed DEM at fluxgate_in and 
  # fluxgate_out (total error is 2*error)
  error = 2.5
  
  # Get flux through upstream gate
  xperp = (-y_in[1]+y_in[0])/np.sqrt((-y_in[1]+y_in[0])**2+(x_in[1]-x_in[0])**2)
  yperp = (x_in[1]-x_in[0])/np.sqrt((-y_in[1]+y_in[0])**2+(x_in[1]-x_in[0])**2)
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
    if len(nans) < 50:
      # If the number of locations without surface velocities is small, let's use the known 
      # surface velocities to interpolate.
      nonnans = np.where(~(np.isnan(vperp_in[i,:])))[0]
      vperp_in[i,nans] = np.interp(l_in[nans],l_in[nonnans],vperp_in[i,nonnans])
      
      # Calculate fluxes
      Q_in[i,:] = (vperp_in[i,:]*(zs_in-zb_in)*dl_in)
      Q_in_max[i,:] = (vperp_in[i,:]*(zs_in+error-zb_in)*dl_in)
      
      # We don't want negative fluxes so let's toss them out.
      Q_in[i,zs_in < zb_in] = 0
      Q_in_max[i,zs_in < zb_in] = 0
      

  # Get flux through downstream gate
  xperp = (-y_out[1]+y_out[0])/np.sqrt((-y_out[1]+y_out[0])**2+(x_out[1]-x_out[0])**2)
  yperp = (x_out[1]-x_out[0])/np.sqrt((-y_out[1]+y_out[0])**2+(x_out[1]-x_out[0])**2)
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
    if len(nans) < 50:
      # If the number of locations without surface velocities is small, let's use the known 
      # surface velocities to interpolate.
      nonnans = np.where(~(np.isnan(vperp_out[i,:])))[0]
      vperp_out[i,nans] = np.interp(l_out[nans],l_out[nonnans],vperp_out[i,nonnans])
      
      # Calculate fluxes
      Q_out[i,:] = (vperp_out[i,:]*(zs_out-zb_out)*dl_out)
      Q_out_max[i,:] = (vperp_out[i,:]*(zs_out-error-zb_out)*dl_out)
    
      # We don't want negative fluxes
      Q_in[i,zs_in < zb_in] = 0
      Q_in_max[i,zs_in < zb_in] = 0
  
  # Get surface area
  shape = shapely.geometry.Polygon(np.column_stack([fluxbox_x,fluxbox_y]))
  A = shape.area
  print A/(1000.0**2)
  
  # Calculate the change in flux between the two gates
  dQ = np.sum(Q_in,1)-np.sum(Q_out,1)
  dQ_max = np.sum(Q_in_max,1)-np.sum(Q_out_max,1)

  dH[:,0] = dQ/A # anticipated dH/dt
  dH[:,1] = (dQ_max-dQ)/A # error in dH/dt
  
  return time, dH

def dem_thinning(glacier,x,y,zs,time,fluxgate_filename):
  
  # Calculate the measured thinning rates in the fluxboxes so that we can compare the measured
  # values to the inferred values from the change in flux (fluxgate_thinning)
  
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
  dH = np.zeros([nt-1,2])
  dH_time = np.zeros([nt-1,2])
  j = 0
  for i in range(1,nt):
    # Calculate time as halfway between the two times
    dH_time[i-1,0] = (time[i]+time[j])/2 # time
    dH_time[i-1,1] = (time[i]-time[j]) # time interval
    
    # Surface elevations inside the fluxbox at times 1 and 2
    zs_inside_t1 = zs[:,:,j].flatten()[inside]
    zs_inside_t2 = zs[:,:,i].flatten()[inside]
    
    # Find locations where we have thinning rates
    nonnan = np.where(~(np.isnan(zs_inside_t2-zs_inside_t1)))[0]
    
    # Only calculate the average thinning rates if the measured thinning rates cover a 
    # sufficient area of the fluxbox and surface elevation profiles are more than a month apart.
    if len(nonnan) > 0.5*ninside and (dH_time[i-1,1] > 1/12.0):
      dH[i-1,0] = np.nanmean(zs_inside_t2-zs_inside_t1)/dH_time[i-1,1] # average thinning rate
      dH[i-1,1] = np.nanstd(zs_inside_t2-zs_inside_t1) # standard deviation of thinning rates in box
      j = j+1
    else:
      dH[i-1,0] = 'NaN'
      dH[i-1,1] = 'NaN'
      if len(np.where(~(np.isnan(zs_inside_t2)))[0]) > len(np.where(~(np.isnan(zs_inside_t1)))[0]):
        j = j+1
    
  return dH_time,dH
  
def fluxbox_geometry(glacier,fluxgate_filename):

  # Find the geometry of the fluxbox by combining the shapefiles for the fluxgates and for the glacier 
  # extent.

  # Get Shapefiles for flux gates and for glacier extent (so we can calculate surface area accurately)
  DIR = os.path.join(os.getenv("HOME"),"Data/ShapeFiles/FluxGates/"+glacier+"/")
  gate_in_sf = shapefile.Reader(DIR+fluxgate_filename+"_in")
  gate_out_sf = shapefile.Reader(DIR+fluxgate_filename+"_out")
  
  # Glacier extent shapefile
  glacier_extent = elmer_mesh.shp_to_xy(os.path.join(os.getenv("HOME"),"Data/ShapeFiles/Glaciers/3D/"+glacier+"/glacier_extent_normal"))
  if glacier == 'Helheim':
    glacier_hole1 = elmer_mesh.shp_to_xy(os.path.join(os.getenv("HOME"),"Data/ShapeFiles/Glaciers/3D/Helheim/glacier_hole1"))
    glacier_hole2 = elmer_mesh.shp_to_xy(os.path.join(os.getenv("HOME"),"Data/ShapeFiles/Glaciers/3D/Helheim/glacier_hole2"))
    glacier_extent = np.row_stack([glacier_extent,glacier_hole1,glacier_hole2])
  
  # Get end points for flux gates
  gate_in_pts = np.array(gate_in_sf.shapes()[0].points)
  sortind = np.argsort(gate_in_pts[:,1])
  gate_in_pts = gate_in_pts[sortind]
  gate_out_pts = np.array(gate_out_sf.shapes()[0].points)
  sortind = np.argsort(gate_out_pts[:,1])
  gate_out_pts = gate_out_pts[sortind]
  
  # Find glacier geometry (i.e., points in glacier extent) that lie between the flux gates, 
  # so we can calculate glacier area.
  
  fluxbox_x = [] # boundary points for area between flux gates 
  fluxbox_x.append(gate_in_pts[:,0].T.tolist())
  fluxbox_y = []
  fluxbox_y.append(gate_in_pts[:,1].T.tolist())
  for i in range(0,len(gate_in_pts)):
    ind_in = np.argmin((gate_in_pts[i,0]-glacier_extent[:,0])**2+(gate_in_pts[i,1]-glacier_extent[:,1])**2) # get closest index
    ind_out = np.argmin((gate_out_pts[i,0]-glacier_extent[:,0])**2+(gate_out_pts[i,1]-glacier_extent[:,1])**2) # get closest index
  
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

        
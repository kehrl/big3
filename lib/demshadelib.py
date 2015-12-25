'''
Module to add shadinging to DEM plots.

LMK, UW, 9/20/2015

'''

import numpy as np
import matplotlib.cm as cm

def set_shade(dem,cmin,cmax,intensity=None,cmap=cm.jet,scale=10.0,azdeg=165.0,altdeg=45.0):
  ''' Set the shading for data array based on an intensity layer or the data's value itself.
  
  Inputs:
  dem: 2-d array or masked array
  intensity: 2-d array (same size as "a") that represents the intensity layer. If none is 
  		given, the data itself is used after getting the hillshade values (see hillshade for 
  		more details).
  cmap: a colormap (e.g matplotlib.colors.LinearSegmentedColormap
              instance)
  scale,azdeg,altdeg: parameters for hillshade function

  Output:
  rgb: an rgb set of the Pegtop soft light composition of the data and intensity '''
  
  # Set up nan mask
  ind = np.where(np.isnan(dem))
  mask = np.zeros_like(dem)
  mask[ind] = 1
  
  dem_mask = np.ma.masked_array(dem,mask)
  
  if intensity is None:
	# hilshading the data
    intensity = hillshade(dem_mask,scale=10.0,azdeg=165.0,altdeg=45.0)
  else:
	# or normalize the intensity
    intensity = (intensity - intensity.min())/(intensity.max() - intensity.min())
  # get rgb of normalized data based on cmap
  rgb = cmap((dem_mask-cmin)/float(cmax-cmin))[:,:,:3]
  # form an rgb eqvivalent of intensity
  d = intensity.repeat(3).reshape(rgb.shape)
  # simulate illumination based on pegtop algorithm.
  rgb = 2*d*rgb+(rgb**2)*(1-2*d)
  
  return rgb

def hillshade(data,scale=10.0,azdeg=165.0,altdeg=45.0):
  ''' 
    Convert data to hillshade based on matplotlib.colors.LightSource class.
    
    Input: 
    data: a 2-d array of data
    scale: scaling value of the data. higher number = lower gradient 
    azdeg: where the light comes from: 0 south ; 90 east ; 180 north ;
    		270 west
    altdeg: where the light comes from: 0 horison ; 90 zenith
    
    Output: 
    intensity: a 2-d array of normalized hilshade 
  '''
  # Set up nan mask
  ind = np.where(np.isnan(data))
  mask = np.zeros_like(data)
  mask[ind] = 1
  
  data_mask = np.ma.masked_array(data,mask)
  
  # convert alt, az to radians
  az = azdeg*np.pi/180.0
  alt = altdeg*np.pi/180.0
  # gradient in x and y directions
  dx, dy = np.gradient(data_mask/float(scale))
  slope = 0.5*np.pi - np.arctan(np.hypot(dx, dy))
  aspect = np.arctan2(dx, dy)
  intensity = np.sin(alt)*np.sin(slope) + np.cos(alt)*np.cos(slope)*np.cos(-az - aspect - 0.5*np.pi)
  intensity = (intensity - intensity.min())/(intensity.max() - intensity.min())
  
  return intensity

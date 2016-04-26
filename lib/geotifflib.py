# This module manipulates geotiff files.
#
# Functions:
# read(filename,xmin=-np.Inf,xmax=np.Inf,ymin=-np.Inf,ymax=np.Inf) - read the geotiff file
#     given minimum and maximum dimensions (so you don't unnecessarily read everything)

from osgeo import gdal
from osgeo import osr
import numpy as np
import math
import matplotlib.mlab as mlab

def extent(filename):

  '''
  xmin,xmax,ymin,ymax = extent(filename)
  
  Get the extent of the geotiff file before we load it. 
  
  Inputs:
  filename: geotiff filename
  
  Outputs:
  xmin,xmax,ymin,ymax: geotiff bounding box
  '''

  geotiffile = gdal.Open(filename,gdal.GA_ReadOnly)

  # Get the size of the dataset
  try: 
    nx = geotiffile.RasterXSize
    ny = geotiffile.RasterYSize
  except ValueError:
    print "The geotiff file does not exist."

  # Get the coordinates of the image
  gt = geotiffile.GetGeoTransform()

  xmin = np.min([gt[0],gt[0]+nx*gt[1]])
  xmax = np.max([gt[0],gt[0]+nx*gt[1]])
  ymin = np.min([gt[3],gt[3]+ny*gt[5]])
  ymax = np.max([gt[3],gt[3]+ny*gt[5]])
     
  return xmin,xmax,ymin,ymax

def read(filename,xmin=-np.Inf,xmax=np.Inf,ymin=-np.Inf,ymax=np.Inf,no_data_value='none'):
  
  '''
  x,y,z = read(filename,xmin=-np.Inf,xmax=np.Inf,ymin=-np.Inf,ymax=np.Inf)
  
  Read geotiff file into python.
  
  Inputs:
  filename: geotiff filename
  xmin,xmax,ymin,ymax: bounding box for the resulting grid that we want
  
  Outputs:
  x,g: coordinates of grid
  z: grid of band values
  '''
  
  # Load geotiff file
  geotiffile = gdal.Open(filename,gdal.GA_ReadOnly)

  # Get the size of the dataset
  try: 
    nx = geotiffile.RasterXSize
    ny = geotiffile.RasterYSize
  except ValueError:
    print "The geotiff file does not exist."

  # Get no data value
  nodatavalue = geotiffile.GetRasterBand(1).GetNoDataValue()
  if no_data_value != 'none':
    nodatavalue = no_data_value
    
  # Get the coordinates of the image
  gt = geotiffile.GetGeoTransform()
  x = np.zeros(nx)
  y = np.zeros(ny)
  dx = math.fabs(gt[1])
  dy = math.fabs(gt[5])
  for i in range(nx):
    x[i] = gt[0]+i*gt[1]
  for i in range(ny):
    y[i] = gt[3]+i*gt[5]

  # Get index extent for the domain we want
  if np.isfinite(xmin):
    j0 = int( np.argmin(abs(xmin-x)/dx) ) 
    j1 = int( np.argmin(abs(xmax-x)/dx) ) 
    if j1 < j0:
      temp = int(j0)
      j0 = int(j1)
      j1 = int(temp)
    i0 = int( np.argmin(abs(ymin-y)/dy) ) 
    i1 = int( np.argmin(abs(ymax-y)/dy) ) 
    if i1 < i0:
      temp = int(i0)
      i0 = int(i1)
      i1 = int(temp)
  else:
    i0 = 0
    j0 = 0
    i1 = ny-1
    j1 = nx-1

  # Sort x,y and crop to desired domain
  x = x[j0:j1+1]
  if x[0] > x[-1]:
    x = x[::-1]
  y = y[i0:i1+1]
  if y[0] > y[-1]:
    y = y[::-1]

  # Load the data from the file
  nx = j1-j0+1
  ny = i1-i0+1
  nbands = geotiffile.RasterCount
  z = np.zeros((ny,nx,nbands))
  for i in range(1,nbands+1):
    z[::-1,:,i-1] = geotiffile.GetRasterBand(1).ReadAsArray(xoff=j0,yoff=i0,\
    				win_xsize=nx,win_ysize=ny)

  # If there is a no data value, we want to save the no data indices to nan's.
  try:
    if not(np.isnan(nodatavalue)):
      z[z==nodatavalue] = float('NaN')
  except:
    pass

  return x,y,z[:,:,0]


def readrgb(filename,xmin=-np.Inf,xmax=np.Inf,ymin=-np.Inf,ymax=np.Inf):
  
  geotiffile = gdal.Open(filename,gdal.GA_ReadOnly)
  
  # Get the size of the dataset
  try: 
    nx = geotiffile.RasterXSize
    ny = geotiffile.RasterYSize
  except ValueError:
    print "The geotiff file does not exist."

  # Load the data from the file
  nbands = geotiffile.RasterCount
  z = np.zeros((ny,nx,nbands))
  for i in range(1,nbands+1):
    z[::-1,:,i-1] = geotiffile.GetRasterBand(i).ReadAsArray()

  # Get no data value
  nodatavalue = geotiffile.GetRasterBand(1).GetNoDataValue()
  
  # Get the coordinates of the image
  gt = geotiffile.GetGeoTransform()
  x = np.zeros(nx)
  y = np.zeros(ny)
  for i in range(nx):
    x[i] = gt[0]+i*gt[1]
  for i in range(ny):
    y[i] = gt[3]+i*gt[5]
  y = y[::-1]

  dx = math.fabs(gt[1])
  dy = math.fabs(gt[5])

  j0 = int( max( 0, (xmin-x[0])/dx ) )
  j1 = int( min( nx, (xmax-x[0])/dx )+1 )
  i0 = int( max( 0, (ymin-y[0])/dy ) )
  i1 = int( min( ny, (ymax-y[0])/dy )+1 )
  
  znormal = np.zeros_like(z)
  for i in range(0,len(z[0,0,:])):
    znormal[:,:,i] = z[:,:,i]/np.max(z[:,:,i])

  # If there is a no data value, we want to save the no data indices to 'nan's.
  try:
    if not(np.isnan(nodatavalue)):
      znormal[znormal==nodatavalue] = float('NaN')
  except:
    pass

  return (x[j0:j1],y[i0:i1],znormal[i0:i1,j0:j1,:])


def write_from_grid(x,y,zgrid,nodatavalue,outfile,proj="+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"):

  nrows,ncols = np.shape(zgrid)
  xres = x[1]-x[0]
  yres = y[1]-y[0]

  geotransform=(x[0],xres,0,y[-1],0, -yres) 
  options= [ 'TILED=YES', 'COMPRESS=LZW']
  output_raster = gdal.GetDriverByName('GTiff').Create(outfile,ncols, nrows, 1 ,gdal.GDT_Float32,options)
  output_raster.SetGeoTransform(geotransform)  
  srs = osr.SpatialReference()                 
  
  srs.ImportFromProj4(proj)               
  output_raster.SetProjection( srs.ExportToWkt() )   
  output_raster.GetRasterBand(1).WriteArray(np.array(zgrid))
  output_raster.GetRasterBand(1).SetNoDataValue(nodatavalue)  
  #output_raster=None

  return 1

def write_from_xyz(x,y,z,outfile,nodatavalue):

  xmin,xmax,ymin,ymax = [min(x),max(x),min(y),max(y)]

  #Grid size
  nx = (xmax - xmin)/100
  ny = (ymax - ymin)/100

  # Generate a regular grid to interpolate the data.
  xi = np.linspace(xmin, xmax, nx)
  yi = np.linspace(ymin, ymax, ny)
  xi, yi = np.meshgrid(xi, yi) 

  # Interpolate the values of z for all points in the rectangular grid
  zgrid = mlab.griddata(x,y,z,xi,yi,interp='linear') #interpolation is 'nn' by default (natural neighbour based on delaunay triangulation) but 'linear' is faster (see http://matplotlib.1069221.n5.nabble.com/speeding-up-griddata-td20906.html)

  #---------------  Write to GeoTIFF ------------------------
  nrows,ncols = np.shape(zgrid)
  xres = (xmax-xmin)/float(ncols)
  yres = (ymax-ymin)/float(nrows)

  geotransform=(xmin,xres,0,ymin,0, yres) 
  output_raster = gdal.GetDriverByName('GTiff').Create(outfile,ncols, nrows, 1 ,gdal.GDT_Int16,['TFW=YES', 'COMPRESS=PACKBITS'])
  output_raster.SetGeoTransform(geotransform)  
  srs = osr.SpatialReference()                 
  srs.ImportFromEPSG(3413)                     
  output_raster.SetProjection( srs.ExportToWkt() )   
  output_raster.GetRasterBand(1).WriteArray(np.array(zgrid))
  output_raster.GetRasterBand(1).SetNoDataValue(nodatavalue)  
  print np.min(zgrid), np.max(zgrid) 
  output_raster=None
  
  return 1
# This script creates mosaics from Dana's TDX DEMs.
#
#
# LMK, UW, 9/29/2015

import os

glacier = "Helheim"
TOPDIR = "/Volumes/insar4/ian/TSX/Helheim/DEM/fromDana/"
OUTDIR =  os.path.join(os.getenv("HOME"),"bigtmp/shean_tdx")


# Get list of directories
dirs = os.listdir(TOPDIR)

# Get files, convert to EPSG 3413, and keep note of the track names so that we can mosaic them later
tracklist = []
for dir in dirs:
  # Check if it is a track directory. If so, convert file to EPSG 3413
  if os.path.isdir(TOPDIR+dir):
    files = os.listdir(TOPDIR+dir)
    for file in files:
      if file.endswith('_DEM.ras'):
        newfile = OUTDIR+file[25:33]+'_'+file[34:38]+'_'+dir[0:7]+'_'+file[41:-4]+'.tif'
        os.system('gdalwarp -t_srs +proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs -et 0 -co TILED=YES -co COMPRESS=LZW -co BIGTIFF=IF_SAFER -r cubic '+TOPDIR+dir+'/'+file+' '+newfile)
        os.system('gdaladdo_lzw.py '+newfile)
    if (dir[0:7] not in tracklist):
      tracklist.append(dir[0:7])

files = os.listdir(OUTDIR)
for track in tracklist:
  filestomosaic=''
  for file in files:
    if track in file:
      filestomosaic = filestomosaic+' '+OUTDIR+file
      fileout = file[0:21]+file[25:]
  os.system('gdalwarp -co TILED=YES -co COMPRESS=LZW -co BIGTIFF=IF_SAFER -r cubic'+filestomosaic+' '+OUTDIR+fileout)

# Use Shean's demtools mask_raster.sh to burn mask out bad sections before aligning the DEM.
# This code takes the original Landsat TIF files and creates 
# RGB files. The output is cropped to keep the file size small. 
# The code requires the glacier name (Kanger, Helheim) as input.
#
# LMK, UW, 6/20/2015

import os
import sys
import shutil
import numpy as np
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"Util/Modules"))
sys.path.append("/usr/local/bin/")
from subprocess import call
import math
import gdal_merge
from geodatlib import *
import glob

# Get arguments
args = sys.argv

DIR=os.path.join(os.getenv("DATA_HOME"),"Imagery/Landsat/")
glacier = args[1][:] # Options: Kanger, Helheim

############################################
# Set extent for cropping based on glacier #
############################################
if glacier == 'Helheim':
  extent='260000 -2600000 325000 -2540000'
elif glacier == 'Kanger':
  extent='440000 -2318000 526000 -2252000'
else:
  sys.exit("Unknown glacier")

###############
# Landsat 4-5 #
###############
files = os.listdir(DIR+glacier+"/OriginalData/L4-5")

#############
# Landsat 7 #
#############
files = os.listdir(DIR+glacier+"/OriginalData/L7")

for file in files:
  print file
  if file.startswith("L"):
    DIRFILE=DIR+glacier+"/OriginalData/L7/"+file+"/"
    os.chdir(DIRFILE)
    metafile = open(DIRFILE+file+"_MTL.txt","r")
    lines = metafile.readlines()
    for line in lines:
      if line.startswith("    DATE_ACQUIRED"):
        year = str(line[20:24])
        month = str(line[25:27])
        day = str(line[28:30])
      if line.startswith("    SCENE_CENTER_TIME"):
        if str(line[23:25])==" \"":
          hour = str(line[25:27])
          min = str(line[28:30])
          sec = str(line[31:33])
        else:
          hour = str(line[24:26])
          min = str(line[27:29])
          sec = str(line[30:32])
    filename=DIR+glacier+"/TIF/"+year+month+day+hour+min+sec+"_"+file+".tif"
    if not(os.path.isfile(filename)):
      print filename
      sys.argv[1:] = ['-separate',file+"_B3.tif",file+"_B2.tif",file+"_B1.tif",'-o','temp1.tif']
      gdal_merge.main()
      os.system('gdal_translate -co PHOTOMETRIC=RGB temp1.tif temp2.tif')
      os.system('otbcli_BundleToPerfectSensor -inp '+file+'_B8.tif  -inxs temp2.tif -out temp3.tif')
      os.system('gdalwarp temp3.tif temp4.tif -t_srs EPSG:3413'+' -te '+extent)
      os.system('gdal_translate -scale 0 60035 0 255 -co TILED=YES -co BIGTIFF=YES -co COMPRESS=LZW temp4.tif '+filename)
      for fname in os.listdir('.'):
        if fname.startswith("temp"):
          os.remove(os.path.join(fname))
  
#############
# Landsat 8 #
#############  
files = os.listdir(DIR+glacier+"/OriginalData/L8")

for file in files:
  if file.startswith("L"):
    DIRFILE=DIR+glacier+"/OriginalData/L8/"+file+"/"
    os.chdir(DIRFILE)
    metafile = open(DIRFILE+file+"_MTL.txt","r")
    lines = metafile.readlines()
    for line in lines:
      if line.startswith("    DATE_ACQUIRED"):
        year = str(line[20:24])
        month = str(line[25:27])
        day = str(line[28:30])
      if line.startswith("    SCENE_CENTER_TIME"):
        if str(line[23:25])==" \"":
          hour = str(line[25:27])
          min = str(line[28:30])
          sec = str(line[31:33])
        else:
          hour = str(line[24:26])
          min = str(line[27:29])
          sec = str(line[30:32])
    filename=DIR+glacier+"/TIF/"+year+month+day+hour+min+sec+"_"+file+".tif"
    if not(os.path.isfile(filename)):
      print filename
      sys.argv[1:] = ['-separate',file+"_B4.tif",file+"_B3.tif",file+"_B2.tif",'-o','temp1.tif']
      gdal_merge.main()
      os.system('gdal_translate -co PHOTOMETRIC=RGB temp1.tif temp2.tif')
      os.system('otbcli_BundleToPerfectSensor -inp '+file+'_B8.tif  -inxs temp2.tif -out temp3.tif')
      os.system('gdalwarp temp3.tif temp4.tif -t_srs EPSG:3413'+' -te '+extent)
      os.system('gdal_translate -scale 0 60035 0 255 -co TILED=YES -co BIGTIFF=YES -co COMPRESS=LZW temp4.tif '+filename)
      for fname in os.listdir('.'):
        if fname.startswith("temp"):
          os.remove(os.path.join(fname))

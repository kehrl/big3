# This code takes the original Landsat TIF files and creates 
# RGB files. The output is cropped to keep the file size small. 
# The code requires the glacier name (Kanger, Helheim) as input.
#
# LMK, UW, 6/20/2015

import os
import sys
import shutil
import numpy as np
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
sys.path.append("/usr/local/bin/")
from subprocess import call
import math
import gdal_merge
from geodat import *
import glob

# Get arguments
args = sys.argv

DIR=os.path.join(os.getenv("HOME"),"Data/Imagery/Landsat/")
glacier = args[1][:] # Options: Kanger, Helheim

############################################
# Set extent for cropping based on glacier #
############################################
if glacier == 'Helheim':
  extent='300000 -2597000 320000 -2560000'
elif glacier == 'Kanger':
  extent='477000 -2308000 516000 -2269000'
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
      hour = str(line[24:26])
      min = str(line[27:29])
      sec = str(line[30:32])
  filename=DIR+glacier+"/TIF/"+year+month+day+hour+min+sec+"_"+file+".tif"
  if not(os.path.isfile(filename)):
    print filename
    sys.argv[1:] = ['-separate',file+"_B3.tif",file+"_B2.tif",file+"_B1.tif",'-o','temp1.tif']
    gdal_merge.main()
    os.system('gdal_translate -co PHOTOMETRIC=RGB temp1.tif temp2.tif')
    os.system('otbcli_BundleToPerfectSensor -inp '+file+'_B8.tif  -inxs temp2.tif -out temp3.tif uint16')
    os.system('gdalwarp temp3.tif '+filename+' -t_srs EPSG:3413'+' -te '+extent)
    os.remove('temp1.tif')
    os.remove('temp2.tif')
    os.remove('temp3.tif')
  
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
        hour = str(line[24:26])
        min = str(line[27:29])
        sec = str(line[30:32])
    filename=DIR+glacier+"/TIF/"+year+month+day+hour+min+sec+"_"+file+".tif"
    if not(os.path.isfile(filename)):
      print filename
      sys.argv[1:] = ['-separate',file+"_B4.tif",file+"_B3.tif",file+"_B2.tif",'-o','temp1.tif']
      gdal_merge.main()
      os.system('gdal_translate -co PHOTOMETRIC=RGB temp1.tif temp2.tif')
      os.system('otbcli_BundleToPerfectSensor -inp '+file+'_B8.tif  -inxs temp2.tif -out temp3.tif uint16')
      os.system('gdalwarp temp3.tif temp4.tif -t_srs EPSG:3413'+' -te '+extent)
      os.system('gdal_translate -ot Byte -scale 0 6535 0 255 -a_nodata "0 0 0" temp4.tif '+filename)
      os.remove('temp1.tif')
      os.remove('temp2.tif')
      os.remove('temp3.tif')
      os.remove('temp4.tif')

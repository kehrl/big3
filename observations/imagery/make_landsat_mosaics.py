# This code takes the original Landsat TIF files and creates 
# RGB files. The output is cropped to keep the file size small. 
# The code requires the glacier name (Kanger, Helheim) as input.
#
# LMK, UW, 6/20/2015

import os
import sys
import shutil
import numpy as np
from subprocess import call
import math
import gdal_merge
from geodatlib import *
import glob

# Get arguments
args = sys.argv
glacier = args[1][:] # Options: Kanger, Helheim

DIR2 = os.path.join(os.getenv("DATA2_HOME"),"Imagery/Landsat/")
DIR = os.path.join(os.getenv("DATA_HOME"),"Imagery/Landsat/")
DIRL45 = DIR2+glacier+"/OriginalData/L4-5/"
DIRL7 = DIR2+glacier+"/OriginalData/L7/"
DIRL8 = DIR2+glacier+"/OriginalData/L8/"
OUTDIR = DIR+glacier+"/TIF/"


############################################
# Set extent for cropping based on glacier #
############################################
if glacier == 'Helheim':
  # Big mosaic
  #extent='260000 -2600000 325000 -2540000'
  extent='270000 -2597000 320000 -2550000'
elif glacier == 'Kanger':
  # Big mosaic
  #extent='440000 -2318000 526000 -2252000'
  extent='450000 -2308000 516000 -2262000'
else:
  sys.exit("Unknown glacier")

###############
# Landsat 4-5 #
###############
files = os.listdir(DIRL45)

#############
# Landsat 7 #
#############
files = os.listdir(DIRL7)
outputfiles = os.listdir(OUTDIR)

for file in files:
  if file.startswith("L"):
    print file[0:21]
    DIRFILE=DIRL7+file[0:21]+"/"
    createfile = 1
    for outputfile in outputfiles:
      if file[0:-7] == outputfile[-25:-4]:
        createfile = 0 
    if createfile == 1:
      inputfile = file[0:21]
      os.chdir(DIRL7)
      os.mkdir(inputfile)
      os.system('tar xf '+file+' -C '+inputfile)
      os.chdir(DIRFILE)
      metafile = open(DIRFILE+inputfile+"_MTL.txt","r")
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
      filename=OUTDIR+year+month+day+hour+min+sec+"_"+inputfile+".tif"

      print filename
      sys.argv[1:] = ['-separate',inputfile+"_B3.tif",inputfile+"_B2.tif",inputfile+"_B1.tif",'-o','temp1.tif']
      gdal_merge.main()
      os.system('gdal_translate -co PHOTOMETRIC=RGB temp1.tif temp2.tif')
      os.system('otbcli_BundleToPerfectSensor -inp '+inputfile+'_B8.tif  -inxs temp2.tif -out temp3.tif')
      os.system('gdalwarp temp3.tif temp4.tif -t_srs EPSG:3413'+' -te '+extent)
      os.system('gdal_translate -scale 0 60035 0 255 -co TILED=YES -co BIGTIFF=YES -co COMPRESS=LZW temp4.tif '+filename)
      for fname in os.listdir('.'):
        if fname.startswith("temp"):
          os.remove(os.path.join(fname))
      os.chdir(DIRL7)
      shutil.rmtree(inputfile)

#############
# Landsat 8 #
#############
files = os.listdir(DIRL8)
outputfiles = os.listdir(OUTDIR)

for file in files:
  if file.startswith("L"):
    print file[0:21]
    DIRFILE=DIRL8+file[0:21]+"/"
    createfile = 1
    for outputfile in outputfiles:
      if file[0:-7] == outputfile[-25:-4]:
        createfile = 0 
    if createfile == 1:
      inputfile = file[0:21]
      os.chdir(DIRL8)
      os.mkdir(inputfile)
      os.system('tar xf '+file+' -C '+inputfile)
      os.chdir(DIRFILE)
      metafile = open(DIRFILE+inputfile+"_MTL.txt","r")
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
      filename=OUTDIR+year+month+day+hour+min+sec+"_"+inputfile+".tif"

      print filename
      sys.argv[1:] = ['-separate',inputfile+"_B3.tif",inputfile+"_B2.tif",inputfile+"_B1.tif",'-o','temp1.tif']
      gdal_merge.main()
      os.system('gdal_translate -co PHOTOMETRIC=RGB temp1.tif temp2.tif')
      os.system('otbcli_BundleToPerfectSensor -inp '+inputfile+'_B8.tif  -inxs temp2.tif -out temp3.tif')
      os.system('gdalwarp temp3.tif temp4.tif -t_srs EPSG:3413'+' -te '+extent)
      os.system('gdal_translate -scale 0 60035 0 255 -co TILED=YES -co BIGTIFF=YES -co COMPRESS=LZW temp4.tif '+filename)
      for fname in os.listdir('.'):
        if fname.startswith("temp"):
          os.remove(os.path.join(fname))
      os.chdir(DIRL8)
      shutil.rmtree(inputfile)

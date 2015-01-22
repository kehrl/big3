import os
import sys
import shutil
import numpy as np
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
sys.path.append("/usr/local/bin/")
import velocity_flowline
from subprocess import call
import math
import gdal_merge

from geodat import *

import glob
import elmer_read

DIR=os.path.join(os.getenv("HOME"),"Data/Imagery/Landsat/")
glacier='Helheim'

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
  print filename
  sys.argv = ['-separate',file+"_B3.tif",file+"_B2.tif",file+"_B1.tif",'-o',filename]
  gdal_merge.main()
  
#############
# Landsat 8 #
#############  
files = os.listdir(DIR+glacier+"/OriginalData/L8")

for file in files:
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
  print filename
  sys.argv = ['-separate',file+"_B4.tif",file+"_B3.tif",file+"_B2.tif",'-o',filename]
  gdal_merge.main()

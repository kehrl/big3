# This is just a small script to move Shean's DEMs to my local machine. I'm leaving it in
# my code so that I remember where the DEMS came from...
#
# LMK, UW, 7/5/2015
#
import os, sys, shutil

glacier = 'Helheim'

# Directory where data is located
if glacier == 'Kanger':
  TOPDIR = "/Volumes/insar3/dshean/helheim_kanger_se/dem/kanger/front/"
else: 
  TOPDIR = "/Volumes/insar3/dshean/helheim_kanger_se/dem/helheim_front/"

# Directory where we want to put the DEMs
OUTDIR = os.path.join(os.getenv("HOME"),"Data/Elevation/Worldview/"+glacier+"/")

# Get contents of directory
try:
  DIRs = os.listdir(TOPDIR)
except:
  sys.exit("Need to load insar3 to volumes.")

# Copy 32 m DEMs to my local machine
for DIR in DIRs:
  if glacier == 'Kanger':
    if DIR.endswith("_align"):
      dir = os.listdir(TOPDIR+DIR)
      files = os.listdir(TOPDIR+DIR)
      for file in files:
        if file.endswith('trans_reference-DEM_32m.tif'):
          print file
          shutil.copy(TOPDIR+DIR+"/"+file,OUTDIR+file)
  if glacier == 'Helheim':
    if DIR.endswith("DEM_32m_trans.tif"):
       shutil.copy(TOPDIR+DIR,OUTDIR+DIR)
  
# This is just a small script to move Shean's DEMs to my local machine. I'm leaving it in
# my code so that I remember where the DEMS came from...
#
# LMK, UW, 7/5/2015
#
import os, sys, shutil

# Directory where data is located
TOPDIR = "/Volumes/insar3/dshean/helheim_kanger_se/dem/"

# Directory where we want to put the DEMs
OUTDIR = os.path.join(os.getenv("HOME"),"Data/Elevation/Worldview/HelheimKanger/")

# Get contents of directory
try:
  DIRs = os.listdir(TOPDIR)
except:
  sys.exit("Need to load insar3 to volumes.")

# Copy 32 m DEMs to my local machine
for DIR in DIRs:
  if DIR.startswith("WV"):
    dir = os.listdir(TOPDIR+DIR)
    files = os.listdir(TOPDIR+DIR+"/"+dir[0])
    for file in files:
      if file.endswith('DEM_32m.tif'):
        shutil.copy(TOPDIR+DIR+"/"+dir[0]+"/"+file,OUTDIR+file)
  
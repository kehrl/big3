# This is just a small script to move Shean's DEMs to my local machine. I'm leaving it in
# my code so that I remember where the DEMS came from...
#
# LMK, UW, 7/5/2015
#
import os, sys, shutil

# Get arguments
args = sys.argv
glacier = args[1][:] # Options: Kanger, Helheim

# Directory where data is located
if glacier == 'Kanger':
  #TOPDIR = "/Volumes/insar3/dshean/helheim_kanger_se/dem/kanger/front/"
  TOPDIR = "/Volumes/insar3/dshean/big3/32m_DEM_movies/kanger/"
elif glacier == 'Helheim':
  #TOPDIR = "/Volumes/insar3/dshean/helheim_kanger_se/dem/helheim_front/"
  TOPDIR = "/Volumes/insar3/dshean/big3/32m_DEM_movies/helheim/"

# Directory where we want to put the DEMs
OUTDIR = os.path.join(os.getenv("DATA_HOME"),"Elevation/Worldview/"+glacier+"/")

# Get contents of directory
try:
  DIRs = os.listdir(TOPDIR)
except:
  sys.exit("Need to load insar3 to volumes.")

os.chdir(OUTDIR)
# Copy 32 m DEMs to my local machine
for DIR in DIRs:
  if DIR.endswith("DEM_32m_trans.tif"):
    shutil.copy(TOPDIR+DIR,OUTDIR+DIR)
    os.system('dem_geoid --geoid EGM2008 '+DIR)

  
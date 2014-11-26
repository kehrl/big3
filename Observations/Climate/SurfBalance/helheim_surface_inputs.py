# Beginning of looking at surface mass balance for the models...

import os
import sys
import numpy as np
sys.path.append(os.path.join(os.getenv("HOME"),"Dropbox/Code/Modules"))
from subprocess import call

##########
# Inputs #
##########

# Jason Box 2004 (polar mm5)
file_lat_in=os.path.join(os.getenv("HOME"),"Data/Surface_balance/Greenland/Box_smb_2004/latitude_101x55_24km.asc")
file_long_in=os.path.join(os.getenv("HOME"),"Data/Surface_balance/Greenland/Box_smb_2004/longitude_101x55_24km.asc")
file_smb_in=os.path.join(os.getenv("HOME"),"Data/Surface_balance/Greenland/Box_smb_2004/2000_smb_101x55_24km.asc")

#############
# Load data #
#############

fid = open(file_lat_in,"r")
lines=fid.readlines()
lat=[]
for line in lines:
  lat.append(float(line))
fid.close()

fid = open(file_long_in,"r")
lines=fid.readlines()
long=[]
for line in lines:
  long.append(float(line))
fid.close()

fid = open(file_smb_in,"r")
lines=fid.readlines()
smb=[]
for line in lines:
  p=line.split()
  for i in range(0,len(p)):
    smb.append(p[i])

fid = open('temp_latlong.dat',"w")    
for i in range(0,len(lat)):
  fid.write('{} {}\n'.format(long[i],lat[i]))
fid.close()

os.system("gdaltransform -s_srs EPSG:4326 -t_srs EPSG:3413 < temp_latlong.dat > temp_xy.dat")
os.system("rm temp_latlong.dat")

fid = open('temp_xy.dat',"r")
lines=fid.readlines()
x=[]
y=[]
for line in lines:
  p=line.split()
  x.append(p[0])
  y.append(p[1])
fid.close()

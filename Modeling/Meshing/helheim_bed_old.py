# This script looks at Helheim's bed, first by making some plots. Eventually we'll
# need to combine several data products to get a full bed profile for Helheim.

# LMK, UW, 10/06/2014

import os
import shutil
import matplotlib.pyplot as plt
import numpy as np
from subprocess import call
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

# Beds
morlighem_file=os.path.join(os.getenv("DATA_HOME"),"Bed/Morlighem_2014/morlighem_bed.grd")
smith_file=os.path.join(os.getenv("DATA_HOME"),"Bed/Ben_smith_feb2014/smith_helheim_bed.grd")
cresis2008_file=os.path.join(os.getenv("DATA_HOME"),"Bed/Cresis/Helheim_2008_2012_Composite/flightlines/Helheim_2008_2012_Composite_Flightlines.txt")
cresis2001_file=os.path.join(os.getenv("DATA_HOME"),"Bed/Cresis/helheim_flightline_05212001_good_nsidc.dat")

try:
  zbed_cresis2008
except:
  # Get the fields we want from the CreSIS flightlines
  fid = open(cresis2008_file,"r")
  x_cresis2008=[];y_cresis2008=[];zbed_cresis2008=[]
  lines=fid.readlines()
  for line in lines[1:-1]:
    fields = line.split()
    x_cresis2008.append(float(fields[0]))
    y_cresis2008.append(float(fields[1]))
    zbed_cresis2008.append(float(fields[3]))
  fid.close()
  zbed_cresis2008=np.array(zbed_cresis2008)
  x_cresis2008 = np.array(x_cresis2008)
  y_cresis2008 = np.array(y_cresis2008)
  del fields, lines, line

  # Write the coordinates to a file so we can use grdtrack to find depths in Morlighem bed map
  fid = open('temp_bed.txt',"w")
  for i in range(0,len(zbed_cresis2008)):
    fid.write('{0} {1}\n'.format(x_cresis2008[i],y_cresis2008[i]))
  fid.close()
  call(["grdtrack","temp_bed.txt","-G"+morlighem_file,">","temp_morlighem.txt"])

  fid = open('temp_morlighem.txt',"r")
  x_mor2008=[];y_mor2008=[];zbed_mor2008=[]
  lines=fid.readlines()
  for line in lines:
    fields = line.split()
    x_mor2008.append(float(fields[0]))
    y_mor2008.append(float(fields[1]))
    zbed_mor2008.append(float(fields[2]))
  fid.close()
  del fields,lines,line

  fid = open(cresis2001_file,"r")
  x_cresis2001=[];y_cresis2001=[];zbed_cresis2001=[];zsur_cresis2001=[]
  lines=fid.readlines()
  for line in lines[1:-1]:
    fields = line.split()
    x_cresis2001.append(float(fields[0]))
    y_cresis2001.append(float(fields[1]))
    zbed_cresis2001.append(float(fields[7]))
    zsur_cresis2001.append(float(fields[6]))
  fid.close()
  zbed_cresis2001=np.array(zbed_cresis2001[3180:3297])
  zsur_cresis2001=np.array(zsur_cresis2001[3180:3297])
  x_cresis2001 = np.array(x_cresis2001)[3180:3297]
  y_cresis2001 = np.array(y_cresis2001)[3180:3297]
  del fields, lines, line
  
  # Write the coordinates to a file so we can use grdtrack to find depths in Morlighem bed map
  fid = open('temp_bed.txt',"w")
  for i in range(0,len(zbed_cresis2001)):
    fid.write('{0} {1}\n'.format(x_cresis2001[i],y_cresis2001[i]))
  fid.close()
  call(["grdtrack","temp_bed.txt","-G"+morlighem_file,">","temp_morlighem.txt"])
  call(["grdtrack","temp_bed.txt","-G"+smith_file,">","temp_smith.txt"])
  
  fid = open('temp_morlighem.txt',"r")
  x_mor2001=[];y_mor2001=[];zbed_mor2001=[]
  lines=fid.readlines()
  for line in lines:
    fields = line.split()
    x_mor2001.append(float(fields[0]))
    y_mor2001.append(float(fields[1]))
    zbed_mor2001.append(float(fields[2]))
  fid.close()
  del fields,lines,line
  
  fid = open('temp_smith.txt',"r")
  x_ben2001=[];y_ben2001=[];zbed_ben2001=[]
  lines=fid.readlines()
  for line in lines:
    fields = line.split()
    x_ben2001.append(float(fields[0]))
    y_ben2001.append(float(fields[1]))
    zbed_ben2001.append(float(fields[2]))
  fid.close()
  del fields,lines,line
  call(["rm","temp_bed.txt","temp_morlighem.txt","temp_smith.txt"])

  # Load Morlighem bed surface
  morfile=Dataset(morlighem_file)
  x=morfile.variables['x_range']
  y=morfile.variables['y_range']
  spacing=morfile.variables['spacing']
  nx = (x[-1]-x[0])/spacing[0]
  ny = (y[-1]-y[0])/spacing[1]
  x_mor=np.linspace(x[0],x[-1],nx)
  y_mor=np.linspace(y[0],y[-1],ny)
  zz = morfile.variables['z']
  zbed_mor = zz[:].reshape(ny, nx)
  del zz,nx,ny,morfile,x,y


# Compare CreSIS flightlines to Morlighem bed
plt.clf()
plt.title('Morlighem vs. Cresis 2008-2012')
plt.imshow(zbed_mor,extent=[x_mor.min(),x_mor.max(),y_mor.min(),y_mor.max()],vmin=-1400,vmax=0)
plt.scatter(x_cresis2008,y_cresis2008,c=zbed_cresis2008,vmin=-1400,vmax=0)
plt.scatter(x_cresis2001,y_cresis2001,c=(zbed_cresis2001),vmin=-1400,vmax=0,marker='d')
plt.xlim([293000,315300])
plt.ylim([-2585200,-2570700])
cbar=plt.colorbar()
cbar.set_label('Bed elevation (m)', rotation=270)
plt.savefig('mor_cresis2008.pdf',format='PDF')

plt.clf()
plt.title('Morlighem vs. Cresis 2008-2012')
plt.imshow(zbed_mor,extent=[x_mor.min(),x_mor.max(),y_mor.min(),y_mor.max()],vmin=-1200,vmax=0)
cbar=plt.colorbar()
cbar.set_label('Bed elevation (m)', rotation=270)
plt.scatter(x_cresis2008,y_cresis2008,c=(zbed_cresis2008-zbed_mor2008),cmap='RdBu_r',vmin=-200,vmax=200)
cbar=plt.colorbar()
cbar.set_label('CreSIS-Morlighem (m)', rotation=270)
plt.xlim([293000,315300])
plt.ylim([-2585200,-2570700])
plt.savefig('mor_cresis2008_dif.pdf',format='PDF')

plt.clf()
dist2001=np.empty([len(zbed_cresis2001),1])
dist2001[1]=0
for i in range(1,len(zbed_cresis2001)):
  dist2001[i]=dist2001[i-1]+np.sqrt((x_cresis2001[i]-x_cresis2001[i-1])**2+(y_cresis2001[i]-y_cresis2001[i-1])**2)

ax=plt.subplot(2,1,1)
plt.imshow(zbed_mor,extent=[x_mor.min(),x_mor.max(),y_mor.min(),y_mor.max()],cmap='RdBu_r',vmin=-1200,vmax=0)
plt.plot(x_cresis2001[0:80],zbed_mor2001[0:80],'k',linewidth=3,label='Morlighem')
plt.plot(x_cresis2001[0:70],zbed_ben2001[0:70],color='0.75',linewidth=3,label='Ben')
plt.plot(x_cresis2001,y_cresis2001,'r.',label='Cresis 2001')
plt.plot(x_cresis2008[273:352],y_cresis2008[273:352],'b.',label='Cresis 2008a')
plt.plot(x_cresis2008[1080:1149],y_cresis2008[1080:1149],'g.',label='Cresis 2008b')
plt.plot(x_cresis2008[21667:21721],y_cresis2008[21667:21721],'c.',label='Cresis 2011')
plt.plot(x_cresis2008[33105:33174],y_cresis2008[33105:33174],'m.',label='Cresis 2012')

plt.ylim([-2581000,-2573000])
plt.xlim([299000,315000])
plt.xticks([])
plt.yticks([])
box = ax.get_position()
ax.set_position([0.12,box.y0,box.width,box.height])
leg = ax.legend(loc = 'center left', bbox_to_anchor = (0.9, 0.5))

plt.subplot(2,1,2)
plt.plot(x_cresis2001,zbed_cresis2001+120,'r.',label='Cresis 2001')
plt.plot(x_cresis2008[273:352],zbed_cresis2008[273:352],'b.',label='Cresis 2008a')
plt.plot(x_cresis2008[1080:1149],zbed_cresis2008[1080:1149],'g.',label='Cresis 2008b')
plt.plot(x_cresis2008[21667:21721],zbed_cresis2008[21667:21721],'c.',label='Cresis 2011')
plt.plot(x_cresis2008[33105:33174],zbed_cresis2008[33105:33174],'m.',label='Cresis 2012')
plt.plot(x_cresis2001[0:80],zbed_mor2001[0:80],'k',linewidth=3,label='Morlighem')
plt.plot(x_cresis2001[0:70],zbed_ben2001[0:70],color='0.75',linewidth=3,label='Ben')
plt.xlabel('Easting (m)')
plt.ylabel('Bed elevation (m)')
plt.xlim([299000,315000])
plt.savefig('cresis2001.pdf',format='PDF')
plt.clf()

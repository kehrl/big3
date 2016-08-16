import sys,os
import elmerreadlib,datelib,glaclib,geotifflib
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import argparse

# Get inputs
parser = argparse.ArgumentParser()
parser.add_argument("-glacier",dest="glacier",required = True, 
  help = "Name of glacier (Kanger or Helheim)")
parser.add_argument("-mesh", dest="meshname", required = True,
  help = "Name of output mesh.")
parser.add_argument("-method", dest="method", required = True,
        help = "Method of inversion (robin or adjoint).")
parser.add_argument("-regpar", dest="regpar", required = True,
        help = "Regularization parameter.")
parser.add_argument("-label", dest="label", required = False,
        help = "subpanel label (a,b,c,etc.).",default = 'none')

# Get arguments
args, _ = parser.parse_known_args(sys.argv)

glacier = args.glacier
meshname = args.meshname
method = args.method
regpar = args.regpar
label = args.label

DIRM = os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/"+meshname+"/")
DIRR = DIRM+"mesh2d/inversion_"+method+"/"

# Find directory that we want to load inversion results from...
dirs = os.listdir(DIRR)
for dir in dirs:
  if dir.startswith("lambda_"+regpar) and os.path.isdir(DIRR+dir):
    try:
      DIRREGPAR
      DIRREGPAR = dir+"/"
      print "More than one result with that regularization parameter, defaulting to the most recent version, "+dir+", rather than "+DIREGPAR
    except:   
      DIRREGPAR = dir+"/"
try:
  DIRREGPAR
except NameError:
  print "Couldn't find a result with regularization parameter "+regpar+" in directory "+meshname
  sys.exit()

# Boundaries
bbed = 4
bsur = 5

# Get date for loading satellite image
fid = open(DIRM+"mesh_info.txt","r")
lines = fid.readlines()
date = lines[1][7:-1]
time = datelib.date_to_fracyear(int(date[0:4]),int(date[4:6]),int(date[6:]))
fid.close()
del fid,lines

# Mesh boundaries
extent = np.loadtxt(DIRM+"inputs/mesh_extent.dat")
try:
  hole1 = np.loadtxt(DIRM+"inputs/mesh_hole1.dat")
  hole2 = np.loadtxt(DIRM+"inputs/mesh_hole2.dat")
  holes=[hole1,hole2]  
except:
  holes = []
  
# Load data for first regularization parameter
bed_3D = elmerreadlib.saveline_boundary(DIRR+DIRREGPAR,method+"_beta.dat",bbed,['velocity','beta'])
surf_3D = elmerreadlib.saveline_boundary(DIRR+DIRREGPAR,method+"_beta.dat",bsur,['velocity','vsurfini'])
taub_3D = elmerreadlib.grid3d(bed_3D,'taub',holes,extent)
vel_3D = elmerreadlib.grid3d(surf_3D,'velocity',holes,extent)
velmes_3D = elmerreadlib.grid3d(surf_3D,'vsurfini',holes,extent)

if glacier == 'Kanger':
  xmin = 468000.
  xmax = 500800.
  ymin = -2300000.
  ymax = -2260200.
elif glacier == 'Helheim':
  xmin = 280200.
  xmax = 313000.
  ymin = -2585100.
  ymax = -2545300.

# Get satellite image for background; we want the satellite image captured on the closest date
#satimage,sattime,sattype = glaclib.load_satimages(glacier, xmin, xmax, ymin, ymax, time1=time, time2=time, data='all')
# Image for plotting
try:
  if glacier == "Helheim":
    xsat,ysat,satimage = geotifflib.readrgb(os.path.join(os.getenv("DATA_HOME"),"Imagery/Landsat/Helheim/TIF/20140704140535_LC82330132014185LGN00.tif"))
  elif glacier == "Kanger":
    xsat,ysat,satimage = geotifflib.readrgb(os.path.join(os.getenv("DATA_HOME"),"Imagery/Landsat/Kanger/TIF/20140706135251_LC82310122014187LGN00.tif"))
  sat_image = bool(True)
except:
  print "No satellite image for background"
  sat_image = bool(False)

fig=plt.figure(figsize=(3.4,4.0))
ax = plt.gca()
if sat_image:
  plt.imshow(satimage[:,:,0],extent=[xsat[0],xsat[-1],ysat[0],ysat[-1]],origin='lower',cmap='gray',clim=[0,0.6])
p = plt.imshow(taub_3D[2]*1e3,extent=[taub_3D[0][0],taub_3D[0][-1],taub_3D[1][0],taub_3D[1][-1]],origin='lower',vmin=0,vmax=500)
plt.xticks([])
plt.yticks([])
ax.axis('equal')
plt.tight_layout()
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])

xmin,xmax = plt.xlim()
ymin,ymax = plt.ylim()

path = matplotlib.path.Path([[0.58*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
  			[0.98*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
  			[0.98*(xmax-xmin)+xmin,0.76*(ymax-ymin)+ymin],
  			[0.58*(xmax-xmin)+xmin,0.76*(ymax-ymin)+ymin],
  			[0.58*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin]])
patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1)
ax.add_patch(patch)
cbaxes = fig.add_axes([0.61, 0.90, 0.28, 0.02]) 
cb = plt.colorbar(p,cax=cbaxes,orientation='horizontal',ticks=[0,200,400]) 
ax.text(xmin+0.60*(xmax-xmin),ymin+0.86*(ymax-ymin),'Basal shear stress',fontsize=8)
ax.text(xmin+0.70*(xmax-xmin),ymin+0.82*(ymax-ymin),'(kPa)',fontsize=8)
if label != 'none':
  ax.text(xmin+0.02*(xmax-xmin),ymin+0.93*(ymax-ymin),label,fontweight='bold',fontsize=10)
cb.ax.tick_params(labelsize=8)
ax.plot([xmin+0.61*(xmax-xmin),xmin+0.61*(xmax-xmin)+5e3],[ymin+0.80*(ymax-ymin),ymin+0.80*(ymax-ymin)],'k',linewidth=1.5)
ax.plot([xmin+0.61*(xmax-xmin),xmin+0.61*(xmax-xmin)],[ymin+0.80*(ymax-ymin),ymin+0.78*(ymax-ymin)],'k',linewidth=1.5)
ax.plot([xmin+0.61*(xmax-xmin)+5e3,xmin+0.61*(xmax-xmin)+5e3],[ymin+0.80*(ymax-ymin),ymin+0.78*(ymax-ymin)],'k',linewidth=1.5)
ax.text(xmin+0.64*(xmax-xmin)+5e3,ymin+0.78*(ymax-ymin),'5 km',fontsize=8)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_taub_"+date+"_"+regpar+".pdf"),format="PDF",dpi=400,transparent=True)
plt.close()


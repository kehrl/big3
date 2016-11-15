# Makes a 3D mesh of Helheim Glacier
# 
# LMK, UW, 06/01/2014
# Last updated 06/14/2014

import os
import shutil
import sys
import vellib, datelib, glaclib, flowparameterlib, meshlib, inverselib, climlib, elmerreadlib
from subprocess import call
from scipy.interpolate import RegularGridInterpolator
import numpy as np
import argparse



##########
# Inputs #
##########

# Get inputs to file
parser = argparse.ArgumentParser()
parser.add_argument("-glacier",dest="glacier",required = True, 
        help = "Name of glacier (Kanger or Helheim)")
parser.add_argument("-output", dest="output", required = True,
        help = "Name of output mesh.")
parser.add_argument("-mesh", dest="meshshp", required = True,
        help = "Name for input shapefile.")
parser.add_argument("-d", dest="date", required = True,
        help = "Date for mesh.")
parser.add_argument("-n", dest="n", required = True,
        help = "Number of partitions.")
parser.add_argument("-bname", dest="bedname", required = False,default='morlighem',
        help = "Name of bed file (smith,morlighem,cresis).")
parser.add_argument("-bmodel", dest="bedmodel", required = False,default='aniso',
        help = "Type of bed (aniso,iso).")
parser.add_argument("-bsmooth", dest="bedsmooth", type=int,required = False,
	default=4,help = "Smoothness of bed (1-8).")
parser.add_argument("-dx", dest="dx", required = False,default='none',
	help = "Grid size for gridded products.")
parser.add_argument("-lc", dest="lc", type=int,required = False,nargs='+',
	default=[1000,1000,3000,5000],\
	help = "Four numbers that define the mesh resolution for grounding-line (1000 m),channels (1000 m),regions near channels (3000 m), and entire mesh (5000 m).")
parser.add_argument("-zb", dest="bottomsurface", default = 'iceshelf',
        help = "Use 'iceshelf' base or 'bed' as bottom surface for mesh.",required = False)
parser.add_argument("-temperature",dest="temp",required=False, default='-10.0',
        help = "Ice temperature in deg C (or model out).")
parser.add_argument("-ssa",dest="ssa",required=False, default=False,
        help = "SSA model.")

#################
# Get arguments #
#################

args, _ = parser.parse_known_args(sys.argv)

date = args.date
partitions = args.n
bedname = args.bedname
bedmodel = args.bedmodel
bedsmoothing = args.bedsmooth
outputmeshname = args.output
meshshp = args.meshshp
glacier = args.glacier
dx = args.dx
bottomsurface = args.bottomsurface
temperature = args.temp
ssa = args.ssa
  
# Mesh refinement
lc3,lc2,lc4,lc1 = args.lc

# Directories
DIRM = os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/"+outputmeshname+"/")
DIRX = os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/Glaciers/3D/"+glacier+"/")
inputs = os.path.join(DIRM+"inputs/")

# Make mesh directories
if not(os.path.isdir(DIRM)):
  os.makedirs(DIRM)
  os.makedirs(DIRM+"/inputs")

# Densities for finding floating ice
rho_i = 917.0
rho_sw = 1025.0
yearinsec = 365.25*24*60*60

# Time
time = datelib.date_to_fracyear(int(date[0:4]),int(date[4:6]),int(date[6:8]))

#################
# Mesh Geometry #
#################

# Mesh exterior
if meshshp.endswith('nofront') or meshshp.endswith('nofront.shp'):
  exterior = glaclib.load_extent(glacier,time,nofront_shapefile=meshshp)
else:
  exterior = meshlib.shp_to_xy(DIRX+meshshp)
np.savetxt(inputs+"mesh_extent.dat",exterior[:,0:2])

# Mesh holes
holes = []
if os.path.isfile(DIRX+"glacier_hole1.shp"):
  hole1 = meshlib.shp_to_xy(DIRX+"glacier_hole1")
  np.savetxt(inputs+"mesh_hole1.dat",hole1[:,0:2])
  holes.append({'xy': hole1})
if os.path.isfile(DIRX+"glacier_hole2.shp"):
  hole2 = meshlib.shp_to_xy(DIRX+"glacier_hole2")
  np.savetxt(inputs+"mesh_hole2.dat",hole2[:,0:2])
  holes.append({'xy': hole2})

# Add locations for refinement
refine = meshlib.shp_to_xy(DIRX+"refine")

#Set output name for gmsh file
file_2d=os.path.join(DIRM+"mesh2d")

##################################################################
# Save file with mesh inputs so we know how the mesh was created #
##################################################################
  
fid = open(DIRM+'mesh_info.txt','w')
fid.write('glacier = {}\n'.format(glacier))
fid.write('date = {}\n'.format(date))
fid.write('meshshapefile = {}\n'.format(meshshp))
fid.write('lc1 = {}\n'.format(lc1))
fid.write('lc2 = {}\n'.format(lc2))
fid.write('lc3 = {}\n'.format(lc3))
fid.write('lc4 = {}\n'.format(lc4))
fid.write('bed = {}\n'.format(bedname))
fid.write('bed model = {}\n'.format(bedmodel))
fid.write('bed smoothness = {}'.format(bedsmoothing))
  
fid.close()
  
#############
# Make mesh #
#############

# Gmsh .geo file
x,y,zbed,zsur,zbot = meshlib.xy_to_gmsh_3d(glacier,date,exterior,holes,refine,DIRM,\
		lc1,lc2,lc3,lc4,bedname=bedname,bedmodel=bedmodel,bedsmoothing=bedsmoothing,\
                rho_i=rho_i,rho_sw=rho_sw,dx=dx,bottomsurface=bottomsurface)

# Create .msh file
call(["gmsh","-1","-2",file_2d+".geo", "-o",os.path.join(os.getenv("HOME"),\
		file_2d+".msh")])

# Create elmer mesh
call(["ElmerGrid","14","2",file_2d+".msh","-autoclean"])

# Partition mesh for parallel processing
os.chdir(DIRM)
call(["ElmerGrid","2","2","mesh2d","dir","-metis",partitions,"0"])

# Output as gmsh file so we can look at it
call(["ElmerGrid","2","4","Elmer"])

##########################################
# Print out velocity data for inversions #
##########################################

# Output files for velocities in x,y directions (u,v)
u,v = vellib.inversion_3D(glacier,x,y,time,inputs,dx=dx)

if ssa:
  fid = open(inputs+"velocity.xyuv","w")
  for i in range(0,len(x)):
    for j in range(0,len(y)):
      fid.write('{0} {1} {2} {3} \n'.format(x[i],y[j],u[j,i],v[j,i]))
  fid.close()

################################################################
# Get climate variables & calculate temperatures at ice divide #
################################################################
  
# Set low resolution mesh (no reason to overkill mesh size for climate variables given 
# spatial resolution of RACMO2.3).
xt2m = np.arange(x[0],x[-1],1e3)
yt2m = np.arange(y[0],y[-1],1e3)
  
# Get average 2-m temperatures and surface mass balance
timet2m,t2m = climlib.racmo_interpolate_to_cartesiangrid(xt2m,yt2m,'t2m',epsg=3413,maskvalues='both',timing='mean')
timesmb,smb = climlib.racmo_interpolate_to_cartesiangrid(xt2m,yt2m,'smb',epsg=3413,maskvalues='both',timing='mean')
#ggrid = bedlib.geothermalflux_grid(xt2m,yt2m,model='davies',method='nearest')

# Set maximum temperature to -1 deg C
ind = np.where(t2m > 272.15)
t2m[ind] = 272.15

fidt2m = open(inputs+"t2m.xy","w")
fidsmb = open(inputs+"smb.xy","w")
fidt2m.write('{}\n{}\n'.format(len(xt2m),len(yt2m)))
fidsmb.write('{}\n{}\n'.format(len(xt2m),len(yt2m)))
#fidgeo = open(inputs+"geothermal.xy","w")
#fidgeo.write('{}\n{}\n'.format(len(xt2m),len(yt2m)))
for i in range(0,len(xt2m)):
  for j in range(0,len(yt2m)):
    fidt2m.write('{0} {1} {2}\n'.format(xt2m[i],yt2m[j],t2m[j,i]))
    fidsmb.write('{0} {1} {2}\n'.format(xt2m[i],yt2m[j],smb[j,i]))
fidt2m.close()
fidsmb.close()
 
# Get surface heights on same grid as temperature and surface mass balance so that
# we can get vertical steady state temperatures.
xgrid,ygrid = np.meshgrid(xt2m,yt2m)
f = RegularGridInterpolator((y,x),zsur-zbed)
Hflat = f((ygrid.flatten(),xgrid.flatten()))
H = np.reshape(Hflat,[len(yt2m),len(xt2m)])
del Hflat,f,xgrid,ygrid
  
# Get 3D grid of temperatures  
T = flowparameterlib.steadystate_vprofile(H,t2m,smb,levels=12)
 
fidT = open(inputs+"tsteady_icedivide.xyz", "w")
fidT.write("{0}\n{1}\n{2}\n".format(len(xt2m), len(yt2m), len(T[0,0,:])))
for j in range(len(xt2m)):
  for i in range(len(yt2m)):
    fidT.write("{0} {1} ".format(xt2m[j], yt2m[i]))
    for k in range(len(T[0,0,:])):
      fidT.write("{0} ".format(T[i, j, k]))
    fidT.write("\n")
fidT.close()

del xt2m,yt2m,timet2m,t2m,fidt2m, H, fidT
# del fidgeo, ggrid

##################################################
# Get ice temperatures from model (if it exists) #
##################################################

print "Trying to pull temperatures from model...\n"

# Set up lower resolution grid and interpolate variables to that grid
xT = np.arange(x[0],x[-1],100)
yT = np.arange(y[0],y[-1],100)
xTgrid,yTgrid = np.meshgrid(xT,yT)
f = RegularGridInterpolator((y,x),zsur)
zsT = np.reshape(f((yTgrid.flatten(),xTgrid.flatten())),[len(yT),len(xT)])
f = RegularGridInterpolator((y,x),zbed)
zbT = np.reshape(f((yTgrid.flatten(),xTgrid.flatten())),[len(yT),len(xT)])
f = RegularGridInterpolator((y,x),u)
uT = np.reshape(f((yTgrid.flatten(),xTgrid.flatten())),[len(yT),len(xT)])
f = RegularGridInterpolator((y,x),v)
vT = np.reshape(f((yTgrid.flatten(),xTgrid.flatten())),[len(yT),len(xT)])  

if temperature == 'model' and not(ssa):
  flowT,flowA,flowU,flowV = flowparameterlib.load_temperature_model(glacier,xT,yT,outputdir=inputs)
if temperature == 'model' and ssa:
  dir = os.path.join(os.getenv("MODEL_HOME"),glacier+"/Outputs/")
  shutil.copy2(dir+"ssa_flowA.xy",inputs+"ssa_flowA.xy")
  
#################################################################
# Calculate basal sliding speed using SIA for inflow boundaries #
#################################################################

print "Calculating basal sliding speed for inflow and ice divide boundaries and guessing a beta...\n"
if temperature == 'model' and ssa:
  xflowA,yflowA,flowA = elmerreadlib.input_file(dir+"ssa_flowA.xy")
  if (len(xflowA) != len(xT)) or (len(yflowA) != len(yT)):
    f = RegularGridInterpolator((yflowA,xflowA),flowA)
    flowA = np.reshape(f((yTgrid.flatten(),xTgrid.flatten())),[len(yT),len(xT)])
    # If there are any nans
    ind = np.where(np.isnan(flowA))
    flowA[ind] = flowparameterlib.arrhenius(263.15)
  del dir 
  ub_all,vb_all,beta_all = inverselib.guess_beta(xT,yT,zsT,zbT,uT,vT,frac=0.5,A=flowA*yearinsec*1e18) 
elif temperature == 'model' and not(ssa):
  # Try to use depth-averaged modeled temperatures for guessing
  ub_all,vb_all,beta_all = inverselib.guess_beta(xT,yT,zsT,zbT,uT,vT,frac=0.5,A=np.mean(flowA,axis=2)*yearinsec*1e18)
else:
  A = flowparameterlib.arrhenius(273.15+float(temperature))*yearinsec*1.0e18
  # If there are no modeled temperatures, then just use a constant flow law parameter
  ub_all,vb_all,beta_all = inverselib.guess_beta(xT,yT,zsT,zbT,uT,vT,frac=0.5,A=A)  

# Write out basal velocities and initial guess for beta
fidub = open(inputs+"ubdem.xy","w")
fidvb = open(inputs+"vbdem.xy","w")
fidbeta = open(inputs+"beta0.xy","w")
fidvb.write('{}\n{}\n'.format(len(xT),len(yT)))
fidub.write('{}\n{}\n'.format(len(xT),len(yT)))
fidbeta.write('{}\n{}\n'.format(len(xT),len(yT)))
for i in range(0,len(xT)):
  for j in range(0,len(yT)):
    fidub.write('{0} {1} {2}\n'.format(xT[i],yT[j],ub_all[j,i]))
    fidvb.write('{0} {1} {2}\n'.format(xT[i],yT[j],vb_all[j,i]))
    fidbeta.write('{0} {1} {2}\n'.format(xT[i],yT[j],beta_all[j,i]))
fidub.close()
fidvb.close()
fidbeta.close()

del xTgrid,yTgrid


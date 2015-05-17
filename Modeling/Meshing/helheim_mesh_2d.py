# Make a flowline mesh of Helheim Glacier using "MshGlacier"
# Currently uses bedrock from Morlighem and others (2014) and surface from GIMP DEM
#
# LMK, UW, 06/01/2014
# Last updated 06/14/2014

import os
import sys
import numpy as np
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Helheim/Tools"))
import helheim_velocity
import elmer_mesh as mesh
import dist
import shapefactor,flowparameter
import elmer_mesh as mesh
import dist
import shapefactor,flowparameter
import subprocess, shutil
import matplotlib.pyplot as plt
from subprocess import call
from scipy import interpolate
from scipy import signal
import geotiff

##########
# Inputs #
##########

MESHNAME='MorlighemNew_SmoothedVelocity'
file_mesh_out="Elmer"

DIRS=os.path.join(os.getenv("HOME"),"Code/Helheim/Modeling/SolverFiles/Flowline/")
DIRM=os.path.join(os.getenv("HOME"),"Models/Helheim/Meshes/Flowline/"+MESHNAME+"/")
DIRR=os.path.join(os.getenv("HOME"),"Models/Helheim/Results/Flowline/"+MESHNAME+"/")

lc=[50,100,500] # characteristic length scales for meshing
lc_d=[0,10000,40000] # transition points between different mesh sizes
width_filt_len = 2000.0 # length along flowline for filtering widths
partitions="4" # Number of partitions
layers=10 # Number of extrusions for mesh

# Shapefiles for flowline
file_flowline_in = os.path.join(os.getenv("HOME"),"Data/ShapeFiles/Glaciers/Flowlines/Helheim/helheim_center_flowline")
file_leftside_in = os.path.join(os.getenv("HOME"),"Data/ShapeFiles/Glaciers/Flowlines/Helheim/helheim_flowline_side2")
file_rightside_in = os.path.join(os.getenv("HOME"),"Data/ShapeFiles/Glaciers/Flowlines/Helheim/helheim_flowline_side1")

# Bed and surface topographies
## Bed inputs
file_bed_in1=os.path.join(os.getenv("HOME"),"Data/Bed/Morlighem_2014/morlighem_bed.tif")
file_bed_in2=os.path.join(os.getenv("HOME"),"Data/Bed/CreSIS/helheim_flightline_05212001_good_nsidc.dat")	

## Surface inputs
file_surf_gimp=os.path.join(os.getenv("HOME"),"Data/Elevation/Gimp/gimpdem3_1.tif")
file_surf_wv=[os.path.join(os.getenv("HOME"),"Data/Elevation/Worldview/Helheim/20120520_1443_102001001BD45E00_102001001C07FB00-DEM_32m_trans.tif"),
			  os.path.join(os.getenv("HOME"),"Data/Elevation/Worldview/Helheim/20120513_1410_102001001B4C6F00_102001001BD7E800-DEM_32m_trans.tif"),
			  os.path.join(os.getenv("HOME"),"Data/Elevation/Worldview/Helheim/20120624_1421_102001001B87EF00_102001001B1FB900-DEM_32m_trans.tif")]

## File to use for terminus position
file_terminus = os.path.join(os.getenv("HOME"),"Data/ShapeFiles/IceFronts/Helheim/2012-174_TSX.shp")

## Velocity for inversion
file_velocity_in=os.path.join(os.getenv("HOME"),"Data/Velocity/TSX/Helheim/track-27794/mosaicOffsets")

###########
# Outputs #
###########

file_shapefactor_out = DIRM+"Inputs/width.dat"
file_bedrock_out = DIRM+"Inputs/roughbed.dat"
file_flowline_out = DIRM+"Inputs/flowline.dat"

# Output files for measured velocities along flowline
file_velocity_out=DIRM+"Inputs/velocity.dat"

###################################################################################
# Make box mesh and then use "MshGlacier" to add surface and bedrock topographies #
###################################################################################

# Make mesh directories
if not(os.path.isdir(DIRM)):
  #os.makedirs(DIRR)
  os.makedirs(DIRM)
  os.makedirs(DIRM+"Inputs")

# Flowline coordinates
flowline = mesh.shp_to_flowline(file_flowline_in)
del file_flowline_in

# Create geo file for flowline
flowline = mesh.xy_to_gmsh_box(flowline,file_terminus,DIRM,file_mesh_out,file_bed_in1,file_surf_gimp,lc,lc_d,layers)

# Adjust bed DEM to 2001 bed profile near terminus
## The present flowline runs off the Morlighem bed DEM, so we need to set the bed near the 
## terminus to the bed elevation measurements from 2001
ind=np.where(flowline[:,1]>309000)
bed2001=np.loadtxt(file_bed_in2,skiprows=1)
bed2001=bed2001[3180:3297,:]
flowline[ind,3]=np.interp(flowline[ind,1],bed2001[:,0],bed2001[:,7]+130)

## Save variable with rest of bed profile for grounded solver
fid = open(DIRM+file_mesh_out+"_bed.dat",'w')
fid.write('{} {} \n'.format(flowline[0,0]-100,flowline[0,3]))
for i in range(0,len(flowline[:,0])):
  fid.write('{} {} \n'.format(flowline[i,0],flowline[i,3]))
fid.write('{} {} \n'.format(flowline[i,0]+100,flowline[i,3]))
fid.close()

# Adjust surface DEM to worldview DEM
flowline_old = np.array(flowline)
for file in file_surf_wv:
  [xs,ys,zs] = geotiff.read(file)
  zs_dem = interpolate.RectBivariateSpline(ys,xs,zs)
  zs_interp = zs_dem.ev(flowline[:,2],flowline[:,1])
  ind=np.where(zs_interp > 100.0)
  flowline[ind,4]=zs_interp[ind]
flowline[ind[0][-1]+1:-1,4] = np.mean(flowline[3370:3390,4])
ind = np.where(abs(np.diff(flowline[:,4])) > 20)
for i in range(0,len(ind[0])):
  satisfied1 = 0; satisfied2 = 0
  ind1=ind[0][i]-1
  ind2=ind[0][i]+1
  while not(satisfied1):
    if ind1 not in ind[0]:
      satisfied1 = 1
    else:
      ind1 = ind1-1
  while not(satisfied2):
    if ind2 not in ind[0]:
      satisfied2 = 1
    else:
      ind2 = ind2+1
  flowline[ind[0][i],4] = np.mean([flowline[ind1,4],flowline[ind2,4]])    
surf_filt_len=float(500)
cutoff=(1/surf_filt_len)/(1/(np.diff(flowline[1:3,0])*2))
b,a=signal.butter(4,cutoff,btype='low')
flowline[:,4]=signal.filtfilt(b,a,flowline[:,4])

fid = open(DIRM+file_mesh_out+"_surf.dat",'w')
fid.write('{} {} \n'.format(flowline[0,0]-100,flowline[0,4]))
for i in range(0,len(flowline[:,0])):
  fid.write('{} {} \n'.format(flowline[i,0],flowline[i,4]))
fid.write('{} {} \n'.format(flowline[i,0]+100,flowline[i,4]))
fid.close()

del lc, lc_d, file_bed_in1, file_bed_in2, file_surf_gimp, bed2001

# Create box mesh for extrusion through MshGlacier
call(["gmsh", "-2",file_mesh_out+".geo", "-o",file_mesh_out+".msh"])

call(["ElmerGrid","14","2",file_mesh_out+".msh"])

# Now we can call MshGlacier
os.chdir(DIRM)
call(["MshGlacier"])

##################
# Partition Mesh #
##################
print "\n## Partitioning mesh into ", partitions, "parts"
os.chdir(DIRM)
call(["ElmerGrid","2","2",file_mesh_out,"dir","-metis",partitions,"0"])

##########################################
# Update the Gmsh file for visualization #
##########################################
call(["ElmerGrid","2","4",file_mesh_out])

###################################################
# Load, find, and write velocities along flowline #
###################################################

# Get velocity along flowline from velocity profiles
helheim_velocity.inversion_2D(flowline[:,1],flowline[:,2],flowline[:,0],file_velocity_in,DIRM+"Inputs/")


#################################
# Print out temperature profile #
#################################
# Get modeled temperatures from Kristin's work
kristin_file=os.path.join(os.getenv("HOME"),"Data/Climate/IceTemperature/Helheim/helheim_TA.xyz")
tempdata=np.genfromtxt(kristin_file,delimiter=',')
tempdata=np.delete(tempdata,(0),axis=0)

# Cut out just one temperature profile for flowline
tempdata=tempdata[135479:141591,:]
xdata,xinds=np.unique(tempdata[:,0],return_index=True)
ydata=tempdata[xinds,1]

# Compute normalized depths for Kristin's temperatures
surf=np.zeros_like(xdata)
bed=np.zeros_like(xdata)
height=np.zeros_like(xdata)
ndata=np.zeros_like(tempdata[:,0])
for i in range(0,len(xdata)):
  colinds=np.array(np.where(tempdata[:,0]==xdata[i]))
  surf[i]=np.max(tempdata[colinds,2]) # Surface elevation
  bed[i]=np.min(tempdata[colinds,2]) # Bed elevation
  height[i]=surf[i]-bed[i] #height
  ndata[colinds]=(tempdata[colinds,2]-bed[i])/height[i] #normalized value
	
# Compute distance along Kristin's temperature profile
ddata=dist.transect(tempdata[:,0],tempdata[:,1])

# Compute distance between our flowline and Kristin's, so that we can adjust
# the distance in Kristin's temperature profile so that it aligns
dist1=np.sqrt((xdata[0]-flowline[0,1])**2+(ydata[0]-flowline[0,2])**2)
dist2=np.sqrt((xdata[1]-flowline[0,1])**2+(ydata[1]-flowline[0,2])**2)
ddata=ddata-ddata[200]*(dist1/(dist1+dist2))

# Write out depths for flowline so that we can use those depths for interpolation
pts=np.zeros([len(flowline[:,0])*layers*2,3])
for i in range(0,len(flowline[:,0])):
  pts[layers*2*i:layers*2*(i+1),0]=flowline[i,0]
  pts[layers*2*i:layers*2*(i+1),1]=np.linspace(0,1,layers*2)
  pts[layers*2*i:layers*2*(i+1),2]=np.linspace(flowline[i,3],flowline[i,4],layers*2)

T=interpolate.griddata(np.column_stack([ddata,ndata]),tempdata[:,3],pts[:,0:2],method='linear')
A=flowparameter.arrhenius(T)

# Write out flow law parameter at each node
nanind=np.where(np.isnan(A))
A=np.delete(A,nanind,axis=0)
pts=np.delete(pts,nanind,axis=0)
fid = open(DIRM+"Inputs/flowparameters.dat","w")
fid.write('{}\n'.format(len(A)))
for i in range(0,len(A)):
  fid.write('{0} {1} {2}\n'.format(pts[i,0],pts[i,2],A[i]))
fid.close() 
del dist1,dist2,xdata,ydata,surf,bed,height,ndata,T,A

#####################################
# Print out bedrock for elmersolver #
#####################################
print "\n## Printing out bedrock from mesh file for grounded solver ##"

meshnode=np.genfromtxt(DIRM+file_mesh_out+"/mesh.nodes")

# Take bed values from mesh for consistency using a Lagrangian mesh
values=np.unique(meshnode[:,2])
flowlinenodes=np.zeros([len(values),2])
for i in range(0,len(values)):
  ind=np.where(meshnode[:,2] == values[i])
  flowlinenodes[i,0]=values[i]
  flowlinenodes[i,1]=np.min(meshnode[ind,3])

# Add bed measurements in front of mesh nodes, so that the glacier can advance
## Interpolate forward using the same grid size as near the terminus
delta = np.min(np.diff(flowlinenodes[:,0]))
xnew = np.arange(flowlinenodes[-1,0],flowline[-1,0],delta)
bedflow=np.zeros([len(xnew)-1,2])
bedflow[:,0]=xnew[1:len(xnew)]
bedflow[:,1]=np.interp(xnew[1:len(xnew)],flowline[:,0],flowline[:,3])

bed=np.row_stack([flowlinenodes,bedflow])

fid = open(file_bedrock_out,'w')
fid.write('{}\n'.format(len(bed[:,0])))
for i in range(0,len(bed[:,0])):
  fid.write("{} {}\n".format(bed[i,0],bed[i,1]))
fid.close()
del fid

####################################################################################
# Print out width and dw/dx for shapefactor and lateral convergence in elmersolver #
####################################################################################

# Calculate width
width = shapefactor.glacierwidth(flowline,file_rightside_in,file_leftside_in,width_filt_len)
width = np.interp(flowlinenodes[:,0],flowline[:,0],width)
thick = np.interp(flowlinenodes[:,0],flowline[:,0],flowline[:,4]-flowline[:,3])

# Calculate shape factor
# f = shapefactor.trapezoid(width,thick)
dwdx = shapefactor.dwdx(flowlinenodes[:,0],width)

# Write out width and dw/dx for elmersolver 
R = len(width)
fid = open(file_shapefactor_out,'w')
fid.write("{0} \n".format(R))
fid.write("{0} {1:.4g} {2}\n".format(-10,width[0],dwdx[0]))
for i in range(0,R):
  fid.write("{0} {1:.4g} {2}\n".format(flowlinenodes[i,0],width[i],dwdx[i]))
fid.close()
del fid

del thick, dwdx, R

#########################################################
# Finally export coordinates of flowline for future use #
#########################################################
print "\n ## Saving flowline coordinates ##"
R = len(flowline[:,0])
fid = open(file_flowline_out,'w')
fid.write("{}\n".format(R))
for i in range(0,R):
  fid.write("{} {} {} {} {}\n".format(flowline[i,0],flowline[i,1],flowline[i,2],flowline[i,3],flowline[i,4]))
fid.close()
del fid 


# Save this file in MESH directory for future reference
shutil.copy('helheim_mesh_2d.py',DIRM+'helheim_mesh_2d.py')
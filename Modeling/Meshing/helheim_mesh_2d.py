# Make a flowline mesh of Helheim Glacier using "MshGlacier"
# Currently uses bedrock from Morlighem and others (2014) and surface from GIMP DEM
#
# LMK, UW, 06/01/2014
# Last updated 06/14/2014

import os
import sys
import numpy as np
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"Util/Modules"))
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"BigThreeGlaciers/Tools"))
import velocity, bed, elevation, fracyear, glacier_flowline, icefronts
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

##########
# Inputs #
##########

glacier = 'Helheim'

## Date for DEM
date = '20120624' # date for terminus DEM
#date = '20110319'

# Densities for finding the bottom of the ice shelf where the ice is floating
rho_i = 917.0
rho_sw = 1020.0

# time for mesh
time = fracyear.date_to_fracyear(int(date[0:4]),int(date[4:6]),int(date[6:8]))

# Mesh name
MESHNAME='DEM'+date+'Low'
file_mesh_out="Elmer"

# Import directories
DIRS=os.path.join(os.getenv("CODE_HOME"),"BigThreeGlaciers/Modeling/SolverFiles/Flowline/")
DIRM=os.path.join(os.getenv("MODEL_HOME"),glacier+"/Meshes/Flowline/"+MESHNAME+"/")
DIRR=os.path.join(os.getenv("MODEL_HOME"),glacier+"/Results/Flowline/"+MESHNAME+"/")

lc=[100,1000,2500] # characteristic length scales for meshing
#lc=[50,100,200] # characteristic length scales for meshing
lc_d=[0,5000,40000] # transition points between different mesh sizes
filt_len = 1000.0 # length along flowline for filtering widths
partitions="4" # Number of partitions
layers=10 # Number of extrusions for mesh

# Shapefiles for flowline
file_flowline_in = 'center_flowline'
file_leftside_in = os.path.join(os.getenv("HOME"),"Data/ShapeFiles/Glaciers/Flowlines/Helheim/helheim_flowline_side2_wall")
file_rightside_in = os.path.join(os.getenv("HOME"),"Data/ShapeFiles/Glaciers/Flowlines/Helheim/helheim_flowline_side1_wall")

###########
# Outputs #
###########

file_shapefactor_out = DIRM+"Inputs/width.dat"
file_bedrock_out = DIRM+"Inputs/roughbed.dat"
file_flowline_out = DIRM+"Inputs/flowline.dat"

#################
# Make box mesh #
#################

# Make mesh directories if mesh directories do not already exist
if not(os.path.isdir(DIRM)):
  #os.makedirs(DIRR)
  os.makedirs(DIRM)
  os.makedirs(DIRM+"Inputs")

# Flowline coordinates
x,y,zb,dists = glacier_flowline.load(glacier,shapefilename=file_flowline_in,filt_len='none')

# Surface elevations along flowline for chosen date
zs = elevation.dem_continuous_flowline(x,y,dists,glacier,date,filt_len=filt_len,verticaldatum='geoid',fillin=True)

# Terminus position for date
terminus = icefronts.position(x,y,dists,glacier,[int(date[0:4]),int(date[4:6]),int(date[6:8])])

# Create geo file for flowline
flowline = mesh.xy_to_gmsh_box(x,y,dists,zb,zs,terminus,glacier,DIRM,file_mesh_out,lc,lc_d,layers,filt_len=filt_len,rho_i=rho_i,rho_sw=rho_sw)

##################################
# Use MshGlacier to extrude mesh #
##################################

os.chdir(DIRM)
# Create box mesh for extrusion through MshGlacier
call(["gmsh", "-2",file_mesh_out+".geo", "-o",file_mesh_out+".msh"])

call(["ElmerGrid","14","2",file_mesh_out+".msh"])

# Now we can call MshGlacier
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
dvel,vel=velocity.inversion_2D(flowline[:,1],flowline[:,2],flowline[:,0],glacier,time,DIRM+"Inputs/",filt_len)

#################################
# Print out temperature profile #
#################################

# Get modeled temperatures from Kristin's work
kristin_file=os.path.join(os.getenv("HOME"),"Data/Climate/IceTemperature/Helheim/helheim_TA.xyz")
tempdata = np.genfromtxt(kristin_file,delimiter=',')
tempdata = np.delete(tempdata,(0),axis=0)

# Cut out just one temperature profile for flowline
tempdata = tempdata[135479:141591,:]
xdata,xinds = np.unique(tempdata[:,0],return_index=True)
ydata = tempdata[xinds,1]

# Compute normalized depths for Kristin's temperatures
surf = np.zeros_like(xdata)
bed = np.zeros_like(xdata)
height = np.zeros_like(xdata)
ndata = np.zeros_like(tempdata[:,0])
elevation = np.zeros_like(tempdata[:,0])
for i in range(0,len(xdata)):
  colinds = np.array(np.where(tempdata[:,0]==xdata[i]))
  surf[i] = np.max(tempdata[colinds,2]) # Surface elevation
  bed[i] = np.min(tempdata[colinds,2]) # Bed elevation
  height[i] = surf[i]-bed[i] #height
  elevation[colinds] = tempdata[colinds,2]
  ndata[colinds] = (tempdata[colinds,2]-bed[i])/height[i] #normalized value
	
# Compute distance along Kristin's temperature profile
ddata = dist.transect(tempdata[:,0],tempdata[:,1])

# Compute distance between our flowline and Kristin's, so that we can adjust
# the distance in Kristin's temperature profile so that it aligns
dist1 = np.sqrt((xdata[0]-flowline[0,1])**2+(ydata[0]-flowline[0,2])**2)
dist2 = np.sqrt((xdata[1]-flowline[0,1])**2+(ydata[1]-flowline[0,2])**2)
ddata = ddata-ddata[200]*(dist1/(dist1+dist2))+flowline[0,0]

# Write out depths for flowline so that we can use those depths for interpolation
pts = np.zeros([len(flowline[:,0])*layers*2,3])
for i in range(0,len(flowline[:,0])):
  pts[layers*2*i:layers*2*(i+1),0]=flowline[i,0]
  pts[layers*2*i:layers*2*(i+1),1]=np.linspace(0,1,layers*2)
  pts[layers*2*i:layers*2*(i+1),2]=np.linspace(flowline[i,3],flowline[i,4],layers*2)

T = interpolate.griddata(np.column_stack([ddata,ndata]),tempdata[:,3],pts[:,0:2],method='linear')
A = flowparameter.arrhenius(T)

# Write out flow law parameter at each node
nanind = np.where(np.isnan(A))
A = np.delete(A,nanind,axis=0)
pts = np.delete(pts,nanind,axis=0)
fid = open(DIRM+"Inputs/flowparameters.dat","w")
fid.write('{}\n'.format(len(A)))
for i in range(0,len(A)):
  fid.write('{0} {1} {2} {3}\n'.format(pts[i,0],pts[i,2],T[i]-273.15,A[i]))
fid.close() 
del dist1,dist2,xdata,ydata,surf,bed,height,ndata,T,A

#####################################
# Print out bedrock for elmersolver #
#####################################

meshnode = np.genfromtxt(DIRM+file_mesh_out+"/mesh.nodes")

# Take bed values from mesh for consistency using a Lagrangian mesh
values=np.unique(meshnode[:,2])
flowlinenodes=np.zeros([len(values),2])
for i in range(0,len(values)):
  ind=np.where(meshnode[:,2] == values[i])
  flowlinenodes[i,0]=values[i]
  flowlinenodes[i,1]=np.min(meshnode[ind,3])

# Add bed measurements in front of mesh nodes, so that the glacier can advance
# Interpolate forward using the same grid size as near the terminus
delta = np.min(np.diff(flowlinenodes[:,0]))
dists_new = np.r_[flowlinenodes[:,0],np.arange(flowlinenodes[-1,0]+delta,flowline[-1,0],delta)]
zb_new = np.interp(dists_new,dists,zb)

fid = open(file_bedrock_out,'w')
fid.write('{}\n'.format(len(dists_new)))
for i in range(0,len(dists_new)):
  fid.write("{0:.6f} {1:.6f}\n".format(dists_new[i],zb_new[i]))
fid.close()
del fid, delta, dists_new,zb_new,ind,values

####################################################################################
# Print out width and dw/dx for shapefactor and lateral convergence in elmersolver #
####################################################################################

# Calculate width
width = shapefactor.glacierwidth(flowline,file_rightside_in,file_leftside_in,filt_len)
width = np.interp(flowlinenodes[:,0],flowline[:,0],width)
thick = np.interp(flowlinenodes[:,0],flowline[:,0],flowline[:,4]-flowline[:,3])

# Calculate dw/dx for mass input due to convergence
dwdx = shapefactor.dwdx(flowlinenodes[:,0],width)

# Write out width and dw/dx for elmersolver 
R = len(width)
fid = open(file_shapefactor_out,'w')
fid.write("{0} \n".format(R))
fid.write("{0} {1:.4g} {2}\n".format(flowlinenodes[0,0]-10,width[0],dwdx[0]))
for i in range(0,R):
  fid.write("{0} {1:.4g} {2}\n".format(flowlinenodes[i,0],width[i],dwdx[i]))
fid.close()

del thick, dwdx, R, fid

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
os.chdir(os.path.join(os.getenv("HOME"),"Code/BigThreeGlaciers/Modeling/Meshing/"))
shutil.copy('helheim_mesh_2d.py',DIRM+'helheim_mesh_2d_'+MESHNAME+'.py')

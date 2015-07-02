# Makes a 3D mesh of Helheim Glacier
# 
# LMK, UW, 06/01/2014
# Last updated 06/14/2014

import os
import shutil
import sys
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
sys.path.append(os.path.join(os.getenv("HOME"),"Code/BigThreeGlaciers/Tools"))
import elmer_mesh as mesh
import elmer_inversion
import helheim_velocity, helheim_bed, helheim_elevation
from subprocess import call
from scipy.interpolate import *
import numpy as np
import flowparameter
import geotiff

##########
# Inputs #
##########

MESHNAME='HighResolution'

# Directories
DIRS=os.path.join(os.getenv("HOME"),"Code/Helheim/Modeling/SolverFiles/3D")
DIRM=os.path.join(os.getenv("HOME"),"Models/Helheim/Meshes/3D/"+MESHNAME+"/")
DIRR=os.path.join(os.getenv("HOME"),"Models/Helheim/Results/3D/")
DIRX=os.path.join(os.getenv("HOME"),"Data/ShapeFiles/Glaciers/3D/Helheim/")
Inputs=os.path.join(DIRM+"Inputs/")

# Make mesh directories
if not(os.path.isdir(DIRM)):
  os.makedirs(DIRM)
  os.makedirs(DIRM+"/Inputs")

# Mesh refinement
lc1=2000 # for entire mesh
lc2=500 # for channels close to the grounding line
lc3=250 # for grounding line
lc4=700 # for regions surrounding channels
lc1=2000 # for entire mesh
lc2=2000 # for channels close to the grounding line
lc3=2000 # for grounding line
lc4=2000 # for regions surrounding channels
levels=5 #levels of extrusion
partitions="4" # Number of partitions

# Bed and surface
file_bed = 'morlighem' 
file_surf = 'gimp'

# Velocity profile for inversion
file_velocity_in = os.path.join(os.getenv("HOME"),"Data/Velocity/TSX/Helheim/track-27794")

#################
# Mesh Geometry #
#################

# Mesh exterior
exterior = mesh.shp_to_xy(DIRX+"glacier_extent_normal")
mesh_file=Inputs+"mesh_extent.dat"
fid = open(mesh_file,"w")
for i in range(0,len(exterior[0])):
  fid.write('{} {}\n'.format(exterior[0][i],exterior[1][i]))
fid.close()

# Mesh holes
hole1 = mesh.shp_to_xy(DIRX+"glacier_hole1")
fid = open(Inputs+"mesh_hole1.dat","w")
for i in range(0,len(hole1[0])):
  fid.write('{} {}\n'.format(hole1[0][i],hole1[1][i]))
fid.close()

hole2 = mesh.shp_to_xy(DIRX+"glacier_hole2")
fid = open(Inputs+"mesh_hole2.dat","w")
for i in range(0,len(hole2[0])):
  fid.write('{} {}\n'.format(hole2[0][i],hole2[1][i]))
fid.close()

holes = []
holes.append({'xy': hole1})
holes.append({'xy': hole2})

# Add locations for refinement
refine = mesh.shp_to_xy(DIRX+"refine")

###########
# Outputs #
###########

#Set output name for gmsh file
file_2d=os.path.join(DIRM+"Planar")
file_3d=os.path.join(DIRM+"Elmer")

#############
# Make mesh #
#############

# Gmsh .geo file
x,y,zbed,zsur = mesh.xy_to_gmsh_3d(exterior,holes,refine,DIRM,lc1,lc2,lc3,lc4,file_bed,file_surf)

# Create .msh file
call(["gmsh","-1","-2",file_2d+".geo", "-o",os.path.join(os.getenv("HOME"),file_2d+".msh")])

# Create elmer mesh
call(["ElmerGrid","14","2",file_2d+".msh","-autoclean"])

# Extrude the mesh with ExtrudeMesh
call(["ExtrudeMesh",file_2d,file_3d,str(levels),"1","1","0","0","0","0",Inputs,"200","2","NaN"])

# Partition mesh for parallel processing
os.chdir(DIRM)
call(["ElmerGrid","2","2","elmer","dir","-metis",partitions])

# Output as gmsh file so we can look at it
call(["ElmerGrid","2","4","elmer"])

##########################################
# Print out velocity data for inversions #
##########################################

# Output files for velocities in x,y directions (u,v)
u,v = velocity.inversion_3D(x,y,file_velocity_in,Inputs,glacier)

#########################################################################
# Import mesh boundary, calculate flow parameter at mesh nodes, and use #
# SIA approximation to get basal sliding speed for the inflow boundary  #
#########################################################################

# Get mesh nodes
nodes_file=DIRM+"/elmer/mesh.nodes"
nodes=np.loadtxt(nodes_file)

# Get modeled temperatures from Kristin's work
kristin_file=os.path.join(os.getenv("HOME"),"Data/Climate/IceTemperature/Helheim/helheim_TA.xyz")
tempdata=np.genfromtxt(kristin_file,delimiter=',')
tempdata=np.delete(tempdata,(0),axis=0)

# Normalize height to 1 in Kristin's temperatures
xdata,xinds=np.unique(tempdata[:,0],return_index=True)
ydata,yinds=np.unique(tempdata[:,1],return_index=True)
inds=np.union1d(xinds,yinds)
del xdata, ydata
xdata=tempdata[inds,0]
ydata=tempdata[inds,1]
tempdata_normalized=tempdata
for i in range(0,len(xdata)):
  colinds=np.where(tempdata[:,0]==xdata[i])
  colinds=np.array(colinds)
  surf=np.max(tempdata[colinds,2])
  bed=np.min(tempdata[colinds,2])
  H=surf-bed
  tempdata_normalized[colinds,2]=(tempdata[colinds,2]-bed)/H
del H,surf,bed,xdata,ydata,xinds,yinds,inds,colinds

# Normalize height to 1 for nodes
junk,xinds=np.unique(nodes[:,2],return_index=True)
junk,yinds=np.unique(nodes[:,3],return_index=True)
inds=np.unique(np.hstack([xinds,yinds]))
xnodes=nodes[inds,2]
ynodes=nodes[inds,3]
surf=np.zeros_like(xnodes)
bed=np.zeros_like(ynodes)
nodes_normalized=np.zeros_like(nodes[:,4])
height=np.zeros_like(xnodes)
for i in range(0,len(inds)):
  xcolinds=np.array(np.where(nodes[:,2]==xnodes[i]))
  ycolinds=np.array(np.where(nodes[:,3]==ynodes[i]))
  colinds=[]
  if len(xcolinds[0]) >= len(ycolinds[0]):
    for j in range(0,len(xcolinds[0])):
      if xcolinds[0,j] in ycolinds[0,:]:
        colinds.append(xcolinds[0,j])
  else:
    for j in range(0,len(ycolinds[0])):
      if ycolinds[0,j] in xcolinds[0,:]:
        colinds.append(ycolinds[0,j])      
  surf[i]=np.max(nodes[colinds,4]) # Surface elevation
  bed[i]=np.min(nodes[colinds,4]) # Bed elevation
  height[i]=surf[i]-bed[i] #height
  nodes_normalized[colinds]=(nodes[colinds,4]-bed[i])/height[i] #normalized value
	
# Now interpolate Kristin's temperatures to the nodes
Temps_lin = griddata(tempdata_normalized[:,0:3],tempdata_normalized[:,3],np.column_stack([nodes[:,2:4], nodes_normalized]),method='linear')
Temps_near = griddata(tempdata_normalized[:,0:3],tempdata_normalized[:,3],np.column_stack([nodes[:,2:4], nodes_normalized]),method='nearest')

nans=np.isnan(Temps_lin)
Temps_lin[nans]=Temps_near[nans]
Anodes=flowparameter.arrhenius(Temps_lin)

# Write out flow law parameter at each node
fid = open(Inputs+"flowparameters.dat","w")
fid.write('{}\n'.format(len(Temps_lin)))
for i in range(0,len(Temps_lin)):
  fid.write('{0} {1} {2} {3} {4}\n'.format(int(nodes[i,0]),nodes[i,2],nodes[i,3],nodes[i,4],Anodes[i]))
fid.close() 
del nans, kristin_file,Temps_near,fid,Temps_lin   

###############################################################
# Print out half of driving stress for initial guess for Beta #
###############################################################


#################################################################
# Calculate basal sliding speed using SIA for inflow boundaries #
#################################################################

frac = 0.5
ub,vb,beta = elmer_inversion.guess_beta(x,y,zsur,zbed,u,v,frac)

# Write out basal velocities and initial guess for beta
fidub = open(Inputs+"ubdem.xy","w")
fidvb = open(Inputs+"vbdem.xy","w")
fidbeta = open(Inputs+"beta0.xy","w")
fidvb.write('{}\n{}\n'.format(len(x),len(y)))
fidub.write('{}\n{}\n'.format(len(x),len(y)))
fidbeta.write('{}\n{}\n'.format(len(x),len(y)))
for i in range(0,len(x)):
  for j in range(0,len(y)):
    fidub.write('{0} {1} {2}\n'.format(x[i],y[j],ub[j,i]))
    fidvb.write('{0} {1} {2}\n'.format(x[i],y[j],vb[j,i]))
    fidbeta.write('{0} {1} {2}\n'.format(x[i],y[j],beta[j,i]))
fidub.close()
fidvb.close()
fidbeta.close()

##############################################################
# Print out bedrock and surface topographies for elmersolver #
##############################################################

fid = open(Inputs+"mesh_bed.dat","w")
fid.write('{}\n'.format(len(nodes)))
for i in range(0,len(nodes[:,1])):
  fid.write("{} {} {} {} \n".format(nodes[i,0],nodes[i,2],nodes[i,3],nodes[i,4]))
fid.close()

# Makes a 3D mesh of Helheim Glacier
# 
# LMK, UW, 06/01/2014
# Last updated 06/14/2014

import os
import shutil
import sys
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Helheim/Tools"))
import elmer_mesh as mesh
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
lc1=1000 # for entire mesh
lc2=250 # for channels close to the grounding line
lc3=100 # for grounding line
lc4=500 # for regions surrounding channels
#lc1=5000 # for entire mesh
#lc2=1000 # for channels close to the grounding line
#lc3=500 # for grounding line
#lc4=2000 # for regions surrounding channels
levels=12 #levels of extrusion
partitions="4" # Number of partitions

# Bed and surface
file_bed = 'morlighem' 
file_surf = 'gimp'

# Velocity profile for inversion
file_velocity_in = os.path.join(os.getenv("HOME"),"Data/Velocity/TSX/Helheim/track-27794")

###########
# Outputs #
###########

#Set output name for gmsh file
file_2d=os.path.join(DIRM+"Planar")
file_3d=os.path.join(DIRM+"Elmer")

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

#############
# Make mesh #
#############

# Set up bedrock and surface topographies for Extrude mesh
if file_bed == 'morlighem':
  x,y,bed = helheim_bed.morlighem_grid(np.min(exterior[0])-5e3,np.max(exterior[0])+5e3,np.min(exterior[1])-5e3,np.max(exterior[1])+5e3,'geoid')
fid = open(Inputs+"bed.xyz","w")
for i in range(0,len(x)):
  for j in range(0,len(y)):
    fid.write('{} {} {}\n'.format(x[i],y[j],bed[j,i]))
fid.close()

if file_surf == 'gimp':
  x,y,surf = helheim_elevation.gimp_grid(np.min(exterior[0])-5e3,np.max(exterior[0])+5e3,np.min(exterior[1])-5e3,np.max(exterior[1])+5e3,'geoid')
fid = open(Inputs+"surf.xyz","w")
for i in range(0,len(x)):
  for j in range(0,len(y)):
    fid.write('{} {} {}\n'.format(x[i],y[j],surf[j,i]))
fid.close()

del x,y,bed,surf
  
# Gmsh .geo file
mesh.xy_to_gmsh_3d(exterior,holes,refine,file_2d+".geo",lc1,lc2,lc3,lc4)

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
helheim_velocity.inversion_3D(file_velocity_in,Inputs)

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

#################################################################
# Calculate basal sliding speed using SIA for inflow boundaries #
#################################################################

## Constants
yearinsec=365.25*24*60*60
g=9.81
rho=917
slope=0.01 # Average slope, slopes from GIMP are so unreliable that this is the best approach for now

# Get mean flow parameter
Amean=np.zeros_like(xnodes)
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
  Amean[i]=np.mean(Anodes[colinds]) 
    
## Get surface velocities
udem=np.loadtxt(Inputs+"UDEM.xy",skiprows=2)
vdem=np.loadtxt(Inputs+"VDEM.xy",skiprows=2)

## Interpolate surface velocities to nodes    
node_u=np.zeros_like(xnodes)
node_v=np.zeros_like(xnodes)
node_surfmag=np.zeros_like(xnodes)
for i in range(0,len(xnodes)):
  ind=((xnodes[i]-udem[:,0])**2+(ynodes[i]-udem[:,1])**2).argmin()
  node_u[i]=udem[ind,2]
  node_v[i]=vdem[ind,2]
  node_surfmag[i]=np.sqrt(node_u[i]**2+node_v[i]**2)

# Calculate basal speed
node_bedmag=np.zeros_like(xnodes)
node_bedmag=node_surfmag-yearinsec*(Amean*(rho*g*slope)**3*(height**4)) 
ind=np.where(node_bedmag < 20)
node_bedmag[ind]=20

# Now put the basal speed in the same direction as the surface speed
node_ub=np.zeros_like(xnodes)
node_vb=np.zeros_like(xnodes)
node_ub=np.sqrt(node_bedmag**2/(1+(node_v**2/node_u**2)))*np.sign(node_u)
node_vb=np.sqrt(node_bedmag**2/(1+(node_u**2/node_v**2)))*np.sign(node_v)

# Write out values  
fid=open(Inputs+"ubdem.xy","w")
fid.write('{}\n'.format(len(xnodes)))
for i in range(0,len(xnodes)):
  fid.write('{0} {1} {2} {3} {4} {5}\n'.format(xnodes[i],ynodes[i],surf[i],bed[i],node_u[i],node_ub[i]))
fid.close()

fid=open(Inputs+"vbdem.xy","w")
fid.write('{}\n'.format(len(xnodes)))
for i in range(0,len(xnodes)):
  fid.write('{0} {1} {2} {3} {4} {5}\n'.format(xnodes[i],ynodes[i],surf[i],bed[i],node_v[i],node_vb[i]))
fid.close()
del xnodes,ynodes,surf,bed,node_v,node_vb,node_u,Amean,height,fid,Anodes,rho,g,yearinsec

##############################################################
# Print out bedrock and surface topographies for elmersolver #
##############################################################

fid = open(Inputs+"mesh_bed.dat","w")
fid.write('{}\n'.format(len(nodes)))
for i in range(0,len(nodes[:,1])):
  fid.write("{} {} {} {} \n".format(nodes[i,0],nodes[i,2],nodes[i,3],nodes[i,4]))
fid.close()

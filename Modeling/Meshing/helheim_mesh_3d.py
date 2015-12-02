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
import velocity, bed, elevation, fracyear, glacier_extent
from subprocess import call
from scipy.interpolate import *
import numpy as np
import flowparameter
import argparse

##########
# Inputs #
##########

# Get inputs to file
parser = argparse.ArgumentParser()
parser.add_argument("-mesh", dest="meshname", required = True,
        help = "Name of mesh.")
parser.add_argument("-d", dest="date", required = True,
            help = "Date for mesh.")
parser.add_argument("-n", dest="npartitions", required = True,
            help = "Number of partitions.")
parser.add_argument("-bname", dest="bedname", required = False,default='smith',
            help = "Name of bed file (smith,morlighem,cresis).")
parser.add_argument("-bmodel", dest="bedmodel", required = False,default='aniso',
            help = "Type of bed (aniso,iso).")
parser.add_argument("-bsmooth", dest="bedsmooth", type=int,required = False,\
			default=4,help = "Smoothness of bed (1-8).")
parser.add_argument("-lc", dest="lc", type=int,required = False,\
			default=[1000,1000,3000,5000],\
			help = "Four numbers that define the mesh resolution for grounding-line (500 m),channels (1000 m),regions near channels (2000 m), and entire mesh (5000 m).")
parser.add_argument("-extrude", dest="extrude", type=int,required = False,\
			default=5,\
			help = "Number of extrusion levels.")

# Get arguments
args, _ = parser.parse_known_args(sys.argv)

date = args.date
partitions = args.npartitions
bedname = args.bedname
bedmodel = args.bedmodel
bedsmoothing = args.bedsmooth
MESHNAME = args.meshname

# Mesh refinement
lc3,lc2,lc4,lc1 = args.lc
levels=args.extrude #levels of extrusion

del args, parser

# File names
glacier = 'Helheim'

# Directories
DIRS = os.path.join(os.getenv("CODE_HOME"),"BigThreeGlaciers/Modeling/SolverFiles/3D")
DIRM = os.path.join(os.getenv("MODEL_HOME"),glacier+"/Meshes/3D/"+MESHNAME+"/")
DIRR = os.path.join(os.getenv("MODEL_HOME"),glacier+"/Results/3D/")
DIRX = os.path.join(os.getenv("DATA_HOME"),"ShapeFiles/Glaciers/3D/"+glacier+"/")
Inputs = os.path.join(DIRM+"Inputs/")

# Make mesh directories
if not(os.path.isdir(DIRM)):
  os.makedirs(DIRM)
  os.makedirs(DIRM+"/Inputs")

# Densities for finding floating ice
rho_i = 917.0
rho_sw = 1020.0

# Time
time = fracyear.date_to_fracyear(int(date[0:4]),int(date[4:6]),int(date[6:8]))

#################
# Mesh Geometry #
#################

# Mesh exterior
exterior = glacier_extent.load(glacier,time)
np.savetxt(Inputs+"mesh_extent.dat",exterior[:,0:2])

# Mesh holes
hole1 = mesh.shp_to_xy(DIRX+"glacier_hole1")
np.savetxt(Inputs+"mesh_hole1.dat",hole1[:,0:2])

hole2 = mesh.shp_to_xy(DIRX+"glacier_hole2")
np.savetxt(Inputs+"mesh_hole2.dat",hole2[:,0:2])

# All holes
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
x,y,zbed,zsur = mesh.xy_to_gmsh_3d(glacier,date,exterior,holes,refine,DIRM,\
		lc1,lc2,lc3,lc4,bedname,bedmodel,bedsmoothing,rho_i,rho_sw)

# Create .msh file
call(["gmsh","-1","-2",file_2d+".geo", "-o",os.path.join(os.getenv("HOME"),\
		file_2d+".msh")])

# Create elmer mesh
call(["ElmerGrid","14","2",file_2d+".msh","-autoclean"])

# Extrude the mesh with ExtrudeMesh
call(["ExtrudeMesh",file_2d,file_3d,str(levels),"1","1","0","0","0","0",Inputs,"200","2","NaN"])

# Partition mesh for parallel processing
os.chdir(DIRM)
call(["ElmerGrid","2","2","Elmer","dir","-metis",partitions])

# Output as gmsh file so we can look at it
call(["ElmerGrid","2","4","Elmer"])

##########################################
# Print out velocity data for inversions #
##########################################

# Output files for velocities in x,y directions (u,v)
u,v = velocity.inversion_3D(glacier,x,y,time,Inputs)

#########################################################################
# Import mesh boundary, calculate flow parameter at mesh nodes, and use #
# SIA approximation to get basal sliding speed for the inflow boundary  #
#########################################################################

# Get mesh nodes
nodes_file=DIRM+"Elmer/mesh.nodes"
nodes=np.loadtxt(nodes_file)

# Get modeled temperatures from Kristin's work
kristin_file=os.path.join(os.getenv("DATA_HOME"),"Climate/IceTemperature/Helheim/helheim_TA.xyz")
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

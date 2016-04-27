# This module uses shapefiles created in QGIS to create Gmsh mesh files. 
#
# Functions:
# shp_to_xy(in_file) - converts a glacier outline shape file to x, y, and boundary number coordinates
# sort_xy(input) - takes an x,y file and sorts it for future use, so neighbors are near each other 
#    (may only work with Helheim data)
# xy_to_gmsh_3d(exterior,holes,refine,filename,lc1,lc2,lc3,file_bed,file_sur) - takes a glacier outline and 
#    holes in the mesh to create a Gmsh .geo file and surface/bed DEMS for Extrude Mesh
# xy_to_gmsh_box(flowline,filename,file_bed,file_surf,lc,lc_d) - creates a box mesh 
#    ([0 1],[0 1]) which can then be extruded using MshGlacier
# xy_to_gmsh_2d(data,filename,lc1) - creates a flowline with the proper surface and bedrock
#    topographies, but doesn't work as well as MshGlacier so I don't presently use this 
#    function (maybe should delete soon)

import os
import shapefile
import numpy as np
import scipy.signal as signal
import math
from scipy import interpolate
import bedlib, floatlib, zslib, coordlib
from shapely.geometry import Polygon,Point
import scipy.interpolate

def shp_to_xy(in_file):

  '''
  exterior = shp_to_xy(in_file)
  
  Converts shapefile "in_file" to a sorted array of x,y points.
  
  inputs:
  in_file: shapefile name
  
  Outputs:
  exterior: three column array with x, y, and boundary number for array (useful for meshing)
  '''

  sf = shapefile.Reader(in_file)
  shapes=sf.shapes()
  try:
    rec = sf.records()
  except:
    pass  
  n = len(shapes)
  x1 = np.zeros(n)
  y1 = np.zeros(n)
  bound1 = np.zeros(n)
  for i in range(0,n):
    x1[i]=float(shapes[i].points[0][0])
    y1[i]=float(shapes[i].points[0][1])
    try:
      bound1[i]=int(rec[i][0])
    except:
      bound1[i]='NaN'
  
  # Sort points
  ind = np.argmin(x1)
  x2=[] 
  x2.append(x1[ind])
  y2=[]
  y2.append(y1[ind])
  bound2=[]
  bound2.append(bound1[ind])
  indices=range(0,len(x1))
  indices.remove(ind)
  n=0
  while len(indices) > 1:
    min_d=1000000000
    d=[]
    for i in indices:
      d=(math.sqrt((x2[n]-x1[i])**2+(y2[n]-y1[i])**2))
      if (d < min_d):
        min_d=d
        ind=i
    x2.append(x1[ind])
    y2.append(y1[ind])
    bound2.append(bound1[ind])
    indices.remove(ind)
    n=n+1
  
  exterior = np.column_stack([x2,y2,bound2])

  return exterior
  
def xy_to_gmsh_3d(glacier,date,exterior,holes,refine,DIRM,lc1,lc2,lc3,lc4,\
		bedname='smith',bedmodel='aniso',bedsmoothing=4,rho_i=917.0,rho_sw=1020.0):
  
  '''
  x,y,zbed,zsur = xy_to_gmsh_3d(glacier,date,exterior,holes,refine,DIRM,lc1,lc2,lc3,lc4,bedname,bedmodel,bedsmoothing)
  
  inputs:
  exterior: list of x,y coordinates for glacier mesh from "shp_to_xy"
  holes: list of holes with x,y coordinates for glacier mesh from "shp_to_xy"
  refine: list of x,y coordinates for refinement, with an index number that states what level of refinment
  DIRM: mesh directory
  lc1, lc2, lc3, lc4: mesh resolution near refinement locations
  file_bed: string for what bed we should use (options: morlighem)
  file_sur: string for what surface we should use (options: gimp)
  
  Output files: surface DEM, bed DEM, Planar.geo
  '''
  
  # Open file "Planar.geo" in the mesh directory
  filename = os.path.join(DIRM+"mesh2d.geo")
  fid = open(filename,'w')
  fid.write("lc1 = %f; \n" % lc1)
  fid.write("lc2 = %f; \n" % lc2)
  fid.write("lc3 = %f; \n" % lc3)
  fid.write("lc4 = %f; \n" % lc4)
  
  n=1 #Set index for counting
  ind_1=[] #Keep track of line indices for the terminus so we can set physical line
  ind_2=[] #Glacier terminus
  ind_3=[] #Keep track of indices for glacier margin near terminus
  ind_4=[] #Keep track of inflow boundary
  
  ############
  # Exterior #
  ############
  
  # Save points for mesh. We check the boundary number for each point and save it with a 
  # different mesh element size (lc*) depending on the boundary number.
  X = len(exterior[:,0])
  d = lc1+1000
  last = 0
  indices = []
  for i in range(0,X):
    if i != 0:
      d = math.sqrt((exterior[last,0]-exterior[i,0])**2+(exterior[last,1]-exterior[i,1])**2)
    if (exterior[i,2] == 2) and (d > lc2/2):
      fid.write('Point({}) = {{{}, {}, 0, lc3}}; \n'.format(n,exterior[i,0],exterior[i,1]))
      last = i
      indices.append(i)
      n = n+1
    elif (exterior[i,2] == 3) and (d > lc3/2):
      fid.write('Point({}) = {{{}, {}, 0, lc3}}; \n'.format(n,exterior[i,0],exterior[i,1]))
      last = i
      indices.append(i)
      n = n+1
    elif (exterior[i,2] == 5) and (d > lc2/2):
      fid.write('Point({}) = {{{}, {}, 0, lc2}}; \n'.format(n,exterior[i,0],exterior[i,1]))
      last=i
      indices.append(i)
      n=n+1
    elif d > lc1/2:
      fid.write('Point({}) = {{{}, {}, 0, lc1}}; \n'.format(n,exterior[i,0],exterior[i,1]))
      last = i
      indices.append(i)
      n = n+1
  pt_end = n-1
  
  # Create lines for gmsh. Keep track of boundary number so that we know which 
  # line segments belong to which boundary number (and thereby boundary condition).    
  for i in range(1,pt_end):
    fid.write('Line({}) = {{{}, {}}}; \n'.format(n,i,i+1))
    # Glacier front
    if (exterior[indices[i]-1,2] == 2) and (exterior[indices[i],2] == 2):
	  ind_2.append(n)
	# Inflow boundary
    elif (exterior[indices[i]-1,2]==4) and (exterior[indices[i],2] == 4):
      ind_4.append(n)
    #Margins
    else:
      ind_1.append(n)
    n = n+1
  
  # Save last line (which completes the line loop and connects the last point to 
  # the first point)
  fid.write('Line({}) = {{{}, {}}}; \n'.format(n,pt_end,1))
  if (exterior[indices[i]-1,2] == 1 or exterior[indices[i],2] == 1):
    ind_1.append(n)
  elif (exterior[indices[i]-1,2] == 4 and exterior[indices[i],2] == 4):
    ind_4.append(n)
  else:
    ind_2.append(n)
  n = n+1
  
  # Save line loop
  fid.write('Line Loop({}) = {{{}'.format(n,pt_end+1))
  for i in range(1,pt_end):
    fid.write(',{}'.format(pt_end+i+1))
  fid.write('}; \n')   
  ext_ind=n # Save index of line loop for later surface calculations
  n = n+1
  
  #########
  # Holes #
  #########
  
  Z = 1
  int_ind = []
  last = 0
  d = lc1+1000
  for j in range(0,len(holes)):
    pt1 = n
    Y = len(holes[j]['xy'][:,0])
    for i in range(0,Y):
      if i != 0:
        d=math.sqrt((holes[j]['xy'][i,0]-holes[j]['xy'][last,0])**2+(holes[j]['xy'][i,1]-holes[j]['xy'][last,1])**2)
      if d > lc2:
        fid.write('Point({}) = {{{}, {}, 0, lc2}}; \n'.format(n,holes[j]['xy'][i,0],holes[j]['xy'][i,1]))
        n = n+1
        last = i
    pt2 = n-1
     
    line1 = n    
    for i in range(pt1,pt2):
      fid.write('Line({}) = {{{}, {}}}; \n'.format(n,i,i+1))
      ind_1.append(n)
      n = n+1   
      
    fid.write('Line({}) = {{{}, {}}}; \n'.format(n,pt2,pt1))
    line2 = n+1
    ind_1.append(n)
    n = n+1
    
    fid.write('Line Loop({}) = {{{}'.format(n,line1))
    int_ind.append(n)
    for i in range(line1+1,line2):
      fid.write(',{}'.format(i))
    fid.write('}; \n')
    n = n+1    
    
  ############################  
  # Add refinement locations #
  ############################
  
  # Add refinement points
  R = len(refine[:,0])
  start = n
  nrefine = 0 
  polyglacier = Polygon(exterior[:,0:2])
  for i in range(0,R):
    polypt = Point(refine[i,0:2])
    if polyglacier.contains(polypt):
      if refine[i,2] == 2:
        fid.write('Point({}) = {{{}, {}, 0, lc3}}; \n'.format(n,refine[i,0],refine[i,1]))
      elif refine[i,2] == 1:
        fid.write('Point({}) = {{{}, {}, 0, lc2}}; \n'.format(n,refine[i,0],refine[i,1]))
      else:
        fid.write('Point({}) = {{{}, {}, 0, lc4}}; \n'.format(n,refine[i,0],refine[i,1]))
      n = n+1
      nrefine = nrefine+1
  
  # Save plane surface
  fid.write('Plane Surface({}) = {{'.format(n))
  for k in range(0,len(int_ind)): 
    fid.write('{},'.format(int_ind[k]))
  fid.write('{}}}; \n'.format(ext_ind))
  surf = n # save surface index for surface
  n = n+1
  
  # Call refinement points in gmsh file
  for j in range(0,nrefine):
    fid.write('Point{{{}}} In Surface{{{}}}; \n'.format(start+j,surf))
  
  # Save physical surfaces and physical lines for calling boundary conditions
  fid.write('Physical Surface({}) = {{{}}}; \n'.format(1,surf))
  fid.write('Physical Line(2) = {{{}'.format(ind_1[0]))
  for m in range(1,len(ind_1)): 
    fid.write(',{}'.format(ind_1[m]))
  fid.write('}; \n')
  fid.write('Physical Line(1) = {{{}'.format(ind_2[0]))
  for m in range(1,len(ind_2)): 
    fid.write(',{}'.format(ind_2[m]))
  fid.write('}; \n') 
  fid.write('Physical Line(3) = {{{}'.format(ind_4[0]))
  for m in range(1,len(ind_4)): 
    fid.write(',{}'.format(ind_4[m]))
  fid.write('}; \n') 
  
  # Close .geo file
  fid.close()
   
  ############################################################ 
  # Set up bedrock and surface topographies for Extrude mesh #
  ############################################################
  
  # Set up grid extent
  xmin = np.min(exterior[:,0])-2.0e3
  xmax = np.max(exterior[:,0])+2.0e3
  ymin = np.min(exterior[:,1])-2.0e3
  ymax = np.max(exterior[:,1])+2.0e3
  
  # Load bed DEM
  if bedname == 'morlighem':
    xbed_grid,ybed_grid,zbed = bedlib.morlighem_grid(xmin,xmax,ymin,ymax,verticaldatum='geoid')
    xbed_grid,ybed_grid = np.meshgrid(xbed,ybed)
  elif bedname == 'smith':
    # irregular triangular grid
    xbed_grid,ybed_grid,zbed_grid = bedlib.smith_grid(glacier,\
    			model=bedmodel,smoothing=bedsmoothing,verticaldatum='geoid')

  # Load surface DEM
  xsur,ysur,zsur_grid = zslib.dem_continuous(glacier,date,\
  			verticaldatum='geoid',fillin=True,blur=False)
  xsur_grid,ysur_grid = np.meshgrid(xsur,ysur)

  # Flatten surface coordinates and interpolate bed elevations to surface coordinates
  # so that we have grids.
  x_flattened = xsur_grid.flatten()
  y_flattened = ysur_grid.flatten()
  zsur_flattened = zsur_grid.flatten()
  zbed_flattened = scipy.interpolate.griddata((ybed_grid.flatten(),xbed_grid.flatten()),zbed_grid.flatten(),\
  					(y_flattened,x_flattened),method='linear')
  
  # Now interpolate surface elevations to bed coordinates, we only do this for finding the 
  # ice bottom.
  f = scipy.interpolate.RegularGridInterpolator((ysur,xsur),zsur_grid,method='linear')
  x_interped = xbed_grid.flatten()
  y_interped = ybed_grid.flatten()
  zsur_interped = f((ybed_grid.flatten(),xbed_grid.flatten()))
  zbed_interped = zbed_grid.flatten()

  zbot_interped = floatlib.icebottom(zbed_interped,zsur_interped,rho_i=rho_i,rho_sw=rho_sw)

  zbed_grid = np.reshape(zbed_flattened,(len(ysur),len(xsur)))
  ind = np.where((zsur_grid-zbed_grid) < 10.)
  zbed_grid[ind] = zsur_grid[ind]-10.
  
  # Print out surface and bed
  fids = open(DIRM+"/inputs/zsdem.xy","w")
  fidb = open(DIRM+"/inputs/zbdem.xy","w")
  fids.write('{}\n{}\n'.format(len(xsur),len(ysur)))
  fidb.write('{}\n{}\n'.format(len(xsur),len(ysur)))
  for i in range(0,len(xsur)):
    for j in range(0,len(ysur)):
      fids.write('{0} {1} {2}\n'.format(xsur[i],ysur[j],zsur_grid[j,i]))
      fidb.write('{0} {1} {2}\n'.format(xsur[i],ysur[j],zbed_grid[j,i]))
  fids.close()
  fidb.close()
   
  # Print out bed
  # np.savetxt(DIRM+"/inputs/roughbed.xyz",np.column_stack((xbed_grid.flatten(),ybed_grid.flatten(),zbed_interped)),fmt='%.6f')
  

  
  return xsur,ysur,zbed_grid,zsur_grid

##########################################################################################
   
def xy_to_gmsh_box(x,y,dists,zb,zs,terminus,glacier,DIRM,filename,lc,lc_d,layers,filt_len='none',rho_i=917.0,rho_sw=1020.0):
  
  '''
  flowline = xy_to_gmsh_box(glacier, flowline, DIRM, filename, date, lc,
  			lc_d, layers, filt_len='none')
  
  Make a box mesh that can be extruded with the Elmer/Ice tool "MshGlacier". 
  The function outputs two files (DIRM+filename+"_surf.dat",DIRM+filename+"_bed.dat") 
  that include the surface elevations and ice bottom elevations, respectively. 
  These files are used by "MshGlacier."

  inputs:
  x,y,dists,zb: output from glacier_flowline
  glacier: glacier name
  date: date for surface DEM
  DIRM: mesh directory
  filename: filename for .geo file
  lc: array for mesh sizes
  lc_d: distance between different mesh sizes
  layers: how many vertical layers for mesh
  filt_len: filter length (in meters) for the surface elevations
  
  Outputs:
  flowline: four column array with distance along flowline, x coordinate, 
  		y coordinate, bed elevation, and surface elevation
  some files in the DIRM directory
  '''
  
  # Change directory to mesh directory
  os.chdir(DIRM)
  
  # Set up variables
  R = len(x) # length of flowline
  d = dists # distance along flowline
  
  # Load terminus file so we can figure out where to cutoff the mesh 
  D = round(terminus-d[0])
  
  # Now start to write out the .geo file
  fid = open(DIRM+filename+".geo",'w')
  fid.write('lc1 = {}; \n'.format(lc[0]/D))
  fid.write('lc2 = {}; \n'.format(lc[0]/D))
  fid.write('lc3 = {}; \n'.format(lc[2]/D))
  
  # Points
  fid.write('Point(1) = {{{}, {}, 0, lc3}}; \n'.format(0,0))
  fid.write('Point(2) = {{{}, {}, 0, lc2}}; \n'.format(1-lc_d[2]/D,0))
  fid.write('Point(3) = {{{}, {}, 0, lc1}}; \n'.format(1-lc_d[1]/D,0))
  fid.write('Point(4) = {{{}, {}, 0, lc1}}; \n'.format(1,0))
  
  # Lines
  fid.write('Line(5) = {1, 2}; \n')
  fid.write('Line(6) = {2, 3}; \n')
  fid.write('Line(7) = {3, 4}; \n')
  
  for i in [5,6,7]:
    fid.write('Extrude {{0, 1, 0}} {{Line{{{0}}}; Layers{{{1}}}; Recombine;}} \n'.format(i,layers))

  
  fid.write('Physical Surface(1) = {11,15,19}; \n');
  
  # Now write out physical lines
  # Glacier bed
  fid.write('Physical Line(1) = {5,6,7}; \n')
  # Glacier terminus
  fid.write('Physical Line(2) = {18}; \n')
  # Glacier surface
  fid.write('Physical Line(3) = {16,12,8}; \n')
  # Upper glacier
  fid.write('Physical Line(4) = {9}; \n')
  fid.close()
    
  # Export the ice surface elevation    
  fid = open(DIRM+filename+"_surf.dat",'w')
  fid.write('{0:.6f} {1:.6f}\n'.format(d[0]-100,zs[0]))
  ind = np.where(d < terminus)[0]
  for i in ind:
    fid.write('{0:.6f} {1:.6f}\n'.format(d[i],zs[i]))
  fid.write('{0:.6f} {1:.6f}'.format(d[i]+100,zs[i]))
  fid.close()
  
  # Calculate and export the ice bottom elevation
  icebottom = floatlib.icebottom(zb,zs,rho_i=rho_i,rho_sw=rho_sw)
  fid = open(DIRM+filename+"_bed.dat",'w')
  fid.write('{0:.6f} {1:.6f}\n'.format(d[0]-100,zb[0]))
  for i in ind:
    fid.write('{0:.6f} {1:.6f}\n'.format(d[i],zb[i]))
  fid.write('{0:.6f} {1:.6f}'.format(d[i]+100,zb[i]))
  fid.close()

  # Now update the "mesh_update.dat" file required for MshGlacier
  num_surf = num_bed = len(ind)

  input_file=os.path.join(DIRM+"../"+"mesh_input_source.dat")
  output_file=os.path.join(DIRM+"mesh_input.dat")
  with open(input_file, "r" ) as source:
    with open(output_file, "w" ) as target:
      data = source.read()
      changed1 = data.replace('BED_PTS',str(num_bed+2))
      changed2 = changed1.replace('SUR_PTS',str(num_surf+2))
      changed3 = changed2.replace('MESH_TERMINUS',str(terminus))
      changed4 = changed3.replace('MESH_START',str(min(d)))
      changed5 = changed4.replace('Elmer',filename)
      target.write( changed4 )
  
  # Save flowline for future interpolation
  flowline=np.column_stack([d,x,y,zb,zs])
  
  return flowline

##############
# xy_gmsh_2d #
##############
# Currently don't use this function, will probably delete soon because
# it's easier to use MshGlacier 
# LMK, UW, 06/01/2014  
def xy_to_gmsh_2d(data,filename,lc1):
  
  #Get rid of the points that have no ice
  R=len(data[0])
  n=[]
  for i in range(0,R):
    if (data[3][i] > 20):
      n.append(i) 
  flowline=[[row[i] for row in data] for i in n] 
      
  #Sort points
  R=len(flowline[:])
  indices = range(0,R)
  pos=0
  ind=[]
  for j in range(0,R):   
    if float(flowline[j][1]) < pos:
      ind=j
      pos=float(flowline[j][1])
  ind2=ind
  x=[]
  y=[]
  z_bed=[]
  z_sur=[]
  n=0
  while len(indices) > 1:
    min_d=1000000000
    d=[]
    x.append(flowline[ind][0])
    y.append(flowline[ind][1])
    z_bedlib.append(flowline[ind][2])
    z_sur.append(flowline[ind][3]) 
    indices.remove(ind)
    for i in indices:
      d=(math.sqrt((x[n]-flowline[i][0])**2+(y[n]-flowline[i][1])**2))
      if (d < min_d):
        min_d=d
        ind=i
    n=n+1
    #print n
  
  d=[]
  d.append(0) 
  for i in range(0,R-2):
  	d_inc=math.sqrt((x[i+1]-x[i])**2+(y[i+1]-y[i])**2)
  	d.append(d[i]+d_inc)
  #Create bottom surface and top surface
  s_bed=[d,z_bed]
  s_sur=[d,z_sur]
  s_sur=[list(reversed(row)) for row in s_sur]
  
  #Now we can actually do some meshing
  fid = open(filename+".geo",'w')
  fid.write("lc1 = %f; \n" % lc1)
  
  n=1 #Set index for counting
  ind_1=[] #Keep track of line indices for the bed
  ind_2=[] #Keep track of line indices for the upper glacier
  ind_3=[] #Keep track of line indices for the surface
  ind_4=[] #Keep track of line indices for glacier terminus
  
  #Points
  f
  for i in range(0,X):
    fid.write('Point({}) = {{{}, {}, 0, lc1}}; \n'.format(n,s_sur[0][i],s_sur[1][i]))
    n=n+1
  
  #Glacier bed    
  fid.write('Spline({}) = {{{}'.format(n,1))
  for i in range(2,X+1):
    fid.write(',{}'.format(i))
  fid.write('}; \n')
  ind_1=n
  n=n+1
  #Upper glacier
  fid.write('Line({}) = {{{}, {}}}; \n'.format(n,X,X+1))
  ind_2=n
  n=n+1
  #Glacier surface
  fid.write('Spline({}) = {{{}'.format(n,X+1))
  for i in range(X+2,2*X+1):
    fid.write(',{}'.format(i))
  fid.write('}; \n')
  ind_3=n
  n=n+1
  #Glacier terminus
  fid.write('Line({}) = {{{}, {}}}; \n'.format(n,2*X,1))
  ind_4=n
  n=n+1

  #Line loop  
  fid.write('Line Loop({}) = {{{}'.format(n,2*X+1))
  for i in range(2*X+2,ind_4+1):
     fid.write(',{}'.format(i))
  fid.write('}; \n')
  n=n+1
  
  #Plane surface
  fid.write('Plane Surface({}) = {{{}}}; \n'.format(n,n-1))
  
  #Physical surface
  fid.write('Physical Surface(1) = {{{}}}; \n'.format(n))
  
  #Physical lines
  #Bed
  fid.write('Physical Line(1) = {{{}}}; \n'.format(ind_1))
  #Glacier terminus
  fid.write('Physical Line(2) = {{{}}}; \n'.format(ind_4))
  #Glacier surface
  fid.write('Physical Line(3) = {{{}}}; \n'.format(ind_3))
  #Upper glacier
  fid.write('Physical Line(4) = {{{}}}; \n'.format(ind_2))
   
  fid.close()
  
  fid2 = open(filename+"_coordlib.txt",'w')
  for i in range(0,R-1):
    fid2.write('{} {} {} \n'.format(x[i],y[i],d[i]))
  fid2.close()  
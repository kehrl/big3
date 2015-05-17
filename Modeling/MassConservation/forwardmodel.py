import meshpy.gmsh_reader
import numpy as np
import os

# Depth-averaged mass continuity with artificial diffusion
# div(Hv)+div(kgrad(H)) = bdot

# Let's use idealized measurements and geometry to start. 
DIRM = os.path.join(os.getenv("HOME"),"Models/Idealized/")

# Mesh directory
#MESHNAME = "MorlighemNew"
#DIRM = os.path.join(os.getenv("HOME"),"Models/Helheim/Meshes/3D/"+MESHNAME+"/")

# Load mesh and get information about elements
mesh = meshpy.gmsh_reader.GmshMeshReceiverNumPy()
meshpy.gmsh_reader.read_gmsh(mesh,DIRM+"Planar.msh",2)
elements = mesh.elements
bound = np.array(mesh.element_markers)

# Find all triangular elements and their nodes
trielem=[]
nodelem=[]
for element in elements:
  if len(element) > 2:
    trielem.append(element)
nodes=mesh.points
x = list(s[0] for s in nodes) 
y = list(s[1] for s in nodes) 
E = len(trielem)
N = len(nodes)
del nodes, element

# Area of elements
A = np.zeros(E)
for i in range(0,E):
  A[i]=(1.0/2)*np.linalg.det([[1, 1, 1],[x[trielem[i][0]], x[trielem[i][1]], x[trielem[i][2]]],[y[trielem[i][0]],y[trielem[i][1]], y[trielem[i][2]]]])

# Interpolate velocity data to nodes
ux = 100*np.ones(N)
uy = 0*np.ones(N)

# Find coefficients for interpolation functions phi for each element
phi=np.zeros([3,3])
alpha=np.zeros([E,3])
beta=np.zeros([E,3])
for i in range(0,E):
  M=np.array([[x[trielem[i][0]], x[trielem[i][1]], x[trielem[i][2]]],[y[trielem[i][0]],y[trielem[i][1]], y[trielem[i][2]]],[1,1,1]]).T
  phi[:,0]=np.dot(np.linalg.pinv(M),[1,0,0])
  phi[:,1]=np.dot(np.linalg.pinv(M),[0,1,0])
  phi[:,2]=np.dot(np.linalg.pinv(M),[0,0,1])
  for j in range(0,3):
    alpha[i,j]=phi[1,j]
    beta[i,j]=phi[2,j]    
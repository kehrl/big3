import meshpy.gmsh_reader
import numpy as np
import os

# Load mesh
MESHNAME = "MorlighemNew"
DIRM = os.path.join(os.getenv("HOME"),"Models/Helheim/Meshes/Flowline/"+MESHNAME+"/")


mesh = meshpy.gmsh_reader.GmshMeshReceiverNumPy()
meshpy.gmsh_reader.read_gmsh(mesh,"Planar.msh",2)

# Forward problem
# div(Hv) = bdot

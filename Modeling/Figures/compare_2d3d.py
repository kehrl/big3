# This is the beginning of an attempt to compare our 2D and 3D results.

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
import elmer_read

##########
# Inputs #
##########

# 3D
meshname='High_Normal'
runname="beta_1e10"
DIRR_3D=os.path.join(os.getenv("HOME"),"Models/Results/3D/Helheim/"+meshname+"/Inversion/")
bbed_3D=3
bsurf_3D=4


# Flowline
meshname='Worldview'
runname = 'beta_1e10'
DIRR_2D=os.path.join(os.getenv("HOME"),"Models/Results/Flowline/Helheim/"+meshname+"/Inversion/")
bbed_2D=1
bsurf_2D=2

#############
# Load Data #
#############

# Load 3D
try:
  bed_3D
except:
  bed_3D = elmer_read.saveline_boundary(DIRR_3D,runname,bbed_3D)
  surf_3D = elmer_read.saveline_boundary(DIRR_3D,runname,bsurf_3D)

# Load 2D
flowline=np.loadtxt(os.path.join(os.getenv("HOME"),"Models/Meshes/Flowline/Helheim/"+meshname+"/Inputs/flowline.dat"),skiprows=1)
try:
  bed_2D
except:
  bed_2D = elmer_read.saveline_boundary(DIRR_2D,runname,bbed_2D)
  surf_2D = elmer_read.saveline_boundary(DIRR_2D,runname,bsurf_2D)
  
###############
# Comparisons #
###############

bed_3D_flow=elmer_read.grid_to_flowline(bed_3D,flowline[:,1],flowline[:,2])
surf_3D_flow=elmer_read.grid_to_flowline(surf_3D,flowline[:,1],flowline[:,2])

plt.plot(flowline[:,0],bed_3D_flow['taub'],'k.',label='3D')
plt.plot(bed_2D['coord1'],bed_2D['taub'],'r.',label='Flowline')
plt.legend()

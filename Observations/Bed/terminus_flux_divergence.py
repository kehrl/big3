# This code considers thinning near the terminus of Helheim Glacier due to 
# longitudinal/transverse strain. Do we expect the terminus to be grounded or floating?

# LMK, UW, 3/2/15

import os
import sys
import numpy as np
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Helheim/Tools"))
import helheim_velocity, helheim_icefronts, helheim_bed, helheim_elevation
import matplotlib.pyplot as plt
import geotiff
import math

#############
# Load data #
#############

# Load Helheim bed
bed = helheim_bed.cresis('2001')

# Load elevations
elevations = helheim_elevation.worldview_at_pts(bed[:,0],bed[:,0],2,[])

# Load velocities
demtimes,dems = helheim_elevation.worldview_at_pts(bed[:,0],bed[:,1],2,[])
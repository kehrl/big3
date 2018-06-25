import os, sys
from shutil import copy2
import numpy as np
import argparse
import elmerreadlib, geotifflib
import matplotlib.pyplot as plt

def get_arguments():

    parser = argparse.ArgumentParser()
    parser.add_argument("-glacier",dest="glacier",required = True,
        help = "Name of glacier (Kanger or Helheim)")
    parser.add_argument("-DIR", dest="DIR", required = True,
        help = "Directory with inversions.")
    parser.add_argument("-SSA",dest="SSA", required = False, default = False,
        help = "If the model is SSA (default is False).")
    parser.add_argument("-beta",dest="beta_suffix", required = True,
        help = "Suffix for beta.")

    args, _ = parser.parse_known_args(sys.argv)
 
    return args

#################
# Get arguments #
#################

args = get_arguments()

glacier = args.glacier
DIRM = args.DIR
SSA = args.SSA
beta_suffix = args.beta_suffix

DIRR = os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/"+DIRM)
DIRO = os.path.join(os.getenv("MODEL_HOME"),glacier+"/Results/"+DIRM)

if not(os.path.isdir(DIRO)):
    os.mkdir(DIRO)

if SSA:
    modelname = '_1e13_SSA_'
    vname = 'ssavelocity'
else:
    modelname = '_1e12_FS_'
    vname = 'velocity'

DIRs = os.listdir(DIRR)
for dir in DIRs:
    if not(dir.endswith('Lcurve')) and not(dir.endswith('steady')):
        data = elmerreadlib.pvtu_file(DIRR+"/"+dir+"/mesh2d/steady"+modelname+beta_suffix+'0001.pvtu',[vname])
        
        # Mesh boundaries
        extent = np.loadtxt(DIRR+"/"+dir+"/inputs/mesh_extent.dat")
        try:
            hole1 = np.loadtxt(DIRR+"/"+dir+"/inputs/mesh_hole1.dat")
            hole2 = np.loadtxt(DIRR+"/"+dir+"/inputs/mesh_hole2.dat")
            holes=[hole1,hole2]
        except:
            holes = []
      
        bed = elmerreadlib.values_in_layer(data,'bed')
        surf = elmerreadlib.values_in_layer(data,'surf')

        # Get basal velocities
        x,y,bed_mod_ub = elmerreadlib.grid3d(bed,vname+' 1',holes,extent)
        x,y,bed_mod_vb = elmerreadlib.grid3d(bed,vname+' 2',holes,extent)
        x,y,bed_mod_mag = elmerreadlib.grid3d(bed,vname,holes,extent)

        # Get surface velocities
        x,y,surf_mod_us = elmerreadlib.grid3d(surf,vname+' 1',holes,extent)
        x,y,surf_mod_vs = elmerreadlib.grid3d(surf,vname+' 2',holes,extent)
        x,y,surf_mod_mag = elmerreadlib.grid3d(surf,vname,holes,extent)

        geotifflib.write_from_grid(x,y,np.flipud(bed_mod_ub),np.float('nan'),DIRO+"/"+dir[0:11]+"_"+beta_suffix+"_bed_mod_ub.tif")
        geotifflib.write_from_grid(x,y,np.flipud(bed_mod_vb),np.float('nan'),DIRO+"/"+dir[0:11]+"_"+beta_suffix+"_bed_mod_vb.tif")
        geotifflib.write_from_grid(x,y,np.flipud(bed_mod_mag),np.float('nan'),DIRO+"/"+dir[0:11]+"_"+beta_suffix+"_bed_mod_umag.tif")

        geotifflib.write_from_grid(x,y,np.flipud(surf_mod_us),float('nan'),DIRO+"/"+dir[0:11]+"_"+beta_suffix+"_surf_mod_us.tif")
        geotifflib.write_from_grid(x,y,np.flipud(surf_mod_vs),float('nan'),DIRO+"/"+dir[0:11]+"_"+beta_suffix+"_surf_mod_vs.tif")
        geotifflib.write_from_grid(x,y,np.flipud(surf_mod_mag),float('nan'),DIRO+"/"+dir[0:11]+"_"+beta_suffix+"_surf_mod_umag.tif")

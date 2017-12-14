import os, sys
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

  args, _ = parser.parse_known_args(sys.argv)
 
  return args

#################
# Get arguments #
#################

args = get_arguments()

glacier = args.glacier
DIRM = args.DIR

DIRR = os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/"+DIRM)
DIRO = os.path.join(os.getenv("MODEL_HOME"),glacier+"/Results/"+DIRM)

if not(os.path.isdir(DIRO)):
  os.mkdir(DIRO)

n = 0
bbed = 4
bsur = 5
DIRs = os.listdir(DIRR)
for dir in DIRs:
  dirs_inv = os.listdir(DIRR+"/"+dir+"/mesh2d/inversion_adjoint/")
  for dir_inv in dirs_inv:
    if dir_inv.startswith('lambda') and not(dir_inv.endswith('.pdf')):
      if n == 0:
        # Mesh boundaries
        extent = np.loadtxt(DIRR+"/"+dir+"/inputs/mesh_extent.dat")
        try:
          hole1 = np.loadtxt(DIRR+"/"+dir+"/inputs/mesh_hole1.dat")
          hole2 = np.loadtxt(DIRR+"/"+dir+"/inputs/mesh_hole2.dat")
          holes=[hole1,hole2]
        except:
         holes = []

      dirvtu = DIRR+"/"+dir+"/mesh2d/inversion_adjoint/"+dir_inv
      #vtufile = elmerreadlib.pvtu_file(dirvtu+"/adjoint_beta0001.pvtu",['vsurfini','velocity','beta'])
      #surf = elmerreadlib.values_in_layer(vtufile,layer='surface')
      #bed = elmerreadlib.values_in_layer(vtufile,layer='bed')

      # Get values at bed and surface
      bed = elmerreadlib.saveline_boundary(dirvtu+"/","adjoint_beta",bbed,['velocity','beta'])
      surf = elmerreadlib.saveline_boundary(dirvtu+"/","adjoint_beta",bsur,['velocity','vsurfini'])

      # Grid bed values
      x,y,bed_mod_taub = elmerreadlib.grid3d(bed,'taub',holes,extent)
      x,y,bed_mod_beta = elmerreadlib.grid3d(bed,'beta',holes,extent)
      x,y,bed_mod_zb = elmerreadlib.grid3d(bed,'z',holes,extent)
      x,y,bed_mod_ub = elmerreadlib.grid3d(bed,'velocity 1',holes,extent)
      x,y,bed_mod_vb = elmerreadlib.grid3d(bed,'velocity 2',holes,extent)
      
      # Grid surface values
      x,y,surf_mea_zs = elmerreadlib.grid3d(surf,'z',holes,extent)
      x,y,surf_mea_us = elmerreadlib.grid3d(surf,'vsurfini 1',holes,extent)      
      x,y,surf_mea_vs = elmerreadlib.grid3d(surf,'vsurfini 2',holes,extent)
      x,y,surf_mod_us = elmerreadlib.grid3d(surf,'velocity 1',holes,extent)
      x,y,surf_mod_vs = elmerreadlib.grid3d(surf,'velocity 2',holes,extent)      

      geotifflib.write_from_grid(x,y,np.flipud(bed_mod_taub),float('nan'),DIRO+"/"+dir[0:11]+"_bed_mod_taub.tif") 
      geotifflib.write_from_grid(x,y,np.flipud(bed_mod_beta),float('nan'),DIRO+"/"+dir[0:11]+"_bed_mod_beta.tif")
      geotifflib.write_from_grid(x,y,np.flipud(bed_mod_ub),float('nan'),DIRO+"/"+dir[0:11]+"_bed_mod_ub.tif")
      geotifflib.write_from_grid(x,y,np.flipud(bed_mod_vb),float('nan'),DIRO+"/"+dir[0:11]+"_bed_mod_vb.tif")
      geotifflib.write_from_grid(x,y,np.flipud(bed_mod_zb),float('nan'),DIRO+"/"+dir[0:11]+"_bed_mod_zb.tif")

      geotifflib.write_from_grid(x,y,np.flipud(surf_mea_zs),float('nan'),DIRO+"/"+dir[0:11]+"_surf_mea_zs.tif")
      geotifflib.write_from_grid(x,y,np.flipud(surf_mea_us),float('nan'),DIRO+"/"+dir[0:11]+"_surf_mea_us.tif")
      geotifflib.write_from_grid(x,y,np.flipud(surf_mea_vs),float('nan'),DIRO+"/"+dir[0:11]+"_surf_mea_vs.tif")
      geotifflib.write_from_grid(x,y,np.flipud(surf_mod_us),float('nan'),DIRO+"/"+dir[0:11]+"_surf_mod_us.tif")
      geotifflib.write_from_grid(x,y,np.flipud(surf_mod_vs),float('nan'),DIRO+"/"+dir[0:11]+"_surf_mod_vs.tif")

  n=n+1

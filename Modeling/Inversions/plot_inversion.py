# Inputs: Runname regpar1 regpar2 regpar3


import sys
import os
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Helheim/Tools"))
import elmer_read
import matplotlib.pyplot as plt
import matplotlib
import pylab
import numpy as np

args = sys.argv

RES=args[1]
regpar = args[2:]
DIRR=os.path.join(os.getenv("HOME"),"Models/Helheim/Results/3D/"+RES+"/Inversion/")
bbed = 3
bsur = 4

# Mesh boundaries
DIRM = os.path.join(os.getenv("HOME"),"Models/Helheim/Meshes/3D/"+RES+"/Inputs/")
extent = np.loadtxt(DIRM+"mesh_extent.dat")
hole1 = np.loadtxt(DIRM+"mesh_hole1.dat")
hole2 = np.loadtxt(DIRM+"mesh_hole2.dat")
holes=[hole1,hole2]  

if len(regpar) == 3:
  # Load data for first regularization parameter
  bed_3D_1 = elmer_read.saveline_boundary(DIRR,"beta_"+regpar[0]+".dat",bbed)
  surf_3D_1 = elmer_read.saveline_boundary(DIRR,"beta_"+regpar[0]+".dat",bsur)
  taub_3D_1=elmer_read.grid3d(bed_3D_1,'taub',holes,extent)
  vel_3D_1=elmer_read.grid3d(surf_3D_1,'vel',holes,extent)
  vel_3D_1_1=elmer_read.grid3d(surf_3D_1,'velmes',holes,extent)

  # Load data for second
  bed_3D_2 = elmer_read.saveline_boundary(DIRR,"beta_"+regpar[1]+".dat",bbed)
  surf_3D_2 = elmer_read.saveline_boundary(DIRR,"beta_"+regpar[1]+".dat",bsur)
  taub_3D_2=elmer_read.grid3d(bed_3D_2,'taub',holes,extent)
  vel_3D_2=elmer_read.grid3d(surf_3D_2,'vel',holes,extent)
  vel_3D_1_2=elmer_read.grid3d(surf_3D_2,'velmes',holes,extent)

  # Load data for third
  bed_3D_3 = elmer_read.saveline_boundary(DIRR,"beta_"+regpar[2]+".dat",bbed)
  surf_3D_3 = elmer_read.saveline_boundary(DIRR,"beta_"+regpar[2]+".dat",bsur)
  taub_3D_3= elmer_read.grid3d(bed_3D_3,'taub',holes,extent)
  vel_3D_3 = elmer_read.grid3d(surf_3D_3,'vel',holes,extent)
  vel_3D_1_3 = elmer_read.grid3d(surf_3D_3,'velmes',holes,extent)

  fig=plt.figure(figsize=(7,4.5))
  
  plt.subplot(231)
  plt.imshow(taub_3D_1[2]*1e3,extent=[np.min(vel_3D_1[0]),np.max(vel_3D_1[0]),np.max(vel_3D_1[1]),np.min(vel_3D_1_1[1])],norm=matplotlib.colors.LogNorm(),clim=[1,1e3])
  plt.title('$\lambda$ = %s'%(regpar[0]))
  plt.gca().invert_yaxis()
  plt.xticks([])
  plt.yticks([])
  plt.xlim([280000,313000])
  plt.ylim([-2586000,-2553000])
  plt.plot(extent[:,0],extent[:,1],'k',hole1[:,0],hole1[:,1],'k',hole2[:,0],hole2[:,1],'k') 
  
  plt.subplot(232)
  plt.imshow(taub_3D_2[2]*1e3,extent=[np.min(vel_3D_1[0]),np.max(vel_3D_1[0]),np.max(vel_3D_1[1]),np.min(vel_3D_1[1])],norm=matplotlib.colors.LogNorm(),clim=[1,1e3])
  plt.title('$\lambda$ = %s'%(regpar[1]))
  plt.gca().invert_yaxis()
  plt.xticks([])
  plt.yticks([])
  plt.xlim([280000,313000])
  plt.ylim([-2586000,-2553000])
  plt.plot(extent[:,0],extent[:,1],'k',hole1[:,0],hole1[:,1],'k',hole2[:,0],hole2[:,1],'k') 

  plt.subplot(233)
  im1=plt.imshow(taub_3D_3[2]*1e3,extent=[np.min(vel_3D_1[0]),np.max(vel_3D_1[0]),np.max(vel_3D_1[1]),np.min(vel_3D_1[1])],norm=matplotlib.colors.LogNorm(),clim=[1,1e3])
  plt.title('$\lambda$ = %s'%(regpar[2]))
  plt.gca().invert_yaxis()
  plt.xticks([])
  plt.yticks([])
  plt.xlim([280000,313000])
  plt.ylim([-2586000,-2553000])
  plt.plot(extent[:,0],extent[:,1],'k',hole1[:,0],hole1[:,1],'k',hole2[:,0],hole2[:,1],'k') 

  plt.subplot(234)
  plt.imshow((vel_3D_1[2]-vel_3D_1_1[2])/1000,extent=[np.min(vel_3D_1[0]),np.max(vel_3D_1[0]),np.max(vel_3D_1[1]),np.min(vel_3D_1[1])],clim=([-1,1]),cmap='RdBu_r')
  plt.gca().invert_yaxis()
  plt.xticks([])
  plt.yticks([])
  plt.xlim([280000,313000])
  plt.ylim([-2586000,-2553000])
  plt.plot(extent[:,0],extent[:,1],'k',hole1[:,0],hole1[:,1],'k',hole2[:,0],hole2[:,1],'k') 

  plt.subplot(235)
  plt.imshow((vel_3D_2[2]-vel_3D_1_2[2])/1000,extent=[np.min(vel_3D_1[0]),np.max(vel_3D_1[0]),np.max(vel_3D_1[1]),np.min(vel_3D_1[1])],clim=([-1,1]),cmap='RdBu_r')
  plt.gca().invert_yaxis()
  plt.xticks([])
  plt.yticks([])
  plt.xlim([280000,313000])
  plt.ylim([-2586000,-2553000])
  plt.plot(extent[:,0],extent[:,1],'k',hole1[:,0],hole1[:,1],'k',hole2[:,0],hole2[:,1],'k') 
  
  plt.subplot(236)
  im2=plt.imshow((vel_3D_3[2]-vel_3D_1_3[2])/1000,extent=[np.min(vel_3D_1[0]),np.max(vel_3D_1[0]),np.max(vel_3D_1[1]),np.min(vel_3D_1[1])],clim=([-1,1]),cmap='RdBu_r')
  plt.gca().invert_yaxis()
  plt.xticks([])
  plt.yticks([])
  plt.xlim([280000,313000])
  plt.ylim([-2586000,-2553000])
  plt.plot(extent[:,0],extent[:,1],'k',hole1[:,0],hole1[:,1],'k',hole2[:,0],hole2[:,1],'k') 

  pylab.tight_layout()
  pylab.subplots_adjust(left=None, bottom=None, right=0.85, top=0.9, wspace=0.0, hspace=0.1)
  cbar_ax1 = fig.add_axes([0.86, 0.495, 0.03, 0.4])
  cf=fig.colorbar(im1, cax=cbar_ax1)
  cf.set_ticks([1e0,1e1,1e2,1e3])
  cf.set_label('kPa',fontsize=14,fontname="Arial")
  cf.ax.tick_params(labelsize=14)
  
  cbar_ax2 = fig.add_axes([0.86, 0.045, 0.03, 0.4])
  cf=fig.colorbar(im2, cax=cbar_ax2)
  cf.set_ticks([-1,0,1])
  cf.set_label('Difference (km/yr)',fontsize=14,fontname="Arial")
  cf.ax.tick_params(labelsize=14)

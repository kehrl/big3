
import sys
import os
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Util/Modules"))
sys.path.append(os.path.join(os.getenv("HOME"),"Code/Helheim/Tools"))
import elmer_read
import matplotlib.pyplot as plt
from pylab import *

args = sys.argv


RES=args[1]
regpars = args[2:]
DIRR=os.path.join(os.getenv("HOME"),"Models/Helheim/Results/3D/"+RES+"/Inversion/")
bbed = 3
bsur = 4

for i in range(0,len(regpars)):
  runname = "beta_"+regpars[i]+".dat"
  # Mesh boundaries
  DIRM = os.path.join(os.getenv("HOME"),"Models/Helheim/Meshes/3D/"+RES+"/Inputs/")
  extent = np.loadtxt(DIRM+"mesh_extent.dat")
  hole1 = np.loadtxt(DIRM+"mesh_hole1.dat")
  hole2 = np.loadtxt(DIRM+"mesh_hole2.dat")
  holes=[hole1,hole2]  

  # Load bed and surface
  bed_3D = elmer_read.saveline_boundary(DIRR,runname,bbed)
  surf_3D = elmer_read.saveline_boundary(DIRR,runname,bsur)

  
  taub_3D=elmer_read.grid3d(bed_3D,'taub',holes,extent)
  vel_3D=elmer_read.grid3d(surf_3D,'vel',holes,extent)
  velmes_3D=elmer_read.grid3d(surf_3D,'velmes',holes,extent)
  beta_3D=elmer_read.grid3d(bed_3D,'beta',holes,extent)
  
  figure(figsize=(16,10))
  ax1 = subplot(141)
  imshow(velmes_3D[2]/1000)
  plt.gca().invert_yaxis()
  plt.clim([0,7])
  cf = plt.colorbar(orientation='horizontal',fraction=0.3,pad=0.03)#,ax=ax1,shrink=0.4)
  setp( ax1.get_xticklabels(), visible=False)
  setp( ax1.get_yticklabels(), visible=False)
  #plt.title('Measured velocity',fontsize=28,fontname="Arial")
  cf.set_ticks([0,2,4,6])
  cf.set_label('km yr$^{-1}$',fontsize=24,fontname="Arial")
  cf.ax.tick_params(labelsize=20) 

  ax2 = subplot(142)
  imshow((vel_3D[2]-velmes_3D[2])/1000)
  plt.gca().invert_yaxis()
  plt.clim([-1,1])
  cf = plt.colorbar(orientation='horizontal',fraction=0.3,pad=0.03)#,ax=ax1,shrink=0.4)
  setp( ax2.get_xticklabels(), visible=False)
  setp( ax2.get_yticklabels(), visible=False)
  #plt.title('Modelled velocity',fontsize=28,fontname="Arial")
  cf.set_ticks([-1,0,1])
  cf.set_label('km yr$^{-1}$',fontsize=24,fontname="Arial")
  cf.ax.tick_params(labelsize=20) 

  ax3 = subplot(143)
  imshow(beta_3D[2]**2,norm=matplotlib.colors.LogNorm())
  plt.gca().invert_yaxis()
  plt.clim([0.0001,0.01])
  cf = plt.colorbar(orientation='horizontal',fraction=0.3,pad=0.03)#,ax=ax1,shrink=0.4)
  setp( ax3.get_xticklabels(), visible=False)
  setp( ax3.get_yticklabels(), visible=False)
  #plt.title('Sliding coefficient',fontsize=28,fontname="ArialMT")
  cf.set_label('MPa m$^{-1}$ yr$^{-1}$',fontsize=24,fontname="Arial")
  cf.set_ticks([0.0001,0.001, 0.01])
  cf.ax.tick_params(labelsize=20) 

  ax4 = subplot(144)
  imshow(taub_3D[2]*10**6,norm=matplotlib.colors.LogNorm())
  plt.gca().invert_yaxis()
  plt.clim([1000,10**5])
  cf = plt.colorbar(orientation='horizontal',fraction=0.3,pad=0.03)#,ax=ax1,shrink=0.4)
  setp( ax4.get_xticklabels(), visible=False)
  setp( ax4.get_yticklabels(), visible=False)
  #plt.title('Basal shear stress',fontsize=28,fontname="Arial")
  cf.set_ticks([1e1,1e3,1e6])
  cf.set_label('Pa',fontsize=24,fontname="Arial")
  cf.ax.tick_params(labelsize=20) 

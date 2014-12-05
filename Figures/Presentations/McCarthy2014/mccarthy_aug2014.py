# This script gathers some basal inversion results and makes some plots for my McCarthy 
# summer school poster.
#
# LMK, UW, 07/30/2014

import os
import sys
import numpy as np
sys.path.append(os.path.join(os.getenv("HOME"),"Dropbox/Code/Modules"))
from geodat import *
from subprocess import call
import math
import glob
import numpy as np
import elmer_read 
import matplotlib as plt
import elmer_mesh as mesh
from pylab import *
import scipy.interpolate
import scipy.signal as signal

##########
# Inputs #
##########

# What do you want to do?
plot_3D=0
plot_Lcurve_3D=1
plot_2D=0
plot_comparison=0

# 3D
DIRR_3D=os.path.join(os.getenv("HOME"),"Results/Elmer/3D/Helheim/Inversion/")
bbed_3D=3
bsurf_3D=4
runname="beta"

# Flowline
DIRR_2D=os.path.join(os.getenv("HOME"),"Results/Elmer/Flowline/Helheim/Inversion/")
DIRS_2D=os.path.join(os.getenv("HOME"),"Dropbox/Code/Solver_Files/Flowline/Helheim/Inputs/")
runname = 'beta'
bbed_2D=1
bsurf_2D=2

# Mesh boundaries
DIRX=os.path.join(os.getenv("HOME"),"Data/Shape_files/Glaciers/3D/Helheim/")
extent = mesh.shp_to_xy(DIRX+"glacier_extent_basin")
hole1 = mesh.shp_to_xy(DIRX+"glacier_hole1")
hole2 = mesh.shp_to_xy(DIRX+"glacier_hole2")
holes = []
holes.append({'xy': np.array(hole1[0:2])})
holes.append({'xy': np.array(hole2[0:2])})
del hole1, hole2

######
# 3D #
######
try:
  bed_3D
except:
  bed_3D = elmer_read.saveline_boundary(DIRR_3D,runname+"_1e10",bbed_3D)
  surf_3D = elmer_read.saveline_boundary(DIRR_3D,runname+"_1e10",bsurf_3D)

if plot_3D:

  taub_3D=elmer_read.grid3d(bed_3D,'taub',holes,extent)
  vel_3D=elmer_read.grid3d(surf_3D,'vel',holes,extent)
  velmes_3D=elmer_read.grid3d(surf_3D,'velmes',holes,extent)
  beta_3D=elmer_read.grid3d(bed_3D,'beta',holes,extent)
  
  
  plt.clf()
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
  plt.clim([0.1,10**6])
  cf = plt.colorbar(orientation='horizontal',fraction=0.3,pad=0.03)#,ax=ax1,shrink=0.4)
  setp( ax4.get_xticklabels(), visible=False)
  setp( ax4.get_yticklabels(), visible=False)
  #plt.title('Basal shear stress',fontsize=28,fontname="Arial")
  cf.set_ticks([1e1,1e3,1e6])
  cf.set_label('Pa',fontsize=24,fontname="Arial")
  cf.ax.tick_params(labelsize=20) 

  #subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.05)

############
# Flowline #
############

try:
  bed_2D
except:
  bed_2D = elmer_read.saveline_boundary(DIRR_2D,runname+"_1e11",bbed_2D)
  surf_2D = elmer_read.saveline_boundary(DIRR_2D,runname+"_1e11",bsurf_2D)
  fid=open(DIRS_2D+"width.dat","r")
  lines = fid.readlines()
  width=np.zeros([2,len(lines)-1])
  for i in range(0,len(lines)-1):
    line = lines[i+1]
    p = line.split()
    width[0,i]=float(p[0])
    width[1,i]=float(p[1])
  fid.close()
  del lines, fid, line


if plot_2D:
  plt.clf()
  plt.close()
  figure(figsize=(8,10))
  ax1 = subplot(411)
  plt.plot(np.array(surf_2D['coord1'])/1000,np.array(surf_2D['vsurfini1'])/1000,'r',linewidth=1.5)
  plt.plot(np.array(surf_2D['coord1'])/1000,np.array(surf_2D['vel1'])/1000,'k',linewidth=1.5)
  plt.legend({"Measured","Modelled"},'lower right')
  plt.xlim([0,41])
  plt.ylim([0,7])
  plt.yticks([0,2,4,6],fontsize=20)
  #plt.ylabel('Velocity \n (km yr$^{-1}$)',fontsize=20,fontname="Arial")
  plt.xticks([0,10,20,30,40])
  setp(ax1.get_xticklabels(), visible=False)
  
  ax2 = subplot(412)
  plt.plot(width[0]/1000,width[1]/1000,'k',linewidth=1.5)
  #plt.ylabel('Width \n (km)',fontsize=20,fontname="Arial")
  plt.ylim([0,25])
  plt.xlim([0,41])
  plt.xticks([0,10,20,30,40])
  plt.yticks([0,10,20],fontsize=20)
  setp(ax2.get_xticklabels(), visible=False)
  
  ax3 = subplot(413)
  plt.plot(np.array(bed_2D['coord1'])/1000,np.array(bed_2D['coord2'])/1000,'k',linewidth=1.5)
  #plt.ylabel('Bed elevation \n (m)',fontsize=20,fontname="Arial")
  plt.xlim([0, 41])
  plt.ylim([-1,0.6])
  plt.yticks([-0.5,0,0.5],fontsize=20)
  plt.xticks([0,10,20,30,40])
  setp(ax3.get_xticklabels(), visible=False)
  
  ax4 = subplot(414)
  plt.plot(np.array(bed_2D['coord1'])/1000,np.array(bed_2D['beta'])**2,'k',linewidth=1.5)
  #plt.ylabel('Sliding coefficient \n (Mpa m$^{-1}$ yr$^{-1}$)',fontsize=20,fontname="Arial")
  plt.ylim([-1e-4,1e-3])
  plt.yticks([0,1e-3],['0','10$^{-2}$'],fontsize=20)
  plt.xlim([0,41])
  plt.xlabel('Distance along flowline (km)',fontsize=24,fontname="Arial")
  plt.xticks([0,10,20,30,40],fontsize=20)
  
  subplots_adjust(left=0.18, bottom=0.07, right=0.98, top=0.98, wspace=-0.1, hspace=None)
  plt.figtext(0.02, 0.92, 'Velocity',fontsize=24,rotation=90,fontname="Arial")
  plt.figtext(0.07, 0.93, '(km yr$^{-1}$)',fontsize=20,rotation=90,fontname="Arial")
  plt.figtext(0.02, 0.25, 'Sliding coefficient',fontsize=24,rotation=90,fontname="Arial")
  plt.figtext(0.07, 0.22, '(MPa m$^{-1}$ yr$^{-1}$)',fontsize=20,rotation=90,fontname="Arial")
  plt.figtext(0.02, 0.49, 'Bed elevation',fontsize=24,rotation=90,fontname="Arial")
  plt.figtext(0.07, 0.43, '(km)',fontsize=20,rotation=90,fontname="Arial")
  plt.figtext(0.02, 0.73, 'Glacier width',fontsize=24,rotation=90,fontname="Arial")
  plt.figtext(0.07, 0.68, '(km)',fontsize=20,rotation=90,fontname="Arial")

# L curve
if plot_Lcurve_3D:
  regpar =[]
  nsim = []
  cost_tot = []
  cost_sur = []
  cost_bed = []
  fid = open(DIRR_3D+"summary.dat","r")
  lines = fid.readlines()
  for line in lines:
    p=line.split()
    if p[0] == 'Lambda':
      pass
    else:
      regpar.append(p[0])
      nsim.append(float(p[1]))
      cost_tot.append(float(p[2]))
      cost_sur.append(float(p[3]))
      cost_bed.append(float(p[4]))
      
  plt.plot(cost_bed,cost_sur,'ko--')
  plt.xlabel('$J_{reg}$',fontsize=24)
  plt.ylabel('$J_o$',fontsize=24)
  for i in range(0,len(nsim)):
    plt.text(cost_bed[i]+0.05,cost_sur[i]+0.25e9,regpar[i],fontsize=20) 
  plt.xticks([0,1,2,3,4],fontsize=16)
  plt.yticks(fontsize=16)
  subplots_adjust(left=None, bottom=0.15, right=0.98, top=0.95, wspace=-0.1, hspace=None)

# Model comparison
if plot_comparison:
  fid = open(DIRS_2D+"flowline.dat","r")
  lines = fid.readlines()
  flowline=np.zeros([5,len(lines)-1])
  for i in range(0,len(lines)-1):
    line = lines[i+1]
    p = line.split()
    flowline[0,i]=float(p[0])
    flowline[1,i]=float(p[1])
    flowline[2,i]=float(p[2])
    flowline[3,i]=float(p[3])
    flowline[4,i]=float(p[4])
  fid.close()
  del lines, fid, line  
  
  sig_d
  taub_flow_3D=np.zeros([1,len(flowline[1,:])])
  bvel_flow_3D=np.zeros([1,len(flowline[1,:])])
  svel_flow_3D=np.zeros([1,len(flowline[1,:])])
  minind=[]
  for i in range(0,len(flowline[1,:])):
    mindist=20000
    for j in range(0,len(bed_3D['coord1'])):
      dist = math.sqrt((flowline[1,i]-bed_3D['coord1'][j])**2+(flowline[2,i]-bed_3D['coord2'][j])**2)
      if dist < mindist:
        minind=j
        mindist=dist
    taub_flow_3D[0,i] = bed_3D['taub'][minind]
    bvel_flow_3D[0,i] = bed_3D['vel1'][minind]
    mindist=20000
    for j in range(0,len(surf_3D['coord1'])):
      dist = math.sqrt((flowline[1,i]-surf_3D['coord1'][j])**2+(flowline[2,i]-surf_3D['coord2'][j])**2)
      if dist < mindist:
        minind=j
        mindist=dist
    svel_flow_3D[0,i] = surf_3D['vel1'][minind]
    
    # Driving stress
    thick=np.array(surf_2D['coord2'])-np.array(bed_2D['coord2'])
    slope=abs(np.diff(np.array(surf_2D['coord2']))/np.diff(np.array(surf_2D['coord1'])))
    b,a=signal.butter(4,0.1,btype='low')
    slope=np.array(signal.filtfilt(b,a,slope))
    sigd=sin(slope)*917*9.81*thick[0:257]
    
    b,a=signal.butter(4,0.1,btype='low')
    smoothed_taub = signal.filtfilt(b,a,np.log10(taub_flow_3D))
    smoothed_taub = 10**np.array(zip(*smoothed_taub))
    smoothed_velb = np.array(zip(*signal.filtfilt(b,a,bvel_flow_3D)))
    smoothed_vels = np.array(zip(*signal.filtfilt(b,a,svel_flow_3D)))
    b,a=signal.butter(4,0.1,btype='low')
    smoothed_sigd=np.array(signal.filtfilt(b,a,sigd))
    
    plt.close()
    figure(figsize=(8,5.5))
    ax1 = plt.gca()
    plt.plot(np.array(surf_2D['coord1'][1:258])/1000,sigd,'k',linewidth=1.5)
    plt.plot(np.array(bed_2D['coord1'])/1000,10**np.array(signal.filtfilt(b,a,np.log10(bed_2D['taub'])))*10**6,'r',linewidth=1.5)
    plt.plot(flowline[0,45:-1]/1000,smoothed_taub[45:-1]*10**6,'b',linewidth=1.5)
    ax1.set_yscale('log')
    plt.xlim([0,41])
    plt.ylabel('Basal shear stress \n (Pa)',fontsize=24,fontname="Arial")
    plt.xlabel('Distance along flowline (km)',fontsize=24,fontname="Arial")
    plt.xticks([0,10,20,30,40],fontsize=20,fontname="Arial")
    plt.yticks([1e-3,1,1e3,1e6],fontsize=20,fontname="Arial")
    plt.ylim([1e-4,1e7])
    legend = plt.legend({"Driving stress","3D","Flowline"},loc='lower left',fontsize=18)
    subplots_adjust(left=0.17, bottom=0.15, right=0.98, top=0.95, wspace=None, hspace=None)
    setp(ax2.get_xticklabels(), visible=False)
    
    #ax2 = subplot(212)
    #plt.plot(np.array(bed_2D['coord1'])/1000,np.array(bed_2D['vel1'])/1000,'r',linewidth=1.5)
    #plt.plot(flowline[0,45:-1]/1000,smoothed_velb[45:-1]/1000,'b',linewidth=1.5)
    #plt.xlim([0,41])
    #plt.ylabel('Sliding speed \n (km yr$^{-1}$)',fontsize=24,fontname="Arial")
    #plt.xlabel('Distance along flowline (km)',fontsize=24,fontname="Arial")
    #plt.xticks([0,10,20,30,40],fontsize=20,fontname="Arial")
    #plt.yticks([0,2,4,6],fontsize=20,fontname="Arial")
    #subplots_adjust(left=0.17, bottom=0.2, right=0.98, top=0.95, wspace=None, hspace=None)
    
    
    
  
  #plt.plot(np.array(bed_2D['coord1'])/1000,np.array(bed_2D['taub']),'k',linewidth=1.5)

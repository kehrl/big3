import os
import sys
import numpy as np
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"Util/Modules"))
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"BigThreeGlaciers/Tools"))
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"Modules/demtools"))
import velocity, icefronts, bed, glacier_flowline, elevation, flotation
import matplotlib.pyplot as plt
import matplotlib, geotiff, fracyear, dem_shading, scipy
from matplotlib.ticker import AutoMinorLocator
import scipy.signal as signal
import gmtColormap, subprocess
cpt_rainbow = gmtColormap.get_rainbow()
plt.register_cmap(cmap=cpt_rainbow)


# Get arguments
args = sys.argv
glacier = args[1][:] # Options: Kanger, Helheim

# What bed to use for thinning gates through fluxgate method
if glacier == 'Helheim':
  bedsource = 'smith'
elif glacier == 'Kanger':
  bedsource = 'cresis'

if glacier == 'Helheim':
  xmin = 305000.0
  xmax = 314000.0
  ymin = -2582000.0
  ymax = -2570000.0
elif glacier == 'Kanger':
  xmin = 480000.0
  xmax = 497000.0
  ymin = -2298000.0
  ymax = -2282000.0

############ 
# Flowline #
############

x,y,zb,dists = glacier_flowline.load(glacier,shapefilename='flowline_flightline',filt_len=2.0e3)

############
# Get DEMs #
############

xdem,ydem,zdem,timedem,errordem = elevation.dem_grid(glacier,xmin,xmax,ymin,ymax,
	years='all',verticaldatum='geoid',return_error=True)

for i in range(0,len(timedem)):

  plt.figure(figsize=(8,4))
  plt.subplot(121)
  plt.imshow(dem_shading.hillshade(zdem[:,:,i]),extent=[xmin,xmax,ymin,ymax],origin='lower',cmap='Greys_r')
  plt.imshow((zdem[:,:,i]),extent=[xmin,xmax,ymin,ymax],origin='lower',clim=[0,400],cmap='cpt_rainbow',alpha=0.6)
  plt.plot(x,y,'w')
  plt.xlim([xmin,xmax])
  plt.ylim([ymin,ymax])
  plt.xticks([])
  plt.yticks([])
  
  
  plt.subplot(122)
  plt.plot(dists/1e3,zb,'k',lw=1.5)
  f=scipy.interpolate.RegularGridInterpolator((ydem,xdem),zdem[:,:,i],bounds_error = False,method='linear',fill_value=float('nan'))
  plt.plot(dists/1e3,f((y,x)),'k')
  plt.plot(dists/1e3,flotation.shelfbase(f(((y,x)))),'k')
  plt.xlim([-5,5])
  if glacier == 'Helheim':
    plt.ylim([-800,200])
  elif glacier == 'Kanger':
    plt.ylim([-1000,200])
  year,month,day=fracyear.fracyear_to_date(timedem[i])
  plt.title(str(year)+'-'+str(month)+'-'+str(day))
  
  plt.tight_layout()
  plt.subplots_adjust(hspace=0.05,wspace=0) 
  
  plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/movie/"+'{0:02g}'.format(i)+'.png'),format='PNG')
  plt.close()

fps=3
ffmpeg_in_opt = "-r %i" % fps
#ffmpeg_out_opt = "-r 5 -y -an -force_fps -s '%ix%i'" % newsize
#-an disables audio
ffmpeg_out_opt = "-y -an -c:v libx264 -pix_fmt yuv420p"
#scale='iw/4:-1'"

outmov = 'movie.mp4'
#cmd = 'ffmpeg {0} -i %04d.jpg {1} {2}'.format(ffmpeg_in_opt, ffmpeg_out_opt, outmov)
#cmd = 'ffmpeg {0} -pattern_type glob -i *_clip.png {1} {2}'.format(ffmpeg_in_opt, ffmpeg_out_opt, outmov)
cmd = 'ffmpeg {0} -i %02d.png {1} {2}'.format(ffmpeg_in_opt, ffmpeg_out_opt, outmov)

print cmd
subprocess.call(cmd, shell=True)


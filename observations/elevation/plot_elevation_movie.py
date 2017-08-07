import os
import sys
import numpy as np
sys.path.append(os.path.join(os.getenv("CODE_HOME"),"Modules/demtools"))
import glaclib, zslib, floatlib, datelib, demshadelib, bedlib
import matplotlib.pyplot as plt
import matplotlib, scipy
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
  xmin = 304700.0
  xmax = 313600.0
  ymin = -2582000.0
  ymax = -2573000.0
elif glacier == 'Kanger':
  xmin = 484400.0
  xmax = 497260.0
  ymin = -2298000.0
  ymax = -2282000.0

############ 
# Flowline #
############

x,y,zb,dists = glaclib.load_flowline(glacier,shapefilename='flowline_flightline',filt_len=2.0e3,verticaldatum='geoid',bedsource='cresis')

############
# Get DEMs #
############

xdem,ydem,zdem,timedem,errordem = zslib.dem_grid(glacier,xmin,xmax,ymin,ymax,
	years='all',verticaldatum='geoid',return_error=True)

####################################
# Get bed elevations near flowline #
####################################

# Get radar thicknesses close to flightline
cresis = bedlib.cresis('all',glacier)
if glacier == 'Helheim':
	cresis2001 = bedlib.cresis('2001',glacier)
	cresis = np.row_stack([cresis,cresis2001])

cutoff = 200.
dcresis = []
zcresis = []
tcresis = []
for i in range(0,len(cresis[:,0])):
	mindist = np.min(np.sqrt((cresis[i,0]-x)**2+(cresis[i,1]-(y))**2))
	if mindist < cutoff:
		minind = np.argmin(np.sqrt((cresis[i,0]-x)**2+(cresis[i,1]-(y))**2))
		dcresis.append(dists[minind])
		zcresis.append(cresis[i,2])
		tcresis.append(cresis[i,4])
dcresis = np.array(dcresis)
zcresis = np.array(zcresis)

###########################
# Make plots for each DEM #
###########################

for i in range(0,len(timedem)):

  plt.figure(figsize=(9,4))
  plt.subplot(121)
  plt.imshow(demshadelib.hillshade(zdem[:,:,i]),extent=[xmin,xmax,ymin,ymax],origin='lower',cmap='Greys_r')
  plt.imshow((zdem[:,:,i]),extent=[xmin,xmax,ymin,ymax],origin='lower',clim=[0,400],cmap='cpt_rainbow',alpha=0.6)
  plt.plot(x,y,'w')
  plt.xlim([xmin,xmax])
  plt.ylim([ymin,ymax])
  plt.xticks([])
  plt.yticks([])
  
  plt.subplot(122)
  plt.plot(dcresis/1e3,zcresis,'.',color='0.5',markersize=3)
  if glacier == 'Helheim':
    plt.plot(dists/1e3,zb,'k',lw=1)
  elif glacier == 'Kanger':
    ind = np.argmin(abs(dists--5e3))
    plt.plot(dists[0:ind]/1e3,zb[0:ind],color='k',linewidth=1.5,label='Smoothed')
  plt.plot([-10,5],[0,0],'0.8',lw=0.5)
  if i != 0:
    plt.plot(dists/1e3,f((y,x)),c='0.7')
    plt.plot(dists/1e3,floatlib.shelfbase(f(((y,x)))),c='0.7')

  f = scipy.interpolate.RegularGridInterpolator((ydem,xdem),zdem[:,:,i],bounds_error = False,method='linear',fill_value=float('nan'))
  plt.plot(dists/1e3,floatlib.height(zb),'k:')
  plt.plot(dists/1e3,f((y,x)),'b',lw=1.5)
  plt.plot(dists/1e3,floatlib.shelfbase(f(((y,x)))),'b',lw=1.5)
  if glacier == 'Helheim':
    plt.ylim([-1200,220])
    plt.xticks(range(-10,5,2),fontsize=10)
    plt.yticks(np.arange(-1200,400,200),fontsize=10)
    plt.xlim([-7,4])
  elif glacier == 'Kanger':
    plt.ylim([-1200,220])
    plt.xticks(np.arange(-10,5,2),fontsize=10)
    plt.yticks(np.arange(-1200,400,200),fontsize=10)
    plt.xlim([-7,4])
  year,month,day=datelib.fracyear_to_date(timedem[i])
  plt.title(str(year)+'-'+str(month)+'-'+str(int(day)))
  
  plt.tight_layout()
  plt.subplots_adjust(hspace=10,wspace=0) 
  
  plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/model_movie/"+'{0:02g}'.format(i)+'.png'),format='PNG',dpi=400)
  plt.close()


##############
# Make movie #
##############

fps=3
ffmpeg_in_opt = "-r %i" % fps
#ffmpeg_out_opt = "-r 5 -y -an -force_fps -s '%ix%i'" % newsize
#-an disables audio
ffmpeg_out_opt = "-y -an -c:v libx264 -pix_fmt yuv420p"
#scale='iw/4:-1'"

os.chdir(os.path.join(os.getenv("HOME"),"Bigtmp/movie/"))
outmov = glacier+'.mp4'
#cmd = 'ffmpeg {0} -i %04d.jpg {1} {2}'.format(ffmpeg_in_opt, ffmpeg_out_opt, outmov)
#cmd = 'ffmpeg {0} -pattern_type glob -i *_clip.png {1} {2}'.format(ffmpeg_in_opt, ffmpeg_out_opt, outmov)
cmd = 'ffmpeg {0} -i %02d.png {1} {2}'.format(ffmpeg_in_opt, ffmpeg_out_opt, outmov)

print cmd
subprocess.call(cmd, shell=True)

cmd = 'rm *.png'
print cmd
subprocess.call(cmd, shell=True)


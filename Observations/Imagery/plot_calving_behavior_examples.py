import os
import sys
import numpy as np
import geotifflib
import matplotlib.pyplot as plt
import matplotlib

LANDSATDIR = os.path.join(os.getenv("DATA_HOME"),"Imagery/Landsat/Helheim/TIF/")

xtab,ytab,ztab = geotifflib.readrgb(LANDSATDIR+"20140307140010_LC82320132014066LGN00.tif")
xnon,ynon,znon = geotifflib.readrgb(LANDSATDIR+"20140501140530_LC82330132014121LGN00.tif")

xmin = 305000.0
xmax = 314000.0
ymin = -2582500.0
ymax = -2572500.0

plt.figure(figsize=(3.5,2.1))
gs = matplotlib.gridspec.GridSpec(1,2)

plt.subplot(gs[0])
plt.imshow(ztab,origin='lower',extent=[np.min(xtab),np.max(xtab),np.min(ytab),np.max(ytab)])

plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
plt.xticks([])
plt.yticks([])
plt.text(xmin+0.05*(xmax-xmin),ymax-0.12*(ymax-ymin),'a',fontsize=9,fontname='arial',fontweight='bold')

plt.subplot(gs[1])
plt.imshow(znon,origin='lower',extent=[np.min(xnon),np.max(xnon),np.min(ynon),np.max(ynon)])
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
plt.xticks([])
plt.yticks([])
plt.text(xmin+0.05*(xmax-xmin),ymax-0.12*(ymax-ymin),'b',fontsize=9,fontname='arial',fontweight='bold')

plt.subplot(gs[0])
ax = plt.gca()
xmin,xmax = plt.xlim()
ymin,ymax = plt.ylim()
path = matplotlib.path.Path([[0.57*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
  			[0.98*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
  			[0.98*(xmax-xmin)+xmin,0.86*(ymax-ymin)+ymin],
  			[0.56*(xmax-xmin)+xmin,0.86*(ymax-ymin)+ymin],
  			[0.56*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin]])
patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1)
ax.add_patch(patch)

ax.plot([xmin+0.60*(xmax-xmin),xmin+0.60*(xmax-xmin)+1e3],[ymin+0.93*(ymax-ymin),ymin+0.93*(ymax-ymin)],'k',linewidth=1.5)
ax.plot([xmin+0.60*(xmax-xmin),xmin+0.60*(xmax-xmin)],[ymin+0.93*(ymax-ymin),ymin+0.91*(ymax-ymin)],'k',linewidth=1.5)
ax.plot([xmin+0.60*(xmax-xmin)+1e3,xmin+0.60*(xmax-xmin)+1e3],[ymin+0.93*(ymax-ymin),ymin+0.91*(ymax-ymin)],'k',linewidth=1.5)
ax.text(xmin+0.62*(xmax-xmin)+1e3,ymin+0.89*(ymax-ymin),'1 km',fontsize=9)

plt.tight_layout()
plt.subplots_adjust(hspace=0.05,wspace=0.05)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/Helheim_calving_example.pdf"),FORMAT='PDF',dpi=600)
plt.close()


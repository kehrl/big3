import elmerreadlib, glaclib, vellib, zslib, datelib, icefrontlib, matplotlib
import matplotlib.pyplot as plt
import numpy as np

glacier = 'Helheim'

# Load flowline
x,y,zb,dist = glaclib.load_flowline(glacier)

# Points for sampling
dists_eul = np.array([-2,-5,-10,-15,-20,-30])*1e3

# Get terminus position along flowline
terminus_val,terminus_time = icefrontlib.distance_along_flowline(x,y,dist,'Helheim',datatypes=['WV','Landsat8','TSX'])

# Load model results
t1 = 1
t2 = 305

model = elmerreadlib.pvtu_timeseries_flowline(x,y,'.','terminusdriven',\
        ['velocity','groundedmask'],layer='surface',debug=True,t1=t1,t2=t2)

# Get model time
model_time = datelib.date_to_fracyear(2008,6,13)+(np.arange(t1,t2+1)/365.25)


vel_model = np.zeros([len(model_time),len(dists_eul)])
zs_model = np.zeros([len(model_time),len(dists_eul)])
for i in range(0,len(dists_eul)):
  ind = np.argmin(abs(dists_eul[i]-dist))
  for j in range(0,len(model_time)):
    vel_model[j,i] = model['velocity'][ind,j]
    zs_model[j,i] = model['z'][ind,j]
    

colors = ['r','b','g','limegreen','gold','k']
fig = plt.figure(figsize=(5,5))
gs = matplotlib.gridspec.GridSpec(5,1)

plt.subplot(gs[0, :])
ax = plt.gca()
plt.plot(terminus_time,terminus_val/1e3,'k.--')
plt.xticks(np.arange(2008,2016,.25))
ax.set_xticklabels([])
plt.xlim([model_time[0],model_time[-1]])
plt.ylim([-2.5,2.5])
plt.yticks(np.arange(-2,3,1),fontsize=8,fontname='Arial')
plt.ylabel('Terminus \n (km)',fontsize=8)


plt.subplot(gs[1:3, :])
ax = plt.gca()
for i in range(0,len(dists_eul)):
  plt.plot(model_time,(vel_model[:,i]-vel_model[0,i])/1e3,color=colors[i],label='H{0:02d}'.format(int(-1*dists_eul[i]/1e3)))
plt.xticks(np.arange(2008,2016,.25))
ax.set_xticklabels([])
plt.xlim([model_time[0],model_time[-1]])
plt.legend(loc=0,borderpad=0.3,fontsize=8,numpoints=1,handlelength=0.2,handletextpad=0.5,labelspacing=0.1,ncol=2,columnspacing=0.8)
plt.yticks(np.arange(-3,1,1),fontsize=8)
plt.ylabel('Velocity \n (km/yr)',fontsize=8,fontname='Arial')

plt.subplot(gs[3:, :])
for i in range(0,len(dists_eul)):
  plt.plot(model_time,zs_model[:,i]-zs_model[0,i],color=colors[i])
plt.xticks(np.arange(2008,2016,.25),fontsize=8)
plt.xlim([model_time[0],model_time[-1]])
plt.yticks(np.arange(-100,25,25),fontsize=8,fontname='Arial')
plt.ylabel('Elevation (m)',fontsize=8)

plt.tight_layout()
plt.subplots_adjust(hspace=0.03,wspace=0.02,top=0.96,right=0.95,left=0.12,bottom=0.07)
plt.savefig('model_relaxation.pdf')
plt.close()


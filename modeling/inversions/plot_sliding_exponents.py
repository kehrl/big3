import os, sys
import elmerreadlib, inverselib, glaclib, datelib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse
import scipy

# Get inputs to file
parser = argparse.ArgumentParser()
parser.add_argument("-glacier", dest = "glacier", required = True,
        help = "Name of glacier.")
parser.add_argument("-temperature", dest = "temperature", required = True,
        help = "Temperature.")
parser.add_argument("-regpar", dest = "regpar", required = False,
        help = "Regularization parameter (default=1e13).", default = '1e13')
parser.add_argument("-SSA", dest = "SSA", required = False, default = True,
        help = "SSA (default=True).")
args, _ = parser.parse_known_args(sys.argv)

temperature = args.temperature
glacier = args.glacier
regpar = args.regpar
SSA = args.SSA

if SSA:
    modelname = 'SSA'
    vname = 'ssavelocity'
else:
    modelname = 'FS'
    vname = 'velocity'

if temperature == 'model':
    temperature_text = 'modelT'
else:
    temperature_text = 'constantT'

if glacier == 'Kanger':
    beta_date = '20120522'
    beta_date_string = '22 May 2012'
    dates = ['20010712','20060505','20100603','20150808']
    date_strings = ['12 July 2001','  5 May 2006',' 3 June 2013','  8 Aug 2015']
elif glacier == 'Helheim':
    beta_date = '20120316'
    beta_date_string = '16 Mar 2012'
    dates = ['20040803','20060825','20120624','20140731']
    date_strings = ['   3 Aug 2004',' 25 Aug 2008','24 June 2012',' 31 July 2014']

# Get flowline
xflow,yflow,zb,dflow = glaclib.load_flowline(glacier)

# Cutoff for calculating mean absolute residual
cutoff = 1000.0
sign = 'over'

# Get directories
if SSA and temperature == 'model':
    maindir = os.path.join(os.getenv("MODEL_HOME"),glacier+'/3D/INV_SSA_ModelT/')
elif SSA and temperature != 'model':
    maindir = os.path.join(os.getenv("MODEL_HOME"),glacier+'/3D/INV_SSA_ConstantT/')
dirs = os.listdir(maindir)

# Get length of model files for setting up variables
n = 0
for dir in dirs:
    if dir.startswith('DEM') and dir.endswith(temperature_text):
        n = n+1

# Get indices where velocity is always greater than the cutoff value
x_cutoff,y_cutoff,vsurfini_cutoff,ind_cutoff_grid, ind_cutoff = inverselib.get_velocity_cutoff(glacier,velocity_cutoff=cutoff,\
        model_dir='INV_SSA_ModelT',sign=sign)

# Sliding exponents
m = [1,2,3,4,5,6,7,8,9,10,20,50,100]

# Set up variables
times = np.zeros([n,])
misfits_ave_mape = np.zeros([len(m),n])
misfits_ave_rmse = np.zeros([len(m),n])
misfits_one_rmse = np.zeros([len(m),n])
misfits_one_mape = np.zeros([len(m),n])
misfits_inv_mape = np.zeros([n,])
misfits_inv_rmse = np.zeros([n,])

velocities_obs = np.zeros([n,])
velocities_inv = np.zeros([n,])
velocities_ave = np.zeros([len(m),n])
velocities_one = np.zeros([len(m),n])

vflow_ave = np.zeros([len(dflow),len(m),len(dates)])
vflow_obs = np.zeros([len(dflow),len(dates)])

n = 0
for dir in dirs:
    if dir.startswith('DEM') and dir.endswith(temperature_text):
        
        # Get time
        date = dir[3:11]
        times[n] = datelib.date_to_fracyear(int(date[0:4]),int(date[4:6]),int(date[6:]))

        # Get files and data
        for i in range(0,len(m)):
            if m[i] == 1:
                slidinglaw = 'linear'
            elif m[i] == 3:
                slidinglaw = 'weertman'
            else:
                slidinglaw = 'm'+str(int(m[i]))
       
            beta_file_ave = maindir+dir+'/mesh2d/steady_'+regpar+\
                '_'+modelname+'_average_2000_2016_'+temperature_text+'_'+slidinglaw+'0001.pvtu'
            beta_file_one = maindir+dir+'/mesh2d/steady_'+regpar+\
                '_'+modelname+'_DEM'+beta_date+'_'+temperature_text+'_'+slidinglaw+'0001.pvtu'
            data_ave = elmerreadlib.pvtu_file(beta_file_ave,['vsurfini',vname])
            surf_ave = elmerreadlib.values_in_layer(data_ave,'surf')
            data_one = elmerreadlib.pvtu_file(beta_file_one,['vsurfini',vname])
            surf_one = elmerreadlib.values_in_layer(data_one,'surf')

            # Get misfits
            misfits_ave_mape[i,n] = np.mean(abs(((surf_ave['vsurfini']-surf_ave[vname])\
                    /surf_ave['vsurfini'])))*100
            misfits_ave_rmse[i,n] = np.sqrt(np.mean((surf_ave['vsurfini']-surf_ave[vname])**2))
            velocities_ave[i,n] = np.mean(surf_ave[vname][ind_cutoff])

            misfits_one_mape[i,n] = np.mean(abs(((surf_one['vsurfini']-surf_one[vname])\
                    /surf_one['vsurfini'])))*100
            misfits_one_rmse[i,n] = np.sqrt(np.mean((surf_one['vsurfini']-surf_one[vname])**2))
            velocities_one[i,n] = np.mean(surf_one[vname][ind_cutoff])

            if date in dates:
                for d in range(0,len(dates)):
                    if dates[d] == date:
                        p = d
                flow = elmerreadlib.grid_to_flowline_surface(surf_ave,xflow,yflow)
                vflow_ave[:,i,p] = flow[vname]
                vflow_obs[:,p] = flow['vsurfini']
                del d

        # Get inversion and modeled velocities
        beta_file_orig = maindir+dir+'/mesh2d/steady_'+regpar+\
                '_'+modelname+'_DEM'+date+'_'+temperature_text+'_linear0001.pvtu'
        data_orig = elmerreadlib.pvtu_file(beta_file_orig,['vsurfini',vname])
        surf_orig = elmerreadlib.values_in_layer(data_orig,'surf')

        misfits_inv_mape[n] = np.mean(abs(((surf_orig['vsurfini']-surf_orig[vname])\
                /surf_orig['vsurfini'])))*100
        misfits_inv_rmse[n] = np.sqrt(np.mean((surf_orig['vsurfini']-surf_orig[vname])**2))
        velocities_obs[n] = np.mean(surf_orig['vsurfini'][ind_cutoff])
        velocities_inv[n] = np.mean(surf_orig[vname][ind_cutoff])

        n = n+1

ind = [[0,0],[0,1],[1,0],[1,1]]
ms = [0,2,4]
iqr_ave = abs(np.mean(misfits_ave_rmse,axis=1)-np.percentile(misfits_ave_rmse,[25,75],axis=1))
iqr_one = abs(np.mean(misfits_one_rmse,axis=1)-np.percentile(misfits_one_rmse,[25,75],axis=1))
iqr_inv = np.percentile(misfits_inv_rmse,[25,75])
fig = plt.figure(figsize=(6.5,3))
matplotlib.rc('font',family='Arial')
gs1 = matplotlib.gridspec.GridSpec(1,1)
gs1.update(left=0.075, right=0.46, bottom=0.13, top=0.98, wspace=0.03, hspace=0.03)
ax1 = plt.subplot(gs1[0])
ax1.tick_params(labelsize=8)
plt.xlabel(r'Sliding exponent $m$',fontsize=8)
plt.ylabel(r'RMSE (m yr$^{-1}$)',fontsize=8)
plt.xticks(np.arange(0,len(m)))
ax1.set_xticklabels(m)
ax1.plot([-1,len(m)+1],[np.mean(misfits_inv_rmse),np.mean(misfits_inv_rmse)],'k--',label='Inversion')
ax1.fill_between([-1,len(m)+1],[iqr_inv[0],iqr_inv[0]],[iqr_inv[1],iqr_inv[1]],color='0.7')
ax1.errorbar(np.arange(0,len(m)),np.mean(misfits_one_rmse,axis=1),yerr=iqr_one,capsize=3,fmt='ko',ecolor='r',\
        markersize=5,markerfacecolor='r',label=beta_date_string)
ax1.errorbar(np.arange(0,len(m)),np.mean(misfits_ave_rmse,axis=1),yerr=iqr_ave,capsize=3,fmt='ks',ecolor='b',\
        markersize=5,markerfacecolor='b',label='Average')
ax1.set_xlim([-1,len(m)])
if glacier == 'Kanger':
    plt.text(0.7,305,'(a)',fontsize=8,fontweight='bold')
elif glacier == 'Helheim':
    plt.text(-0.6,335,'(a)',fontsize=8,fontweight='bold')
plt.legend(loc=0,labelspacing=0.4,handletextpad=1,handlelength=1.5,fontsize=8)
gs2 = matplotlib.gridspec.GridSpec(2,2)
gs2.update(left=0.55, right=0.99, bottom=0.13, top=0.98, wspace=0.03, hspace=0.03)
for i in range(0,len(dates)):
    ax2 = plt.subplot(gs2[ind[i][0],ind[i][1]])
    plt.plot(dflow/1e3,vflow_obs[:,i]/1e3,'k',label='Observed',lw=3)
    if glacier == 'Kanger':
        plt.text(-16,1.5,date_strings[i],fontsize=8)
    elif glacier == 'Helheim':
        plt.text(-21,1.5,date_strings[i],fontsize=8)
    for j in ms:
        plt.plot(dflow/1e3,vflow_ave[:,j,i]/1e3,'--',label=r'$m=$'+'{0}'.format(int(m[j])))
    ax2.tick_params(labelsize=8)
    if ind[i][1] != 0:
        ax2.set_yticklabels([])
        if glacier == 'Kanger':
            plt.text(-27,10.3,['(c)','(e)'][ind[i][0]],fontsize=8,fontweight='bold')
        elif glacier == 'Helheim':
            plt.text(-39.2,9,['(c)','(e)'][ind[i][0]],fontsize=8,fontweight='bold')
    elif ind[i][1] == 0 and ind[i][0] == 1:
        plt.ylabel('                              Velocity (km yr$^{-1}$)',fontsize=8)
    else:
        plt.legend(loc=2,labelspacing=0.25,handletextpad=1,handlelength=1.5,fontsize=8)
        if glacier == 'Kanger':
            plt.text(-27,5,'(b)',fontsize=8,fontweight='bold')
        elif glacier == 'Helheim':
            plt.text(-39.2,4.5,'(b)',fontsize=8,fontweight='bold')
    if ind[i][0] == 0:
        ax2.set_xticklabels([])
    elif ind[i][1] == 0 and ind[i][0] == 1:
        plt.xlabel('                                                  Distance along flowline (km)',fontsize=8)
        if glacier == 'Kanger':
            plt.text(-27,10.3,'(d)',fontsize=8,fontweight='bold')
        elif glacier == 'Helheim':
            plt.text(-39.2,9,'(d)',fontsize=8,fontweight='bold')
    if glacier == 'Kanger':
        plt.xlim([-28.1,-4]); plt.ylim([1.2,11.5])
    elif glacier == 'Helheim':
        plt.ylim([1.2,10])
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+'_sliding_exponent_'+modelname+'_'+\
        temperature_text+'_RMSE.pdf'),FORMAT='PDF',dpi=400)
print "Saving as "+os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+'_sliding_exponent_'+modelname+'_'+\
        temperature_text+'_RMSE.pdf')
plt.close()



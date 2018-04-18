import os, sys
import elmerreadlib, inverselib, glaclib, datelib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse

# Get inputs to file
parser = argparse.ArgumentParser()
parser.add_argument("-glacier",dest="glacier",required = True,
        help = "Name of glacier.")
parser.add_argument("-temperature", dest="temperature", required = True,
        help = "Temperature.")
parser.add_argument("-regpar",dest="regpar",required = False,
        help = "Regularization parameter", default = '1e13')
args, _ = parser.parse_known_args(sys.argv)

vname = 'ssavelocity'
temperature = args.temperature
glacier = args.glacier
regpar = args.regpar

if temperature == 'model':
    temperature_text = 'modelT'
else:
    temperature_text = 'constantT'

# Cutoff for calculating mean absolute residual
cutoff = 1000

# Get directories
maindir = os.path.join(os.getenv("MODEL_HOME"),glacier+'/3D/')
dirs = os.listdir(maindir)

# Get indices where velocity is always greater than the cutoff value
x_cutoff,y_cutoff,vsurfini_cutoff,ind_cutoff_grid, ind_cutoff  = inverselib.get_velocity_cutoff(glacier,velocity_cutoff=cutoff)

# Set up variables
misfit_inversion_mape = []
misfit_inversion_rmse = []
misfit_m1_ave_mape = []
misfit_m1_ave_rmse = []
misfit_m3_ave_mape = []
misfit_m3_ave_rmse = []
velocities_obs = []
velocities_inv = []
velocities_m1_ave = []
velocities_m3_ave = []
times = []

n = 0
for dir in dirs:
    if dir.startswith('DEM') and dir.endswith(temperature_text):
        beta_date = dir[3:11]
        beta_file_m1 = maindir+dir+'/mesh2d/steady_'+regpar+\
                '_SSA_DEM'+beta_date+'_'+temperature_text+'_linear0001.pvtu'
        beta_file_m3 = maindir+dir+'/mesh2d/steady_'+regpar+\
                '_SSA_DEM'+beta_date+'_'+temperature_text+'_weertman0001.pvtu'
        beta_file_m1_ave = maindir+dir+'/mesh2d/steady_'+regpar+\
                '_SSA_average_2000_2016_'+temperature_text+'_linear0001.pvtu'
        beta_file_m3_ave = maindir+dir+'/mesh2d/steady_'+regpar+\
                '_SSA_average_2000_2016_'+temperature_text+'_weertman0001.pvtu'

        data_m1 = elmerreadlib.pvtu_file(beta_file_m1,['vsurfini',vname])
        surf_m1 = elmerreadlib.values_in_layer(data_m1,'surf')

        data_m3 = elmerreadlib.pvtu_file(beta_file_m3,['vsurfini',vname])
        surf_m3 = elmerreadlib.values_in_layer(data_m3,'surf')

        data_m1_ave = elmerreadlib.pvtu_file(beta_file_m1_ave,['vsurfini',vname])
        surf_m1_ave = elmerreadlib.values_in_layer(data_m1_ave,'surf')

        data_m3_ave = elmerreadlib.pvtu_file(beta_file_m3_ave,['vsurfini',vname])
        surf_m3_ave = elmerreadlib.values_in_layer(data_m3_ave,'surf')

        # Calculate misfits
        times.append(datelib.date_to_fracyear(int(beta_date[0:4]),int(beta_date[4:6]),int(beta_date[6:])))
        velocities_obs.append(np.mean(surf_m1['vsurfini'][ind_cutoff]))
        misfit_inversion_mape.append(np.mean(abs(surf_m1['vsurfini'][ind_cutoff]-\
                surf_m1[vname][ind_cutoff])/surf_m1['vsurfini'][ind_cutoff])*100)
        misfit_inversion_rmse.append(np.sqrt(np.mean((surf_m1['vsurfini'][ind_cutoff]-\
                surf_m1[vname][ind_cutoff])**2)))
        velocities_inv.append(np.mean(surf_m1[vname][ind_cutoff]))
        misfit_m1_ave_mape.append(np.mean(abs(surf_m1_ave['vsurfini'][ind_cutoff]-\
                surf_m1_ave[vname][ind_cutoff])/surf_m1_ave['vsurfini'][ind_cutoff])*100)
        misfit_m1_ave_rmse.append(np.sqrt(np.mean((surf_m1_ave['vsurfini'][ind_cutoff]-\
                surf_m1_ave[vname][ind_cutoff])**2)))
        velocities_m1_ave.append(np.mean(surf_m1_ave[vname][ind_cutoff]))
        misfit_m3_ave_mape.append(np.mean(abs(surf_m3_ave['vsurfini'][ind_cutoff]-\
                surf_m3_ave[vname][ind_cutoff])/surf_m3_ave['vsurfini'][ind_cutoff])*100)
        misfit_m3_ave_rmse.append(np.sqrt(np.mean((surf_m3_ave['vsurfini'][ind_cutoff]-\
                surf_m3_ave[vname][ind_cutoff])**2)))
        velocities_m3_ave.append(np.mean(surf_m3_ave[vname][ind_cutoff]))

times = np.array(times)
velocities_obs = np.array(velocities_obs)
velocities_inv = np.array(velocities_inv)
velocities_m1_ave = np.array(velocities_m1_ave)
velocities_m3_ave = np.array(velocities_m3_ave)
misfit_inversion_mape = np.array(misfit_inversion_mape)
misfit_inversion_rmse = np.array(misfit_inversion_rmse)
misfit_m1_ave_mape = np.array(misfit_m1_ave_mape)
misfit_m1_ave_rmse = np.array(misfit_m1_ave_rmse)
misfit_m3_ave_mape = np.array(misfit_m3_ave_mape)
misfit_m3_ave_rmse = np.array(misfit_m3_ave_rmse)

fig = plt.figure(figsize=(7,6))
gs = matplotlib.gridspec.GridSpec(3,1)
ax = plt.subplot(gs[0])
plt.plot(times,velocities_obs,'ks',label='Observed',markersize=10,markerfacecolor='w')
plt.plot(times,velocities_inv,'k^',label='Inversion',markerfacecolor='g')
plt.plot(times,velocities_m1_ave,'ko',markerfacecolor='r',label=r'Average $m=1$')
plt.plot(times,velocities_m3_ave,'ko',markerfacecolor='b',label=r'Average $m=3$')
plt.legend(loc=2,labelspacing=0.1,handletextpad=1,handlelength=0.1)
plt.xlim([2001,2016])
if glacier == 'Kanger':
    plt.ylim([2200,4700])
elif glacier == 'Helheim':
    plt.ylim([2900,4100])
ax.set_xticklabels([])
plt.ylabel(r'$\bar{u}$ (m / yr)')
ax = plt.subplot(gs[1])
plt.plot(times,misfit_inversion_rmse,'k^',label='Inversion',markerfacecolor='g')
plt.plot(times,misfit_m1_ave_rmse,'ko',label=r'Average $m=1$',markerfacecolor='r')
plt.plot(times,misfit_m3_ave_rmse,'ko',label=r'Average $m=3$',markerfacecolor='b')
plt.xlim([2001,2016])
plt.ylabel('RMSE (m / yr)')
ax.set_xticklabels([])
plt.subplot(gs[2])
plt.plot(times,misfit_inversion_mape,'k^',label='Inversion',markerfacecolor='g')
plt.plot(times,misfit_m1_ave_mape,'ko',label=r'Average $m=1$',markerfacecolor='r')
plt.plot(times,misfit_m3_ave_mape,'ko',label=r'Average $m=3$',markerfacecolor='b')
plt.ylabel('MAR (%)')
plt.xlim([2001,2016])
plt.tight_layout()
plt.subplots_adjust(hspace=0.03,wspace=0.03,top=0.99,right=0.96,left=0.1,bottom=0.05)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+'_beta_'+temperature_text+'_MAR.pdf'),FORMAT='PDF',dpi=400)
print "Saving as "+os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+'_beta_'+temperature_text+'_MAR.pdf')
plt.close()

fig = plt.figure(figsize=(3.5,4))
gs = matplotlib.gridspec.GridSpec(2,1)
ind = np.where(times > 2011)
ax = plt.subplot(gs[0])
plt.plot(12*(times[ind]-np.floor(times[ind])),misfit_inversion_rmse[ind],'k^',label='Inversion',markerfacecolor='g')
plt.plot(12*(times[ind]-np.floor(times[ind])),misfit_m1_ave_rmse[ind],'ko',label=r'Average $m=1$',markerfacecolor='b')
plt.plot(12*(times[ind]-np.floor(times[ind])),misfit_m3_ave_rmse[ind],'ko',label=r'Average $m=3$',markerfacecolor='r')
plt.xlim([0,12])
plt.ylim([50,350])
plt.legend(labelspacing=0.1,handletextpad=1,handlelength=0.1,loc=2)
ax.set_xticks(np.arange(0,13))
ax.set_xticklabels([])
plt.yticks(np.arange(100,400,100))
plt.ylabel('RMSE (m / yr)')
ax = plt.subplot(gs[1])
plt.plot(12*(times[ind]-np.floor(times[ind])),misfit_inversion_mape[ind],'k^',label='Inversion',markerfacecolor='g')
plt.plot(12*(times[ind]-np.floor(times[ind])),misfit_m1_ave_mape[ind],'ko',label=r'Average $m=1$',markerfacecolor='r')
plt.plot(12*(times[ind]-np.floor(times[ind])),misfit_m3_ave_mape[ind],'ko',label=r'Average $m=3$',markerfacecolor='b')
plt.ylabel('MAR (%)')
plt.xlim([0,12])
plt.ylim([1.5,10])
ax.set_xticks(np.arange(0,13))
ax.set_xticklabels(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan'],\
        rotation=45,ha='right')
plt.tight_layout()
plt.subplots_adjust(hspace=0.03,wspace=0.03,top=0.98,right=0.98,left=0.17,bottom=0.09)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+'_beta_'+temperature_text+'_MAR_seasonal.pdf'),FORMAT='PDF',dpi=400)
print "Saving as "+os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+'_beta_'+temperature_text+'_MAR_seasonal.pdf')
plt.close()

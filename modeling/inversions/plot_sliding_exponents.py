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
    date = '20120522'
elif glacier == 'Helheim':
    date = '20120316'

# Cutoff for calculating mean absolute residual
cutoff = 1000.0
sign = 'over'

# Get directories
maindir = os.path.join(os.getenv("MODEL_HOME"),glacier+'/3D/')
dirs = os.listdir(maindir)

# Get length of model files for setting up variables
n = 0
for dir in dirs:
    if dir.startswith('DEM') and dir.endswith(temperature_text):
        n = n+1

# Get indices where velocity is always greater than the cutoff value
x_cutoff,y_cutoff,vsurfini_cutoff,ind_cutoff_grid, ind_cutoff = inverselib.get_velocity_cutoff(glacier,velocity_cutoff=cutoff,sign=sign)

# Sliding exponents
m = 10

# Set up variables
times = np.zeros([n,])
misfits_ave_mape = np.zeros([m,n])
misfits_ave_rmse = np.zeros([m,n])
misfits_one_rmse = np.zeros([m,n])
misfits_one_mape = np.zeros([m,n])
misfits_inv_mape = np.zeros([n,])
misfits_inv_rmse = np.zeros([n,])

velocities_obs = np.zeros([n,])
velocities_inv = np.zeros([n,])
velocities_ave = np.zeros([m,n])
velocities_one = np.zeros([m,n])

n = 0
for dir in dirs:
    if dir.startswith('DEM') and dir.endswith(temperature_text):
        
        # Get time
        beta_date = dir[3:11]
        times[n] = datelib.date_to_fracyear(int(beta_date[0:4]),int(beta_date[4:6]),int(beta_date[6:]))

        # Get files and data
        for i in range(0,m):
            if i == 0:
                slidinglaw = 'linear'
            elif i == 2:
                slidinglaw = 'weertman'
            else:
                slidinglaw = 'm'+str(int(i+1))
       
            beta_file_ave = maindir+dir+'/mesh2d/steady_'+regpar+\
                '_'+modelname+'_average_2000_2016_'+temperature_text+'_'+slidinglaw+'0001.pvtu'
            beta_file_one = maindir+dir+'/mesh2d/steady_'+regpar+\
                '_'+modelname+'_DEM'+date+'_'+temperature_text+'_'+slidinglaw+'0001.pvtu'
            data_ave = elmerreadlib.pvtu_file(beta_file_ave,['vsurfini',vname])
            surf_ave = elmerreadlib.values_in_layer(data_ave,'surf')
            data_one = elmerreadlib.pvtu_file(beta_file_one,['vsurfini',vname])
            surf_one = elmerreadlib.values_in_layer(data_one,'surf')

            # Get misfits
            misfits_ave_mape[i,n] = np.mean(abs(((surf_ave['vsurfini'][ind_cutoff]-surf_ave[vname][ind_cutoff])\
                    /surf_ave['vsurfini'][ind_cutoff])))*100
            misfits_ave_rmse[i,n] = np.sqrt(np.mean((surf_ave['vsurfini'][ind_cutoff]-surf_ave[vname][ind_cutoff])**2))
            velocities_ave[i,n] = np.mean(surf_ave['vsurfini'][ind_cutoff])

            misfits_one_mape[i,n] = np.mean(abs(((surf_one['vsurfini'][ind_cutoff]-surf_one[vname][ind_cutoff])\
                    /surf_one['vsurfini'][ind_cutoff])))*100
            misfits_one_rmse[i,n] = np.sqrt(np.mean((surf_one['vsurfini'][ind_cutoff]-surf_one[vname][ind_cutoff])**2))
            velocities_one[i,n] = np.mean(surf_one['vsurfini'][ind_cutoff])

        # Get inversion and modeled velocities
        beta_file_orig = maindir+dir+'/mesh2d/steady_'+regpar+\
                '_'+modelname+'_DEM'+beta_date+'_'+temperature_text+'_linear0001.pvtu'
        data_orig = elmerreadlib.pvtu_file(beta_file_orig,['vsurfini',vname])
        surf_orig = elmerreadlib.values_in_layer(data_orig,'surf')

        misfits_inv_mape[n] = np.mean(abs(((surf_orig['vsurfini'][ind_cutoff]-surf_orig[vname][ind_cutoff])\
                /surf_orig['vsurfini'][ind_cutoff])))*100
        misfits_inv_rmse[n] = np.sqrt(np.mean((surf_orig['vsurfini'][ind_cutoff]-surf_orig[vname][ind_cutoff])**2))

        velocities_obs[n] = np.mean(surf_orig['vsurfini'][ind_cutoff])
        velocities_inv[n] = np.mean(surf_orig[vname][ind_cutoff])
        
        n = n+1

fig = plt.figure(figsize=(3.5,3.0))
ax1 = plt.gca()
iqr_ave = abs(np.mean(misfits_ave_rmse,axis=1)-np.percentile(misfits_ave_rmse,[25,75],axis=1))
iqr_one = abs(np.mean(misfits_one_rmse,axis=1)-np.percentile(misfits_one_rmse,[25,75],axis=1))
iqr_inv = np.percentile(misfits_inv_rmse,[25,75])
plt.xlabel(r'Sliding exponent $m$')
plt.ylabel('RMSE (m / yr)')
plt.xticks(np.arange(1,m+1,1))
ax1.plot([0,m+1],[np.mean(misfits_inv_rmse),np.mean(misfits_inv_rmse)],'k--',label='Inversions')
ax1.fill_between([0,m+1],[iqr_inv[0],iqr_inv[0]],[iqr_inv[1],iqr_inv[1]],color='0.7')
ax1.errorbar(np.arange(1,m+1),np.mean(misfits_one_rmse,axis=1),yerr=iqr_one,capsize=3,fmt='ko',ecolor='r',markerfacecolor='r',label='20120522')
ax1.errorbar(np.arange(1,m+1),np.mean(misfits_ave_rmse,axis=1),yerr=iqr_ave,capsize=3,fmt='ks',ecolor='b',markerfacecolor='b',label='Average')
ax1.set_xlim([0,m+1])
plt.legend(labelspacing=0.1,handletextpad=1,handlelength=1)

plt.tight_layout()
plt.subplots_adjust(hspace=0.03,wspace=0.03,top=0.98,right=0.98,left=0.15,bottom=0.15)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+'_sliding_exponent_'+modelname+'_'+temperature_text+'_'+sign+str(int(cutoff))+'_RMSE.pdf'),FORMAT='PDF',dpi=400)
print "Saving as "+os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+'_sliding_exponent_'+modelname+'_'+temperature_text+'_'+sign+str(int(cutoff))+'_RMSE.pdf')
plt.close()

fig = plt.figure(figsize=(3.5,3.0))
ax1 = plt.gca()
iqr_ave = abs(np.mean(misfits_ave_mape,axis=1)-np.percentile(misfits_ave_mape,[25,75],axis=1))
iqr_one = abs(np.mean(misfits_one_mape,axis=1)-np.percentile(misfits_one_mape,[25,75],axis=1))
iqr_inv = np.percentile(misfits_inv_mape,[25,75])
plt.xlabel(r'Sliding exponent $m$')
plt.ylabel('MAR (%)')
plt.xticks(np.arange(1,m+1,1))
ax1.plot([0,m+1],[np.mean(misfits_inv_mape),np.mean(misfits_inv_mape)],'k--',label='Inversions')
ax1.fill_between([0,m+1],[iqr_inv[0],iqr_inv[0]],[iqr_inv[1],iqr_inv[1]],color='0.7')
ax1.errorbar(np.arange(1,m+1),np.mean(misfits_one_mape,axis=1),yerr=iqr_one,capsize=3,fmt='ko',ecolor='r',markerfacecolor='r',label='20120522')
ax1.errorbar(np.arange(1,m+1),np.mean(misfits_ave_mape,axis=1),yerr=iqr_ave,capsize=3,fmt='ks',ecolor='b',markerfacecolor='b',label='Average')
ax1.set_xlim([0,m+1])
plt.legend(labelspacing=0.1,handletextpad=1,handlelength=1)

plt.tight_layout()
plt.subplots_adjust(hspace=0.03,wspace=0.03,top=0.98,right=0.98,left=0.15,bottom=0.15)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+'_sliding_exponent_'+modelname+'_'+temperature_text+'_'+sign+str(int(cutoff))+'_MAR.pdf'),FORMAT='PDF',dpi=400)
print "Saving as "+os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_sliding_exponent_"+modelname+'_'+temperature_text+'_'+sign+str(int(cutoff))+"_MAR.pdf")
plt.close()


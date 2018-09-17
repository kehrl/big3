import os, sys
import geotifflib, vellib, floatlib, climlib, icefrontlib, elmerreadlib, inverselib, glaclib, datelib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib, scipy
import argparse
import scipy.stats

# Get inputs to file
parser = argparse.ArgumentParser()
parser.add_argument("-glacier",dest="glacier",required = True,
        help = "Name of glacier.")
args, _ = parser.parse_known_args(sys.argv)

glacier = args.glacier

# Cutoff for calculating mean absolute residual
cutoff = 1000.0
if glacier == 'Kanger':
    beta_suffix = 'DEM20120522'
elif glacier == 'Helheim':
    beta_suffix = 'DEM20120316'

# Get directories
DIRG = os.path.join(os.getenv("MODEL_HOME"),glacier+'/Results/')
DIR_SSA_MT = DIRG+"INV_SSA_ModelT/"
DIR_SSA_CT = DIRG+"INV_SSA_ConstantT/"
DIR_FS_MT = DIRG+"INV_FS_ModelT/"
DIR_FS_CT = DIRG+"INV_FS_ConstantT/"

# Get indices where velocity is always greater than the cutoff value
clip = 0.0
x_cutoff_SSA,y_cutoff_SSA,vsurfini_cutoff_SSA,ind_cutoff_grid,ind_cutoff_SSA  = inverselib.get_velocity_cutoff(glacier,\
        velocity_cutoff=cutoff,model_dir='INV_SSA_ModelT',SSA=True,sign='over',clip_boundary=clip)
x_cutoff_SSA_fast,y_cutoff_SSA_fast,vsurfini_cutoff_SSA_fast,ind_cutoff_grid_fast,ind_cutoff_SSA_fast  = inverselib.get_velocity_cutoff(glacier,\
        velocity_cutoff=4000,model_dir='INV_SSA_ModelT',SSA=True,sign='over',clip_boundary=clip)
x_cutoff_SSA_slow,y_cutoff_SSA_slow,vsurfini_cutoff_SSA_slow,ind_cutoff_grid_slow,ind_cutoff_SSA_slow  = inverselib.get_velocity_cutoff(glacier,\
        velocity_cutoff=1000,model_dir='INV_SSA_ModelT',SSA=True,sign='under',clip_boundary=clip)

# Get ice front positions
xflow,yflow,zflow,dflow = glaclib.load_flowline(glacier)
term_pos,term_time = icefrontlib.distance_along_flowline(xflow,yflow,dflow,glacier)

# Some constants
rho_i = 917.
rho_sw = 1025.
g = 9.8

# Get runoff
ind = np.argmin(abs(dflow+5e3))
xrac,yrac,runoff,runoff_time = climlib.racmo_at_pts(xflow[ind],yflow[ind],'runoff',filt_len=14)
year_runoff,day1_runoff,day2_runoff,meltlength_runoff,total_runoff = climlib.seasonlength(runoff_time,runoff,'runoff')

# Get dates
dirs = os.listdir(DIR_FS_CT)
dates = []
for dir in dirs:
    if dir.startswith('DEM') and not(dir[3:11] in dates):
        dates.append(dir[3:11])

# Set up variables
modelnames = ['FS-CT','FS-MT','SSA-CT','SSA-MT']
temperature_texts = ['constantT','modelT','constantT','modelT']
DIRs = [DIR_FS_CT, DIR_FS_MT, DIR_SSA_CT, DIR_SSA_MT]
misfit_inv = np.zeros([len(dates),4])
misfit_linear = np.zeros([len(dates),4])
us_obs = np.zeros([len(dates),])
us_inv = np.zeros([len(dates),4])
ub_inv = np.zeros([len(dates),4])
us_obs_fast = np.zeros([len(dates),])
us_inv_fast = np.zeros([len(dates),4])
ub_inv_fast = np.zeros([len(dates),4])
us_obs_med = np.zeros([len(dates),])
us_inv_med = np.zeros([len(dates),4])
ub_inv_med = np.zeros([len(dates),4])
us_obs_slow = np.zeros([len(dates),])
us_inv_slow = np.zeros([len(dates),4])
us_linear = np.zeros([len(dates),4])
ub_linear = np.zeros([len(dates),4])
beta = np.zeros([len(dates),4])
beta_fast = np.zeros([len(dates),4])
beta_med = np.zeros([len(dates),4])
beta_slow = np.zeros([len(dates),4])
taub = np.zeros([len(dates),4])
taub_fast = np.zeros([len(dates),4])
taub_med = np.zeros([len(dates),4])
taub_slow = np.zeros([len(dates),4])
taud = np.zeros([len(dates),])
taud_fast = np.zeros([len(dates),])
taud_med = np.zeros([len(dates),])
taud_slow = np.zeros([len(dates),])
H = np.zeros([len(dates),])
P_i = np.zeros([len(dates),])
P_w = np.zeros([len(dates),])
P_i_med = np.zeros([len(dates),])
P_i_fast = np.zeros([len(dates),])
P_w_med = np.zeros([len(dates),])
P_w_fast = np.zeros([len(dates),])
H_fast = np.zeros([len(dates),])
H_med = np.zeros([len(dates),])
Hf = np.zeros([len(dates),])
Hf_fast = np.zeros([len(dates),])
Hf_med = np.zeros([len(dates),])
times = np.zeros([len(dates),])

for i in range(0,len(dates)):
    for j in range(0,len(DIRs)):
        times[i] = datelib.date_to_fracyear(int(dates[i][0:4]),int(dates[i][4:6]),int(dates[i][6:]))

        if 'FS' in DIRs[j] and 'ConstantT' in DIRs[j]:
            # Get heights
            x,y,zs_grid = geotifflib.read(DIRs[j]+'DEM'+dates[i]+'_surf_mea_zs.tif')
            x,y,zb_grid = geotifflib.read(DIRs[j]+'DEM'+dates[i]+'_bed_mod_zb.tif')
            grid = zs_grid - zb_grid
            grid[ind_cutoff_grid_fast] = np.float('nan')
            H_fast[i] = np.nanmean(grid)
            P_i_fast[i] = np.nanmean(rho_i*g*grid)
            P_w_fast[i] = np.nanmean(-rho_sw*g*zb_grid)
            grid = zs_grid - zb_grid
            grid[ind_cutoff_grid] = np.float('nan')
            zb_grid[ind_cutoff_grid] = np.float('nan')
            H_med[i] = np.nanmean(grid[vsurfini_cutoff_SSA < 4000])
            P_i_med[i] = np.nanmean(rho_i*g*grid[vsurfini_cutoff_SSA < 4000])
            P_w_med[i] = np.nanmean(-rho_sw*g*zb_grid[vsurfini_cutoff_SSA < 4000])
            H[i] = np.nanmean(grid)
            P_i[i] = np.nanmean(rho_i*g*grid)
            P_w[i] = np.nanmean(-rho_sw*g*zb_grid)

            #Get height above flotation
            grid = zs_grid - floatlib.height(zb_grid)
            grid[ind_cutoff_grid_fast] = np.float('nan')
            Hf_fast[i] = np.nanmean(grid)
            grid = zs_grid - floatlib.height(zb_grid)
            grid[ind_cutoff_grid] = np.float('nan')
            Hf_med[i] = np.nanmean(grid[vsurfini_cutoff_SSA < 4000])
            Hf[i] = np.nanmean(grid)

            # Get observed velocity
            x,y,us_obs_grid = geotifflib.read(DIRs[j]+'DEM'+dates[i]+'_surf_mea_us.tif')
            x,y,vs_obs_grid = geotifflib.read(DIRs[j]+'DEM'+dates[i]+'_surf_mea_vs.tif')
            grid = np.sqrt(us_obs_grid**2+vs_obs_grid**2)
            grid[ind_cutoff_grid_fast] = np.float('nan')
            us_obs_fast[i] = np.nanmean(grid)
            grid = np.sqrt(us_obs_grid**2+vs_obs_grid**2)
            grid[ind_cutoff_grid_slow] = np.float('nan')
            us_obs_slow[i] = np.nanmean(grid)
            grid = np.sqrt(us_obs_grid**2+vs_obs_grid**2)
            grid[ind_cutoff_grid] = np.float('nan')
            us_obs_med[i] = np.nanmean(grid[vsurfini_cutoff_SSA < 4000])
            us_obs[i] = np.nanmean(grid)

            # Get driving stress
            xt,yt,taud_grid = elmerreadlib.input_file(os.path.join(os.getenv("MODEL_HOME"),glacier+'/3D/INV_FS_ConstantT/DEM'+dates[i]+'_constantT/inputs/taud.xy'))
            f = scipy.interpolate.RegularGridInterpolator((yt,xt),taud_grid)
            xmesh,ymesh = np.meshgrid(x,y)
            grid = f((ymesh.flatten(),xmesh.flatten())).reshape([len(y),len(x)])
            grid[ind_cutoff_grid_fast] = np.float('nan')
            taud_fast[i] = np.nanmean(grid)
            grid = f((ymesh.flatten(),xmesh.flatten())).reshape([len(y),len(x)])
            grid[ind_cutoff_grid_slow] = np.float('nan')
            taud_slow[i] = np.nanmean(grid)
            grid = f((ymesh.flatten(),xmesh.flatten())).reshape([len(y),len(x)])
            grid[ind_cutoff_grid] = np.float('nan')
            taud_med[i] = np.nanmean(grid[vsurfini_cutoff_SSA < 4000])
            taud[i] = np.nanmean(grid)

        # Get surface velocity from inversion
        x,y,us_inv_grid = geotifflib.read(DIRs[j]+'DEM'+dates[i]+'_surf_mod_us.tif')
        x,y,vs_inv_grid = geotifflib.read(DIRs[j]+'DEM'+dates[i]+'_surf_mod_vs.tif')
        grid = np.sqrt(us_inv_grid**2+vs_inv_grid**2)
        grid[ind_cutoff_grid_fast] = np.float('nan')
        us_inv_fast[i,j] = np.nanmean(grid)
        grid = np.sqrt(us_inv_grid**2+vs_inv_grid**2)
        grid[ind_cutoff_grid_slow] = np.float('nan')
        us_inv_slow[i,j] = np.nanmean(grid)
        grid = np.sqrt(us_inv_grid**2+vs_inv_grid**2)
        grid[ind_cutoff_grid] = np.float('nan')
        us_inv_med[i,j] = np.nanmean(grid[vsurfini_cutoff_SSA < 4000])
        us_inv[i,j] = np.nanmean(grid)

        # Calculate misfit between observed and modeled surface velocity for inversion
        grid = np.sqrt((us_obs_grid-us_inv_grid)**2+(vs_obs_grid-vs_inv_grid)**2)
        grid[ind_cutoff_grid] = np.float('nan')
        misfit_inv[i,j] = np.nanmean(grid)

        # Get taub from inversion
        x,y,taub_grid = geotifflib.read(DIRs[j]+'DEM'+dates[i]+'_bed_mod_taub.tif')
        grid = np.array(taub_grid) 
        grid[ind_cutoff_grid_fast] = np.float('nan')
        taub_fast[i,j] = np.nanmean(grid)
        grid = np.array(taub_grid)
        grid[ind_cutoff_grid_slow] = np.float('nan')
        taub_slow[i,j] = np.nanmean(grid)
        grid = np.array(taub_grid)
        grid[ind_cutoff_grid] = np.float('nan')
        taub_med[i,j] = np.nanmean(grid[vsurfini_cutoff_SSA < 4000])
        taub[i,j] = np.nanmean(grid)

        # Get beta from inversion
        x,y,beta_grid = geotifflib.read(DIRs[j]+'DEM'+dates[i]+'_bed_mod_beta.tif')
        grid = np.array(beta_grid)
        grid[ind_cutoff_grid_fast] = np.float('nan')
        beta_fast[i,j] = np.nanmean(grid)
        grid = np.array(beta_grid)
        grid[ind_cutoff_grid_slow] = np.float('nan')
        beta_slow[i,j] = np.nanmean(grid)
        grid = np.array(beta_grid)
        grid[ind_cutoff_grid] = np.float('nan')
        beta_med[i,j] = np.nanmean(grid[vsurfini_cutoff_SSA < 4000])
        beta[i,j] = np.nanmean(grid)

        # Get basal velocity from inversion
        x,y,ub_inv_grid = geotifflib.read(DIRs[j]+'DEM'+dates[i]+'_bed_mod_ub.tif')
        x,y,vb_inv_grid = geotifflib.read(DIRs[j]+'DEM'+dates[i]+'_bed_mod_vb.tif')
        grid = np.sqrt(ub_inv_grid**2+vb_inv_grid**2)
        grid[ind_cutoff_grid_fast] = np.float('nan')
        ub_inv_fast[i,j] = np.nanmean(grid)
        grid = np.sqrt(ub_inv_grid**2+vb_inv_grid**2)
        ub_inv_med[i,j] = np.nanmean(grid[vsurfini_cutoff_SSA < 4000])
        grid[ind_cutoff_grid] = np.float('nan')
        ub_inv[i,j] = np.nanmean(grid)

        # Get basal velocity from linear sliding law
        x,y,ub_linear_grid = geotifflib.read(DIRs[j]+'DEM'+dates[i]+'_'+beta_suffix+'_'+temperature_texts[j]+'_linear_bed_mod_ub.tif')
        x,y,vb_linear_grid = geotifflib.read(DIRs[j]+'DEM'+dates[i]+'_'+beta_suffix+'_'+temperature_texts[j]+'_linear_bed_mod_vb.tif')
        grid = np.sqrt(ub_linear_grid**2+vb_linear_grid**2)
        grid[ind_cutoff_grid] = np.float('nan')
        ub_linear[i,j] = np.nanmean(grid)

        # Get surface velocity from linear sliding law
        x,y,us_linear_grid = geotifflib.read(DIRs[j]+'DEM'+dates[i]+'_'+beta_suffix+'_'+temperature_texts[j]+'_linear_surf_mod_us.tif')
        x,y,vs_linear_grid = geotifflib.read(DIRs[j]+'DEM'+dates[i]+'_'+beta_suffix+'_'+temperature_texts[j]+'_linear_surf_mod_vs.tif')
        grid = np.sqrt(us_linear_grid**2+vs_linear_grid**2)
        grid[ind_cutoff_grid] = np.float('nan')
        us_linear[i,j] = np.nanmean(grid)
        
        # Calculate misfit between observed and modeled surface velocity for linear sliding law
        grid = np.sqrt((us_obs_grid-us_linear_grid)**2+(vs_obs_grid-vs_linear_grid)**2)
        grid[ind_cutoff_grid] = np.float('nan')
        misfit_linear[i,j] = np.nanmean(grid)

if glacier == 'Helheim':
    time1 = 2004; time2 = 2015.2
elif glacier == 'Kanger':
    time1 = 2001; time2 = 2016.2
colors=['c','b','orange','r']
fig = plt.figure(figsize=(3.5,4.75))
gs = matplotlib.gridspec.GridSpec(11,1)
matplotlib.rc('font',family='Arial')

ax1 = plt.subplot(gs[8:]); ax1.tick_params(labelsize=8)
ax1.set_ylabel(r'${\tau}$ (kPa)',fontsize=8,fontname='Arial')
for i in range(0,len(modelnames)):
    ax1.plot(times,taub_med[:,i]*1e3,'ko--',markerfacecolor=colors[i],markersize=4,label=modelnames[i])
plt.xlim([time1,time2])
#ax2 = ax1.twinx(); ax2.tick_params(labelsize=8)
ax1.plot(times,taud_med*1e3,'ko--',markerfacecolor='k',markersize=4,label=r'$\tau_d$')
#ax2.set_ylabel(r'$\bar{{\tau}}_d$ (kPa)',fontsize=8,fontname='Arial')
if glacier == 'Kanger':
    ax1.text(2011.9,290,r'Upper glacier',fontname='Arial',fontsize=8)
    plt.xticks(np.arange(2002,2017,4))
    ax1.set_ylim([135,310]); #ax2.set_ylim([160,310])
elif glacier == 'Helheim':
    plt.xticks(np.arange(2004,2016,2))
    ax1.text(2011.95,280,r'Upper glacier',fontname='Arial',fontsize=8)
    ax1.set_ylim([150,330]); #ax1.set_ylim([180,330])
ymin,ymax = plt.ylim()
plt.text(time1+0.01*(time2-time1),ymin+0.03*(ymax-ymin),'(d)',fontsize=8,fontweight='bold')
ax1.legend(loc=4,ncol=3,labelspacing=0.1,columnspacing=0.8,borderpad=0.1,handletextpad=0.5,handlelength=0.1,fontsize=8)

ax1 = plt.subplot(gs[5:8]); ax1.tick_params(labelsize=8)
#ax2 = ax1.twinx(); ax2.tick_params(labelsize=8)
for i in range(0,len(modelnames)):
    ax1.plot(times,taub_fast[:,i]*1e3,'ko--',markerfacecolor=colors[i],markersize=4,label=modelnames[i])
ax1.plot(times,taud_fast*1e3,'ko--',markerfacecolor='k',markersize=4,label=r'$\bar{\tau}_d$')
plt.xlim([time1,time2])
ax1.set_ylabel(r'${\tau}$ (kPa)',fontsize=8,fontname='Arial')
#ax2.set_ylabel(r'$\bar{{\tau}}_d$ (kPa)',fontsize=8,fontname='Arial')
if glacier == 'Kanger':
    plt.xticks(np.arange(2002,2017,4))
    ax1.set_ylim([40,390]); #ax2.set_ylim([90,390])
    ax1.set_yticks(np.arange(50,400,50))
    ax1.text(2011.9,350,r'Lower glacier',fontsize=8,fontname='Arial')
elif glacier == 'Helheim':
    plt.xticks(np.arange(2004,2016,2),fontsize=8)
    ax1.set_ylim([50,330])
    ax1.set_yticks(np.arange(50,350,50))
    #ax2.set_ylim([180,330])
    ax1.text(2011.95,270,r'Lower glacier',fontname='Arial',fontsize=8)
ymin,ymax = plt.ylim()
plt.text(time1+0.01*(time2-time1),ymin+0.03*(ymax-ymin),'(c)',fontsize=8,fontweight='bold')
ax1.set_xticklabels([])

ax1 = plt.subplot(gs[3:5]); ax1.tick_params(labelsize=8)
plt.xlim([time1,time2])
ax1.set_xticklabels([])
plt.ylabel(r'Upper $H_f$ (m)',fontsize=8)
ax2 = ax1.twinx(); ax2.tick_params(labelsize=8)
ax2.set_ylabel(r'Lower $H_f$ (m)',fontsize=8)
ax2.plot(times,Hf_fast,'ko--',markerfacecolor='w',markersize=4,label='Lower glacier')
ax1.plot(times,Hf_med,'ko--',markerfacecolor='k',markersize=4,label='Upper glacier')
ax1.plot(0,0,'ko--',markerfacecolor='w',markersize=4,label='Lower glacier')
if glacier == 'Kanger':
    ax1.set_xticks(np.arange(2002,2017,4));
    ax1.set_yticks(np.arange(500,720,100))
    ax2.set_yticks(np.arange(100,360,100))
    ax1.set_ylim([450,680])
    ax2.set_ylim([130,350])
elif glacier == 'Helheim':
    ax1.set_xticks(np.arange(2004,2016,2))
    ax1.set_ylim([625,675])
    ax1.set_yticks(np.arange(630,680,20))
    ax2.set_yticks(np.arange(320,370,20))
    ax2.set_ylim([320,370])
ymin,ymax = plt.ylim()
plt.text(time1+0.01*(time2-time1),ymin+0.045*(ymax-ymin),'(b)',fontsize=8,fontweight='bold')
ax1.legend(loc=1,labelspacing=0.1,columnspacing=0.8,borderpad=0.1,handletextpad=0.5,handlelength=0.1,fontsize=8)

ax1 = plt.subplot(gs[0:3]); ax1.tick_params(labelsize=8)
plt.plot(times,us_obs_fast/1e3,'ks',markersize=7,markerfacecolor='w')
markers = ['o','^','x','+']
for i in range(0,len(modelnames)):
    plt.plot(times,us_inv_fast[:,i]/1e3,'x',markeredgecolor=colors[i], markerfacecolor='None', marker=markers[i], label=modelnames[i])
    #plt.plot(times,us_linear[:,i]/1e3,'ko',markerfacecolor=colors[i],markersize=4)
plt.plot(2000,2.000,'ks',markerfacecolor='w',label='Observed',markersize=7)
plt.xlim([time1,time2])
if glacier == 'Kanger':
    plt.xticks(np.arange(2002,2017,4))
    plt.ylim([3.7,9.5])
    plt.yticks(np.arange(4,10,1.00))
    plt.text(2011.9,8.85,'Lower glacier',fontname='Arial',fontsize=8)
    plt.legend(loc=4,ncol=3,labelspacing=0.2,columnspacing=0.4,handletextpad=0.3,handlelength=1,fontsize=8)
elif glacier == 'Helheim':
    plt.xticks(np.arange(2004,2016,2)); plt.yticks(np.arange(5,7.5,1))
    plt.ylim([4.5,7.5])
    plt.text(2011.95,6.3,'Lower glacier',fontname='Arial',fontsize=8)
    plt.legend(loc=1,ncol=3,labelspacing=0.2,columnspacing=0.4,handletextpad=0.3,handlelength=1,fontsize=8)
ymin,ymax = plt.ylim()
plt.text(time1+0.01*(time2-time1),ymin+0.03*(ymax-ymin),'(a)',fontsize=8,fontweight='bold')
ax1.set_xticklabels([])
plt.ylabel(r'$u$ (km yr$^{-1}$)',fontsize=8,fontname='Arial')

plt.subplots_adjust(hspace=0.0,wspace=0.0,top=0.99,right=0.86,left=0.14,bottom=0.05)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+'_velocity_taub_timeseries.pdf'),FORMAT='PDF',dpi=300)
print "Saving as "+os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+'_velocity_taub_timeseries.pdf')
plt.close()

# Print stats for long-term behavior for paper
# First print ratios
for j in range(0,len(modelnames)):
    print "upper tau_d/tau_b, "+modelnames[j]+" :", np.round(np.mean(taub_med[:,j]/taud_med),2)
for j in range(0,len(modelnames)):
    print "lower tau_d/tau_b, "+modelnames[j]+" :", np.round(np.mean(taub_fast[:,j]/taud_fast),2)

# Now print changes
if glacier == 'Kanger':
    inds = [[1,2],[7,8]]#[1,3],[2,7],[7,8]]
if glacier == 'Helheim':
    inds = [[0,1],[1,2],[16,17]]
for i in range(0,len(inds)):
    print "Change from "+dates[inds[i][0]]+" to "+dates[inds[i][1]]
    print r'upper u_s : ', np.round(us_obs_med[inds[i][0]]), 'm/yr to', np.round(us_obs_med[inds[i][1]]), \
            'm/yr, change of', np.round((us_obs_med[inds[i][1]]-us_obs_med[inds[i][0]])), 'm/yr'
    print r'upper H_f : ', np.round(Hf_med[inds[i][0]]), 'm to', np.round(Hf_med[inds[i][1]]), \
            'm, change of', np.round((Hf_med[inds[i][1]]-Hf_med[inds[i][0]])), 'm'
    print r'upper tau_d : ', np.round(taud_med[inds[i][0]]*1e3), 'kPa to', np.round(taud_med[inds[i][1]]*1e3), \
           'kPa, change of', np.round((taud_med[inds[i][1]]-taud_med[inds[i][0]])*1e3), 'kPa' 
    for j in range(0,len(modelnames)):
        print r'upper tau_b, '+modelnames[j]+' : ',np.round(taub_med[inds[i][0],j]*1e3), 'kPa to', np.round(taub_med[inds[i][1],j]*1e3), \
            'kPa, change of', np.round((taub_med[inds[i][1],j]-taub_med[inds[i][0],j])*1e3), 'kPa'

    print r'lower u_s : ', np.round(us_obs_fast[inds[i][0]]), 'm/yr to', np.round(us_obs_fast[inds[i][1]]), \
            'm/yr, change of', np.round((us_obs_fast[inds[i][1]]-us_obs_fast[inds[i][0]])), 'm/yr'
    print r'lower H_f : ', np.round(Hf_fast[inds[i][0]]), 'm to', np.round(Hf_fast[inds[i][1]]), \
            'm, change of', np.round((Hf_fast[inds[i][1]]-Hf_fast[inds[i][0]])), 'm'
    print r'lower tau_d : ', np.round(taud_fast[inds[i][0]]*1e3), 'kPa to', np.round(taud_fast[inds[i][1]]*1e3), \
            'kPa, change of', np.round((taud_fast[inds[i][1]]-taud_fast[inds[i][0]])*1e3), 'kPa'
    for j in range(0,len(modelnames)):
        print r'lower tau_b, '+modelnames[j]+' : ',np.round(taub_fast[inds[i][0],j]*1e3), 'kPa to', np.round(taub_fast[inds[i][1],j]*1e3), \
            'kPa, change of', np.round((taub_fast[inds[i][1],j]-taub_fast[inds[i][0],j])*1e3), 'kPa'
    print "\n"

m = 3.0
if glacier == 'Kanger':
    inds = [[1,2],[2,7],[7,8]]
elif glacier == 'Helheim':
    inds = [[0,1],[1,2]]
for i in range(0,len(inds)):
    print dates[inds[i][0]], "to", dates[inds[i][1]], "with thicknesss change of ", H[inds[i][1]] - H[inds[i][0]], "m"
    for j in range(0,4):
        ub_1 = ub_inv[inds[i][0],j]; ub_2 = ub_inv[inds[i][1],j]
        tb_1 = taub[inds[i][0],j]; tb_2 = taub[inds[i][1],j]
        H_1 = H[inds[i][0]]; H_2 = H[inds[i][1]]
        P_i_1 = P_i[inds[i][0]]; P_i_2 = P_i[inds[i][1]]
        P_w_1 = P_w[inds[i][0]]; P_w_2 = P_w[inds[i][1]]
        f_N_weertman = ((tb_1/tb_2)*(ub_2/ub_1)**(1/m))**(1/m)
        f_N_coulomb = tb_2/tb_1
        print modelnames[j]+":","Effective pressure: weertman", np.round((f_N_weertman-1)*100,1),"%, coulomb:",\
                np.round((f_N_coulomb-1)*100,1),"%"
        f_P_weertman = (P_i_1-(P_i_1-P_w_1)*f_N_weertman)/P_w_1
        f_P_coulomb = (P_i_1-(P_i_1-P_w_1)*f_N_coulomb)/P_w_1
        print modelnames[j]+":","P_w: weertman", np.round((f_P_weertman-1)*100,1),"%, coulomb:",\
                np.round((f_P_coulomb-1)*100,1),"%"
        f_H_weertman = ((rho_i*g*H_1-P_w_1)*f_N_weertman + P_w_1)/(rho_i*g)/H_1
        f_H_coulomb = ((rho_i*g*H_1-P_w_1)*f_N_coulomb + P_w_1)/(rho_i*g)/H_1
        print modelnames[j]+":","H: weertman", np.round(f_H_weertman*H_1-H_1), "m, coulomb:",\
                np.round(f_H_coulomb*H_1-H_1), "m"

# Get u_obs for all TSX velocity estimates
x_tsx,y_tsx,grids_tsx,time_tsx = vellib.velocity_grid(glacier,x_cutoff_SSA[0],\
        x_cutoff_SSA[-1],y_cutoff_SSA[0],y_cutoff_SSA[-1],x_cutoff_SSA[1]-x_cutoff_SSA[0])
grids_tsx[ind_cutoff_grid] = np.float('nan')
us_obs_tsx = np.zeros([len(time_tsx),])
for i in range(0,len(time_tsx)):
    us_obs_tsx[i] = np.nanmean(grids_tsx[:,:,i])

if glacier == 'Kanger':
    time1 = 2011.0; time2 = 2014.5
elif glacier == 'Helheim':
    time1 = 2011.0; time2 = 2015.
fig = plt.figure(figsize=(3.5,5))
gs = matplotlib.gridspec.GridSpec(9,1)
matplotlib.rc('font',family='Arial')

ax = plt.subplot(gs[-4:])
ax.tick_params(labelsize=8)
for i in range(0,len(year_runoff)):
    path = matplotlib.path.Path([[day1_runoff[i],100],[day1_runoff[i],350],[day2_runoff[i],350],[day2_runoff[i],100]])
    patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
    ax.add_patch(patch)
for i in range(0,len(modelnames)):
    plt.plot(times,taub[:,i]*1e3,'ko--',\
            markersize=4,markerfacecolor=colors[i],label=modelnames[i])
plt.plot(times,taud*1e3,'ko--',markerfacecolor='k',markersize=4, label=r'$\tau_d$')
plt.ylabel(r'${\tau}$ (kPa)',fontsize=8)
#plt.legend(labelspacing=0.1,handletextpad=1,handlelength=0.1,fontsize=8,ncol=4,loc=2)
plt.ylim([162,315]); plt.yticks(np.arange(175,325,25))
if glacier == 'Kanger':
    plt.legend(labelspacing=0.1,handletextpad=0.5,columnspacing=0.8, handlelength=0.1,fontsize=8,ncol=3)
else:
    plt.legend(labelspacing=0.1,handletextpad=0.5,columnspacing=0.8, handlelength=0.1,fontsize=8,ncol=3, loc='upper left', bbox_to_anchor=(0.0, 0.8))
ax.set_xticks(np.arange(np.floor(time1),1+np.ceil(time2)))
ymin,ymax = plt.ylim()
plt.text(time1+0.01*(time2-time1),ymin+0.03*(ymax-ymin),'(e)',fontsize=8,fontweight='bold')
plt.xlim([time1,time2])

#ax = plt.subplot(gs[-3])
#ax.tick_params(labelsize=8)
#for i in range(0,len(year_runoff)):
#    path = matplotlib.path.Path([[day1_runoff[i],250],[day1_runoff[i],350],[day2_runoff[i],350],[day2_runoff[i],250]])
#    patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
#    ax.add_patch(patch)
#plt.plot(times,taud*1e3,'ko--',markerfacecolor='k',markersize=4)
#plt.ylabel(r'${\tau}_d$ (kPa)',fontsize=8)
#if glacier == 'Kanger':
#    plt.ylim([260,300])
#    plt.yticks(np.arange(275,305,25))
#elif glacier == 'Helheim':
#    plt.ylim([290,329])
#    plt.yticks(np.arange(300,330,25))
#ymin,ymax = plt.ylim()
#plt.text(time1+0.01*(time2-time1),ymin+0.06*(ymax-ymin),'(e)',fontsize=8,fontweight='bold')
#ax.set_xticks(np.arange(np.floor(time1),1+np.ceil(time2))); ax.set_xticklabels([])
#plt.xlim([time1,time2])

ax = plt.subplot(gs[3:5])
ax.tick_params(labelsize=8)
for i in range(0,len(year_runoff)):
    path = matplotlib.path.Path([[day1_runoff[i],1],[day1_runoff[i],5],[day2_runoff[i],5],[day2_runoff[i],1]])
    patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
    ax.add_patch(patch)
for i in range(0,len(modelnames)):
    plt.plot(times,ub_inv[:,i]/1e3,'ko--',\
            markersize=4,markerfacecolor=colors[i],label=modelnames[i])
plt.ylabel(r'${u}_b$ (km yr$^{-1}$)',fontsize=8)
#plt.legend(labelspacing=0.1,handletextpad=0.5,columnspacing=0.8, handlelength=0.1,fontsize=8,ncol=4,loc=0)
if glacier == 'Kanger':
    plt.ylim([1.9,3.9])
    plt.yticks(np.arange(2.5,4,0.5))
elif glacier == 'Helheim':
    plt.ylim([1.8,3.8])
ymin,ymax = plt.ylim()
plt.text(time1+0.01*(time2-time1),ymin+0.03*(ymax-ymin),'(d)',fontsize=8,fontweight='bold')
ax.set_xticks(np.arange(np.floor(time1),1+np.ceil(time2))); ax.set_xticklabels([])
plt.xlim([time1,time2])

ax = plt.subplot(gs[2])
ax.tick_params(labelsize=8)
for i in range(0,len(year_runoff)):
    path = matplotlib.path.Path([[day1_runoff[i],0],[day1_runoff[i],5],[day2_runoff[i],5],[day2_runoff[i],0]])
    patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
    ax.add_patch(patch)
#plt.plot(time_tsx,us_obs_tsx/1e3,'ks--',markerfacecolor='w',markersize=4,label='Observed')
plt.plot(0,0,'ko',markersize=4,label='Inversion')
plt.plot(times,us_obs/1e3,'ko--',markerfacecolor='w',markersize=4)
#for i in range(0,len(modelnames)):
#    plt.plot(times,us_inv/1e3,'ko',markerfacecolor=colors[i],markersize=4)
plt.ylabel(r'${u}$ (km yr$^{-1}$)',fontsize=8)
if glacier == 'Kanger':
    plt.ylim([2.9,3.9])
    plt.yticks(np.arange(3.0,4.0,.5))
elif glacier == 'Helheim':
    plt.ylim([2.7,3.7])
ymin,ymax = plt.ylim()
plt.text(time1+0.01*(time2-time1),ymin+0.06*(ymax-ymin),'(c)',fontsize=8,fontweight='bold')
ax.set_xticks(np.arange(np.floor(time1),1+np.ceil(time2))); ax.set_xticklabels([])
plt.xlim([time1,time2])

ax = plt.subplot(gs[1])
ax.tick_params(labelsize=8)
for i in range(0,len(year_runoff)):
    path = matplotlib.path.Path([[day1_runoff[i],0],[day1_runoff[i],2000],[day2_runoff[i],2000],[day2_runoff[i],0]])
    patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
    ax.add_patch(patch)
plt.plot(times, Hf,'ko--',markersize=4,markerfacecolor='w')
plt.ylabel(r'${H}_f$ (m)',fontsize=8)
if glacier == 'Kanger':
    plt.ylim([440,485]); plt.yticks(np.arange(450,490,25))
elif glacier == 'Helheim':
    plt.ylim([485,505])
ymin,ymax = plt.ylim()
plt.text(time1+0.01*(time2-time1),ymin+0.06*(ymax-ymin),'(b)',fontsize=8,fontweight='bold')
ax.set_xticks(np.arange(np.floor(time1),1+np.ceil(time2))); ax.set_xticklabels([])
plt.xlim([time1,time2])

ax = plt.subplot(gs[0])
ax.tick_params(labelsize=8)
for i in range(0,len(year_runoff)):
    path = matplotlib.path.Path([[day1_runoff[i],-5],[day1_runoff[i],5],[day2_runoff[i],5],[day2_runoff[i],-5]])
    patch = matplotlib.patches.PathPatch(path,facecolor='0.85',edgecolor='none',lw=0)
    ax.add_patch(patch)
plt.plot(term_time,term_pos/1e3,'k.',markersize=5)
plt.ylabel('Terminus\n(km)',fontsize=8, multialignment='center')
ax.set_xticks(np.arange(np.floor(time1),1+np.ceil(time2))); ax.set_xticklabels([])
if glacier == 'Kanger':
    plt.ylim([-4.5,4.5]); plt.yticks(np.arange(-3,4,3))
elif glacier == 'Helheim':
    plt.ylim([-2.75,2.75]); plt.yticks(np.arange(-2,3,2))
ymin,ymax = plt.ylim()
plt.text(time1+0.01*(time2-time1),ymin+0.06*(ymax-ymin),'(a)',fontsize=8,fontweight='bold')
plt.xlim([time1,time2])

plt.subplots_adjust(hspace=0.0,wspace=0.0,top=0.98,right=0.96,left=0.15,bottom=0.05)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+'_seasonal.pdf'),FORMAT='PDF',dpi=300)
print "Saving as "+os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+'_seasonal.pdf')
plt.close()

# Print out seasonal results for paper
ind = np.where((times > 2011) & (times < 2015))[0]
if glacier == 'Kanger':
    inds = [[1,3],[5,6]]
    for i in range(0,len(inds)):
        for j in range(0,len(modelnames)):
            print dates[ind[inds[i][0]]], "to", dates[ind[inds[i][1]]]+',', modelnames[j]+':', \
                    "increase of", np.round((taub[ind[inds[i][1]],j]-taub[ind[inds[i][0]],j])*1e3), "kPa"
            print dates[ind[inds[i][0]]], "to", dates[ind[inds[i][1]]]+',', modelnames[j]+':', \
                    "decrease of", np.round(ub_inv[ind[inds[i][0]],j]-ub_inv[ind[inds[i][1]],j]), "m/yr"
            print dates[ind[inds[i][0]]], "to", dates[ind[inds[i][1]]]+',', modelnames[j]+':', \
                    "u_b/u_s", np.round((ub_inv[ind[inds[i][0]],j]-ub_inv[ind[inds[i][1]],j])/(us_inv[ind[inds[i][0]],j]-us_inv[ind[inds[i][1]],j]),2)

m = 3.0
if glacier == 'Kanger':
    inds = [[1,3],[5,6]]
elif glacier == 'Helheim':
    inds = [[1,2],[5,6]]
for i in range(0,len(inds)):
    print dates[ind[inds[i][0]]],'to',dates[ind[inds[i][1]]]
    for j in range(0,4):
        ub_1 = ub_inv[ind[inds[i][0]],j]; ub_2 = ub_inv[ind[inds[i][1]],j]
        tb_1 = taub[ind[inds[i][0]],j]; tb_2 = taub[ind[inds[i][1]],j]
        H_1 = H[ind[inds[i][0]]]; H_2 = H[ind[inds[i][1]]]
        P_i_1 = P_i[ind[inds[i][0]]]; P_i_2 = P_i[ind[inds[i][1]]]
        P_w_1 = P_w[ind[inds[i][0]]]; P_w_2 = P_w[ind[inds[i][1]]]
        f_N_weertman = ((tb_1/tb_2)*(ub_2/ub_1)**(1/m))**(1/m)
        f_N_coulomb = tb_2/tb_1
        print modelnames[j]+":","Effective pressure: weertman", np.round((f_N_weertman-1)*100,1),"%, coulomb:",\
                np.round((f_N_coulomb-1)*100,1),"%"
        f_P_weertman = (P_i_2-(P_i_1-P_w_1)*f_N_weertman)/P_w_1
        f_P_coulomb = (P_i_2-(P_i_1-P_w_1)*f_N_coulomb)/P_w_1
        print modelnames[j]+":","P_w: weertman", np.round((f_P_weertman-1)*100,1),"%, coulomb:",\
                np.round((f_P_coulomb-1)*100,1),"%"
        f_H_weertman = ((rho_i*g*H_1-P_w_1)*f_N_weertman + P_w_2)/(rho_i*g)/H_1
        f_H_coulomb = ((rho_i*g*H_1-P_w_1)*f_N_coulomb + P_w_2)/(rho_i*g)/H_1
        print modelnames[j]+":","H: weertman", np.round(f_H_weertman*H_1-H_1), "m, coulomb:",\
                np.round(f_H_coulomb*H_1-H_1), "m"

# Plot u_b vs tau_b
if glacier == 'Kanger':
    subpanels = ['(a)','(b)','(c)']
    labels = ['Lower KG','Upper KG','Combined KG']
elif glacier == 'Helheim':
    subpanels = ['(d)','(e)','(f)']
    labels = ['Lower HG','Upper HG','Combined HG']    
fig = plt.figure(figsize=(6.5,2.25))
gs = matplotlib.gridspec.GridSpec(1,3)
matplotlib.rc('font',family='Arial')

ax = plt.subplot(gs[0]); ax.tick_params(labelsize=8); plt.grid()
for i in range(0,len(modelnames)):
    slope,intercept,rvalue,pvalue,stderr = scipy.stats.linregress((ub_inv_fast[:,i]-np.mean(ub_inv_fast[:,i])),(taub_fast[:,i]-np.mean(taub_fast[:,i]))*1e3)
    if (pvalue < 0.05):
        plt.plot((ub_inv_fast[:,i]-np.mean(ub_inv_fast[:,i])),slope*(ub_inv_fast[:,i]-np.mean(ub_inv_fast[:,i]))+intercept,color=colors[i])
    else:
        plt.plot((ub_inv_fast[:,i]-np.mean(ub_inv_fast[:,i])),slope*(ub_inv_fast[:,i]-np.mean(ub_inv_fast[:,i]))+intercept,':',color=colors[i])
    plt.plot((ub_inv_fast[:,i]-np.mean(ub_inv_fast[:,i])),(taub_fast[:,i]-np.mean(taub_fast[:,i]))*1e3,'ko',markerfacecolor=colors[i],markersize=5,zorder=10,label=modelnames[i])
    print modelnames[i],': lower R, pvalue = ',np.round(rvalue,2), np.round(pvalue, 4)
if glacier == 'Kanger':
    plt.legend(loc=0,labelspacing=0.05,columnspacing=0.2,borderpad=0.1,handletextpad=0.3,handlelength=1,fontsize=8)
plt.ylabel(r'$\tau_b$ (kPa)',fontsize=8,fontname='arial')
xmin,xmax = plt.xlim(); ymin,ymax = plt.ylim()
plt.text(xmin+0.02*(xmax-xmin),ymin+0.025*(ymax-ymin),subpanels[0],fontsize=8,fontweight='bold',fontname='Arial')
plt.text(xmin+0.12*(xmax-xmin),ymin+0.025*(ymax-ymin),labels[0],fontsize=8,fontname='arial')

ax = plt.subplot(gs[1]); ax.tick_params(labelsize=8); plt.grid()
for i in range(0,len(modelnames)):
    slope,intercept,rvalue,pvalue,stderr = scipy.stats.linregress((ub_inv_med[:,i]-np.mean(ub_inv_med[:,i])),(taub_med[:,i]-np.mean(taub_med[:,i]))*1e3)
    if (pvalue < 0.05):
        plt.plot((ub_inv_med[:,i]-np.mean(ub_inv_med[:,i])),slope*(ub_inv_med[:,i]-np.mean(ub_inv_med[:,i]))+intercept,colors[i])
    else:
        plt.plot((ub_inv_med[:,i]-np.mean(ub_inv_med[:,i])),slope*(ub_inv_med[:,i]-np.mean(ub_inv_med[:,i]))+intercept,':',color=colors[i])
    plt.plot((ub_inv_med[:,i]-np.mean(ub_inv_med[:,i])),(taub_med[:,i]-np.mean(taub_med[:,i]))*1e3,'ko',markerfacecolor=colors[i],markersize=5,zorder=10)
    print modelnames[i],': upper R, pvalue = ',np.round(rvalue,2), np.round(pvalue,4)
plt.xlabel(r'$u_b$ (m yr$^{-1}$)',fontsize=8,fontname='arial')
xmin,xmax = plt.xlim(); ymin,ymax = plt.ylim()
plt.text(xmin+0.02*(xmax-xmin),ymin+0.025*(ymax-ymin),subpanels[1],fontsize=8,fontweight='bold',fontname='Arial')
plt.text(xmin+0.12*(xmax-xmin),ymin+0.025*(ymax-ymin),labels[1],fontsize=8,fontname='arial')

ax = plt.subplot(gs[2]); ax.tick_params(labelsize=8); plt.grid()
for i in range(0,len(modelnames)):
    slope,intercept,rvalue,pvalue,stderr = scipy.stats.linregress((ub_inv[:,i]-np.mean(ub_inv[:,i])),(taub[:,i]-np.mean(taub[:,i]))*1e3)
    if (pvalue < 0.05):
        plt.plot((ub_inv[:,i]-np.mean(ub_inv[:,i])),slope*(ub_inv[:,i]-np.mean(ub_inv[:,i]))+intercept,color=colors[i])
    else:
        plt.plot((ub_inv[:,i]-np.mean(ub_inv[:,i])),slope*(ub_inv[:,i]-np.mean(ub_inv[:,i]))+intercept,':',color=colors[i])
    plt.plot((ub_inv[:,i]-np.mean(ub_inv[:,i])),(taub[:,i]-np.mean(taub[:,i]))*1e3,'ko',markerfacecolor=colors[i],markersize=5,zorder=10)
    print modelnames[i],': combined R, pvalue = ',np.round(rvalue,2), np.round(pvalue,4)
xmin,xmax = plt.xlim(); ymin,ymax = plt.ylim()
plt.text(xmin+0.02*(xmax-xmin),ymin+0.025*(ymax-ymin),subpanels[2],fontsize=8,fontweight='bold',fontname='Arial')
plt.text(xmin+0.12*(xmax-xmin),ymin+0.025*(ymax-ymin),labels[2],fontsize=8,fontname='arial')

plt.subplots_adjust(hspace=0.2,wspace=0.2,top=0.97,right=0.98,left=0.075,bottom=0.19)
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+'_ub_taub.pdf'),FORMAT='PDF',dpi=300)
print "Saving as "+os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+'_ub_taub.pdf')
plt.close()

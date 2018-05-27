import sys, os, argparse
import scipy, matplotlib
import numpy as np
import matplotlib.pyplot as plt
import elmerreadlib, inverselib, geotifflib

##########
# Inputs #
##########

args = sys.argv

# Get inputs to script
parser = argparse.ArgumentParser()
parser.add_argument("-glacier", dest="glacier",required = True,
                        help = "Glacier name.")

args, _ = parser.parse_known_args(sys.argv)
glacier = args.glacier

# Get velocity cutoff for MAR/MAPE calculations
cutoff = 1000.0

# Get date for steady simulations
if glacier == "Helheim":
    date = '20120316'
elif glacier == "Kanger":
    date = '20120522'

# Model directories
DIR_FS_ModelT = os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/INV_FS_ModelT/DEM"+date+"_modelT_steady/")
DIR_FS_ConstantT = os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/INV_FS_ConstantT/DEM"+date+"_constantT_steady/")
DIR_SSA_ModelT = os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/INV_SSA_ModelT/DEM"+date+"_modelT_steady/")
DIR_SSA_ConstantT = os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/INV_SSA_ConstantT/DEM"+date+"_constantT_steady/")

modnames = ['FS-CT','FS-MT','SSA-CT','SSA-MT']
DIRs = [DIR_FS_ConstantT,DIR_FS_ModelT,DIR_SSA_ConstantT,DIR_SSA_ModelT]

# Inversion filenames for each model
if glacier == 'Kanger':
    adjointfiles = [DIR_FS_ConstantT+"mesh2d/inversion_adjoint/lambda_1e12_20180415/adjoint_beta0001.pvtu",\
                    DIR_FS_ModelT+"mesh2d/inversion_adjoint/lambda_1e12_20180413/adjoint_beta0001.pvtu",\
                    DIR_SSA_ConstantT+"mesh2d/inversion_adjoint/lambda_1e13_20180413/adjoint_beta_ssa0001.pvtu",\
                    DIR_SSA_ModelT+"mesh2d/inversion_adjoint/lambda_1e13_20180413/adjoint_beta_ssa0001.pvtu"]
elif glacier == 'Helheim':
    adjointfiles = [DIR_FS_ConstantT+"mesh2d/inversion_adjoint/lambda_1e12_20180319/adjoint_beta0001.pvtu",\
                    DIR_FS_ModelT+"mesh2d/inversion_adjoint/lambda_1e12_20180316/adjoint_beta0001.pvtu",\
                    DIR_SSA_ConstantT+"mesh2d/inversion_adjoint/lambda_1e13_20180226/adjoint_beta_ssa0001.pvtu",\
                    DIR_SSA_ModelT+"mesh2d/inversion_adjoint/lambda_1e13_20180226/adjoint_beta_ssa0001.pvtu"]

# Get mesh boundary
extent = np.loadtxt(DIR_FS_ModelT+"inputs/mesh_extent.dat")
try:
    hole1 = np.loadtxt(DIR_FS_ModelT+"inputs/mesh_hole1.dat")
    hole2 = np.loadtxt(DIR_FS_ModelT+"inputs/mesh_hole2.dat")
    holes=[hole1,hole2]
except:
    holes = []

# Plot options
vmin = -75 # min residual
vmax = 75 # max residual
cmap = matplotlib.colors.LinearSegmentedColormap.from_list('map',['#400000','#800000','#a50026','#d73027','#f46d43','#fdae61',\
    '#ffffff','#abd9e9','#74add1','#4575b4','#313695','#191970','#000059'])

# Get indices where velocity is always greater than the cutoff value
# For SSA-CT and SSA-MT
x_cutoff,y_cutoff,vsurfini_cutoff,ind_cutoff_grid,ind_cutoff_SSA  = inverselib.get_velocity_cutoff(glacier,velocity_cutoff=cutoff,\
        SSA=True,model_dir='INV_SSA_ModelT')
# For FS-CT and FS-MT (no indices may vary from SSA-CT and SSA-MT)
x_cutoff,y_cutoff,vsurfini_cutoff,ind_cutoff_grid,ind_cutoff_FS  = inverselib.get_velocity_cutoff(glacier,velocity_cutoff=cutoff,\
        SSA=False,model_dir='INV_FS_ModelT')
ind_cutoffs = [ind_cutoff_FS,ind_cutoff_FS,ind_cutoff_SSA,ind_cutoff_SSA]

# Set up figure
if glacier == 'Kanger':
    height = 6.0
    cbar_z = 0.065
    cbar_height = 0.02
    top_z = 0.96
    bot_z = 0.095
elif glacier == 'Helheim':
    height = 8.2
    cbar_z = 0.055
    cbar_height = 0.014
    top_z = 0.97
    bot_z = 0.08

fig = plt.figure(figsize=(6.5,height))
matplotlib.rc('font',family='Arial')
gs = matplotlib.gridspec.GridSpec(4,5)

# Plot mesh extent and holes as block onlines
for i in range(0,4):
    for j in range(0,5):
        plt.subplot(gs[i,j])
        plt.xticks([])
        plt.yticks([])
        plt.plot(np.r_[extent[:,0],extent[0,0]],np.r_[extent[:,1],extent[0,1]],lw=1.5,color='k',zorder=1)
        if len(holes) > 0:
            for hole in holes:
                plt.plot(np.r_[hole[:,0],hole[0,0]],np.r_[hole[:,1],hole[0,1]],lw=1.5,color='k',zorder=1)

# Now cycle through model directories
column_labels = ['a','b','c','d','e']
for i in range(0,len(DIRs)):
    DIR = DIRs[i]
    ind_cutoff = ind_cutoffs[i]

    # Get filenames for pvtu files
    files = os.listdir(DIR+"mesh2d")
    for file in files:
        if file.endswith('.pvtu'):
            if ("FS" in file) and ("modelT" in file):
                file_fs_mt = DIR+"mesh2d/"+file
            elif ("FS" in file) and ("constantT" in file):
                file_fs_ct = DIR+"mesh2d/"+file
            elif ("SSA" in file) and ("modelT" in file):
                file_ssa_mt = DIR+"mesh2d/"+file
            elif ("SSA" in file) and ("constantT" in file):
                file_ssa_ct = DIR+"mesh2d/"+file

    # Get pvtu files
    if 'FS' in DIR:
        velocityname = 'velocity'
    else:
        velocityname = 'ssavelocity'
    fs_mt = elmerreadlib.pvtu_file(file_fs_mt,[velocityname,'vsurfini'])
    fs_ct = elmerreadlib.pvtu_file(file_fs_ct,[velocityname,'vsurfini'])
    ssa_mt = elmerreadlib.pvtu_file(file_ssa_mt,[velocityname,'vsurfini'])
    ssa_ct = elmerreadlib.pvtu_file(file_ssa_ct,[velocityname,'vsurfini'])

    # Get values at the surface for this forward model
    surf_fs_mt = elmerreadlib.values_in_layer(fs_mt,'surface')
    surf_fs_ct = elmerreadlib.values_in_layer(fs_ct,'surface')
    surf_ssa_mt = elmerreadlib.values_in_layer(ssa_mt,'surface')
    surf_ssa_ct = elmerreadlib.values_in_layer(ssa_ct,'surface')

    # Get values at the bed for this inversion
    ax = plt.subplot(gs[i,0])
    adjoint = elmerreadlib.pvtu_file(adjointfiles[i],[velocityname,'beta'])
    bed = elmerreadlib.values_in_layer(adjoint,'bed')
    x,y,taub = elmerreadlib.grid3d(bed,'taub',holes,extent)

    grid = np.array(taub.data)
    grid[ind_cutoff_grid] = np.float('nan')
    ind = np.where(grid*1e3 <= 10)
    nonnan = np.where(~(np.isnan(grid)))
    print DIR
    print np.round(float(len(ind[0]))/len(nonnan[0])*100)

    # Plot taub
    im1 = plt.imshow(taub*1e3,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',vmin=0,vmax=400)
    #plt.contour(x,y,taub*1e3,[20],colors='w',linewidth=1)#,linestyles='dashed')
    plt.ylabel(r'$\bf{('+str(i+1)+') }$ '+modnames[i],fontname='Arial',fontsize=8)
    if glacier == 'Kanger':
        bot = 0.76
        label_bot = 0.11
    if glacier == 'Helheim':
        bot = 0.82
        label_bot = 0.08
    xmin,xmax = plt.xlim()
    ymin,ymax = plt.ylim()
    if i == 0:
        path = matplotlib.path.Path([[0.67*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                [0.97*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                [0.97*(xmax-xmin)+xmin,bot*(ymax-ymin)+ymin],
                [0.67*(xmax-xmin)+xmin,bot*(ymax-ymin)+ymin],
                [0.67*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin]])
        patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1,zorder=2)
        ax.add_patch(patch)
        if glacier == 'Kanger':
            ax.plot([xmin+0.71*(xmax-xmin),xmin+0.71*(xmax-xmin)+5e3],\
                [ymax-0.07*(ymax-ymin),ymax-0.07*(ymax-ymin)],'k')
            ax.plot([xmin+0.71*(xmax-xmin),xmin+0.71*(xmax-xmin)],\
                [ymax-0.07*(ymax-ymin),ymax-0.11*(ymax-ymin)],'k')
            ax.plot([xmin+0.71*(xmax-xmin)+5e3,xmin+0.71*(xmax-xmin)+5e3],\
                [ymax-0.07*(ymax-ymin),ymax-0.11*(ymax-ymin)],'k')
            ax.text(xmin+0.7*(xmax-xmin),ymax-0.21*(ymax-ymin),'5 km',fontsize=8)
        elif glacier == 'Helheim':
            ax.plot([xmin+0.71*(xmax-xmin),xmin+0.71*(xmax-xmin)+5e3],\
                [ymax-0.06*(ymax-ymin),ymax-0.06*(ymax-ymin)],'k')
            ax.plot([xmin+0.71*(xmax-xmin),xmin+0.71*(xmax-xmin)],\
                [ymax-0.06*(ymax-ymin),ymax-0.09*(ymax-ymin)],'k')
            ax.plot([xmin+0.71*(xmax-xmin)+5e3,xmin+0.71*(xmax-xmin)+5e3],\
                [ymax-0.06*(ymax-ymin),ymax-0.08*(ymax-ymin)],'k')
            ax.text(xmin+0.7*(xmax-xmin),ymax-0.16*(ymax-ymin),'5 km',fontsize=8)
    ax.text(xmin+0.04*(xmax-xmin),ymax-label_bot*(ymax-ymin),'(a'+str(i+1)+')',fontsize=8,\
            fontweight='bold')
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    
    # Plot residual for FS-CT friction coefficient
    ax = plt.subplot(gs[0,i+1])
    x,y,velocity = elmerreadlib.grid3d(surf_fs_ct,velocityname,holes,extent)
    x,y,vsurfini = elmerreadlib.grid3d(surf_fs_ct,'vsurfini',holes,extent)
    im2 = plt.imshow(100*(velocity-vsurfini)/vsurfini,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',\
            cmap=cmap.reversed(),vmin=vmin,vmax=vmax)
    plt.contour(x_cutoff,y_cutoff,vsurfini_cutoff,[cutoff],colors='k',linestyles='dashed',linewidths=1.5)
    MAPE = 100*np.mean(abs((surf_fs_ct[velocityname][ind_cutoff]-surf_fs_ct['vsurfini'][ind_cutoff])\
            /surf_fs_ct['vsurfini'][ind_cutoff]))
    xmin,xmax = plt.xlim()
    ymin,ymax = plt.ylim()
    path = matplotlib.path.Path([[0.67*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.97*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.97*(xmax-xmin)+xmin,bot*(ymax-ymin)+ymin],
                        [0.67*(xmax-xmin)+xmin,bot*(ymax-ymin)+ymin],
                        [0.67*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin]])
    patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1,zorder=2)
    ax.add_patch(patch)
    if glacier == 'Kanger':
        plt.text(xmax-0.28*(xmax-xmin),ymax-0.11*(ymax-ymin),'MAR',fontsize=8,zorder=2)
        plt.text(xmax-0.31*(xmax-xmin),ymax-0.21*(ymax-ymin),'{0:>4.1f}%'.format(MAPE),fontsize=8,zorder=2)
    elif glacier == 'Helheim':
        plt.text(xmax-0.28*(xmax-xmin),ymax-0.09*(ymax-ymin),'MAR',fontsize=8,zorder=2)
        plt.text(xmax-0.31*(xmax-xmin),ymax-0.16*(ymax-ymin),'{0:>4.1f}%'.format(MAPE),fontsize=8,zorder=2)
    ax.text(xmin+0.04*(xmax-xmin),ymax-label_bot*(ymax-ymin),'(b'+str(i+1)+')',fontsize=8,\
            fontweight='bold')

    # Plot residual for FS-MT friction coefficient
    ax = plt.subplot(gs[1,i+1])
    x,y,velocity = elmerreadlib.grid3d(surf_fs_mt,velocityname,holes,extent)
    x,y,vsurfini = elmerreadlib.grid3d(surf_fs_mt,'vsurfini',holes,extent)
    plt.imshow(100*(velocity-vsurfini)/vsurfini,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',\
            cmap=cmap.reversed(),vmin=vmin,vmax=vmax)
    plt.contour(x_cutoff,y_cutoff,vsurfini_cutoff,[cutoff],colors='k',linestyles='dashed',linewidths=1.5)
    MAPE = 100*np.mean(abs((surf_fs_mt[velocityname][ind_cutoff]-surf_fs_mt['vsurfini'][ind_cutoff])\
            /surf_fs_mt['vsurfini'][ind_cutoff]))
    path = matplotlib.path.Path([[0.67*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.97*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.97*(xmax-xmin)+xmin,bot*(ymax-ymin)+ymin],
                        [0.67*(xmax-xmin)+xmin,bot*(ymax-ymin)+ymin],
                        [0.67*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin]])
    patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1,zorder=2)
    ax.add_patch(patch)
    if glacier == 'Kanger':
        plt.text(xmax-0.28*(xmax-xmin),ymax-0.11*(ymax-ymin),'MAR',fontsize=8,zorder=2)
        plt.text(xmax-0.31*(xmax-xmin),ymax-0.21*(ymax-ymin),'{0:>4.1f}%'.format(MAPE),fontsize=8,zorder=2)
    elif glacier == 'Helheim':
        plt.text(xmax-0.28*(xmax-xmin),ymax-0.09*(ymax-ymin),'MAR',fontsize=8,zorder=2)
        plt.text(xmax-0.31*(xmax-xmin),ymax-0.16*(ymax-ymin),'{0:>4.1f}%'.format(MAPE),fontsize=8,zorder=2)
    ax.text(xmin+0.04*(xmax-xmin),ymax-label_bot*(ymax-ymin),'(c'+str(i+1)+')',fontsize=8,\
            fontweight='bold')

    # Plot residual for SSA-CT friction coefficient
    ax = plt.subplot(gs[2,i+1])
    x,y,velocity = elmerreadlib.grid3d(surf_ssa_ct,velocityname,holes,extent)
    x,y,vsurfini = elmerreadlib.grid3d(surf_ssa_ct,'vsurfini',holes,extent)
    plt.imshow(100*(velocity-vsurfini)/vsurfini,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',\
            cmap=cmap.reversed(),vmin=vmin,vmax=vmax)
    plt.contour(x_cutoff,y_cutoff,vsurfini_cutoff,[cutoff],colors='k',linestyles='dashed',linewidths=1.5)
    MAPE = 100*np.mean(abs((surf_ssa_ct[velocityname][ind_cutoff]-surf_ssa_ct['vsurfini'][ind_cutoff])\
            /surf_ssa_ct['vsurfini'][ind_cutoff]))
    path = matplotlib.path.Path([[0.67*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.97*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.97*(xmax-xmin)+xmin,bot*(ymax-ymin)+ymin],
                        [0.67*(xmax-xmin)+xmin,bot*(ymax-ymin)+ymin],
                        [0.67*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin]])
    patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1,zorder=2)
    ax.add_patch(patch)
    if glacier == 'Kanger':
        plt.text(xmax-0.28*(xmax-xmin),ymax-0.11*(ymax-ymin),'MAR',fontsize=8,zorder=2)
        plt.text(xmax-0.31*(xmax-xmin),ymax-0.21*(ymax-ymin),'{0:>4.1f}%'.format(MAPE),fontsize=8,zorder=2)
    elif glacier == 'Helheim':
        plt.text(xmax-0.28*(xmax-xmin),ymax-0.09*(ymax-ymin),'MAR',fontsize=8,zorder=2)
        plt.text(xmax-0.31*(xmax-xmin),ymax-0.16*(ymax-ymin),'{0:>4.1f}%'.format(MAPE),fontsize=8,zorder=2)
    ax.text(xmin+0.04*(xmax-xmin),ymax-label_bot*(ymax-ymin),'(d'+str(i+1)+')',fontsize=8,\
            fontweight='bold')

    # Plot residual for SSA-MT friction coefficient
    ax = plt.subplot(gs[3,i+1])
    x,y,velocity = elmerreadlib.grid3d(surf_ssa_mt,velocityname,holes,extent)
    x,y,vsurfini = elmerreadlib.grid3d(surf_ssa_mt,'vsurfini',holes,extent)
    plt.imshow(100*(velocity-vsurfini)/vsurfini,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',\
            cmap=cmap.reversed(),vmin=vmin,vmax=vmax)
    plt.contour(x_cutoff,y_cutoff,vsurfini_cutoff,[cutoff],colors='k',linestyles='dashed',linewidths=1.5)
    MAPE = 100*np.mean(abs((surf_ssa_mt[velocityname][ind_cutoff]-surf_ssa_mt['vsurfini'][ind_cutoff])\
            /surf_ssa_mt['vsurfini'][ind_cutoff]))
    path = matplotlib.path.Path([[0.67*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.97*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin],
                        [0.97*(xmax-xmin)+xmin,bot*(ymax-ymin)+ymin],
                        [0.67*(xmax-xmin)+xmin,bot*(ymax-ymin)+ymin],
                        [0.67*(xmax-xmin)+xmin,0.98*(ymax-ymin)+ymin]])
    patch = matplotlib.patches.PathPatch(path,edgecolor='k',facecolor='w',lw=1,zorder=2)
    ax.add_patch(patch)
    if glacier == 'Kanger':
        plt.text(xmax-0.28*(xmax-xmin),ymax-0.11*(ymax-ymin),'MAR',fontsize=8,zorder=2)
        plt.text(xmax-0.31*(xmax-xmin),ymax-0.21*(ymax-ymin),'{0:>4.1f}%'.format(MAPE),fontsize=8,zorder=2)
    elif glacier == 'Helheim':
        plt.text(xmax-0.28*(xmax-xmin),ymax-0.09*(ymax-ymin),'MAR',fontsize=8,zorder=2)
        plt.text(xmax-0.31*(xmax-xmin),ymax-0.16*(ymax-ymin),'{0:>4.1f}%'.format(MAPE),fontsize=8,zorder=2)
    ax.text(xmin+0.04*(xmax-xmin),ymax-label_bot*(ymax-ymin),'(e'+str(i+1)+')',fontsize=8,\
            fontweight='bold')

    if i == 0:
        plt.subplot(gs[i,0])
        plt.title(r'$\bf{(a)}$ ${\tau_b}$',fontsize=8,fontname='Arial')
        for j in range(1,5):
            plt.subplot(gs[i,j])
            plt.title(r'$\bf{('+column_labels[j]+') }$ '+modnames[j-1],fontsize=8,fontname='Arial')

# Set up colorbars
cbar_ax1 = fig.add_axes([0.03,cbar_z,0.19,cbar_height])
cb1 = fig.colorbar(im1, cax=cbar_ax1,orientation='horizontal',extend='max',ticks=np.arange(0,410,100))
cb1.set_label(r'$\tau_b$ (kPa)',fontsize=8,fontname='Arial')
cbar_ax2 = fig.add_axes([0.415,cbar_z,0.385,cbar_height])
cb2 = fig.colorbar(im2, cax=cbar_ax2,orientation='horizontal',extend='both')
cb2.set_label('Residual (%)',fontname='Arial',fontsize=8)

plt.tight_layout()
plt.subplots_adjust(hspace=0.0,wspace=0.0,top=top_z,right=0.99,left=0.03,bottom=bot_z)

plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_beta_simulations.pdf"),format="PDF",dpi=300)
plt.close()

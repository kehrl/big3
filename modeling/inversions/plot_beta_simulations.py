import sys, os, argparse
import scipy, matplotlib
import numpy as np
import matplotlib.pyplot as plt
import elmerreadlib, geotifflib

##########
# Inputs #
##########

args = sys.argv

# Get inputs to file
parser = argparse.ArgumentParser()
parser.add_argument("-glacier", dest="glacier",required = True,
                        help = "Glacier name.")

args, _ = parser.parse_known_args(sys.argv)
glacier = args.glacier

# Get date for steady simulations

if glacier == "Helheim":
    date = '20120316'
    ximage,yimage,image = geotifflib.readrgb(os.path.join(os.getenv("DATA_HOME"),"Imagery/Landsat/Helheim/"+\
        "TIF/20140704140535_LC82330132014185LGN00.tif"))
elif glacier == "Kanger":
    date = '20120213'
    ximage,yimage,image = geotifflib.readrgb(os.path.join(os.getenv("DATA_HOME"),"Imagery/Landsat/Kanger/"+\
        "TIF/20140706135251_LC82310122014187LGN00.tif"))

# Directories
DIR_FS_ModelT = os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/INV_FS_ModelT/DEM"+date+"_modelT_steady/")
DIR_FS_ConstantT = os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/INV_FS_ConstantT/DEM"+date+"_constantT_steady/")
DIR_SSA_ModelT = os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/INV_SSA_ModelT/DEM"+date+"_modelT_steady/")
DIR_SSA_ConstantT = os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/INV_SSA_ConstantT/DEM"+date+"_constantT_steady/")

modnames = ['FS-CT','FS-MT','SSA-CT','SSA-MT']
DIRs = [DIR_FS_ConstantT,DIR_FS_ModelT,DIR_SSA_ConstantT,DIR_SSA_ModelT]
if glacier == 'Kanger':
    adjointfiles = [DIR_FS_ConstantT+"mesh2d/inversion_adjoint/lambda_1e11_20180208/adjoint_beta0001.pvtu",\
                    DIR_FS_ModelT+"mesh2d/inversion_adjoint/lambda_1e11_20180208/adjoint_beta0001.pvtu",\
                    DIR_SSA_ConstantT+"mesh2d/inversion_adjoint/lambda_1e9_20180209/adjoint_beta_ssa0001.pvtu",\
                    DIR_SSA_ModelT+"mesh2d/inversion_adjoint/lambda_1e9_20180209/adjoint_beta_ssa0001.pvtu"]
elif glacier == 'Helheim':
    adjointfiles = [DIR_FS_ConstantT+"mesh2d/inversion_adjoint/lambda_1e11_20170915/adjoint_beta0001.pvtu",\
                    DIR_FS_ModelT+"mesh2d/inversion_adjoint/lambda_1e11_20170912/adjoint_beta0001.pvtu",\
                    DIR_SSA_ConstantT+"mesh2d/inversion_adjoint/lambda_1e9_20170917/adjoint_beta_ssa0001.pvtu",\
                    DIR_SSA_ModelT+"mesh2d/inversion_adjoint/lambda_1e9_20180112/adjoint_beta_ssa0001.pvtu"]

# Get mesh boundary
extent = np.loadtxt(DIR_FS_ModelT+"inputs/mesh_extent.dat")
try:
    hole1 = np.loadtxt(DIR_FS_ModelT+"inputs/mesh_hole1.dat")
    hole2 = np.loadtxt(DIR_FS_ModelT+"inputs/mesh_hole2.dat")
    holes=[hole1,hole2]
except:
    holes = []


# Plot options
vmin = -2000
vmax = 2000
cmap = 'RdYlBu_r'

# Set up figure
if glacier == 'Kanger':
    height = 7
    cbar_z = 0.07
    cbar_height = 0.03
    top_z = 0.96
    bot_z = 0.11
elif glacier == 'Helheim':
    height = 9.55
    cbar_z = 0.05
    cbar_height = 0.022
    top_z = 0.97
    bot_z = 0.08
fig = plt.figure(figsize=(7.5,height))
matplotlib.rc('font',family='Arial')
gs = matplotlib.gridspec.GridSpec(4,5)

for i in range(0,4):
    for j in range(0,5):
        plt.subplot(gs[i,j])
        plt.xticks([])
        plt.yticks([])
        plt.imshow(image[:,:,0],extent=[ximage[0],ximage[-1],yimage[0],yimage[-1]],\
                cmap='Greys_r',origin='lower',clim=[0,0.6])
        plt.plot(np.r_[extent[:,0],extent[0,0]],np.r_[extent[:,1],extent[0,1]],lw=1.5,color='k')

for i in range(0,len(DIRs)):
    DIR = DIRs[i]

    # Get filenames for pvtu files
    files = os.listdir(DIR+"mesh2d")
    n_fs_mt = 0; n_fs_ct = 0; n_ssa_mt = 0; n_ssa_ct = 0;
    for file in files:
        if file.endswith('.pvtu'):
            if ("fs_mt" in file) and (int(file[-9:-5]) > n_fs_mt):
                n_fs_mt = int(file[-9:-5])
                file_fs_mt = DIR+"mesh2d/"+file
            elif ("fs_ct" in file) and (int(file[-9:-5]) > n_fs_ct):
                n_fs_ct = int(file[-9:-5])
                file_fs_ct = DIR+"mesh2d/"+file
            elif ("ssa_mt" in file) and (int(file[-9:-5]) > n_ssa_mt):
                n_ssa_mt = int(file[-9:-5])
                file_ssa_mt = DIR+"mesh2d/"+file
            elif ("ssa_ct" in file) and (int(file[-9:-5]) > n_ssa_ct):
                n_ssa_ct = int(file[-9:-5])
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

    # Get values at the surface
    surf_fs_mt = elmerreadlib.values_in_layer(fs_mt,'surface')
    surf_fs_ct = elmerreadlib.values_in_layer(fs_ct,'surface')
    surf_ssa_mt = elmerreadlib.values_in_layer(ssa_mt,'surface')
    surf_ssa_ct = elmerreadlib.values_in_layer(ssa_ct,'surface')

    # Get taub for this model
    
    plt.subplot(gs[i,0])
    adjoint = elmerreadlib.pvtu_file(adjointfiles[i],[velocityname,'beta'])
    bed = elmerreadlib.values_in_layer(adjoint,'bed')
    x,y,taub = elmerreadlib.grid3d(bed,'taub',holes,extent)
    im1 = plt.imshow(taub*1e3,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',vmin=0,vmax=500)
    plt.ylabel(modnames[i],fontname='Arial',fontsize=10)

    plt.subplot(gs[i,1])
    x,y,velocity = elmerreadlib.grid3d(surf_fs_ct,velocityname,holes,extent)
    x,y,vsurfini = elmerreadlib.grid3d(surf_fs_ct,'vsurfini',holes,extent)
    im2 = plt.imshow(velocity-vsurfini,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',cmap=cmap,\
            vmin=vmin,vmax=vmax,norm=matplotlib.colors.SymLogNorm(100))

    plt.subplot(gs[i,2])
    x,y,velocity = elmerreadlib.grid3d(surf_fs_mt,velocityname,holes,extent)
    x,y,vsurfini = elmerreadlib.grid3d(surf_fs_mt,'vsurfini',holes,extent)
    plt.imshow(velocity-vsurfini,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',cmap=cmap,\
            vmin=vmin,vmax=vmax,norm=matplotlib.colors.SymLogNorm(100))

    plt.subplot(gs[i,3])
    x,y,velocity = elmerreadlib.grid3d(surf_ssa_ct,velocityname,holes,extent)
    x,y,vsurfini = elmerreadlib.grid3d(surf_ssa_ct,'vsurfini',holes,extent)
    plt.imshow(velocity-vsurfini,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',cmap=cmap,\
            vmin=vmin,vmax=vmax,norm=matplotlib.colors.SymLogNorm(100))

    plt.subplot(gs[i,4])
    x,y,velocity = elmerreadlib.grid3d(surf_ssa_mt,velocityname,holes,extent)
    x,y,vsurfini = elmerreadlib.grid3d(surf_ssa_mt,'vsurfini',holes,extent)
    plt.imshow(velocity-vsurfini,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',cmap=cmap,\
            vmin=vmin,vmax=vmax,norm=matplotlib.colors.SymLogNorm(100))

    if i == 0:
        for j in range(1,5):
            plt.subplot(gs[i,j])
            plt.title(modnames[j-1],fontsize=10,fontname='Arial')

cbar_ax1 = fig.add_axes([0.03,cbar_z,0.19,cbar_height])
cb1 = fig.colorbar(im1, cax=cbar_ax1,orientation='horizontal')
cb1.set_label(r'$\tau_b$ (kPa)',fontsize=10,fontname='Arial')
cbar_ax2 = fig.add_axes([0.415,cbar_z,0.385,cbar_height])
cb2 = fig.colorbar(im2, cax=cbar_ax2,orientation='horizontal')
cb2.set_label('Modeled-Measured (m/yr)',fontname='Arial',fontsize=10)

plt.tight_layout()
plt.subplots_adjust(hspace=0.01,wspace=0.01,top=top_z,right=0.99,left=0.03,bottom=bot_z)

plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_beta_simulations.pdf"),format="PDF",dpi=600)
plt.close()

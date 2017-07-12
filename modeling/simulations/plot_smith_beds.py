import elmerreadlib, glaclib, vellib, zslib, datelib, icefrontlib, matplotlib
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import numpy as np
import os
import matplotlib.path
from mpl_toolkits.axes_grid1 import AxesGrid

glacier = 'Helheim'

# Directories
DIR_bed3 = os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/TD_20110615_20120615_bsmooth3/")
DIR_bed4 = os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/TD_20110615_20120615_bsmooth4/")
DIR_bed5 = os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/TD_20110615_20120615_bsmooth5/")
DIR_bed8 = os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/TD_20110615_20120615_bsmooth8/")

# Load beds
x,y,zbed3 = elmerreadlib.input_file(DIR_bed3+'inputs/bedrock.xy')
x,y,zbed4 = elmerreadlib.input_file(DIR_bed4+'inputs/bedrock.xy')
x,y,zbed5 = elmerreadlib.input_file(DIR_bed5+'inputs/bedrock.xy')
x,y,zbed8 = elmerreadlib.input_file(DIR_bed8+'inputs/bedrock.xy')

xgrid,ygrid = np.meshgrid(x,y)

# Load mesh extents for cropping
mesh_extent = np.loadtxt(DIR_bed3+'inputs/mesh_extent.dat')
mesh_hole1 = np.loadtxt(DIR_bed3+'inputs/mesh_hole1.dat')
mesh_hole2 = np.loadtxt(DIR_bed3+'inputs/mesh_hole2.dat')

path = matplotlib.path.Path(np.column_stack([mesh_extent[:,0],mesh_extent[:,1]]))
inmesh = path.contains_points(np.column_stack([xgrid.flatten(),ygrid.flatten()]))
inmesh = inmesh.reshape(len(y),len(x))

path = matplotlib.path.Path(np.column_stack([mesh_hole1[:,0],mesh_hole1[:,1]]))
inhole1 = path.contains_points(np.column_stack([xgrid.flatten(),ygrid.flatten()]))
inhole1 = inhole1.reshape(len(y),len(x))
path = matplotlib.path.Path(np.column_stack([mesh_hole2[:,0],mesh_hole2[:,1]]))
inhole2 = path.contains_points(np.column_stack([xgrid.flatten(),ygrid.flatten()]))
inhole2 = inhole2.reshape(len(y),len(x))

del path

zbed3[~inmesh] = float('nan')
zbed3[inhole1] = float('nan')
zbed3[inhole2] = float('nan')
zbed4[~inmesh] = float('nan')
zbed4[inhole1] = float('nan')
zbed4[inhole2] = float('nan')
zbed5[~inmesh] = float('nan')
zbed5[inhole1] = float('nan')
zbed5[inhole2] = float('nan')
zbed8[~inmesh] = float('nan')
zbed8[inhole1] = float('nan')
zbed8[inhole2] = float('nan')

# Load bsmooth3
bed3_0 = elmerreadlib.pvtu_file(DIR_bed3+'mesh2d/terminusdriven0001.pvtu',['velocity','zs top'])
bed3_0 = bed3_0[bed3_0['zs top'] != 0]
bed3_1 = elmerreadlib.pvtu_file(DIR_bed3+'mesh2d/terminusdriven0367.pvtu',['velocity','zs top'])
bed3_1 = bed3_1[bed3_1['zs top'] != 0]

bed3_z0 = griddata((bed3_0['x'],bed3_0['y']),bed3_0['z'],(xgrid,ygrid))
bed3_z1 = griddata((bed3_1['x'],bed3_1['y']),bed3_1['z'],(xgrid,ygrid))

bed3_z0[~inmesh] = float('nan')
bed3_z0[inhole1] = float('nan')
bed3_z0[inhole2] = float('nan')
bed3_z1[~inmesh] = float('nan')
bed3_z1[inhole1] = float('nan')
bed3_z1[inhole2] = float('nan')

# Load bsmooth4
bed4_0 = elmerreadlib.pvtu_file(DIR_bed4+'mesh2d/terminusdriven0001.pvtu',['velocity','zs top'])
bed4_0 = bed4_0[bed4_0['zs top'] != 0]
bed4_1 = elmerreadlib.pvtu_file(DIR_bed4+'mesh2d/terminusdriven0367.pvtu',['velocity','zs top'])
bed4_1 = bed4_1[bed4_1['zs top'] != 0]

bed4_z0 = griddata((bed4_0['x'],bed4_0['y']),bed4_0['z'],(xgrid,ygrid))
bed4_z1 = griddata((bed4_1['x'],bed4_1['y']),bed4_1['z'],(xgrid,ygrid))

bed4_z0[~inmesh] = float('nan')
bed4_z0[inhole1] = float('nan')
bed4_z0[inhole2] = float('nan')
bed4_z1[~inmesh] = float('nan')
bed4_z1[inhole1] = float('nan')
bed4_z1[inhole2] = float('nan')

# Load bsmooth5
bed5_0 = elmerreadlib.pvtu_file(DIR_bed5+'mesh2d/terminusdriven0001.pvtu',['velocity','zs top'])
bed5_0 = bed5_0[bed5_0['zs top'] != 0]
bed5_1 = elmerreadlib.pvtu_file(DIR_bed5+'mesh2d/terminusdriven0367.pvtu',['velocity','zs top'])
bed5_1 = bed5_1[bed5_1['zs top'] != 0]

bed5_z0 = griddata((bed5_0['x'],bed5_0['y']),bed5_0['z'],(xgrid,ygrid))
bed5_z1 = griddata((bed5_1['x'],bed5_1['y']),bed5_1['z'],(xgrid,ygrid))

bed5_z0[~inmesh] = float('nan')
bed5_z0[inhole1] = float('nan')
bed5_z0[inhole2] = float('nan')
bed5_z1[~inmesh] = float('nan')
bed5_z1[inhole1] = float('nan')
bed5_z1[inhole2] = float('nan')

# Load bsmooth8
bed8_0 = elmerreadlib.pvtu_file(DIR_bed8+'mesh2d/terminusdriven0001.pvtu',['velocity','zs top'])
bed8_0 = bed8_0[bed8_0['zs top'] != 0]
bed8_1 = elmerreadlib.pvtu_file(DIR_bed8+'mesh2d/terminusdriven0367.pvtu',['velocity','zs top'])
bed8_1 = bed8_1[bed8_1['zs top'] != 0]

bed8_z0 = griddata((bed8_0['x'],bed8_0['y']),bed8_0['z'],(xgrid,ygrid))
bed8_z1 = griddata((bed8_1['x'],bed8_1['y']),bed8_1['z'],(xgrid,ygrid))

bed8_z0[~inmesh] = float('nan')
bed8_z0[inhole1] = float('nan')
bed8_z0[inhole2] = float('nan')
bed8_z1[~inmesh] = float('nan')
bed8_z1[inhole1] = float('nan')
bed8_z1[inhole2] = float('nan')


# Make plot

fig = plt.figure(figsize=(7.5,4))
grid = AxesGrid(fig,[0,0.05,0.95,0.92],  # similar to subplot(122)
                    nrows_ncols=(2, 4),
                    axes_pad=0.10,
                    label_mode="1",
                    share_all=True,
                    cbar_location="right",
                    cbar_mode="edge",
                    cbar_size="7%",
                    cbar_pad="2%",
                    )


xmin = np.min(bed3_0['x'])
xmax = np.max(bed3_0['x'])
ymin = np.min(bed3_0['y'])
ymax = np.max(bed3_0['y'])

cmaps = [plt.get_cmap("RdBu_r"),plt.get_cmap("RdBu")]

grid[0].imshow(zbed3,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',cmap=cmaps[0],clim=[-800,800])

grid[4].imshow(bed3_z1-bed3_z0,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',cmap=cmaps[1],clim=[-25,25])

grid[1].imshow(zbed4,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',cmap=cmaps[0],clim=[-800,800])

grid[5].imshow(bed4_z1-bed4_z0,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',cmap=cmaps[1],clim=[-25,25])

grid[2].imshow(zbed5,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',cmap=cmaps[0],clim=[-800,800])

grid[6].imshow(bed5_z1-bed5_z0,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',cmap=cmaps[1],clim=[-25,25])

im=grid[3].imshow(zbed8,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',cmap=cmaps[0],clim=[-800,800])
grid.cbar_axes[0].colorbar(im)

im=grid[7].imshow(bed8_z1-bed8_z0,extent=[x[0],x[-1],y[0],y[-1]],origin='lower',cmap=cmaps[1],clim=[-25,25])
grid.cbar_axes[1].colorbar(im)

grid.axes_llc.set_xticks([])
grid.axes_llc.set_yticks([])
grid.axes_llc.set_xlim([xmin,xmax])
grid.axes_llc.set_ylim([ymin,ymax])

cax= grid.cbar_axes[0]
cax.toggle_label(True)
cax.axis[cax.orientation].set_label('Bed elevation (m)')
cax= grid.cbar_axes[1]
cax.toggle_label(True)
cax.axis[cax.orientation].set_label(r'$z_o-z_f$ (m)')

plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/"+glacier+"_bsmooth_diff.pdf"),FORMAT='PDF',dpi=600)
plt.close()


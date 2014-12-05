#Plot bed map in polar coordinates

gmtset FONT_LABEL 9
gmtset FONT_ANNOT_PRIMARY 10
gmtset PS_MEDIA letter
gmtset MAP_ANNOT_ORTHO ver_text
gmtset MAP_FRAME_TYPE plain
gmtset COLOR_BACKGROUND 255
gmtset COLOR_FOREGROUND 255/0/0

#R=-R287950/314050/-2583000/-2559000
R=-R226000/319000/-2595000/-2515000
#R=-R270000/315000/-2590000/-2540000 
J=-Jx1:1000000
B=-B100000000m100000000m100000000

FILE=mesh_polar.ps

# Radar images
IMAGE=$HOME/Data/Mosaics/Greenland/mosaic100m.00-01geo.grd

# Outlines
COAST=$HOME/Data/Shape_files/greenland_coast_polar.txt
MESH=$HOME/Data/Meshes/3D/Helheim/planar/mesh.nodes
GLACIER=$HOME/Dropbox/Code/Solver_files/3D/Helheim/Inputs/mesh_extent.dat
FLOWLINE=$HOME/Data/Shape_files/Glaciers/Flowlines/Helheim/helheim_flowline.dat
HOLE1=$HOME/Dropbox/Code/Solver_files/3D/Helheim/Inputs/mesh_hole1.dat
HOLE2=$HOME/Dropbox/Code/Solver_files/3D/Helheim/Inputs/mesh_hole2.dat

makecpt -Cgray -T0/255/1 -Z > grayscale.cpt

awk '{print $3,$4}' $MESH > nodes.txt

psbasemap $J $R $B -P -K -V > $FILE
grdimage $IMAGE -Cgray $J $R $B -K -O -V >> $FILE
#grdimage bed1.grd -Ccolorscale.cpt $J $R $B -K -Q -O -V >> $FILE
#psxy $FLOWLINE $R $J -V -O -K >> $FILE
psxy $GLACIER $HOLE1 $HOLE2 -W1p,blue $R $J -L -O -K -V >> $FILE
psxy nodes.txt -Sc2p -Wblue -Gblue $R $J -O -K -V >> $FILE

#Scale bars
psxy $R $J -O -K -W1p <<END>> $FILE 
230000	-2588000 
230000	-2586000 
240000	-2586000
240000	-2588000
240000	-2586000
250000	-2586000 
250000	-2588000 
END

psxy $R $J -O -K -W1p -G255/255/255 <<END>> $FILE 
309200 -2562900
309200 -2561800
END

pstext $R $J -O <<END>> $FILE
230000 -2520000 12 1 1 TL a
229000  -2589000 11 0 0 TL 0
239000  -2589000 11 0 0 TL 10
249000  -2589000 11 0 0 TL 20 km
END

ps2raster -Tf mesh_polar.ps
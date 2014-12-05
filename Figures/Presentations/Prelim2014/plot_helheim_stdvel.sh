#Plot bed map in polar coordinates

gmtset FONT_LABEL 9
gmtset FONT_ANNOT_PRIMARY 10
gmtset PS_MEDIA letter
gmtset MAP_ANNOT_ORTHO ver_text
gmtset MAP_FRAME_TYPE plain
gmtset COLOR_BACKGROUND 255
gmtset COLOR_FOREGROUND 255/0/0

#R=-R287950/314050/-2583000/-2559000
R=-R283000/314000/-2586000/-2553000
#R=-R270000/315000/-2590000/-2540000 
J=-Jx1:400000
B=-B100000000m100000000m100000000

FILE=stdvel_polar.ps

# Radar images
IMAGE=$HOME/Data/Mosaics/Greenland/mosaic100m.00-01geo.grd

# Outlines
COAST=$HOME/Data/ShapeFiles/greenland_coast_polar.txt
GLACIER=$HOME/Models/Helheim/Meshes/3D/High_Normal/Inputs/mesh_extent.dat
FLOWLINE=$HOME/Data/ShapeFiles/Glaciers/Flowlines/Helheim/helheim_flowline.dat
HOLE1=$HOME/Models/Helheim/Meshes/3D/High_Normal/Inputs/mesh_hole1.dat
HOLE2=$HOME/Models/Helheim/Meshes/3D/High_Normal/Inputs/mesh_hole2.dat

xyz2grd velocities_all_std.dat -Gtaub1.grd $R -I100
#triangulate taub_1e10.dat -Gtaub1.grd $R -I100
grdmask $HOLE1 -Ghole1.grd $R -I100 -N1/1/NaN
grdmask $HOLE2 -Ghole2.grd $R -I100 -N1/1/NaN
grdmask $GLACIER -Gextent.grd $R -I100 -NNan/1/1
grdmath hole1.grd hole2.grd extent.grd MUL MUL = mask.grd
grdmath taub1.grd mask.grd MUL = taub2.grd
grdmath taub2.grd 1000 DIV = taub2.grd

makecpt -Cgray -T0/255/1 -Z > grayscale.cpt

psbasemap $J $R $B -P -K -V > $FILE
grdimage $IMAGE -Cgray $J $R $B -K -O -V >> $FILE
grdimage taub2.grd -Cvelocity_change.cpt $J $R $B -K -Q -O -V >> $FILE
#psxy $FLOWLINE $R $J -V -O -K >> $FILE
psxy $GLACIER $HOLE1 $HOLE2 -W1p $R $J -L -O -K -V >> $FILE

psxy $R $J -O -K -L -W0.5p -G255/255/255 <<END>> $FILE
313400 -2553500
313400 -2561000
302500 -2561000
302500 -2553500
END

psxy $R $J -O -K -L -W0.5p -G255/255/255 <<END>> $FILE
313400 -2566700
313400 -2563300
303700 -2563300
303700 -2566700
END

#Scale bars
psscale -D2.45i/2.95i/0.9i/0.15ih -Cvelocity_change.cpt -Qo -B0::/:: -O -K -V >> $FILE
psxy $R $J -O -K -W1p <<END>> $FILE 
305000	-2564500 
305000	-2564000 
310000	-2564000 
310000	-2564500 
END

psxy $R $J -O -K -W1p -G255/255/255 <<END>> $FILE 
308700 -2556000
308700 -2557500
END

psxy $R $J -O -K -W1p -G255/255/255 <<END>> $FILE 
311200 -2556000
311200 -2557500
END

psxy $R $J -O -K -W1p -G255/255/255 <<END>> $FILE 
306000 -2556000
306000 -2557500
END

pstext $R $J -O <<END>> $FILE
284100 -2555000 12 1 1 TL b
303000  -2554000 10 0 0 TL Velocity change
304500  -2565000 11 0 0 TL 0
309500  -2565000 11 0 0 TL 5 km
#308500  -2567200 11 0 0 TL km
308000  -2558000 11 0 0 TL 100
#312400  -2558000 11 0 0 TL 1000
306000  -2559200 11 0 0 TL m/yr
305200  -2558000 11 0 0 TL 10
#305073.776937 -2578289.20561 9 0 0 TL 0 km
#289749.968209 -2566725.38141 9 0 0 TL 20 km
END

ps2raster -Tf stdvel_polar.ps

rm taub2.grd taub1.grd mask.grd extent.grd hole1.grd hole2.grd grayscale.cpt gmt.history gmt.conf
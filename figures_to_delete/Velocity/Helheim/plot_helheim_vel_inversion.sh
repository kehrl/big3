#Plot a location map for Helheim Glacier
#Last edited: 06/05/2014, LMK, UW

gmtset FONT_LABEL 9
gmtset FONT_ANNOT_PRIMARY 9
gmtset PS_MEDIA letter
gmtset MAP_ANNOT_ORTHO ver_text
gmtset MAP_FRAME_TYPE plain
gmtset COLOR_BACKGROUND 127.5
gmtset COLOR_FOREGROUND 255/0/0


FILE=helheim_vel.ps

R=-R290000/315000/-2585000/-2567000
J=-Jx1:400000
B=-B10000000m10000000m10000000

#IMAGE=$HOME/Data/Mosaics/Helheim/mosaicHelheim.2013-062.148.31711_1-20mgeo.tif
IMAGE=$HOME/Data/Mosaics/Helheim/mosaicHelheim.2012-032.148.25699_1-20mgeo.tif
gdal_translate -b 1 -of GMT $IMAGE radar_image.grd
#xyz2grd radar_image.xyz -Gradar_image.grd $R -I20 

COAST=$HOME/Data/Shape_files/greenland_coast_polar.txt
GLACIER=$HOME/Data/Shape_files/Glaciers/3D/Helheim/glacier_extent.txt
FLOWLINE=$HOME/Data/Shape_files/Glaciers/Flowlines/Helheim/helheim_flowline.txt
HOLE1=$HOME/Data/Shape_files/Glaciers/3D/Helheim/glacier_hole1.txt
HOLE2=$HOME/Data/Shape_files/Glaciers/3D/Helheim/glacier_hole2.txt

#Velocity
DIR=$HOME/Data/Velocity/TSX/Helheim/Outputs/
#DATA=$DIR/vel_track-30542.txt #Dec. 21 2012
#DATA=$DIR/vel_track-31544.txt #Feb. 26, 2013
#DATA=$DIR/vel_track-32713.txt #May 14, 2013
DATA=$DIR/vel_track-17850.txt #Sept 2, 2010
VEL1=vel1.grd
VEL2=vel2.grd
VEL3=vel3.grd

#Get flowline coordinates
awk '{print $2,$3}' $FLOWLINE > temp2.txt

#Grid file, if needed
awk '{print $1,$2,$3}' $DATA > $DIR/temp.txt
xyz2grd $DIR/temp.txt -G$VEL1 -NNaN -I100 $R
rm $DIR/temp.txt
#grdmask $HOLE1 -Ghole1.grd $R -I100 -N1/1/NaN
#grdmask $HOLE2 -Ghole2.grd $R -I100 -N1/1/NaN
grdmask $GLACIER -Gextent.grd $R -I100 -NNan/1/1
grdmath extent.grd MUL = mask.grd
grdmath vel1.grd mask.grd MUL = $VEL2
grdmath $VEL2 1000 DIV = $VEL3
grd2xyz $VEL3 > velocities.txt


makecpt -Cgray -T0/255/1 -Z > grayscale.cpt
makecpt -Cseis -T0/11/0.1 > velocity.cpt

psbasemap $J $R $B -P -K -V > $FILE
grdimage radar_image.grd  -Cgrayscale.cpt $R $J -O -K -V >> $FILE
grdimage $VEL3 -Cvelocity.cpt $R $J -Q -O -K -V >> $FILE
psxy $GLACIER $HOLE1 $HOLE2 $R $J -L -O -K -V >> $FILE
psxy $FLOWLINE $R $J -W1p -K -O  -V >> $FILE


psxy $R $J -O -K -L -W0.5p -G255/255/255 <<END>> $FILE
290500 -2584500
296500 -2584500
296500 -2582000
290500 -2582000
END

psxy $R $J -O -K -L -W0.5p -G255/255/255 <<END>> $FILE
314500 -2567500
314500 -2571100
307000 -2571100
307000 -2567500
END

#Scale bars
psscale -D2.05i/1.66i/0.6i/0.1ih -Cvelocity.cpt -B0::/:: -O -K -V >> $FILE
psxy $R $J -O -K -W1p <<END>> $FILE 
291000	-2583000 
291000	-2582600 
296000	-2582600 
296000	-2583000 
END

psxy $R $J -O -K -Sc0.15c -G255/255/255 -W <<END>> $FILE 
309960  -2577501
305028  -2576891
300296  -2575325
297477  -2571340
294940  -2567064
END

pstext $R $J -O <<END>> $FILE
#290500 -2567800 10 1 1 TL a
290800  -2583200 9 0 0 TL 0
293400  -2583200 9 0 0 TL 5 km
309700 -2569500 9 0 0 TL 10 km/yr
307500 -2569500 9 0 0 TL 0
310500  -2577100 9 0 0 TL 0 km
291500  -2567464 9 0 0 TL 20 km
END

#ps2pdf helheim_vel.ps helheim_vel.pdf
#ps2raster -A -Tf -E300 helheim_map.ps
#rm *.grd *.cpt gmt* temp2.txt

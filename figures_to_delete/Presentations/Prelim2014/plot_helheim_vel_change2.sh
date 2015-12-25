#Plot a location map for Helheim Glacier
#Last edited: 06/05/2014, LMK, UW

gmtset FONT_LABEL 9
gmtset FONT_ANNOT_PRIMARY 9
gmtset PS_MEDIA letter
gmtset MAP_ANNOT_ORTHO ver_text
gmtset MAP_FRAME_TYPE plain
gmtset COLOR_BACKGROUND 127.5
gmtset COLOR_FOREGROUND 255/0/0


FILE=velocity_change.ps

R=-R288000/316000/-2583000/-2559000
#R=-R278000/315000/-2600000/-2550000
J=-Jx1:400000
B=-B10000000m10000000m10000000

#IMAGE=$HOME/Data/Mosaics/Helheim/mosaicHelheim.2013-062.148.31711_1-20mgeo.tif
IMAGE=$HOME/Data/Mosaics/Greenland/mosaic100m.05-06geo.grd

COAST=$HOME/Data/Shape_files/greenland_coast_polar.txt
GLACIER=$HOME/Dropbox/Code/Solver_files/3D/Helheim/Inputs/mesh_extent.dat
FLOWLINE=$HOME/Data/Shape_files/Glaciers/Flowlines/Helheim/helheim_flowline.dat
HOLE1=$HOME/Dropbox/Code/Solver_files/3D/Helheim/Inputs/mesh_hole1.dat
HOLE2=$HOME/Dropbox/Code/Solver_files/3D/Helheim/Inputs/mesh_hole2.dat

#Velocity

#2 km retreat during stability
#DIR=$HOME/Data/Velocity/TSX/Helheim/Outputs/
#DATA=$DIR/vel_track-31544.txt 
#DATA=$DIR/vel_track-32713.txt 
#RETREAT=$HOME/Data/Shape_files/Retreat_area/retreat_2013day51_to_2013day128.dat

#5 l, retreat during 2001-2006
DIR=$HOME/Data/Velocity/RADARSAT/Helheim/Outputs/
#DATA=$DIR/vel_winter00-01.txt 
DATA=$DIR/vel_winter05-06.txt 
RETREAT=$HOME/Data/Shape_files/Retreat_area/retreat_wint01_to_wint05.dat

#Grid file, if needed
awk '{print $1,$2,$3}' $DATA > $DIR/temp.txt
xyz2grd $DIR/temp.txt -Gvel1.grd -NNaN -I100 $R
rm $DIR/temp.txt


grdmask $HOLE1 -Ghole1.grd $R -I100 -N1/1/NaN
grdmask $HOLE2 -Ghole2.grd $R -I100 -N1/1/NaN
grdmask $GLACIER -Gextent.grd $R -I100 -NNan/1/1
grdmath hole2.grd hole1.grd extent.grd MUL MUL = mask.grd
grdmath vel1.grd mask.grd MUL = vel2.grd
grdmath vel2.grd 1000 DIV = vel3.grd


makecpt -Cgray -T0/255/1 -Z > grayscale.cpt
makecpt -Cdrywet -T0/10/0.1 -Z > velocity.cpt

psbasemap $J $R $B -P -K -V > $FILE
grdimage $IMAGE -Cgrayscale.cpt $R $J -O -K -V >> $FILE
grdimage vel3.grd -Cvelocity.cpt $R $J -Q -O -K -V >> $FILE
psxy $GLACIER $HOLE1 $HOLE2 $R $J -L -O -K -V >> $FILE
#psxy $RETREAT $R $J -W1p -L -Gp600/8 -K -O -V >> $FILE


#psxy $R $J -O -K -L -W0.5p -G255/255/255 <<END>> $FILE
#300500 -2574500
#306500 -2574500
#306500 -2572000
#300500 -2572000
#END

psxy $R $J -O -K -L -W0.5p -G255/255/255 <<END>> $FILE
315500 -2559500
315500 -2568700
305000 -2568700
305000 -2559500
END

#Scale bars
psscale -D2.2i/2.08i/0.9i/0.1ih -Cvelocity.cpt -B0::/:: -O -K -V >> $FILE
psxy $R $J -O -K -W1p <<END>> $FILE 
306800	-2567100 
306800	-2566600 
311800	-2566600 
311800	-2567100 
END

psxy $R $J -O -K -W1p -G255/255/255 <<END>> $FILE 
314900 -2563300
314900 -2561800
END

psxy $R $J -O -K -W1p -G255/255/255 <<END>> $FILE 
310350 -2563300
310350 -2561800
END

psxy $R $J -O -K -W1p -G255/255/255 <<END>> $FILE 
305700 -2563300
305700 -2561800
END

pstext $R $J -O <<END>> $FILE
#289000 -2560500 10 1 1 TL b
307400  -2560000 10 0 0 TL Velocity
306500  -2567300 10 0 0 TL 0
311500  -2567300 10 0 0 TL 5 km
305500  -2563500 10 0 0 TL 0
310000  -2563500 10 0 0 TL 5
313700  -2563500 10 0 0 TL 10
308500  -2564700 10 0 0 TL km/yr
#305500  -2563200 10 0 0 TL 0
END

#ps2pdf helheim_vel.ps helheim_vel.pdf
#ps2raster -A -Tf -E300 helheim_map.ps
rm vel1.grd gmt* temp2.txt

ps2raster -Tf velocity_change.ps
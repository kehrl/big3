#Plot a location map for Helheim Glacier
#Last edited: 06/05/2014, LMK, UW

gmtset FONT_LABEL 9
gmtset FONT_ANNOT_PRIMARY 10
gmtset PS_MEDIA letter
gmtset MAP_ANNOT_ORTHO ver_text
gmtset MAP_FRAME_TYPE plain
gmtset COLOR_BACKGROUND 127.5
gmtset COLOR_FOREGROUND 255/0/0

FILE=helheim_map.ps

R=-R290000/315000/-2585000/-2567000
J=-Jx1:400000
B=-B10000000m10000000m10000000

IMAGE=$HOME/Data/Mosaics/Helheim/mosaicHelheim.2012-032.148.25699_1-20mgeo.tif
gdal_translate -b 1 -of GMT $IMAGE radar_image.grd
#xyz2grd radar_image.xyz -Gradar_image.grd $R -I20 

COAST=$HOME/Data/Shape_files/greenland_coast_polar.txt
GLACIER=$HOME/Data/Shape_files/Glaciers/Helheim/glacier_extent.txt
FLOWLINE=$HOME/Data/Shape_files/Glaciers/Flowlines/helheim_flowline.txt
HOLE1=$HOME/Data/Shape_files/Glaciers/Helheim/glacier_hole1.txt
HOLE2=$HOME/Data/Shape_files/Glaciers/Helheim/glacier_hole2.txt
DIR=$HOME/Data/Velocity/TSX/Helheim/Outputs/
DATA=$DIR/vel_track-10335.txt


VEL=$DIR/vel_track-10335.grd

#Flowline coordinates
awk '{print $2,$3}' $FLOWLINE > flowline_coords.txt
#Grid file, if needed
#awk '{print $1,$2,$3}' $DATA > $DIR/temp.txt
#xyz2grd $DIR/temp.txt -G$VEL -N-1000 -I100 $R
#rm $DIR/temp.txt

makecpt -Cgray -T0/255/1 -Z > grayscale.cpt

psbasemap $J $R $B -P -K -V > $FILE
#psxy $COAST $R $J -O -K -V >> $FILE
#psxy $GLACIER $R $J -L -O -K -V >> $FILE
grdimage radar_image.grd  -Cgrayscale.cpt $R $J -O -K -V >> $FILE
psxy flowline_coords.txt $R $J -W1p -K -O  -V >> $FILE
#psscale -D5.5i/0.8i/2i/0.2ih -Ccolorscale.cpt -A -B5000::/:: -O -K -V >> $FILE

psxy $R $J -O -K -L -W0.5p -G255/255/255 <<END>> $FILE
290500 -2584500
296500 -2584500
296500 -2582000
290500 -2582000
END

#Scale bar
psxy $R $J -O -K -W1p <<END>> $FILE 
291000 -2583000
291000 -2582600 
296000 -2582600
296000 -2583000
END

pstext $R $J -O <<END>> $FILE
290800 -2583200 9 0 0 TL 0
293400 -2583200 9 0 0 TL 5 km
END

ps2pdf helheim_map.ps helheim_map.pdf
#ps2raster -A -Tf -E300 helheim_map.ps
#rm *.ps
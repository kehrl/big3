#Plot MEASURES velocity data in lat,long

gmtset FONT_LABEL 10
gmtset FONT_ANNOT_PRIMARY 10
gmtset PS_MEDIA letter
gmtset MAP_FRAME_TYPE plain
gmtset MAP_TICK_LENGTH_SECONDARY 0.0
gmtset COLOR_BACKGROUND 127.5
gmtset COLOR_FOREGROUND 255/0/0


COAST=$HOME/Data/Shape_files/Coastlines/greenland_coast_wgs84.txt
VELDATA=$HOME/Data/Velocity/Random/Greenland/track-07to10/vel_07to10_v.dat

#R=-R-76/-10/59.3/83.9
R=-R-60/58/15/80r
#JM2i
J=-JA-45/70/2i
B=-B0

FILE=velocity.ps

# Convert velocity tif file 
#gdalwarp -s_srs EPSG:3413 -t_srs EPSG:4326 $VELDATA vel.tif
#gdal_translate vel.tif -of XYZ vel.dat
#xyz2grd vel.dat -Gvel.grd -I0.1 $R -NNaN

#grdmask $COAST -Gmask.grd $R -I0.1 -NNan/1/1
#grdmath mask.grd vel.grd MUL = vel2.grd

#makecpt -Cdrywet -T0.1/6000/3 -Qo -Z > colorscale.cpt

psbasemap $J $R $B -L-62/70/-50/100k -P -K -V > $FILE
grdimage vel2.grd -Ccolorscale.cpt $J $R -Q -K -V -O >> $FILE

#psxy flowline.txt $R $J -O -V -K >> $FILE

psxy $COAST -W1p $R $J $B -O -K -V >> $FILE
psscale -D0.32i/0.43i/0.55i/0.1ih -Ccolorscale.cpt -B0 -Q -O -K -V >> $FILE


psxy $R $J -Sc0.2 -W2p -Gred -V -O -K -P <<END>> $FILE
#-33.0000 68.63333
-38.2000 66.3500
#-49.8333 69.1667
END

psxy $R $J -A -W1p -V -O -K -P <<END>> $FILE
-57.1 62
-56.5 60.7
END

psxy $R $J -A -W1p -V -O -K -P <<END>> $FILE
-52.9 62.5
-52.4 61
END

#psxy $R $J -W1p -V -O -K -P <<END>> $FILE
#-72.8 65
#-72.8 63.5
#END


pstext $R $J -O <<END>> $FILE
-62 63 11 0 0 TL Velocity 
-58 59.5 10 0 0 TL m/yr
-58.2 60.4 10 0 0 TL 10
-54 60.8 10 0 0 TL 1000
-90.5 80.5 14 1 1 TL a
#-36 71 10 0 0 TL KL
#-48 67 10 0 0 TL HH
-37 65.5 12 0 0 TL Helheim
#-58 70 10 0 0 TL JK
END

rm vel.tif
ps2raster -Tf velocity.ps
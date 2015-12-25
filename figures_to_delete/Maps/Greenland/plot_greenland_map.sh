#Plot Greenland in lat long
gmtset FONT_LABEL 11
gmtset FONT_ANNOT_PRIMARY 14
gmtset MAP_FRAME_TYPE plain

COAST=$HOME/Data/ShapeFiles/Coastlines/greenland_coast_wgs84.txt

R=-R-76.5/-9.5/58.5/84.5
J=-JM2i
B=-B10m10m10

FILE=greenland_map.ps

#pscoast $J $R $B -W1 -P -V > $FILE
psbasemap $R $J -L-64.5/65.3/-50/500k -K -V > $FILE
psxy $COAST -W0.5p $R $J -O -K -V >> $FILE
#grdtrack flowline.txt -Gvel.grd -S > flowline_vel.txt

psxy $R $J -Sc0.4 -W1p -Gred -V -O -P -K <<END>> $FILE
#-33.0000 68.63333
-38.2000 66.3500
#-49.8333 69.1667
END


pstext $R $J -O <<END>> $FILE
#-34 66.5 10 0 0 TL Velocity (m/yr)
#-36 71 10 0 0 TL KL
-37 65.5 16 0 0 TL Helheim
#-58 70 10 0 0 TL JK
END

ps2raster -A -Tf greenland_map.ps
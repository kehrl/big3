#Plot MEASURES velocity data in lat,long

gmtset LABEL_FONT_SIZE 9
gmtset ANNOT_FONT_SIZE_PRIMARY 10
gmtset PAPER_MEDIA letter+
gmtset Y_AXIS_TYPE ver_text
gmtset BASEMAP_TYPE plain
gmtset COLOR_FOREGROUND 255/0/0
gmtset COLOR_BACKGROUND 127.5

data=kanger_bed.txt

#R=-R-2589750/-2510250/250250/329750
#J=-Jx1:500000
#B=-B1000000m1000000m1000000
R=-R-33.77/-32.5/68.42/68.95
#R=-R-76/-10/59/84
J=-JM2i
B=-B10m10m10

FILE=bed.ps

#Grid bed data
xyz2grd $data -Gbed.grd -I600e $R


makecpt -Cpolar -T-1200/1200/50 > colorscale.cpt

psbasemap $J $R $B -P -K -V > $FILE
grdimage bed.grd -Ccolorscale.cpt $R $J -O -K >> $FILE
#psxy greenland_coast.txt $R -W1.5p $J -O -V -K >> $FILE

psxy $R $J $J -G255/255/255 -W1p0/0/0 -O -K <<END>> $FILE
-33.19 68.79
-32.53 68.79
-32.53 68.94
-33.19 68.94
END

psscale -D1.47i/1.92i/0.6i/0.1ih -Ccolorscale.cpt -A -B1200::/:: -O -K -V >> $FILE

#grdtrack flowline.txt -Gvel.grd -S > flowline_vel.txt

pstext $R $J -O <<END>> $FILE
-33.05 68.83 9 0 0 TL Elevation (m)
END

ps2raster -A -Tf bed.ps
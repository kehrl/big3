#Plot MEASURES velocity data in lat,long

gmtset LABEL_FONT_SIZE 9
gmtset ANNOT_FONT_SIZE_PRIMARY 10
gmtset PAPER_MEDIA letter+
gmtset Y_AXIS_TYPE ver_text
gmtset BASEMAP_TYPE plain
gmtset COLOR_BACKGROUND 127.5

data=jakobshavn_bed.txt

#R=-R-2589750/-2510250/250250/329750
#J=-Jx1:500000
#B=-B1000000m1000000m1000000
R=-R-51.1/-45/68.6/70
#R=-R-76/-10/59/84
J=-JM2i
B=-B400m400m400

FILE=bed.ps

#Grid bed data
xyz2grd $data -Gbed.grd -I600e $R
grdmask greenland_coast.txt -Gmask.grd -I600e $R -N-10000/1/1
grdmath bed.grd mask.grd MUL = bed_masked.grd


makecpt -Cpolar -T-1200/1200/50 > colorscale.cpt

psbasemap $J $R $B -P -K -V > $FILE
grdimage bed.grd -Ccolorscale.cpt $R $J -O -K >> $FILE
#psxy greenland_coast.txt $R -W1.5p $J -O -V -K >> $FILE
#psbasemap -J -R -L-45.75/69.45/-45.75/50k -P -K -O -V >> $FILE


#psxy $R $J $J -G255/255/255 -W1p0/0/0 -O -K <<END>> $FILE
#-47.7 69.33
#-45.1 69.33
#-45.1 69.95
#-47.7 69.95
#END

#psscale -D1.55i/0.95i/0.5i/0.1ih -Ccolorscale.cpt -A -B1500::/:: -O -K -V >> $FILE

#grdtrack flowline.txt -Gvel.grd -S > flowline_vel.txt

pstext $R $J -O <<END>> $FILE
#-47.5 69.49 9 0 0 TL Elevation (m)
END

ps2raster -A -Tf bed.ps
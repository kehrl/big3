#Plot bed map in lat,long

gmtset FONT_LABEL 9
gmtset FONT_ANNOT_PRIMARY 10
gmtset PS_MEDIA letter
gmtset MAP_ANNOT_ORTHO ver_text
gmtset MAP_FRAME_TYPE plain
gmtset COLOR_BACKGROUND 127.5
gmtset COLOR_FOREGROUND 255/0/0

data=/Users/laurakehrl/Data/Bed/Cresis/Helheim_2008_2012_Composite/grids/helheim_bed_epsg4326.txt

#R=-R-2589750/-2510250/250250/329750
#J=-Jx1:500000
#B=-B1000000m1000000m1000000
R=-R-39.5/-37.5/66.2/67.02
#R=-R-76/-10/59/84
J=-JM2i
B=-B10m10m10

FILE=bed_wgs84.ps

#Grid bed data
xyz2grd $data -Gbed.grd -I600e $R


makecpt -Cpolar -T-1200/1200/50 -V -Z > colorscale.cpt

psbasemap $J $R $B -P -K -V > $FILE
grdimage bed.grd -Ccolorscale.cpt $R $J -O -V -K >> $FILE
#psxy greenland_coast.txt $R -W1.5p $J -O -V -K >> $FILE

#psxy $R $J $J -G255/255/255 -W1p0/0/0 -O -K <<END>> $FILE
#-38.6 66.75
#-37.55 66.75
#-37.55 67
#-38.6 67
#END

#psscale -D1.45i/1.72i/0.6i/0.1ih -Ccolorscale.cpt -A -B1500::/:: -O -K -V >> $FILE

#grdtrack flowline.txt -Gvel.grd -S > flowline_vel.txt

pstext $R $J -O <<END>> $FILE
#-38.4 66.82 9 0 0 TL Elevation (m)
-35 70.5 10 0 0 TL KL
-46.5 67 10 0 0 TL HH
-57 70 10 0 0 TL JK

END

ps2raster -A -Tf bed_wgs84.ps

#Plot MEASURES velocity data in lat,long

gmtset LABEL_FONT_SIZE 10
gmtset ANNOT_FONT_SIZE_PRIMARY 10
gmtset PAPER_MEDIA letter+
gmtset Y_AXIS_TYPE ver_text
gmtset BASEMAP_TYPE plain

#DIRECTORY=/Users/kehrl/data/velocity/measures/regions/Wcoast-69.95N/TSX_Mar-27-2011_Apr-07-2011_09-56-36
#DIRECTORY=/Users/kehrl/data/velocity/measures/2008
DIRECTORY=/Users/kehrl/data/velocity/twila_composite
data=$DIRECTORY/xyv.txt

#R=-R-2589750/-2510250/250250/329750
#J=-Jx1:500000
B=-B1000000m1000000m1000000
#R=-R-51.2/-41.1/68.3/70.4
R=-R-51.1/-45/68.6/70
J=-JM2i
#B=-B10m10m10

FILE=vel.ps

#Grid bed data
#xyz2grd $data -Gvel.grd -I1000e $R


makecpt -Crainbow -I -T0/9000/50 -V -Z > colorscale.cpt

psbasemap $J $R $B -P -K -V > $FILE
#grdimage vel.grd -Ccolorscale.cpt $R $J -O -V -K >> $FILE
grdimage $DIRECTORY/mosaic07to10_red.grd $DIRECTORY/mosaic07to10_green.grd $DIRECTORY/mosaic07to10_blue.grd $R $J -O -V -K >> $FILE

#psxy greenland_coast.txt $R -W1.5p $J -O -V -K >> $FILE

psxy $R $J $J -G255/255/255 -W1p0/0/0 -O -K <<END>> $FILE
-46.45 69.57
-45.1 69.57
-45.1 69.95
-46.45 69.95
END
psbasemap -J -R -L-45.75/69.9/-45.75/25k -P -K -O -V >> $FILE

#psscale -D2.4i/1.55i/0.8i/0.1ih -Ccolorscale.cpt -A -B9000::/:: -O -K -V >> $FILE

#grdtrack flowline.txt -Gvel.grd -S > flowline_vel.txt

pstext $R $J -O <<END>> $FILE
#-47.1 69.6 10 0 0 TL Velocity (m/yr)
END

ps2raster -A -Tf vel.ps
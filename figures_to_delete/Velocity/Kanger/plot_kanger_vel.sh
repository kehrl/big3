#Plot MEASURES velocity data in lat,long

gmtset LABEL_FONT_SIZE 10
gmtset ANNOT_FONT_SIZE_PRIMARY 10
gmtset PAPER_MEDIA letter+
gmtset Y_AXIS_TYPE ver_text
gmtset BASEMAP_TYPE plain

#DIRECTORY=/Users/kehrl/data/velocity/measures/regions/Ecoast-68.80N/TSX_Mar-09-2011_Mar-20-2011_19-40-41
DIRECTORY=/Users/kehrl/data/velocity/measures/2008
data=$DIRECTORY/xyv.txt
DIRECTORY=/Users/kehrl/data/velocity/twila_composite


#R=-R-2589750/-2510250/250250/329750
#J=-Jx1:500000
#B=-B1000000m1000000m1000000
R=-R-33.77/-32.5/68.42/68.95
#R=-R-40/-36/66/67
J=-JM2i
B=-B10m10m10

FILE=vel.ps

#Grid bed data
xyz2grd $data -Gvel.grd -I1000e $R


makecpt -Crainbow -I -T0/10000/50 -V -Z > colorscale.cpt

psbasemap $J $R $B -P -K -V > $FILE
#grdimage vel.grd -Ccolorscale.cpt $R $J -O -V -K >> $FILE
grdimage $DIRECTORY/mosaic07to10_red.grd $DIRECTORY/mosaic07to10_green.grd $DIRECTORY/mosaic07to10_blue.grd $R $J -O -V -K >> $FILE
#psxy greenland_coast.txt $R -W1.5p $J -O -V -K >> $FILE

psxy $R $J $J -G255/255/255 -W1p0/0/0 -O -K <<END>> $FILE
-33.17 68.79
-32.53 68.79
-32.53 68.94
-33.17 68.94
END

psxy $R $J $J -G255/255/255 -W1p0/0/0 -O -K <<END>> $FILE
-32.9 68.78
-32.53 68.78
-32.53 68.69
-32.9 68.69
END
psbasemap -J -R -L-32.72/68.76/-33.52/25k -P -K -O -V >> $FILE

#psscale -D0.65i/3.1i/0.8i/0.1ih -Ccolorscale.cpt -A -B10000::/:: -O -K -V >> $FILE

#grdtrack flowline.txt -Gvel.grd -S > flowline_vel.txt

pstext $R $J -O <<END>> $FILE
-33.1 68.83 9 0 0 TL Velocity (m/yr)
END

ps2raster -A -Tf vel.ps
#Plot bed map in polar coordinates

gmtset FONT_LABEL 9
gmtset FONT_ANNOT_PRIMARY 10
gmtset PS_MEDIA letter
gmtset MAP_ANNOT_ORTHO ver_text
gmtset MAP_FRAME_TYPE plain
gmtset COLOR_BACKGROUND 255
gmtset COLOR_FOREGROUND 255/0/0

#R=-R287950/314050/-2583000/-2559000
R=-R226000/319000/-2595000/-2515000
#R=-R270000/315000/-2590000/-2540000 
J=-Jx1:1000000
B=-B100000000m100000000m100000000

FILE=bed_polar.ps

# Beds
BED=$HOME/Data/Bed/Morlighem_2014/morlighem_helheim_bed.txt
FLOWLINE_BED=$HOME/Data/Bed/Morlighem_2014/morlighem_helheim_bed_flowline.txt
#BED=$HOME/Data/Bed/Ben_smith_feb2014/smith_helheim_bed.grd
#BED=$HOME/Data/Bed/Bamber_2001/bamber_5km_bed.grd

xyz2grd $BED -Gbed1.grd $R -I300 -NNaN

# Radar images
IMAGE=$HOME/Data/Mosaics/Greenland/mosaic100m.00-01geo.grd
#gdal_translate -b 1 -of GMT $IMAGE radar_image.grd
#xyz2grd radar_image.xyz -Gradar_image.grd $R -I20 

# Outlines
COAST=$HOME/Data/Shape_files/greenland_coast_polar.txt
GLACIER=$HOME/Dropbox/Code/Solver_files/3D/Helheim/Inputs/mesh_extent.dat
FLOWLINE=$HOME/Data/Shape_files/Glaciers/Flowlines/Helheim/helheim_flowline.dat
HOLE1=$HOME/Dropbox/Code/Solver_files/3D/Helheim/Inputs/mesh_hole1.dat
HOLE2=$HOME/Dropbox/Code/Solver_files/3D/Helheim/Inputs/mesh_hole2.dat

#Grid bed data
grdmask $HOLE1 -Ghole1.grd $R -I300 -N1/1/NaN
grdmask $HOLE2 -Ghole2.grd $R -I300 -N1/1/NaN
grdmask $GLACIER -Gextent.grd $R -I300 -NNan/1/1
grdmath hole1.grd hole2.grd extent.grd MUL MUL = mask.grd
grdmath bed1.grd mask.grd MUL = bed2.grd


makecpt -Cgray -T0/255/1 -Z > grayscale.cpt
makecpt -Cpolar -T-1000/1000/50 -V -Z > colorscale.cpt

psbasemap $J $R $B -P -K -V > $FILE
grdimage $IMAGE -Cgray $J $R $B -K -O -V >> $FILE
grdimage bed2.grd -Ccolorscale.cpt $J $R $B -K -Q -O -V >> $FILE
#psxy $FLOWLINE $R $J -V -O -K >> $FILE
psxy $GLACIER $HOLE1 $HOLE2 -W1p $R $J -L -O -K -V >> $FILE

#Calculate bed along flowline
#awk '{print $2,$3}' $FLOWLINE > junk.txt
#grdtrack junk.txt -G$BED > junk2.txt
#paste $FLOWLINE junk2.txt | awk '{print $1,$2,$3,$8}' > $FLOWLINE_BED
#rm junk2.txt junk.txt

psxy $R $J -O -K -L -W0.5p -G255/255/255 <<END>> $FILE
288800 -2516500
288800 -2536000
317500 -2536000
317500 -2516500
END

#Scale bars
psscale -D3.04i/2.85i/1i/0.15ih -Ccolorscale.cpt -B0::/:: -O -K -V >> $FILE
psxy $R $J -O -K -W1p <<END>> $FILE 
230000	-2588000 
230000	-2586000 
240000	-2586000
240000	-2588000
240000	-2586000
250000	-2586000 
250000	-2588000 
END

#psxy $R $J -O -K -Sc0.15c -G255/255/255 -W <<END>> $FILE 
#308068.140459 -2577500.18262
#303071.264922 -2577411.6871
#298820.865145 -2575260.75625
#296554.766007 -2570828.61875
#293749.968209 -2566725.38141
#290206.375958 -2563197.34722
#END

psxy $R $J -O -K -W1p -G255/255/255 <<END>> $FILE 
303800 -2523000
303800 -2526000
END

pstext $R $J -O <<END>> $FILE
230000 -2520000 12 1 1 TL c
290000  -2518000 12 0 0 TL Bed Elevation
229000  -2589000 11 0 0 TL 0
239000  -2589000 11 0 0 TL 10
249000  -2589000 11 0 0 TL 20 km
302500  -2528000 12 0 0 TL 0
314300  -2528000 12 0 0 TL 1
301800  -2532100 12 0 0 TL km
290500  -2528000 12 0 0 TL -1
#305073.776937 -2578289.20561 9 0 0 TL 0 km
#289749.968209 -2566725.38141 9 0 0 TL 20 km
END

ps2raster -Tf bed_polar.ps
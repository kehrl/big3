#Plot bed map in polar coordinates

gmtset FONT_LABEL 9
gmtset FONT_ANNOT_PRIMARY 10
gmtset PS_MEDIA letter
gmtset MAP_ANNOT_ORTHO ver_text
gmtset MAP_FRAME_TYPE plain
#gmtset COLOR_BACKGROUND 255
#gmtset COLOR_FOREGROUND 255/0/0

#R=-R287950/314050/-2583000/-2559000
R=-R286000/317000/-2583000/-2555000
#R=-R270000/315000/-2590000/-2540000 
J=-Jx1:400000
B=-B100000000m100000000m100000000

FILE=glacier_polar.ps

# Radar images
IMAGE1=$HOME/Data/Mosaics/Helheim/mosaicHelheim.2011-067.148.3959_1-20mgeo.tif
gdal_translate -b 1 -of GMT $IMAGE1 radar_image.grd
grdclip radar_image.grd -Gradar_image2.grd $R -Sb1/NaN
IMAGE2=$HOME/Data/Mosaics/Greenland/mosaic100m.00-01geo.grd


# Outlines
COAST=$HOME/Data/Shape_files/greenland_coast_polar.txt
GLACIER=$HOME/Dropbox/Code/Solver_files/3D/Helheim/Inputs/mesh_extent.dat
FLOWLINE=$HOME/Data/Shape_files/Glaciers/Flowlines/Helheim/helheim_flowline.dat
HOLE1=$HOME/Dropbox/Code/Solver_files/3D/Helheim/Inputs/mesh_hole1.dat
HOLE2=$HOME/Dropbox/Code/Solver_files/3D/Helheim/Inputs/mesh_hole2.dat

#cp $FLOWLINE flowline.dat

makecpt -Cgray -T0/255/1 -Z > grayscale.cpt


psbasemap $J $R $B -P -K -V > $FILE
grdimage $IMAGE2  -Cgrayscale.cpt $R $J -O -K -V >> $FILE
grdimage radar_image2.grd  -Cgrayscale.cpt $R $J -O -K -V >> $FILE
psxy flowline.dat $R $J -W2p -V -O -K >> $FILE
#psxy $GLACIER $HOLE1 $HOLE2 -W1p $R $J -L -O -K -V >> $FILE

#psxy $R $J -O -K -L -W0.5p -G255/255/255 <<END>> $FILE
#313400 -2553500
#313400 -2561000
#302500 -2561000
#302500 -2553500
#END

psxy $R $J -O -K -L -W0.5p -G255/255/255 <<END>> $FILE
316400 -2559200
316400 -2555800
306400 -2555800
306400 -2559200
END

#Scale bars
#psscale -D2.45i/2.95i/0.9i/0.15ih -Ccolorscale.cpt -B0::/:: -Qo -O -K -V >> $FILE
psxy $R $J -O -K -W1p <<END>> $FILE 
308000	-2557000 
308000	-2556500 
313000	-2556500 
313000	-2557000 
END

#psxy $R $J -O -K -W1p -G255/255/255 <<END>> $FILE 
#309500 -2556000
#309500 -2557500
#END

#psxy $R $J -O -K -W1p -G255/255/255 <<END>> $FILE 
#306400 -2556000
#306400 -2557500
#END

pstext $R $J -O <<END>> $FILE
287100 -2556500 12 1 1 TL b
#302800 -2554200 9 0 0 TL Basal shear stress
307500  -2557500 11 0 0 TL 0
312500  -2557500 11 0 0 TL 5 km
#308500  -2567200 11 0 0 TL km
#306000  -2558000 11 0 0 TL 10
#308800  -2558000 11 0 0 TL 100
#306000  -2559400 11 0 0 TL kPa
#303200  -2558000 11 0 0 TL 1
#305073.776937 -2578289.20561 9 0 0 TL 0 km
#289749.968209 -2566725.38141 9 0 0 TL 20 km
END

ps2raster -Tf glacier_polar.ps
#Plot bed map in polar coordinates

gmtset FONT_LABEL 9
gmtset FONT_ANNOT_PRIMARY 10
gmtset PS_MEDIA letter
gmtset MAP_ANNOT_ORTHO ver_text
gmtset MAP_FRAME_TYPE plain
gmtset COLOR_BACKGROUND 255
gmtset COLOR_FOREGROUND 255/0/0

#R=-R287950/314050/-2583000/-2559000
R=-R214000/319000/-2600000/-2510000
#R=-R270000/315000/-2590000/-2540000 
J=-Jx1:700000
B=-B100000000m100000000m100000000

FILE=mesh_polar.ps

# Radar images
IMAGE=$HOME/Data/Mosaics/Greenland/mosaic100m.05-06geo.grd

# Outlines
VEL=VMAG.xy
COAST=$HOME/Data/ShapeFiles/greenland_coast_polar.txt
FLOWLINE=flowline.dat
GLACIER=$HOME/Models/Helheim/Meshes/3D/InverseClass/Inputs/mesh_extent.dat
HOLE1=$HOME/Models/Helheim/Meshes/3D/InverseClass/Inputs/mesh_hole1.dat
HOLE2=$HOME/Models/Helheim/Meshes/3D/InverseClass/Inputs/mesh_hole2.dat

xyz2grd $VEL -Gvel1.grd $R -I300 -NNaN

grdmask $HOLE1 -Ghole1.grd $R -I300 -N1/1/NaN
grdmask $HOLE2 -Ghole2.grd $R -I300 -N1/1/NaN
grdmask velocity_extent.dat -Gextent.grd $R -I300 -NNan/1/1
grdmath hole1.grd hole2.grd extent.grd MUL MUL = mask.grd
grdmath vel1.grd mask.grd MUL = vel2.grd

#makecpt -Cdrywet -T1/100000/3 -Z -Qo > velscale.cpt
makecpt -Cgray -T0/255/1 -Z > grayscale.cpt

psbasemap $J $R $B -P -K -V > $FILE
grdimage $IMAGE -Cgray $J $R $B -K -O -V >> $FILE
grdimage vel2.grd -Cvelscale.cpt $J $R $B -K -O -V -Q >> $FILE
#grdimage bed1.grd -Ccolorscale.cpt $J $R $B -Q -O -V -K >> $FILE
#psxy $FLOWLINE $R $J -V -O -K >> $FILE
psxy $GLACIER $HOLE1 $HOLE2 -W2p,black $R $J -L -O -K -V >> $FILE
psxy $FLOWLINE -W2p,black,- $R $J -O -K -V >> $FILE
#psxy nodes.txt -Sc2p -Wblue -Gblue $R $J -O -K -V >> $FILE

#Scale bars
psscale -D0.8i/4.5i/1.2i/0.2ih -Cvelscale.cpt -Qo -B0 -O -K -V >> $FILE
psxy $R $J -O -K -W1p -G255/255/255 <<END>> $FILE 
223000 -2520000
223000 -2525000
END

psxy $R $J -O -K -W1p -G255/255/255 <<END>> $FILE 
233500 -2520000
233500 -2525000
END

psxy $R $J -O -K -W1p <<END>> $FILE 
218000	-2533000 
218000	-2531000 
228000	-2531000 
228000	-2533000 
END

psxy $R $J -O -K -W1p -G255/255/255 <<END>> $FILE 
309200 -2562900
309200 -2561800
END

#psxy $R $J -O -K -W1p,- -L <<END>> $FILE 
#283000 -2586000
#314000 -2586000
#314000 -2553000
#283000 -2553000
#END

# Plot colored points
#psxy $R $J -O -K -Sc0.4 -W1p,black -G0/0/255 <<END>> $FILE
#307603.127391 -2576603.36739
#END
#psxy $R $J -O -K -Sc0.4 -W1p,black -G0/255/255 <<END>> $FILE
#302606.528319 -2576451.94854
#END
#psxy $R $J -O -K -Sc0.4 -W1p,black -G0/255/0 <<END>> $FILE
#298279.549287 -2574311.03185
#END
#psxy $R $J -O -K -Sc0.4 -W1p,black -G255/255/0 <<END>> $FILE
#293715.243861 -2565508.40342
#END
#psxy $R $J -O -K -Sc0.4 -W1p,black -G255/0/0 <<END>> $FILE
#286903.450133 -2559818.03087
#END

psxy $R $J -O -K -Sc0.2 -W1p,black -Gwhite <<END>> $FILE
305635.040792 -2576526.65039
END

psxy $R $J -O -K -Sc0.2 -W1p,black -Gwhite <<END>> $FILE
#305635.040792 -2576526.65039
297418.041796 -2572551.87064
292253.770448 -2564250.8552
284595.101792 -2557834.64926
277156.208252 -2551158.76415
269579.909803 -2544632.88213
262387.862518 -2537683.87666
255579.00153  -2530361.96072
249726.980328 -2522292.22662
END

pstext $R $J -O -K <<END>> $FILE
306635.040792 -2576526.65039 14 0 0 TR 80 km
297918.041796 -2571851.87064 14 0 0 BL 70
292253.770448 -2563250.8552  14 0 0 BL 60
285595.101792 -2557834.64926 14 0 0 BL 50
278156.208252 -2551158.76415 14 0 0 BL 40
270579.909803 -2544632.88213 14 0 0 BL 30
263387.862518 -2537683.87666 14 0 0 BL 20
256579.00153  -2530361.96072 14 0 0 BL 10
251756.980328 -2522292.22662 14 0 0 TL 0 km
END

pstext $R $J -O <<END>> $FILE
217000 -2513000 24 0 0 TL Velocity (m/yr)
221000 -2525000 16 0 0 TL 10
231000 -2525000 16 0 0 TL 1000
217000  -2534000 16 0 0 TL 0
227000  -2534000 16 0 0 TL 10 
220500  -2537000 16 0 0 TL km
END

ps2raster -Tf mesh_polar.ps
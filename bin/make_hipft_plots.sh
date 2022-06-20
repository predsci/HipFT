#!/bin/bash

# This script can be used to plot the output of HipFT.
# It requires the PSI tool "plot2d", and ffmpeg with the H.264 codec.
# It assumes a linux machine.

DATADIR=$1
OUTFILE=$2

mkdir $DATADIR/plots
cd $DATADIR/plots

for file in $(ls ../*.h5)
do
  plot2d -cmin -25 -cmax 25 -dpi 92 -tp -ll -finegrid -unit_label Gauss $file -o "$(basename $file .h5).png"
done

mkdir tmp
cd tmp

i=1

for file in $(ls ../*.png)
do

  ln -s $file  movie$i.png

  i=$((i+1))

done

ffmpeg -framerate 15 -i "movie%d.png" -pix_fmt yuv420p -c:a copy -crf 20 -r 15 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -codec:v libx264 "movie.mov"

cp movie.mov ../${OUTFILE}.mov

cd ..
rm -fr tmp







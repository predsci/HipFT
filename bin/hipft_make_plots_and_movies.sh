#!/bin/bash

# This script can be used to plot the output of HipFT.
# It requires the PSI tool "plot2d", and ffmpeg with the H.264 codec.
# It assumes a linux machine.

# Need to add UT times!!!
# To do that, check for a ut version fo the hipft dump text file.
# If not there, than do not use UT, otherwise use that list to plot files,
# Adding the UT as a title on all plots.

#Directory of output files
DATADIR=$1
#Base filename for output mov movies.
OUTFILE=$2
#Slice index (if not specified all slices wil be plotted)
Slice=$3

DPI=128
LABEL="Gauss"
CMIN="-25"
CMAX="25"

CURRDIR=$PWD

mkdir $DATADIR/plots
cd $DATADIR/plots

idx=0
for file in ../*.h5 
do

  idx=$(($idx+1))
  
#   ADD UT TIME STUFF HERE.
#   If no UT, use hours.

# If the UT time dump file exists, use index to extract UT time for title.
# If not, use original dump file to add  hours to title.

  TITLE=" "

  TWOD=$(h5dump -d dim3 -H $file 2>/dev/null |tr '"' '\n' | grep DATASPACE\ \ SIMPLE | head -1 | tr -dc '^[0-9]/' | tr '/' '\n'| head -1)

  if [[ $TWOD -eq "" ]]; then
    plot2d -title ${TITLE} -cmin ${CMIN} -cmax ${CMAX} -dpi $DPI -tp -ll -finegrid -unit_label $LABEL $file -o "$(basename $file .h5).png"
  elif [[ "$Slice" = "all" ]]; then
    for i in $(seq 1 $TWOD)
    do
      hipft_extract_realization.py $file -r $i -o tmp_file.h5
      fslice=$(printf %06d\\n $i)
      plot2d -title ${TITLE} -cmin ${CMIN} -cmax ${CMAX} -dpi $DPI -tp -ll -finegrid -unit_label $LABEL tmp_file.h5 -o "$(basename $file .h5)_r$fslice.png"
    done
  elif [[ $Slice -ne "" ]]; then
    if [[ $Slice -gt $TWOD  ]]; then
      echo "Slice requested is outside range defaulting to last slice : $TWOD"
      hipft_extract_realization.py $file -r $TWOD -o tmp_file.h5
      fslice=$(printf %06d\\n $TWOD)
      plot2d -title ${TITLE} -cmin ${CMIN} -cmax ${CMAX} -dpi $DPI -tp -ll -finegrid -unit_label $LABEL tmp_file.h5 -o "$(basename $file .h5)_r$fslice.png"
    else
      hipft_extract_realization.py $file -r $Slice -o tmp_file.h5
      fslice=$(printf %06d\\n $Slice)
      plot2d -title ${TITLE} -cmin ${CMIN} -cmax ${CMAX} -dpi $DPI -tp -ll -finegrid -unit_label $LABEL tmp_file.h5 -o "$(basename $file .h5)_r$fslice.png"
    fi
  else
    hipft_extract_realization.py $file -r 1 -o tmp_file.h5
    plot2d -title ${TITLE} -cmin ${CMIN} -cmax ${CMAX} -dpi $DPI -tp -ll -finegrid -unit_label $LABEL tmp_file.h5 -o "$(basename $file .h5)_r000001.png"
  fi
done

rm tmp_file.h5 2>/dev/null

mkdir tmp
cd tmp

i=1

if [[ $TWOD -eq "" ]]; then
  for file in ../*.png
  do
    ln -s $file  movie$i.png
    i=$((i+1))
  done
  ffmpeg -framerate 15 -i "movie%d.png" -pix_fmt yuv420p -c:a copy -crf 20 -r 15 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -codec:v libx264 "movie.mov"
  cp movie.mov ${CURRDIR}/${OUTFILE}.mov
elif [[ "$Slice" = "all" ]]; then
  for j in $(seq 1 $TWOD)
  do
    i=1
    fslice=$(printf %06d\\n ${j})
    for file in ../*_r$fslice.png
    do
      ln -s $file  movie$i.png
      i=$((i+1))
    done
    ffmpeg -framerate 15 -i "movie%d.png" -pix_fmt yuv420p -c:a copy -crf 20 -r 15 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -codec:v libx264 "movie.mov"
    cp movie.mov ${CURRDIR}/${OUTFILE}_r$fslice.mov
    rm movie.mov
  done
elif [[ $Slice -ne "" ]]; then
  if [[ $Slice -gt $TWOD  ]]; then
    fslice=$(printf %06d\\n $TWOD)
    for file in ../*_r$fslice.png
    do
      ln -s $file  movie$i.png
      i=$((i+1))
    done
    ffmpeg -framerate 15 -i "movie%d.png" -pix_fmt yuv420p -c:a copy -crf 20 -r 15 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -codec:v libx264 "movie.mov"
    cp movie.mov ${CURRDIR}/${OUTFILE}_r$fslice.mov
  else
    fslice=$(printf %06d\\n $Slice)
    for file in ../*_r$fslice.png
    do
      ln -s $file  movie$i.png
      i=$((i+1))
    done
    ffmpeg -framerate 15 -i "movie%d.png" -pix_fmt yuv420p -c:a copy -crf 20 -r 15 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -codec:v libx264 "movie.mov"
    fslice=$(printf %06d\\n $Slice)
    cp movie.mov ${CURRDIR}/${OUTFILE}_r$fslice.mov
  fi
else
  for file in ../*idx*.png
  do
    ln -s $file  movie$i.png
    i=$((i+1))
  done
  ffmpeg -framerate 15 -i "movie%d.png" -pix_fmt yuv420p -c:a copy -crf 20 -r 15 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -codec:v libx264 "movie.mov"
  cp movie.mov ${CURRDIR}/${OUTFILE}_r000001.mov
fi

cd ..
rm -fr tmp


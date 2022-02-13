#!/bin/bash

cse=1
width=0.0
for mean in 0.005 0.01 0.015 0.02 0.025 0.03 0.035 0.04 0.045 0.05
do
   
   cd src/ #go to src; in src folder
   sed -i s/"double   hchi_mean         =  .*;"/"double   hchi_mean         =  -$mean;"/g global.c
   sed -i s/"double   hchi_width        =  .*;"/"double   hchi_width        =  $width;"/g global.c
   cd .. #in MD_hchi folder
   make
   make install
   mkdir -p results_in_here/case1/width_$width-mean_$mean/
   cp results/* results_in_here/case$cse/width_$width-mean_$mean/
   cd results_in_here/case$cse/width_$width-mean_$mean/ #in case/width... folder
   nohup ./rolat.out &
   cd ../../../ # in MD_hchi folder
done
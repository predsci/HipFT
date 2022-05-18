#!/bin/bash
set -x

HIPFT_DIR="/home/sumseq/Dropbox/PSI/SWQU/OFT/codes/hipft/git"
INPUT_TEMPLATE_DIR="/home/sumseq/Dropbox/PSI/SWQU/OFT/codes/hipft/tests/soccer_ball"
FLOW_METHOD="2"

runname="soccer_advect_rk3_uw"

#Time-steps:
dt_vec=(0.3 0.15 0.075 0.0375 0.01875 0.009375 0.0046875 0.00234375)
#Minimum time-step (stable for minimum dx):
dt_min=0.035
#Default time-step (stable for default dx):
dt_default=0.30

#Resolutions (inverse as these are number of theta points):
dx_inv_vec=(16 32 64 128 256 512 1024 2048)
#Maximum resoluton (inverse) to use to isolate time error:
dx_max_inv=4096
#Default resolution (inverse) to use for time-step scaling runs:
dx_default_inv=512

omega_km_s=1.8076624395655556
# 1.5 deg/day (max|| of DR) ~ 0.21 km/s
# ~10 m/s for MR  

#First do maximum resolution run with stable time-step:

#  Stable DT: 0.038970, so use 0.035

if [ 1 == 1 ] 
then
dx_inv=${dx_max_inv}
p=$((${dx_inv}*2))
dt_inv=${dt_min}
output_tag="DT${dt_inv}_DX${dx_inv}"
input_file=${INPUT_TEMPLATE_DIR}/soccer_ball_scipy_pt_t${dx_inv}_p${p}.h5
${HIPFT_DIR}/trunk/src/hipft $input_file -vrun -time 672 -flow -fm ${FLOW_METHOD} -vpomega ${omega_km_s} -dtmax ${dt_inv} ${runname} 1>${runname}_${output_tag}_hipft.log 2>${runname}_${output_tag}_hipft.err
${HIPFT_DIR}/trunk/bin/hipft_plot_histories.py
mv history_flux.png ${runname}_${output_tag}_history_flux.png
mv history_br.png ${runname}_${output_tag}_history_br.png
mv history_val.png ${runname}_${output_tag}_history_val.png
mv history_num.dat ${runname}_${output_tag}_history_num.dat
mv history_sol.dat ${runname}_${output_tag}_history_sol.dat
mv ${runname}_initial.h5 ${runname}_${output_tag}_initial.h5
mv ${runname}_final.h5 ${runname}_${output_tag}_final.h5
rm ${runname}_final_analytic.h5 
rm ${runname}_initial_0.h5
diffh ${runname}_${output_tag}_initial.h5 ${runname}_${output_tag}_final.h5 > ${runname}_${output_tag}_diffh.out
plot2d -tp -ll -finegrid -cmin -1000 -cmax 1000 ${runname}_${output_tag}_initial.h5        -o ${runname}_${output_tag}_initial.png
plot2d -tp -ll -finegrid -cmin -1000 -cmax 1000 ${runname}_${output_tag}_final.h5          -o ${runname}_${output_tag}_final.png
fi


#Now do minimum time-step run with default resolution:
if [ 1 == 1 ] 
then
dx_inv=${dx_default_inv}
p=$((${dx_inv}*2))
dt_inv=${dt_min}
output_tag="DT${dt_inv}_DX${dx_inv}"
input_file=${INPUT_TEMPLATE_DIR}/soccer_ball_scipy_pt_t${dx_inv}_p${p}.h5
${HIPFT_DIR}/trunk/src/hipft $input_file -vrun -time 672 -flow -fm ${FLOW_METHOD} -vpomega ${omega_km_s} -dtmax ${dt_inv} ${runname} 1>${runname}_${output_tag}_hipft.log 2>${runname}_${output_tag}_hipft.err
${HIPFT_DIR}/trunk/bin/hipft_plot_histories.py
mv history_flux.png ${runname}_${output_tag}_history_flux.png
mv history_br.png ${runname}_${output_tag}_history_br.png
mv history_val.png ${runname}_${output_tag}_history_val.png
mv history_num.dat ${runname}_${output_tag}_history_num.dat
mv history_sol.dat ${runname}_${output_tag}_history_sol.dat
mv ${runname}_initial.h5 ${runname}_${output_tag}_initial.h5
mv ${runname}_final.h5 ${runname}_${output_tag}_final.h5
rm ${runname}_final_analytic.h5 
rm ${runname}_initial_0.h5
diffh ${runname}_${output_tag}_initial.h5 ${runname}_${output_tag}_final.h5 > ${runname}_${output_tag}_diffh.out
plot2d -tp -ll -finegrid -cmin -1000 -cmax 1000 ${runname}_${output_tag}_initial.h5        -o ${runname}_${output_tag}_initial.png
plot2d -tp -ll -finegrid -cmin -1000 -cmax 1000 ${runname}_${output_tag}_final.h5          -o ${runname}_${output_tag}_final.png
fi

#Now do runs over all resolutions at dt_min:
if [ 1 == 1 ] 
then
for i in ${!dx_inv_vec[@]}; do
  dx_inv=${dx_inv_vec[$i]}
  p=$((${dx_inv}*2))
  dt_inv=${dt_min}
  output_tag="DT${dt_inv}_DX${dx_inv}"
  echo ${output_tag}
  input_file=${INPUT_TEMPLATE_DIR}/soccer_ball_scipy_pt_t${dx_inv}_p${p}.h5
  ${HIPFT_DIR}/trunk/src/hipft $input_file -vrun -time 672 -flow -fm ${FLOW_METHOD} -vpomega ${omega_km_s} -dtmax ${dt_inv} ${runname} 1>${runname}_${output_tag}_hipft.log 2>${runname}_${output_tag}_hipft.err
  ${HIPFT_DIR}/trunk/bin/hipft_plot_histories.py
  mv history_flux.png ${runname}_${output_tag}_history_flux.png
  mv history_br.png ${runname}_${output_tag}_history_br.png
  mv history_val.png ${runname}_${output_tag}_history_val.png  
  mv history_num.dat ${runname}_${output_tag}_history_num.dat
  mv history_sol.dat ${runname}_${output_tag}_history_sol.dat
  mv ${runname}_initial.h5 ${runname}_${output_tag}_initial.h5
  mv ${runname}_final.h5 ${runname}_${output_tag}_final.h5
  rm ${runname}_final_analytic.h5 
  rm ${runname}_initial_0.h5
  diffh ${runname}_${output_tag}_initial.h5 ${runname}_${output_tag}_final.h5 > ${runname}_${output_tag}_diffh.out
  plot2d -tp -ll -finegrid -cmin -1000 -cmax 1000 ${runname}_${output_tag}_initial.h5        -o ${runname}_${output_tag}_initial.png
  plot2d -tp -ll -finegrid -cmin -1000 -cmax 1000 ${runname}_${output_tag}_final.h5          -o ${runname}_${output_tag}_final.png
done
fi

#Now do runs over all time-steps:

if [ 1 == 1 ] 
then
for i in ${!dt_vec[@]}; do
  dt_inv=${dt_vec[$i]}
  dx_inv=${dx_default_inv}
  p=$((${dx_inv}*2))
  output_tag="DT${dt_inv}_DX${dx_inv}"
  echo ${output_tag}
  input_file=${INPUT_TEMPLATE_DIR}/soccer_ball_scipy_pt_t${dx_inv}_p${p}.h5
  ${HIPFT_DIR}/trunk/src/hipft $input_file -vrun -time 672 -flow -fm ${FLOW_METHOD} -vpomega ${omega_km_s} -dtmax ${dt_inv} ${runname} 1>${runname}_${output_tag}_hipft.log 2>${runname}_${output_tag}_hipft.err
  ${HIPFT_DIR}/trunk/bin/hipft_plot_histories.py
  mv history_flux.png ${runname}_${output_tag}_history_flux.png
  mv history_br.png ${runname}_${output_tag}_history_br.png
  mv history_val.png ${runname}_${output_tag}_history_val.png  
  mv history_num.dat ${runname}_${output_tag}_history_num.dat
  mv history_sol.dat ${runname}_${output_tag}_history_sol.dat
  mv ${runname}_initial.h5 ${runname}_${output_tag}_initial.h5
  mv ${runname}_final.h5 ${runname}_${output_tag}_final.h5
  rm ${runname}_final_analytic.h5 
  rm ${runname}_initial_0.h5
  diffh ${runname}_${output_tag}_initial.h5 ${runname}_${output_tag}_final.h5 > ${runname}_${output_tag}_diffh.out
  plot2d -tp -ll -finegrid -cmin -1000 -cmax 1000 ${runname}_${output_tag}_initial.h5        -o ${runname}_${output_tag}_initial.png
  plot2d -tp -ll -finegrid -cmin -1000 -cmax 1000 ${runname}_${output_tag}_final.h5          -o ${runname}_${output_tag}_final.png
done
fi

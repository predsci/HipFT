#!/bin/bash

set -x

TEST=test2_1cr_v3_strang_rk3tvd
VISC=200.0
OUTCAD=2.0
INITMAP=/home/sumseq/Dropbox/PSI/SWQU/OFT/codes/hipft/tests/data_assimilation_dev/test1_1day/init_map/br_final_pt.h5
INDATACSV=/home/sumseq/Dropbox/PSI/SHARED-INTERNS/OFT/data/download_results/download_results_2021-01-02T00_00_00_28day.csv
INDATADIR=/home/sumseq/Dropbox/PSI/SHARED-INTERNS/OFT/data/hmi_map
OUTDATADIR=/data/OFT/runs/$TEST
TMAX=672.0

./hipft $INITMAP -flow -dr 1 -mf 1 -va -diff -visc $VISC -diffusion_subcycles 1 -fm 2 -omct $OUTCAD -assimilate_data -assimilate_data_map_list_filename $INDATACSV -assimilate_data_map_root_dir $INDATADIR -output_map_directory $OUTDATADIR -time $TMAX $TEST 1>${TEST}.log 2>${TEST}.err 

~/Dropbox/PSI/SWQU/OFT/codes/hipft/git/trunk/bin/hipft_plot_histories.py

~/Dropbox/PSI/SWQU/OFT/codes/hipft/git/trunk/bin/make_hipft_plots.sh $OUTDATADIR $TEST

cp ${OUTDATADIR}/plots/${TEST}.mov .






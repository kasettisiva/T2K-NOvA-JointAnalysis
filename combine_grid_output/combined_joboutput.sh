#!/bin/bash

#output_dir=${1}
osc_par=${1}
asimov=${2}
date=${3}

expt=("T2K" "NOvA" "T2KNOvA")
#asimov=("asimov1" "asimov0" "asimov4")

output_dir="margTemp_${osc_par}_wRC_t2knova${asimov}_${date}/"
root -l -b -q "MergeToyMargTemplates.cxx(\"${output_dir}/root_files/*\",\"${output_dir}/allfits.root\")"

for iexpt in 0 1 2; do

root -l -b -q "PlotSensiAsimov_1D_parallel_modified.C(\"${output_dir}\",\"${expt[${iexpt}]}\")"
root -l -b -q "PlotSensiAsimov_1D_combined.C(\"${output_dir}\",\"${expt[${iexpt}]}\",\"100k\",\"${output_dir}\")"

done

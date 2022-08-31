#!/bin/bash

#source /cvmfs/juno.ihep.ac.cn/centos7_amd64_gcc830/Pre-Release/J20v1r1-branch/bashrc.sh
source /cvmfs/sft.cern.ch/lcg/views/LCG_96/x86_64-slc6-gcc8-opt/setup.sh


#let bin_start="${start_bin}*100"
let start=0
nbins=102

path=`pwd`

let bin_start="${start_bin}*1000"

#toyxp="Validation_ToyXP_AsimovA_2018_500_Asimovs_01Feb2022.root"
toyxp="Validation_ToyXP_AsimovA_2018_16Feb2022.root"
templates="$(printf "subset%04d.root" "${bin_start}")"
outbase="$(printf "bins%04d" "${bin_start}")"
osc_throw="$(printf "Osc_wRC_%04d.root" "${bin_start}")"

#templates="$(printf "subset%04d.root" "${bin_start}")"

echo "Templates: $templates"

#osc_throw="Osc_wRC_20k.root"

oscparfile="osc_OA2019_dCP.yaml"

#outbase="$(printf "bins%04d" "${bin_start}")"

cd $path/P-theta/Minimal
cp $path/nova-sl7-novat2k_v4_fixcosmicsrock.sif ./

ln -sf $path/P-theta/Minimal/NOvA_files/jf_data_asimov1_fake/data_nova2020ana_fake_asimov1.root $path/P-theta/Minimal/NOvA_files/jf_data_asimov1_fake/data.container.root

$path/P-theta_install/bin/PTThrowOsc \
  -c $path/P-theta/Minimal/scripts/script2_OA2019_temp.yaml \
  --osc-param-file $path/P-theta/Minimal/scripts/osc_OA2019_throwosc_RC.yaml \
  -o "${path}/${osc_throw}" \
  -n 1000 \

$path/P-theta_install/bin/PTMargTemplates \
 -c $path/P-theta/Minimal/scripts/script2_OA2019_temp.yaml \
 --osc-param-file "$path/P-theta/Minimal/scripts/$oscparfile" \
 --osc-file "${path}/$osc_throw" \
 --templates-file "${path}/$templates" \
 --data-file "${path}/$toyxp" \
 --start-bin "$start" \
 --number-bins "$nbins" \
 --t2k-nova-joint-fit 1 \
 --number-toy-xp 1 \
 -o "$path/${outbase}.root" \

cd -

#          --stat-only 1 \
echo "Done making sensitivity.!"
        

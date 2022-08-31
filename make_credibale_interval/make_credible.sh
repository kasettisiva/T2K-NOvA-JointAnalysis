#!/bin/bash

path=${1}
nbins=${2}
for expt in "T2K" "NOvA" "T2KNOvA"; do
    [ -d "${path}/credible_${expt}" ] || mkdir "${path}/credible_${expt}"
    for mh in "" "_IH"; do
        [ -d "${path}/credible_${expt}${mh}" ] || mkdir "${path}/credible_${expt}${mh}"
	./Bicubic_1D  -f ${path}/hist_${expt}${mh}.root -o ${path}/smoothed_${expt}${mh}.root -n ${nbins} -h cont
        ./Credible_1D -f ${path}/smoothed_${expt}${mh}.root  -o ${path}/credible_${expt}${mh} -t chi2 -g 1 
    done
    ./Credible_1D_both -fn ${path}/smoothed_${expt}.root -fi ${path}/smoothed_${expt}_IH.root -o ${path}/credible_${expt} -t chi2 -g 1
done 

#!/bin/bash

#echo "Getting jobid .."
i=0
while read line;
do
i=$((i+1))
#echo $i
 if [[ $(( $i % 2 )) -ne 0 ]]; then
   continue
 fi
#echo $line
#file=`echo $line | cut -d'/' -f9`
#echo $file
#jobid=`echo $line | cut -d':' -f15 | cut -d'}' -f1`
#jobid=`echo $line | cut -d':' -f16 | cut -d'}' -f1`
jobid=`echo $line | cut -d':' -f17 | cut -d'}' -f1`
echo $jobid
#python monitor_dirac_job.py ${jobid}
python output_dirac_job.py ${jobid}
done < jobid_dCP_1D_28June2022.txt

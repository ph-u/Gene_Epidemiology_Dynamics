#!/bin/env bash
# author: ph-u
# script: sangerHPC.sh
# desc: HPC script for transposon model
# in: bash sangerHPC.sh
# out: NA
# arg: 0
# date: 20260702, 20260824

PATH="/software/isg/private/wrappers/apptainer/1.4.0:$PATH"

mkdir -p ../data

[[ -f gene_epidemiology_dynamics_latest.sif ]] && rm gene_epidemiology_dynamics_latest.sif
apptainer pull docker://ghcr.io/ph-u/gene_epidemiology_dynamics:latest


if [[ $USER!="ubuntu" ]];then
  sed -e "s/mAx/$(( `wc -l < ../upload/seed.csv` ))/" sangerHPC_child.sh > sC.sh
  bsub < sC.sh
else
  lIm=3
  sed -e "s/LSB_JOBINDEX/1/" -e "s/ --rm//" sangerHPC_child.sh > sC.sh
  i0=`wc -l < ../data/seed.csv`
  for i in `seq 1 ${i0}`;do
    while (( $(jobs -rp | wc -l) >= lIm ));do
      wait -n
    done
    bash sC.sh ${i} 1> ../work/pj01-${i}.txt &
    echo -e "`date`: Submitted ${i}"
  done
fi

exit

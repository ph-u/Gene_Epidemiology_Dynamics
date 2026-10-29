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

sed -e "s/mAx/$(( `wc -l < ../upload/seed.csv` ))/" sangerHPC_child.sh > sC.sh

bsub < sC.sh

exit

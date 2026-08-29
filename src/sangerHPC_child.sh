#!/bin/env bash
# author: ph-u
# script: sangerHPC_child.sh
# desc: HPC script for NFDS individual-based model
# in: bsub sangerHPC_child.sh
# out: NA
# arg: 0
# date: 20260702, 20260824

#BSUB -G team377f
#BSUB -o ../work/pj01-%J-%I.o
#BSUB -e ../work/pj01-%J-%I.e
#BSUB -q basement
#BSUB -n 10
#BSUB -M 6000
#BSUB -R "select[mem>6000] rusage[mem=6000] span[hosts=1]"
#BSUB -J "gdy[1-mAx]"

PATH="/software/isg/private/wrappers/apptainer/1.4.0:$PATH"
export NFDS_SRC=${PWD}/src

apptainer run --bind ${PWD}/../data:/data --pwd /src gene_epidemiology_dynamics_latest.sif bash run.sh ${LSB_JOBINDEX}

exit
## any BSUB -J with the correct format, no matter at any line, later line will replace earlier lines # BSUB -J "tdy[1-360]"

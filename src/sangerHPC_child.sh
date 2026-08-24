#!/bin/env bash
# author: ph-u
# script: sangerHPC_child.sh
# desc: HPC script for NFDS individual-based model
# in: bsub sangerHPC_child.sh
# out: NA
# arg: 0
# date: 20260702, 20260824

#BSUB -G team377f
#BSUB -o ../work/pj02-%J-%I.o
#BSUB -e ../work/pj02-%J-%I.e
#BSUB -q normal
#BSUB -n 10
#BSUB -M 4000
#BSUB -R "select[mem>4000] rusage[mem=4000] span[hosts=1]"
#BSUB -J "gdy[1-mAx]"

PATH="/software/isg/private/wrappers/apptainer/1.4.0:$PATH"

apptainer run --bind ${PWD}/../data:/data --pwd /src gene_epidemiology_dynamics_latest.sif bash run.sh ${LSB_JOBINDEX}

exit
## any BSUB -J with the correct format, no matter at any line, later line will replace earlier lines # BSUB -J "tdy[1-360]"

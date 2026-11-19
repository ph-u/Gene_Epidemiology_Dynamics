#!/bin/env bash
# author: ph-u
# script: run.sh
# desc: run ABC-SMC.NFDS individual-based model parameter estimations: pf, fSel, wSel, vSel, m
# in: bash run.sh [seed number]
# out: NA
# arg: 1
# date: 20260824

export NFDS_SRC=/src
#Rscript model.r ../data/input.csv $1
julia --project=$NFDS_SRC $NFDS_SRC/ABC_SMC_NFDS.jl $1

exit

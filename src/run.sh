#!/bin/env bash
# author: ph-u
# script: run.sh
# desc: run ABC-SMC.NFDS individual-based model parameter estimations: pf, fSel, wSel, vSel, m
# in: bash run.sh [seed number]
# out: NA
# arg: 1
# date: 20260824

export NFDS_SRC=/src
Rscript model.r ../raw/input.csv $1

exit


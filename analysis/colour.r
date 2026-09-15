#!/bin/env Rscript
# author: ph-u
# script: colour.r
# desc: customized colours
# in: source("colour.r")
# out: NA
# arg: 0
# date: 20260903

##### colour #####
cBp = c(); set.seed(123); for(i in c("Okabe-Ito", "alphabet", "polychrome 36", "dark 2", "set 1", "classic tableau")){
  if(i == "Okabe-Ito"){
    cBp = c(cBp, rev(palette.colors(palette = i, alpha=1, recycle = F)))
  }else{
    cBp.t = palette.colors(palette = i, alpha=1, recycle = F)
    cBp = c(cBp, sample(cBp.t, length(cBp.t)))
  }};rm(i, cBp.t);cBp = unique(cBp)


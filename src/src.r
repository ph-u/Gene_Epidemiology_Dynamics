#!/bin/env Rscript
# author: ph-u
# script: src.r
# desc: technical functions of model.r
# in: source("src.r")
# out: NA
# arg: 0
# date: 20260819

##### get values from input.csv #####
g = function(x, f = f.in){return(f$Value[f$Type==x])}

##### barcoding #####
bcod = function(df = unique(d0[,-(1:2)])){
  nC = ceiling(log(nrow(df)) / log(length(LETTERS)))
  b = data.frame(a = LETTERS, b = rep(LETTERS, each = length(LETTERS)))
  if(nC > 2){
    for(i in 1:(nC-2)){
      b = cbind(b, rep(LETTERS, each = nrow(b)))
    }
  }
  return(apply(b,1,paste, collapse = "")[1:nrow(df)])
}

##### VT*SC distribution vector #####
vtsc = function(idx, cls = vtsc.idx, nLev = length(vtsc.lev)){
  return(tabulate(cls[idx], nbins = nLev))
}

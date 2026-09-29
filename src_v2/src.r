#!/usr/bin/env Rscript
# author: ph-u
# script: src.r
# desc: helper functions for the NFDS ABCSMC model v2
# in: source("src.r")
# out: NA
# arg: 0
# date: 20260919

##### Read a value from the input CSV #####
## Like looking up a word in a dictionary: give it the name, get back the value.
g = function(x, f = f.in){ return(f$Value[f$Type == x]) }

##### Give each unique genotype a short letter code #####
## Like giving each student a unique student-ID card made of letters.
bcod = function(df = unique(d0[, -(1:2)])){
  nC = ceiling(log(nrow(df)) / log(length(LETTERS)))
  b  = data.frame(a = LETTERS, b = rep(LETTERS, each = length(LETTERS)))
  if(nC > 2) for(i in seq_len(nC - 2)) b = cbind(b, rep(LETTERS, each = nrow(b)))
  return(apply(b, 1, paste, collapse = "")[1:nrow(df)])
}

##### Count bacteria by (vaccine-type profile x SC x deme) #####
## Like tallying students by (subject x year group x school) all at once.
## idx = genotype row index of each bacterium
## dm  = deme index of each bacterium
vtsc = function(idx, dm, nLev = length(vtsc.lev)){
  tabulate((dm - 1L) * nLev + vtsc.idx[idx], nbins = nLev * nD)
}

##### Summarise an infection history table #####
## Reads the history book of cleared infections and prints headline statistics:
## how many were recorded, how long infections lasted, and how often hosts
## carried more than one bacterial type at once (co-infection).
hist_stats = function(hH){
  if(!is.data.frame(hH) || nrow(hH) == 0)
    return(data.frame(n_inf = 0L, n_host = 0L, mean_dur = NA_real_,
                      median_dur = NA_real_, pct_coInf = NA_real_,
                      mean_coInf_size = NA_real_))
  hH$dur = hH$tClr - hH$tInf
  data.frame(
    n_inf           = nrow(hH),
    n_host          = length(unique(hH$hID)),
    mean_dur        = mean(hH$dur,   na.rm = TRUE),
    median_dur      = median(hH$dur, na.rm = TRUE),
    pct_coInf       = 100 * mean(hH$nCoInf > 1L, na.rm = TRUE),
    mean_coInf_size = { cx = hH$nCoInf[hH$nCoInf > 1L]
                        if(length(cx)) mean(cx) else NA_real_ })
}


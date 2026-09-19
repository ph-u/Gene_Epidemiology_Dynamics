#!/bin/env Rscript
# author: ph-u
# script: dyRecap.r
# desc: Plot check if dynamics recap published dynamics
# in: Rscript dyRecap.r
# out: res/dyRecap--*
# arg: 0
# date: 20260901

dAte = "20260919"

##### env #####
source("colour.r")
f.in = list.files(paste0("../data/",dAte), pattern = "all_", full.names = T)
sEed = read.table(text = gsub("_","-",basename(f.in)), sep = "-")[,2]
for(i in seq_len(length(f.in))){
  d0 = read.csv(f.in[i], header = T)
  d0$f.in = sEed[i]
  if(i>1){d.all = rbind(d.all, d0)}else{d.all = d0}
};rm(i, f.in, d0)
d.all$plot = d.all$dist1 < quantile(d.all$dist1, probs = .1)
# d.all$plot = d.all$dist1 %in% d.all$dist1[order(d.all$dist1)[seq_len(9000)]]

## Re-transformation to parameter values in paper
# d.all$fSelected = exp(d.all$fSelected)
# d.all$wSelected = exp(d.all$wSelected)
d.all$vSelected = d.all$vSelected * .5
d.all$migration = d.all$migration * .2

##### Mark parameter values in paper #####
tRuth = c(pf = 0.2483, sigma_f = 0.1363, sigma_w = 0.0023, sigma_v = 0.0812, m = 0.0044) # 2017 publication
t1 = c(pf = 0.1197, sigma_f = 0.0213, sigma_w = 0.0010, sigma_v = 0.0491, m = 0.0015)
t2 = c(pf = 0.5448, sigma_f = 0.2113, sigma_w = 0.0514, sigma_v = 0.1254, m = 0.0165)
# tRuth[3] = log(tRuth[3])
# t1[3] = log(t1[3])
# t2[3] = log(t2[3])

jpeg("../res/dyRecap--parHist.jpeg", height = 1200, width = 1600, res = 300)
par(mfrow = c(3,2), mar = c(5,4,1,0)+.1)
for(i in seq_len(length(tRuth))){
  for(i0 in seq_len(length(sEed))){
    if(i %in% 3){
      xLim = c(0,.05)
    }else{
      xLim = range(quantile(d.all[,i+2], probs = c(0,1)), tRuth[i], t1[i], t2[i])
      print(xLim)
    }
    hist(d.all[which(d.all$f.in==sEed[i0] & d.all$plot == T),i+2], xlim = xLim, breaks = 100, freq = T, xlab = paste0(ifelse(i %in% 0,"log( ", ""), colnames(d.all)[i+2], ifelse(i %in% 0," )", "")), main = "", border = NA, col = paste0(substr(cBp[i0],1,7),"33"), add = !(i0==1))
  };rm(i0)
  abline(v = tRuth[i], lwd = 2, col = "#ff00ffff")
  abline(v = t1[i], lwd = 2, lty = 2, col = "#ff00ffff")
  abline(v = t2[i], lwd = 2, lty = 2, col = "#ff00ffff")
};rm(i, xLim)
invisible(dev.off())

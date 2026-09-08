#!/bin/env Rscript
# author: ph-u
# script: dyRecap.r
# desc: Plot check if dynamics recap published dynamics
# in: Rscript dyRecap.r
# out: res/dyRecap--*
# arg: 0
# date: 20260901

d0 = read.csv("../data/all_accepted_particles.csv", header = T)

tRuth = c(pf = 0.2483, sigma_f = 0.1363, sigma_w = 0.0023, sigma_v = 0.0812, m = 0.0044) # 2017 publication
t1 = c(pf = 0.1197, sigma_f = 0.0213, sigma_w = 0.0010, sigma_v = 0.0491, m = 0.0015)
t2 = c(pf = 0.5448, sigma_f = 0.2113, sigma_w = 0.0514, sigma_v = 0.1254, m = 0.0165)

jpeg("../res/dyRecap--parHist.jpeg", height = 1200, width = 1600, res = 300)
par(mfrow = c(3,2), mar = c(5,4,1,0)+.1)
for(i in seq_len(length(tRuth))){
  hist(d0[,i+2], xlim = range(d0[,i+2], tRuth[i], t1[i], t2[i]), breaks = 25, freq = T, xlab = colnames(d0)[i+2], main = "", col = "#00000000")
  abline(v = tRuth[i], lwd = 2, col = "#ff00ffff")
  abline(v = t1[i], lwd = 2, lty = 2, col = "#ff7700bb")
  abline(v = t2[i], lwd = 2, lty = 2, col = "#ff7700bb")
};rm(i)
invisible(dev.off())

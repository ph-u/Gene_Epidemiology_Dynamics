#!/bin/env Rscript
# author: ph-u
# script: convergenceCheck.r
# desc: check for SMC convergence
# in: Rscript convergenceCheck.r
# out: NA
# arg: 0
# date: 20260909

dAte = "20260919"

##### env #####
source("colour.r")
f.in0 = list.files(paste0("../data/",dAte), pattern = "thres", full.names = T)
f.in1 = list.files(paste0("../data/",dAte), pattern = "all_", full.names = T)
f.in2 = list.files(paste0("../data/",dAte), pattern = "last_", full.names = T)
sEed = read.table(text = gsub("_","-",basename(f.in0)), sep = "-")[,2]
for(i in seq_len(length(f.in0))){
  d0 = read.csv(f.in0[i], header = T)
  d0$f.in = sEed[i]
  d0$drop = c(NA, -diff(d0$dist1))
  d0$dropRatio = c(NA, d0$drop[-1]/d0$dist1[-nrow(d0)])
  d1 = read.csv(f.in1[i], header = T)
  if(i==1){pArams = colnames(d1)[3:7]; par.Nam = c(colnames(d1)[1], paste0(rep(pArams, each = 2), ".", c("median", "sd")))}
  d1.agg = as.data.frame(matrix(unlist(aggregate(d1[pArams], list(gen = d1$gen), function(v) c(med = median(v, na.rm = T), sd = sd(v)))), ncol = length(par.Nam)))
  colnames(d1.agg) = par.Nam
  d1$f.in = d1.agg$f.in = sEed[i]
  if(i>1){
    d0.all = rbind(d0.all, d0)
    d1.all = rbind(d1.all, d1)
    d1.allAgg = rbind(d1.allAgg, d1.agg)
  }else{
    d0.all = d0
    d1.all = d1
    d1.allAgg = d1.agg
  }
};rm(i, f.in0, d0, f.in1, d1, d1.agg)

##### Check plateau #####
d0 = reshape(d0.all[,c(1,4,ncol(d0.all))], direction = "wide", timevar = "f.in", idvar = "gen")
pdf("../res/convergenceCheck--dropRatio.pdf")
par(mar = c(4,4,0,0)+.1)
matplot(x = d0[,1], y = d0[,-1], type = "l", lty = 1, xlab = "SMC generation", ylab = "Data-simulation discrepancy drop ratio")
invisible(dev.off())

##### Posterior stability #####
d1 = reshape(d1.allAgg, direction = "wide", timevar = "f.in", idvar = "gen")
pdf("../res/convergenceCheck--params.pdf")
par(mfrow = c(5,2), mar = c(4,4,1,0)+.1)
for(i in seq_len(length(par.Nam)-1)){
  i0 = d1[,c(1,grep(par.Nam[i+1], colnames(d1)))]
  matplot(x = i0[,1], y = i0[,-1], type = "l", lty = 1, xlab = "SMC generation", ylab = par.Nam[i+1])
};rm(i, i0)
invisible(dev.off())

##### Between-seed agreement #####
d2 = do.call(rbind, lapply(seq_along(f.in2), function(i) transform(read.csv(f.in2[i], header = T), seed = sEed[i])))
pdf("../res/convergenceCheck--paramsLastAccept.pdf")
par(mfrow = c(3,2), mar = c(4,4,1,0)+.1)
for(i in seq_len(length(pArams))){
  boxplot(d2[,which(colnames(d2)==pArams[i])] ~ d2$seed, col = "#00000000", xlab = "Simulations", ylab = pArams[i], xaxt = "n")
};rm(i)
invisible(dev.off())

sapply(pArams, function(p){
  m = tapply(d2[[p]], d2$seed, median)
  c(between_sd = sd(m), within_sd = mean(tapply(d2[[p]], d2$seed, sd)))
})

##### Effective sample size ##### aim: >10% of particle count
(sum(d2$pWeight)^2 / sum(d2$pWeight^2))/nrow(d2)

##### Recovery on synthetic data #####
sapply(pArams, function(p){
  c(p.025 = tapply(d2[[p]], d2$seed, quantile, probs = .025), p.500 = tapply(d2[[p]], d2$seed, quantile, probs = .5), p.975 = tapply(d2[[p]], d2$seed, quantile, probs = .975))
})
pdf("../res/convergenceCheck--paramsPairwiseLastAccept.pdf")
# plot(d1.all, col = "#00000022") # file size = 187.9 MB
plot(d2[,3:7], cex = .1)
invisible(dev.off())

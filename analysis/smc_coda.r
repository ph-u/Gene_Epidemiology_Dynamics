#!/bin/env Rscript
# author: ph-u
# script: smc_coda.r
# desc: coda package MC chain test
# in: Rscript smc_coda.r
# out: diagnostic tests
# arg: 0
# date: 20260925

argv = (commandArgs(T))
dAte = argv[1]

library(coda)
pr <- c("propStrong","fSelected","wSelected","vSelected","migration")

ff <- Sys.glob(paste0("../data/",dAte,"/run_*-last_accepted_particles.csv"))
mc <- mcmc.list(lapply(ff, function(f) {
  z <- read.csv(f)
  mcmc(as.matrix(z[, pr]))          # one "chain" = one seed's final particles
}))

print(gelman.diag(mc, autoburnin = FALSE, multivariate = FALSE))

print(summary(mc))                 # per-parameter quantiles; SE column is unreliable (assumes a chain)

png(paste0("../res/mcChainAnalysis_v1_",dAte,".png"), width = 1500, height = 1500, res = 300)
par(mar = c(5,3,3,0)+.1)
plot(mc, trace = FALSE)     # densities only -- trace plots are meaningless for particles
invisible(dev.off())

print(crosscorr(mc))               # the propStrong/wSelected correlation, properly across seeds
print(HPDinterval(mc, prob = 0.95))

print(effectiveSize(mc))     # assumes autocorrelation; particles are independent
print(autocorr.diag(mc))     # same
print(raftery.diag(mc))      # quantile accuracy for a chain
print(heidel.diag(mc))       # stationarity of a time series
print(geweke.diag(mc))       # compares chain segments -- generations aren't segments

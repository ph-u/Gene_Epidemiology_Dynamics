#!/bin/env Rscript
# author: ph-u
# script: nfds.r
# desc: NFDS model
# in: source("nfds.r")
# out: NA
# arg: 0
# date: 20260820

##### NFDS simulation #####
m.nfds = function(propStrong, fSelected, wSelected, vSelected, migration, meanStandardize = F, keepGenotypes = F){

  ##### Runtime acceleration #####
  hEad = 1024
  G     = rbind(G0, matrix(0, hEad, ncol(G0))); storage.mode(G) = "double"
  VTu   = c(vt0,  rep(0, hEad));  SCu   = c(sc0,  rep(NA, hEad))
  mProb = c(mP0,  rep(0, hEad));  tagU  = c(tag0, rep(NA, hEad))
  clsU  = c(vtsc.idx, rep(NA, hEad))
  bOrn  = c(rep(0, nrow(G0)), rep(NA, hEad)) # generation of origin
  pAra  = c(rep(NA_integer_, nrow(G0) + hEad)) # parent genotype row
  nU    = nrow(G0) # active rows

  ##### selection pressure per gene #####
  sEl = rep(log1p(wSelected), length(gNam))
  sEl[order(selMode$strength)[seq_len(floor(length(gNam) * propStrong))]] = log1p(fSelected)
  
  ##### Initial generation #####
  num.Infect = ceiling(k * as.numeric(g("percentage initial infected")) / 100)
  gen1 = sample(tag.pre, num.Infect, replace = T)
  rec.eQm = matrix(0, nrow = length(vtsc.lev), ncol = length(eQm.date)); j = 1
  
  ##### Later generations #####
  tIme = (min(eQm$Month) + 1):nGen
  for(i in seq_len(length(tIme))){
    t2r = drop(tabulate(gen1, nbins = nrow(G)) %*% G) / length(gen1) # https://www.geeksforgeeks.org/r-language/r-operators/; https://www.mathsisfun.com/algebra/matrix-multiplying.html
    fIt.u = exp(drop(G %*% ((eqm.pre - t2r) * sEl))) * (1 - vSelected * vCurve[i+1] * VTu) * (1 - migration) # fitness for unique tags
    pi.Omega = fIt.u[gen1]
    fIt.adj = k / num.Infect
    if(meanStandardize){ fIt.adj = fIt.adj / mean(pi.Omega) }
    oFf = rpois(num.Infect, pi.Omega * fIt.adj)
    if(anyNA(oFf) || sum(oFf) >= popRunaway){return(NULL)} # defend against runaway population size
    gen1 = rep(gen1, oFf) # new generation offspring
    
  ##### Recombination (future) #####
    
  ##### Migrations #####
    if(sum(oFf) < 1){return(NULL)} # if whole bacterial population wiped out
    if(migration > 0){
      nMig = rbinom(1, k, min(1, migration * k / num.Infect))
      gen1 = c(gen1, sample(migIdx, nMig, replace = T, prob = mP0[migIdx]))
    }
    num.Infect = length(gen1)
        
  ##### Simulation records #####
    if(tIme[i] %in% eQm.date){
      rec.eQm[,j] = vtsc(gen1[sample.int(num.Infect, as.numeric(nObs[as.character(tIme[i])]), replace = T)]) # convert index for tags to index for VT|SC types, and sample
      j = j + 1
    }
  };rm(i)
  
  if(!keepGenotypes){ return(rec.eQm) }
  kEep = seq_len(nU)
  return(list(
    ss = rec.eQm,
    genotypes = data.frame(row = kEep, tag = tagU[kEep], VT = VTu[kEep], SC = SCu[kEep],
                           born = bOrn[kEep], parent = pAra[kEep],
                           finalCount = tabulate(gen1, nbins = nrow(G))[kEep]),
    G = G[kEep, , drop = FALSE]))
}

##### Jensen-Shannon divergence (first translation from cpp model by Claude.ai) #####
jsd = function(p, q){ # p, q = phenotype counts
  p = p / sum(p); q = q / sum(q) # rescale counts into population fraction
  m = (p + q) / 2 # midpoints
  s = 0
  i = m > 0 & p > 0; s = s + 0.5 * sum(p[i] * log(p[i] / m[i]))
  i = m > 0 & q > 0; s = s + 0.5 * sum(q[i] * log(q[i] / m[i]))
  return(s)
}

##### Model using JSD index #####
nfds_jsd = function(x, ss_obs) {
  sim <- m.nfds(x[["propStrong"]], x[["fSelected"]], x[["wSelected"]],
                x[["vSelected"]], x[["migration"]])
  if(is.null(sim) || any(is.na(sim))){ return(ncol(ss_obs) * log(2)) }  # extinction
  return( sum(vapply(seq_len(ncol(ss_obs)),
                     function(j) jsd(sim[,j], ss_obs[,j]), 0)) )
}


#!/bin/env Rscript
# author: ph-u
# script: nfds.r
# desc: NFDS model
# in: source("nfds.r")
# out: NA
# arg: 0
# date: 20260820

pAr = list(
  propStrong = .25,
  fSelected = .1,
  wSelected = .0001,
  vSelected = .2,
  migration = .05
)

##### NFDS simulation #####
m.nfds = function(propStrong, fSelected, wSelected, vSelected, migration, meanStandardize = F){
  k = as.numeric(g("population")); nGen = eQm[nrow(eQm), 1] ## simulation constants
  
  ##### selection pressure per gene #####
  sEl = rep(log1p(wSelected), length(gNam))
  sEl[order(selMode$strength)[seq_len(floor(length(gNam) * propStrong))]] = log1p(fSelected)
  
  ##### Initial generation #####
  num.Infect = ceiling(k * as.numeric(g("percentage initial infected")) / 100)
  gen1 = sample(match(d0$tag[d0$Time <= 0], d0.u$tag), num.Infect, replace = T)
  rec.eQm = c()
  
  ##### Later generations #####
  for(i in seq_len(nGen)){
    G = as.matrix(d0.u[,gNam]); storage.mode(G) = "double" # Matrix for speed
    t2r = drop(tabulate(gen1, nbins = nrow(G)) %*% G) / length(gen1) # https://www.geeksforgeeks.org/r-language/r-operators/; https://www.mathsisfun.com/algebra/matrix-multiplying.html
    fIt.u = exp(drop(G %*% ((eqm.pre - t2r) * sEl))) * (1 - vSelected * vCurve[i] * d0.u$VT) * (1 - migration) # fitness for unique tags
    pi.Omega = fIt.u[gen1]
    fIt.adj = k / num.Infect
    if(meanStandardize){ fIt.adj = fIt.adj / mean(pi.Omega) }
    gen1 = rep(gen1, rpois(num.Infect, pi.Omega * fIt.adj)) # new generation offspring
    
  ##### Migrations #####
    if(migration > 0){
      nMig = rbinom(1, k, min(1, migration * k / num.Infect))
      gen1 = c(gen1, sample(migIdx, nMig, replace = T, prob = d0.u$migProb[migIdx]))
    }
    num.Infect = length(gen1)
    
  ##### Simulation records #####
    if(num.Infect < 1){return(NULL)} # if whole bacterial population wiped out
    if(i %in% eQm.date){
      vtsc.lev = sort(unique(d0.u$paste))
      rec.eQm = c(rec.eQm, 
        tabulate(match(d0.u$paste, vtsc.lev)[gen1[sample.int(num.Infect, as.numeric(nObs[as.character(i)]), replace = T)]], nbins = length(vtsc.lev)) # convert index for tags to index for VT|SC types, and sample
        )
    }
  };rm(i)
  return(matrix(rec.eQm, ncol = length(eQm.date)))

  ##### Group genes that face strong / weak selection #####
#  selMode$category[order(selMode$strength)[1:floor(nrow(selMode)*propStrong)]] = "strong"
#  pStrong = rep(log1p(wSelected), length(gNam))
#  pStrong[which(selMode$category=="strong")] = log1p(fSelected)

  ##### Initial population #####
#  num.Infect = ceiling(as.numeric(g("population"))*as.numeric(g("percentage initial infected"))/100)
#  gen0 = sample(d0$tag[d0$Time <= 0], num.Infect, replace = T)
#  gen1 = match(gen0, d0.u$tag)
  # rec.gen0 = rep(NA, eQm[nrow(eQm),1] + 1)
  # rec.gen0[1] = gsub(" ","",paste(apply(as.data.frame(table(gen0)), 1, paste, collapse = "!"), collapse = ";")) # record offspring population structure

  ##### Simulation #####
#  i = 1; rec.eQm = c(); repeat{
#    cat(date(),": gen",i,"\n")

  ##### Matrix for rapid gene matrix calculations #####
#    G = as.matrix(d0.u[,gNam])
#    storage.mode(G) = "double"

  ##### Calculate selection pressure strength #####
#    for(i0 in 1:nrow(d0.u)){
#      pO = piOmega(d0.u$data[i0], eqm.pre - t2r(gen1), pStrong)
#      d0.u[i0,c("d.omega", "d.pi")] = pO[c(2,1)]
#    };rm(i0, pO) # runtime bottleneck
#    pi.omega = exp(d0.u$d.pi[gen1] * log1p(fSelected)) * exp(d0.u$d.omega[gen1] * log1p(wSelected)) * (1 - vSelected * d0.u$VT[gen1] * vCurve[i]) * (1 - migration) ## raw NFDS selection coefficient for each individual
#    offspring = rpois(num.Infect, pi.omega * as.numeric(g("population")) / num.Infect / ifelse(meanStandardize, mean(pi.omega), 1))
    # rec.gen0[i+1] = gsub(" ","",paste(gen0,offspring, sep = "!", collapse = ";")) # record offspring population structure

  ##### Count VT*SC type #####
#    gen0 = rep(gen0, offspring); gen1 = match(gen0, d0.u$tag)

  ##### Migration #####
#    if(migration > 0){
#      gen0 = c(gen0, sample(d0.u$tag[grep("[[:upper:]]", d0.u$tag)], rbinom(1, as.numeric(g("population")), min(1, migration * as.numeric(g("population")) / num.Infect)), replace = T, prob = d0.u$migProb[grep("[[:upper:]]", d0.u$tag)])) # 2017 publication
#      gen1 = match(gen0, d0.u$tag)
#    }
    
#    if(length(gen0) < 1){return(NULL)}
#    if(i %in% eQm.date){rec.eQm = c(rec.eQm, vtsc(sample(gen0, sum(d0$Time == i), replace = T)))}
#    if(i >= eQm[nrow(eQm),1]){break}else{i = i + 1; num.Infect = length(gen0)}
#  }
#  return(matrix(rec.eQm, ncol = length(eQm.date)))
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


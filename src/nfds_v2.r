#!/usr/bin/env Rscript
# author: ph-u
# script: nfds.r
# desc: NFDS model -- metapopulation, multi-vaccine, graded roll-out, co-infection
# in: source("nfds.r")
# out: NA
# arg: 0
# date: 20260917

##### NFDS simulation #####
m.nfds = function(propStrong, fSelected, wSelected, vSel, migration,
                  coInf, vSelMod, mWithin,
                  meanStandardize = F, keepGenotypes = F){

  ##### Parameter transformation (all priors are unif(0,1)) #####
  vSel      = vSel * .5      # per-vaccine efficacy, 0 - 0.5 / month (vector, length nVac)
  migration = migration * .2 # external immigration, 0 - 0.2
  mWithin   = mWithin * .1   # between-deme movement, 0 - 0.1 / month

  U = nrow(G0); L = length(gNam)

  ##### Per-locus NFDS weight (one vector; the C++ cogWeights) #####
  sEl = rep(log1p(wSelected), L)
  sEl[order(selMode$strength)[seq_len(floor(L * propStrong))]] = log1p(fSelected)

  ##### Roll-out coverage, nVac x nD x nMonth, built once per call #####
  ## vSelMod widens the step into a cohort ramp; vSelMod -> 0 recovers a step function
  sLope = max(vSelMod * 24, 1e-6)
  ## linear cohort ramp: exactly 0 before roll-out, 1 after sLope months
  vCov  = vapply(vMonth, function(tt) pmin(1, pmax(0, (tt - vStart) / sLope)), matrix(0, nVac, nD))

  ##### Initial generation: resample pre-vaccination isolates within each deme #####
  nInit = pmax(1, round(kD * as.numeric(g("percentage initial infected")) / 100))
  dEme  = rep(seq_len(nD), nInit)
  gen1  = unlist(lapply(seq_len(nD), function(dd) sample(tag.pre[[dd]], nInit[dd], replace = T)))

  rec.eQm = matrix(0, nrow = length(vtsc.lev) * nD, ncol = length(eQm.date)); j = 1

  ##### Generations #####
  tIme = (min(eQm$Month) + 1):nGen
  for(i in seq_along(tIme)){

    ## per-deme locus frequencies -- one matrix product for the whole metapopulation
    cnt  = vapply(seq_len(nD), function(dd) tabulate(gen1[dEme == dd], nbins = U), numeric(U))
    nPer = colSums(cnt)
    if(any(nPer < 1)){ return(NULL) }                 # a deme went extinct
    t2r  = t(crossprod(cnt, G0) / nPer)               # L x nD

    ## vaccine pressure: product over vaccines, deme-specific timing
    vNow   = vSel * matrix(vCov[, , i + 1], nVac, nD) # nVac x nD
    vPress = exp(VTmat %*% log1p(-vNow))              # U x nD

    ## fitness, per genotype per deme
    fIt.u    = exp(G0 %*% ((eqm.pre - t2r) * sEl)) * vPress * (1 - migration)   # U x nD
    pi.Omega = fIt.u[cbind(gen1, dEme)]
    fIt.adj  = (kD / nPer)[dEme]
    if(meanStandardize){ fIt.adj = fIt.adj / mean(pi.Omega) }

    ## reproduction -- co-infection widens the transmission bottleneck
    ## (coInf -> 0 gives Poisson exactly, so the single-infection model is nested)
    oFf = if(coInf < 1e-5) rpois(length(gen1), pi.Omega * fIt.adj)
          else rnbinom(length(gen1), mu = pi.Omega * fIt.adj, size = 1 / coInf)
    if(anyNA(oFf) || sum(oFf) < 1 || sum(oFf) >= popRunaway){ return(NULL) }
    kEep = rep(seq_along(gen1), oFf)
    gen1 = gen1[kEep]; dEme = dEme[kEep]              # offspring stay in the parent's deme

    ## between-deme movement (metapopulation mixing)
    if(mWithin > 0){
      sWap = runif(length(gen1)) < mWithin
      if(any(sWap)) dEme[sWap] = sample.int(nD, sum(sWap), replace = T, prob = dProp)
    }

    ## external immigration: destination deme ~ deme size; source pool SC-balanced within it
    if(migration > 0){
      nMig = rbinom(1, round(k), min(1, migration * k / length(gen1)))
      if(nMig > 0){
        nPerD = tabulate(sample.int(nD, nMig, replace = T, prob = dProp), nD)
        gen1  = c(gen1, unlist(lapply(seq_len(nD), function(dd)
                  if(nPerD[dd] > 0) sample.int(U, nPerD[dd], replace = T, prob = mP0[, dd])
                  else integer(0))))
        dEme  = c(dEme, rep(seq_len(nD), nPerD))
      }
    }

    ## sample, stratified by deme to match observed effort
    if(tIme[i] %in% eQm.date){
      cc = numeric(length(vtsc.lev) * nD)
      for(dd in seq_len(nD)){
        nTd = nObs[as.character(tIme[i]), dd]
        if(nTd > 0){
          iD = which(dEme == dd)
          if(length(iD) < 1){ return(NULL) }
          s  = iD[sample.int(length(iD), nTd, replace = T)]
          cc = cc + vtsc(gen1[s], dEme[s])
        }
      }
      rec.eQm[, j] = cc; j = j + 1
    }
  };rm(i)

  if(!keepGenotypes){ return(rec.eQm) }
  return(list(
    ss = rec.eQm,
    genotypes = data.frame(row = seq_len(U), tag = tag0, SC = sc0,
                           as.data.frame(VTmat), finalCount = tabulate(gen1, nbins = U)),
    demeCount = vapply(seq_len(nD), function(dd) tabulate(gen1[dEme == dd], nbins = U), numeric(U)),
    G = G0))
}

##### Jensen-Shannon divergence (translated from the cpp model) #####
jsd = function(p, q){
  p = p / sum(p); q = q / sum(q)
  m = (p + q) / 2
  s = 0
  i = m > 0 & p > 0; s = s + 0.5 * sum(p[i] * log(p[i] / m[i]))
  i = m > 0 & q > 0; s = s + 0.5 * sum(q[i] * log(q[i] / m[i]))
  return(s)
}

##### Model using JSD index #####
nfds_jsd = function(x, ss_obs){
  ## workers are fresh R processes: globalenv() does not travel with this function
  if(!exists("m.nfds", envir = globalenv(), inherits = FALSE)){
    d = Sys.getenv("NFDS_SRC", unset = getwd())
    owd = setwd(d); on.exit(setwd(owd), add = TRUE)
    sT = Sys.getenv("NFDS_STATE", unset = "../data/setupState.RData")
    if(file.exists(sT)){ load(sT, envir = globalenv()) }
    else{ sys.source("setup.r", envir = globalenv()) }
  }
  vSel = vapply(seq_len(nVac), function(v) x[[paste0("vSelected", v)]], 0)
  sim  = m.nfds(x[["propStrong"]], x[["fSelected"]], x[["wSelected"]], vSel, x[["migration"]],
                x[["coInf"]], x[["vSelMod"]], x[["mWithin"]])
  if(is.null(sim) || any(is.na(sim))){ return(ncol(ss_obs) * log(2)) }
  d = sum(vapply(seq_len(ncol(ss_obs)), function(j) jsd(sim[, j], ss_obs[, j]), 0))
  if(!is.finite(d)){ return(ncol(ss_obs) * log(2)) }
  return(d)
}


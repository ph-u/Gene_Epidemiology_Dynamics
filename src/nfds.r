#!/usr/bin/env Rscript
# author: ph-u
# script: nfds.r
# desc: NFDS ABCSMC v2 -- individual-based simulation
#
# ============================================================
# VERSION 2 FEATURES
#
# METAPOPULATION
#   Bacteria live across nD demes (cities / study sites).
#   HETEROGENEOUS MIGRATION: mWithin is an nD x nD matrix.
#   mWithin[d1, d2] = probability per generation that a HOST
#   in deme d1 travels to deme d2.  All bacteria carried by
#   that host travel with them (they stay in the same person).
#   Migration is resolved at the HOST level to keep co-infected
#   bacteria together.
#
# MULTIPLE VACCINES
#   Each vaccine has its own efficacy (vSel, a vector) and its
#   own roll-out schedule (vStart, aCoh, uPt -- all nVac x nD).
#
# HOST AGE STRUCTURE
#   Each bacterium lives in a person of a certain age (hAge).
#   Age is redrawn every generation because transmission moves
#   the bacterium into a DIFFERENT person.  Vaccine coverage is
#   computed from the person's age and when the vaccine started
#   in their deme -- not fitted.
#
# CO-INFECTION
#   coInf = 0: each person carries one strain (Poisson offspring).
#   coInf > 0: burst transmission -- a person can receive multiple
#   bacteria from the same parent (NB offspring), AND bacteria from
#   different parents can share the same person (cross-genotype
#   co-infection, controlled by the NB variance).
#
# INFECTION HISTORY (colonisation history)
#   Every person (host) gets a unique ID (hID).  When the infection
#   clears (the bacterium fails to produce offspring), we write the
#   episode into the history book (hHist): who cleared, which
#   bacterial type, which deme, at what age, start and end month,
#   and how many bacteria were sharing that person at the time.
# ============================================================
#
# in: source("nfds.r")
# out: NA
# arg: 0
# date: 20260919 (infection/colonisation history + heterogeneous mWithin matrix),
#       20260918 (age stratification, metapopulation, multi-vaccine), 20260820

##### NFDS individual-based simulation #####
m.nfds = function(propStrong, fSelected, wSelected, vSel, migration,
                  coInf, mWithin,       # mWithin: nD x nD matrix, zero diagonal
                  meanStandardize = FALSE, keepGenotypes = FALSE){

  ## --- Parameter transformation (all priors unif(0,1); rescale here) ---
  vSel      = vSel * .5        # per-vaccine efficacy: 0-0.5 per month
  migration = migration * .2   # external immigration: 0-0.2 per month

  ## === HETEROGENEOUS MIGRATION MATRIX ===
  ## mWithin[d1, d2] = probability a host in deme d1 moves to deme d2 per generation.
  ## Element-wise rescaling keeps all rates in 0-0.1 per month.
  mWithin = mWithin * .1
  stopifnot(is.matrix(mWithin),
            nrow(mWithin) == nD, ncol(mWithin) == nD,
            all(diag(mWithin) == 0),
            all(rowSums(mWithin) < 1 + 1e-9))   # total emigration < 1 per generation

  U = nrow(G0); L = length(gNam)

  ## --- Per-locus NFDS weights ---
  sEl = rep(log1p(wSelected), L)
  sEl[order(selMode$strength)[seq_len(floor(L * propStrong))]] = log1p(fSelected)

  ## Pre-compute log vaccine fitness cost (reused every generation)
  lVs = log1p(-vSel)           # length nVac, <= 0

  ## --- Initial population ---
  nInit = pmax(1L, round(kD * as.numeric(g("percentage initial infected")) / 100))
  dEme  = rep(seq_len(nD), nInit)
  gen1  = unlist(lapply(seq_len(nD), function(dd)
            sample(tag.pre[[dd]], nInit[dd], replace = TRUE)))

  ## === INFECTION HISTORY: initialise host tracking ===
  ## hID  = unique integer ID of the person currently carrying each bacterium.
  ## tInf = month when the current infection started.
  ## hCtr = counter; increments every time a new person enters the simulation.
  hCtr = length(gen1)
  hID  = seq_along(gen1)                          # initial: one unique host per bacterium
  tInf = rep(min(eQm$Month), length(gen1))        # all pre-vaccine isolates start at t_min

  ## History book: one row per cleared infection episode.
  ## Columns: hID | geno | deme | age | tInf | tClr | nCoInf
  maxHist = 500L   # cap per generation (memory bound); increase if needed
  if(keepGenotypes){
    hHist = data.frame(hID    = integer(0), geno   = integer(0),
                       deme   = integer(0), age    = numeric(0),
                       tInf   = integer(0), tClr   = integer(0),
                       nCoInf = integer(0))
  }

  rec.eQm = matrix(0L, nrow = length(vtsc.lev) * nD, ncol = length(eQm.date)); j = 1L

  ## --- Main generation loop ---
  tIme = (min(eQm$Month) + 1L):nGen
  for(i in seq_along(tIme)){
    tNow = tIme[i]

    ## Per-deme locus frequencies (one matrix-vector product for the whole metapopulation)
    cnt  = vapply(seq_len(nD),
                  function(dd) tabulate(gen1[dEme == dd], nbins = U),
                  numeric(U))
    nPer = colSums(cnt)
    if(any(nPer < 1L)){ return(NULL) }
    t2r  = t(crossprod(cnt, G0) / nPer)    # L x nD

    ## === HOST AGE STRUCTURE ===
    ## Each transmission gives the bacterium a NEW host, so age is redrawn.
    hAge = numeric(length(gen1))
    for(dd in seq_len(nD)){
      s = which(dEme == dd)
      if(length(s))
        hAge[s] = aLev[sample.int(length(aLev), length(s),
                                  replace = TRUE, prob = aProb[, dd])]
    }

    ## === MULTIPLE VACCINES ===
    ## A person is covered by vaccine v in their deme if: the vaccine had started
    ## AND they were young enough when it began.  Coverage is DATA, not fitted.
    lvp = numeric(length(gen1))
    for(v in seq_len(nVac)){
      sV  = vStart[v, dEme]
      cOv = ((tNow >= sV) & (hAge - (tNow - sV) < aCoh[v, dEme])) * uPt[v, dEme]
      lvp = lvp + cOv * VTmat[cbind(gen1, v)] * lVs[v]
    }

    ## Fitness = NFDS term x vaccine pressure x density regulation
    nfdsU    = exp(G0 %*% ((eqm.pre - t2r) * sEl))    # U x nD
    pi.Omega = nfdsU[cbind(gen1, dEme)] * exp(lvp) * (1 - migration)
    fIt.adj  = (kD / nPer)[dEme]
    if(meanStandardize){ fIt.adj = fIt.adj / mean(pi.Omega) }

    ## === CO-INFECTION: overdispersed offspring ===
    oFf = if(coInf < 1e-5) rpois(length(gen1), pi.Omega * fIt.adj)
          else rnbinom(length(gen1), mu = pi.Omega * fIt.adj, size = 1 / coInf)
    if(anyNA(oFf) || sum(oFf) < 1L || sum(oFf) >= popRunaway){ return(NULL) }

    ## === INFECTION HISTORY: record clearances ===
    ## A bacterium that produced zero offspring has cleared from its host.
    if(keepGenotypes){
      clr = which(oFf == 0L)
      if(length(clr)){
        uid        = unique(hID)
        nCoInf_all = tabulate(match(hID, uid))
        nCoInf_clr = nCoInf_all[match(hID[clr], uid)]
        smp = if(length(clr) <= maxHist) clr else sample(clr, maxHist)
        hHist = rbind(hHist, data.frame(
          hID    = hID[smp],    geno   = gen1[smp], deme = dEme[smp],
          age    = hAge[smp],   tInf   = tInf[smp], tClr = tNow,
          nCoInf = nCoInf_clr[match(smp, clr)]))
      }
    }

    ## === INFECTION HISTORY: assign new host IDs to offspring ===
    ## Each reproducing parent creates a new "person" for its offspring.
    ## All offspring from the same parent share that person (same-genotype burst).
    hasOff  = oFf > 0L
    nNew    = sum(hasOff)
    newHIDs = seq.int(hCtr + 1L, hCtr + nNew)
    hCtr    = hCtr + nNew
    parHID           = integer(length(gen1))
    parHID[hasOff]   = newHIDs

    ## === CO-INFECTION: cross-genotype sharing ===
    ## Bacteria from DIFFERENT parents sometimes land in the same person.
    ## The NB variance drives same-genotype bursts; this block merges
    ## distinct-genotype hosts, making the sharing cross-genotype.
    if(coInf > 0 && nNew >= 2L){
      pMrg = coInf / (1 + coInf)
      nMrg = rbinom(1L, nNew %/% 2L, pMrg)
      if(nMrg > 0L && nMrg * 2L <= nNew){
        pm = sample.int(nNew)
        for(m in seq_len(nMrg)){
          a = newHIDs[pm[m]]; b = newHIDs[pm[nMrg + m]]
          parHID[parHID == b] = a    # all of B's offspring move into A's host
        }
      }
    }

    ## Expand all vectors to offspring generation
    kEep = rep(seq_along(gen1), oFf)
    gen1 = gen1[kEep]; dEme = dEme[kEep]
    hID  = parHID[kEep]              # each offspring inherits its parent's new host ID
    tInf = rep(tNow, length(gen1))   # infection starts NOW for all offspring

    ##### Recombination (future) #####

    ## ==============================================================
    ## HETEROGENEOUS METAPOPULATION MIGRATION (nD x nD matrix)
    ## mWithin[d1, d2] = probability per generation that a HOST in
    ## deme d1 travels to deme d2.  Diagonal = 0; row sums < 1.
    ##
    ## Migration is resolved at the HOST level so that co-infected
    ## bacteria (sharing the same hID) always move together.
    ##
    ## Algorithm:
    ##   1. Freeze every host's current deme (uHdeme0).
    ##   2. For each source deme d1, decide independently for each
    ##      HOST (not bacterium) whether it travels, and if so where.
    ##   3. Record the new deme in newDeme[].
    ##   4. After processing all source demes, update dEme for every
    ##      bacterium in one vectorised step.
    ## ==============================================================
    if(any(mWithin > 0) && nD > 1L){
      uH       = unique(hID)
      uHdeme0  = dEme[match(uH, hID)]   # frozen deme per unique host
      newDeme  = uHdeme0                 # will be updated for traveling hosts

      for(dd in seq_len(nD)){
        hHere = uH[uHdeme0 == dd]       # hosts currently in source deme dd
        if(!length(hHere)) next
        rOut = mWithin[dd, seq_len(nD)[-dd]]  # emigration rates: mWithin[from, to]
        pLv  = sum(rOut); if(pLv <= 0) next
        ## Each HOST decides independently whether to travel
        trav = hHere[runif(length(hHere)) < pLv]
        if(!length(trav)) next
        dest = seq_len(nD)[-dd]
        newD = dest[sample.int(nD - 1L, length(trav), replace = TRUE, prob = rOut / pLv)]
        ## Record destination for each traveling host
        newDeme[match(trav, uH)] = newD
      }

      ## Apply all movements in one step (vectorised via match)
      hostIdx  = match(hID, uH)          # for each bacterium: index into uH
      changed  = newDeme[hostIdx] != uHdeme0[hostIdx]
      dEme[changed] = newDeme[hostIdx[changed]]
    }

    ## External immigration (bacteria from outside the whole system)
    ## Destination deme ~ deme size; source is SC-balanced within that deme.
    ## Immigrants receive fresh host IDs and tInf = now.
    if(migration > 0){
      nMig = rbinom(1L, round(k), min(1, migration * k / length(gen1)))
      if(nMig > 0L){
        nPerD = tabulate(sample.int(nD, nMig, replace = TRUE, prob = dProp), nD)
        nImm  = sum(nPerD)
        newG  = unlist(lapply(seq_len(nD), function(dd)
                  if(nPerD[dd] > 0L) sample.int(U, nPerD[dd],
                                                replace = TRUE, prob = mP0[, dd])
                  else integer(0)))
        gen1 = c(gen1, newG)
        dEme = c(dEme, rep(seq_len(nD), nPerD))
        hID  = c(hID,  seq.int(hCtr + 1L, hCtr + nImm))
        tInf = c(tInf, rep(tNow, nImm))
        hCtr = hCtr + nImm
      }
    }

    ## Sample the simulated population to match observed genomic effort
    if(tNow %in% eQm.date){
      cc = numeric(length(vtsc.lev) * nD)
      for(dd in seq_len(nD)){
        nTd = nObs[as.character(tNow), dd]
        if(nTd > 0L){
          iD = which(dEme == dd)
          if(!length(iD)){ return(NULL) }
          s  = iD[sample.int(length(iD), nTd, replace = TRUE)]
          cc = cc + vtsc(gen1[s], dEme[s])
        }
      }
      rec.eQm[, j] = cc; j = j + 1L
    }
  };rm(i)

  if(!keepGenotypes){ return(rec.eQm) }

  ## === INFECTION HISTORY: snapshot of current hosts (final generation) ===
  uid_fin    = unique(hID)
  nCoInf_fin = tabulate(match(hID, uid_fin))[match(hID, uid_fin)]
  hostState  = data.frame(
    hID          = hID,
    geno         = gen1,
    tag          = tag0[gen1],
    SC           = sc0[gen1],
    deme         = dEme,
    tInf         = tInf,
    duration_now = nGen - tInf,
    nCoInf       = nCoInf_fin,
    coInfected   = nCoInf_fin > 1L)

  list(ss        = rec.eQm,
       genotypes = data.frame(row        = seq_len(U), tag = tag0, SC = sc0,
                              as.data.frame(VTmat),
                              finalCount = tabulate(gen1, nbins = U)),
       demeCount = vapply(seq_len(nD),
                          function(dd) tabulate(gen1[dEme == dd], nbins = U),
                          numeric(U)),
       hostState  = hostState,
       hHist      = if(keepGenotypes) hHist else NULL,
       hHistStats = hist_stats(if(keepGenotypes) hHist else data.frame()),
       G          = G0)
}

##### Jensen-Shannon divergence (translated from functions.cpp) #####
jsd = function(p, q){
  p = p / sum(p); q = q / sum(q); m = (p + q) / 2; s = 0
  i = m > 0 & p > 0; s = s + 0.5 * sum(p[i] * log(p[i] / m[i]))
  i = m > 0 & q > 0; s = s + 0.5 * sum(q[i] * log(q[i] / m[i]))
  return(s)
}

##### Distance function for BRREWABC (called by each parallel worker) #####
nfds_jsd = function(x, ss_obs){
  ## Workers are fresh R processes -- globalenv() was not serialised.
  ## The bootstrap block restores all setup.r objects from the state file.
  if(!exists("m.nfds", envir = globalenv(), inherits = FALSE)){
    d  = Sys.getenv("NFDS_SRC", unset = getwd())
    ow = setwd(d); on.exit(setwd(ow), add = TRUE)
    sT = Sys.getenv("NFDS_STATE", unset = "../data/setupState.RData")
    if(file.exists(sT)) load(sT, envir = globalenv())
    else sys.source("setup.r", envir = globalenv())
  }
  ## Assemble per-vaccine vector from named parameters vSelected1, vSelected2, ...
  vSel = vapply(seq_len(nVac), function(v) x[[paste0("vSelected", v)]], 0)
  ## Assemble nD x nD migration matrix from named parameters mWithin_d1_d2
  mW = matrix(0, nD, nD)
  for(d1 in seq_len(nD)) for(d2 in seq_len(nD))
    if(d1 != d2) mW[d1, d2] = x[[paste0("mWithin_", d1, "_", d2)]]
  sim = m.nfds(x[["propStrong"]], x[["fSelected"]], x[["wSelected"]], vSel,
               x[["migration"]], x[["coInf"]], mW)
  if(is.null(sim) || any(is.na(sim))){ return(ncol(ss_obs) * log(2)) }
  d = sum(vapply(seq_len(ncol(ss_obs)), function(j) jsd(sim[,j], ss_obs[,j]), 0))
  if(!is.finite(d)){ return(ncol(ss_obs) * log(2)) }
  return(d)
}


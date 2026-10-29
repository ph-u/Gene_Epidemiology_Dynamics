#!/usr/bin/env Rscript
# author: ph-u, nickjcroucher
# script: setup.r
# desc: NFDS ABCSMC v2 -- data preparation
# Features: metapopulation | multi-vaccine | host age structure | co-infection |
#           heterogeneous migration matrix | infection history
# in: Rscript setup.r [../raw/input.csv] [seed entry]  (usually sourced by model.r)
# out: NA
# arg: 1
# date: 20260919

##### env #####
source("src.r"); source("nfds.r")
if(!exists("argv")){
  argv = commandArgs(T)
  if(length(argv) != 2){ argv = c("../raw/input.csv", 1) }
}
f.in = read.csv(argv[1], header = T)

k          = as.numeric(g("population"))
popRunaway = 10 * k

##### Data #####
## Column order: Time | Deme | Age | VT1..VTn | SC | <gene columns>
## Time  = month (integer; 0 = earliest vaccine introduction anywhere)
## Deme  = site label, e.g. "UK", "NL"
## Age   = host age in MONTHS at time of sampling
## SC    = sequence cluster -- must be the LAST metadata column
d0   = read.table(g("data"), header = T, sep = "\t")
mEnd = which(colnames(d0) == "SC")
stopifnot(all(c("Time","Deme","Age","SC") %in% colnames(d0)[1:mEnd]))

## === FEATURE: METAPOPULATION ===
## Encode deme labels as integers 1..nD so they index into matrices.
dLev    = sort(unique(d0$Deme))
d0$Deme = match(d0$Deme, dLev)
nD      = length(dLev)

## === FEATURE: MULTIPLE VACCINES ===
vtCols = grep("^VT[0-9]+$", colnames(d0), value = TRUE)
nVac   = length(vtCols)
stopifnot(nVac > 0, nD > 0)

##### Intermediate-frequency filter (pre-vaccine window only) #####
## Over 30 years a gene that swept during the study would be removed by a
## whole-series filter -- that is the wrong reason.  Restrict to Time <= 0.
pre = d0$Time <= 0
stopifnot(sum(pre) > 0)
fPre    = colMeans(d0[pre, -(1:mEnd), drop = FALSE])
d0.keep = fPre > .05 & fPre < .95
d0      = d0[, c(colnames(d0)[1:mEnd], names(d0.keep)[which(d0.keep)])]
gNam    = colnames(d0)[-(1:mEnd)]

##### Barcoding unique genotypes #####
d0$data = apply(d0[, gNam, drop = FALSE], 1, paste, collapse = "")
uNq     = unique(d0[, c(gNam, "data")])
d0.u    = cbind(tag = bcod(df = uNq), uNq)
d0$tag  = d0.u$tag[match(d0$data, d0.u$data)]
d0$data = d0.u$data = NULL

##### Genotype attributes #####
## Deme and Age are HOST properties -- the same genotype appears in many demes
## and at many ages.  They are NOT stored in d0.u.
d0.u[vtCols] = d0[match(d0.u$tag, d0$tag), vtCols]
d0.u$SC      = d0$SC[match(d0.u$tag, d0$tag)]
stopifnot(!anyDuplicated(d0.u$tag),
          all(tapply(d0$SC, d0$tag, function(z) length(unique(z))) == 1))
for(v in vtCols) stopifnot(all(tapply(d0[[v]], d0$tag, function(z) length(unique(z))) == 1))

VTmat = as.matrix(d0.u[, vtCols, drop = FALSE]); storage.mode(VTmat) = "double"
G0    = as.matrix(d0.u[, gNam]);                 storage.mode(G0)    = "double"
tag0  = d0.u$tag
sc0   = d0.u$SC

## Summary-statistic classes: VT profile x SC (deme added at sampling time)
d0.u$paste = do.call(paste, c(d0.u[vtCols], list(d0.u$SC), sep = ";"))
vtsc.lev   = sort(unique(d0.u$paste))
vtsc.idx   = match(d0.u$paste, vtsc.lev)

##### Equilibrium gene frequencies #####
eQm       = data.frame(Month = sort(unique(d0$Time)))
eQm[gNam] = t(sapply(split(d0[, gNam, drop = FALSE], d0$Time), colMeans))
nGen      = max(eQm$Month, 1)
eQm.date  = eQm$Month[eQm$Month > 0]
stopifnot(all(eQm.date %in% seq_len(nGen)), min(eQm$Month) < nGen)

## Per-DEME pre-vaccination equilibrium (L x nD, isolate-weighted).
## NFDS is local: each deme's bacteria are selected against that deme's equilibrium.
eqm.pre = vapply(seq_len(nD), function(dd){
  s = pre & d0$Deme == dd
  if(!any(s)) stop("deme ", dLev[dd], " has no pre-vaccination isolates")
  as.numeric(colMeans(d0[s, gNam, drop = FALSE]))
}, numeric(length(gNam)))
stopifnot(nrow(eqm.pre) == length(gNam), ncol(eqm.pre) == nD)

## Global pre-vaccination equilibrium -- used ONLY to rank loci into strong/weak classes.
eqm.glb = as.numeric(colMeans(d0[pre, gNam, drop = FALSE]))

##### Strong / weak NFDS locus ranking #####
selMode = (colMeans(eQm[eQm$Month > 0, gNam, drop = FALSE]) - eqm.glb)^2 /
          (1 - eqm.glb * (1 - eqm.glb))
selMode = data.frame(gene = gNam, strength = as.numeric(selMode), category = "weak")

##### === FEATURE: MULTIPLE VACCINES === #####
## Roll-out table: Deme | Vaccine | StartMonth | AgeWindow | Uptake
rOut      = read.table(g("rollout"), header = T, sep = "\t")
rOut$Deme = match(rOut$Deme, dLev)
stopifnot(all(rOut$Vaccine %in% vtCols), !anyNA(rOut$Deme),
          all(c("StartMonth","AgeWindow") %in% colnames(rOut)), all(rOut$AgeWindow > 0))
vIx    = cbind(match(rOut$Vaccine, vtCols), rOut$Deme)
vStart = matrix(Inf, nVac, nD, dimnames = list(vtCols, dLev))
aCoh   = matrix(0,   nVac, nD, dimnames = list(vtCols, dLev))
uPt    = matrix(0,   nVac, nD, dimnames = list(vtCols, dLev))
vStart[vIx] = rOut$StartMonth
aCoh[vIx]   = rOut$AgeWindow
uPt[vIx]    = if("Uptake" %in% colnames(rOut)) rOut$Uptake else 1
stopifnot(all(uPt >= 0), all(uPt <= 1))
vMonth = seq(min(eQm$Month), max(eQm$Month))    # month sequence (bounds reference only)

##### === FEATURE: HOST AGE STRUCTURE === #####
## Probability that a bacterium is in a person of a given age, per deme.
## Derived from the sampled isolates; if surveillance is age-biased, supply
## an independent distribution here.
aLev  = sort(unique(d0$Age))
aProb = vapply(seq_len(nD), function(dd)
          as.numeric(table(factor(d0$Age[d0$Deme == dd], levels = aLev))),
          numeric(length(aLev)))
aProb = t(t(aProb) / colSums(aProb))            # nA x nD, columns sum to 1
stopifnot(all(abs(colSums(aProb) - 1) < 1e-8), all(is.finite(aProb)))

##### === FEATURE: METAPOPULATION === #####
## Partition carrying capacity in proportion to deme size.
dProp = as.numeric(table(factor(d0$Deme, levels = seq_len(nD)))) / nrow(d0)
kD    = k * dProp
stopifnot(all(kD > 0))

## External immigrant pool: SC-balanced WITHIN each deme (U x nD, columns sum to 1).
mP0 = vapply(seq_len(nD), function(dd){
  s     = d0$Deme == dd
  nGeno = as.numeric(table(factor(d0$tag[s], levels = d0.u$tag)))
  nSC   = table(factor(d0$SC[s], levels = sort(unique(d0$SC))))
  nSCu  = as.numeric(nSC[as.character(d0.u$SC)])
  Sd    = sum(nSC > 0)
  ifelse(nGeno > 0, nGeno / (Sd * nSCu), 0)
}, numeric(nrow(d0.u)))
stopifnot(all(abs(colSums(mP0) - 1) < 1e-8))
migIdx = seq_len(nrow(d0.u))

##### Starting pools and sampling effort #####
tag.pre = lapply(seq_len(nD), function(dd) match(d0$tag[pre & d0$Deme == dd], d0.u$tag))
stopifnot(all(lengths(tag.pre) > 0))
nObs = table(factor(d0$Time, levels = eQm$Month),
             factor(d0$Deme, levels = seq_len(nD)))

##### Observed summary statistic #####
mIg0 = vapply(eQm.date, function(tt){
  s = which(d0$Time == tt)
  vtsc(match(d0$tag[s], d0.u$tag), d0$Deme[s])
}, numeric(length(vtsc.lev) * nD))

message(sprintf(paste("classes = %d (%d VT|SC x %d demes) | zero cells = %.1f%%",
                      "| median isolates/timepoint = %.0f | ages = %d | params = %d"),
                nrow(mIg0), length(vtsc.lev), nD, 100 * mean(mIg0 == 0),
                median(colSums(mIg0)), length(aLev),
                5L + nD * (nD - 1L) + nVac))


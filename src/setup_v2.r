#!/usr/bin/env Rscript
# author: ph-u, nickjcroucher
# script: setup.r
# desc: NFDS ABCSMC -- metapopulation, multi-vaccine, graded roll-out, co-infection
# in: Rscript setup.r [../raw/input.csv] [seed entry]
# out: NA
# arg: 1
# date: 20260918

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
## column order: Time | Deme | VT1..VTn | SC | <gene columns>   (SC = LAST metadata column)
d0   = read.table(g("data"), header = T, sep = "\t")
mEnd = which(colnames(d0) == "SC")
stopifnot(all(c("Time", "Deme", "SC") %in% colnames(d0)[1:mEnd]))

dLev    = sort(unique(d0$Deme))
d0$Deme = match(d0$Deme, dLev)                       # demes as 1..nD
nD      = length(dLev)
vtCols  = grep("^VT[0-9]+$", colnames(d0), value = TRUE)
nVac    = length(vtCols)
stopifnot(nVac > 0, nD > 0)

## intermediate-frequency filter on the PRE-vaccine window only:
## over 30 years a gene that swept during the study would otherwise be dropped
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
d0.u[vtCols] = d0[match(d0.u$tag, d0$tag), vtCols]
d0.u$SC      = d0$SC[match(d0.u$tag, d0$tag)]
## migProb and the class definition both assume genotype -> (SC, VT profile) is one-to-one
stopifnot(!anyDuplicated(d0.u$tag),
          all(tapply(d0$SC, d0$tag, function(z) length(unique(z))) == 1))

VTmat = as.matrix(d0.u[, vtCols, drop = FALSE]); storage.mode(VTmat) = "double"
G0    = as.matrix(d0.u[, gNam]);                 storage.mode(G0)    = "double"
tag0  = d0.u$tag
sc0   = d0.u$SC

## summary-statistic classes: VT profile x SC   (deme is added at sampling time)
d0.u$paste = do.call(paste, c(d0.u[vtCols], list(d0.u$SC), sep = ";"))
vtsc.lev   = sort(unique(d0.u$paste))
vtsc.idx   = match(d0.u$paste, vtsc.lev)

##### Equilibrium #####
eQm       = data.frame(Month = sort(unique(d0$Time)))
eQm[gNam] = t(sapply(split(d0[, gNam, drop = FALSE], d0$Time), colMeans))
nGen      = max(eQm$Month, 1)
eQm.date  = eQm$Month[eQm$Month > 0]
stopifnot(all(eQm.date %in% seq_len(nGen)))

## PER-DEME pre-vaccination equilibrium frequencies (L x nD), isolate-weighted.
## NFDS is local: each deme's genes are selected against that deme's own equilibrium.
eqm.pre = vapply(seq_len(nD), function(dd){
  s = pre & d0$Deme == dd
  if(!any(s)) stop("deme ", dLev[dd], " has no pre-vaccination isolates")
  as.numeric(colMeans(d0[s, gNam, drop = FALSE]))
}, numeric(length(gNam)))
stopifnot(nrow(eqm.pre) == length(gNam), ncol(eqm.pre) == nD)

## GLOBAL pre-vaccination frequencies -- used only to rank loci into strong/weak
## classes, because that is a property of the gene, not of a deme.
eqm.glb = as.numeric(colMeans(d0[pre, gNam, drop = FALSE]))

##### Strong / weak NFDS classes (global ranking) #####
selMode = (colMeans(eQm[eQm$Month > 0, gNam, drop = FALSE]) - eqm.glb)^2 /
          (1 - eqm.glb * (1 - eqm.glb))
selMode = data.frame(gene = gNam, strength = as.numeric(selMode), category = "weak")

##### Vaccine roll-out: deme- and vaccine-specific #####
rOut      = read.table(g("rollout"), header = T, sep = "\t")   # Deme | Vaccine | StartMonth
rOut$Deme = match(rOut$Deme, dLev)
stopifnot(all(rOut$Vaccine %in% vtCols), !anyNA(rOut$Deme))
vStart = matrix(Inf, nVac, nD, dimnames = list(vtCols, dLev))  # Inf = never introduced there
vStart[cbind(match(rOut$Vaccine, vtCols), rOut$Deme)] = rOut$StartMonth
vMonth = seq(min(eQm$Month), max(eQm$Month))

##### Deme sizes: kappa allocation and external-immigrant weighting #####
dProp = as.numeric(table(factor(d0$Deme, levels = seq_len(nD)))) / nrow(d0)
kD    = k * dProp                                              # per-deme carrying capacity
stopifnot(all(kD > 0))

##### External immigrant pool: SC-balanced WITHIN each deme (U x nD, columns sum to 1) #####
mP0 = vapply(seq_len(nD), function(dd){
  s     = d0$Deme == dd
  nGeno = as.numeric(table(factor(d0$tag[s], levels = d0.u$tag)))   # n_{d,u}
  nSC   = table(factor(d0$SC[s], levels = sort(unique(d0$SC))))     # n_{d,s}
  nSCu  = as.numeric(nSC[as.character(d0.u$SC)])                    # n_{d,s(u)}
  Sd    = sum(nSC > 0)                                              # SCs present in deme d
  ifelse(nGeno > 0, nGeno / (Sd * nSCu), 0)
}, numeric(nrow(d0.u)))
stopifnot(all(abs(colSums(mP0) - 1) < 1e-8))
migIdx = seq_len(nrow(d0.u))

##### Starting pools and sampling effort #####
tag.pre = lapply(seq_len(nD), function(dd) match(d0$tag[pre & d0$Deme == dd], d0.u$tag))
stopifnot(all(lengths(tag.pre) > 0))
nObs = table(factor(d0$Time, levels = eQm$Month), factor(d0$Deme, levels = seq_len(nD)))

##### Observed summary statistic: (VT profile x SC x deme) rows, timepoints as columns #####
mIg0 = vapply(eQm.date, function(tt){
  s = which(d0$Time == tt)
  vtsc(match(d0$tag[s], d0.u$tag), d0$Deme[s])
}, numeric(length(vtsc.lev) * nD))

message(sprintf("classes = %d (%d VT|SC x %d demes) | zero cells = %.1f%% | median isolates/timepoint = %.0f",
                nrow(mIg0), length(vtsc.lev), nD, 100 * mean(mIg0 == 0), median(colSums(mIg0))))


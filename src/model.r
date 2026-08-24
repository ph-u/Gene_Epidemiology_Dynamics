#!/bin/env Rscript
# author: ph-u, nickjcroucher
# script: model.r
# desc: NFDS ABCSMC model
# in: Rscript model.r [../raw/input.csv] [seed entry]
# out: NA
# arg: 1
# date: 20260819

argv = (commandArgs(T))
if(length(argv)!=2){argv = c("../raw/input.csv", 1)}

##### env #####
nCPU = as.integer(Sys.getenv("LSB_DJOB_NUMPROC", unset = "1"))
sEed = read.csv("../raw/seed.csv", header = F)[[1]][as.numeric(argv[2])]
set.seed(sEed)
library(BRREWABC) # https://github.com/GaelBn/BRREWABC
source("src.r");source("nfds.r")
f.in = read.csv(argv[1], header = T)

popRunaway = 10 * as.numeric(g("population"))
d0 = read.table(g("data"), header = T, sep = "\t")
mEnd = which(colnames(d0)=="SC")
d0.keep = colMeans(d0[,-(1:mEnd)]) > 0 & colMeans(d0[,-(1:mEnd)]) < 1
d0 = d0[,c(colnames(d0)[1:mEnd], names(d0.keep)[which(d0.keep)])]

##### Barcoding unique genotypes & construct difference dictionary #####
gNam = colnames(d0)[-(1:mEnd)]
eqm.pre = as.numeric(colMeans(d0[d0$Time <= 0, gNam]))
d0$data = apply(d0[-(1:mEnd)],1,paste, collapse = "")
d0.u = cbind(bcod(df = unique(d0[,-(1:mEnd)])), unique(d0[,-(1:mEnd)]))
colnames(d0.u)[1] = "tag"
d0$tag = d0.u[match(d0$data, d0.u$data),1]
d0$data = d0.u$data = NULL

##### Relevant data #####
d0.u$VT = d0$VT[match(d0.u$tag, d0$tag)]
d0.u$SC = d0$SC[match(d0.u$tag, d0$tag)]
d0.u$paste = paste(as.character(d0.u$VT),as.character(d0.u$SC), sep = ";")
d0.u$migProb = as.numeric(table(factor(d0$tag, levels = d0.u$tag))) / (length(unique(d0$SC)) * as.numeric(table(d0$SC)[as.character(d0.u$SC)]))

##### Equilibrium calculations #####
eQm = as.data.frame(matrix(nrow = length(unique(d0$Time)), ncol = ncol(d0)-mEnd))
colnames(eQm) = c("Month", colnames(d0)[-c(1:mEnd,ncol(d0))])
eQm$Month = sort(unique(d0$Time))
eQm[,-1] = t(sapply(split(d0[,-c(1:mEnd,ncol(d0))], d0$Time), colMeans))
k = as.numeric(g("population")); nGen = max(eQm[nrow(eQm), 1],1)

##### Vaccine duration #####
stopifnot(length(as.numeric(g("vaccine start month"))) == 1)
vCurve = seq(min(eQm[,1]), max(eQm[,1]))
vCurve = ifelse(vCurve < as.numeric(g("vaccine start month")), 0, 1)

##### Group genes that face strong / weak selection (constant dataframe) #####
selMode = (colMeans(eQm[eQm[,1] > 0,-1]) - eqm.pre)^2/(1 - eqm.pre * (1 - eqm.pre))
selMode = data.frame(gene = names(selMode), strength = as.numeric(selMode), category = "weak")

##### Runtime acceleration #####
nObs = table(d0$Time) # number of observations in each time points
G0 = as.matrix(d0.u[,gNam]); storage.mode(G0) = "double" # Matrix for speed

tag0 = d0.u$tag
vt0 = d0.u$VT
sc0 = d0.u$SC
mP0 = d0.u$migProb
vtsc.lev = sort(unique(d0.u$paste))
vtsc.idx = match(d0.u$paste, vtsc.lev)

migIdx = seq_len(nrow(d0.u))
stopifnot(sum(d0$Time <= 0) > 0, length(eqm.pre) == length(gNam))
tag.pre = match(d0$tag[d0$Time <= 0], d0.u$tag)
eQm.date = eQm[which(eQm[,1] > 0), 1]
stopifnot(all(eQm.date %in% seq_len(nGen)))

##### Migration probability (VT & SC balanced, 2017 paper) #####
mIg0 = c(); for(i in 1:length(eQm.date)){ mIg0 = c(mIg0, vtsc(match(d0$tag[which(d0$Time == eQm.date[i])], d0.u$tag))) };rm(i)
mIg0 = matrix(mIg0, ncol = length(eQm.date))

##### Parameters #####
#tRuth = c(sigma_f = 0.15, sigma_v = 0.10, m = 0.02, pf = 0.25, sigma_w = 0.003) # 2017 publication

prior_dist <- list(nfds = list(c("propStrong", "unif", 0, 1),
                               c("fSelected", "unif", 1e-6, .22),
                               c("wSelected", "unif", 1e-6, .15),
                               c("vSelected", "unif", 0, .5),
                               c("migration", "unif", 0, .2)))
##### ABCSMC-NFDS (first draft by Claude.ai) #####
res <- abcsmc(
  model_list             = list(nfds = nfds_jsd),
  prior_dist             = prior_dist,
  ss_obs                 = mIg0,
  nb_threshold           = 1,
  nb_acc_prtcl_per_gen   = 3,
  max_number_of_gen      = 10,
  new_threshold_quantile = 0.8,
  distance_threshold_min = .14, # Claude Opus 5 suggested
  acceptance_rate_min    = 0.01,
  experiment_folderpath  = "../data",
  max_concurrent_jobs    = nCPU,
  use_lhs_for_first_iter = TRUE,
  verbose                = TRUE,
  progressbar            = FALSE
)

##### Export #####
pOst = res$particles[res$particles$gen == max(res$particles$gen), ]
oUt = vector("list", nrow(pOst))
for(p in seq_len(nrow(pOst))){
  sP = sEed + p; set.seed(sP) # recorded, so the run is reproducible
  r = m.nfds(pOst$propStrong[p], pOst$fSelected[p], pOst$wSelected[p],
             pOst$vSelected[p], pOst$migration[p], keepGenotypes = T)
  gT = r$genotypes; gT$particle = p; gT$seed = sP
  oUt[[p]] = gT
  saveRDS(list(G = r$G, meta = gT), paste0("../data/nfdsG_", sP, ".rds"))
}
write.csv(do.call(rbind, oUt), paste0("../data/nfdsGenotypes_", gsub(" ", "-", date()), "_", sEed, ".csv"), row.names = F, quote = F)


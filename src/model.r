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
sEed = read.csv("../raw/seed.csv", header = F)[[1]][as.numeric(argv[2])]
set,seed(sEed)
#library(BRREWABC) # https://github.com/GaelBn/BRREWABC
source("src.r");source("nfds.r")
f.in = read.csv(argv[1], header = T)
d0 = read.table(g("data"), header = T, sep = "\t")
mEnd = which(colnames(d0)=="SC")
d0.keep = colMeans(d0[,-(1:mEnd)]) > 0 & colMeans(d0[,-(1:mEnd)]) < 1
d0 = d0[,c(colnames(d0)[1:mEnd], names(d0.keep)[which(d0.keep)])]

##### Genotype hex reformat #####
# i = (which(colnames(d0)=="SC") + 1):ncol(d0)
# if(length(i)%%31 > 0){i0 = rep(NA,31-length(i)%%31)}else{i0 = c()}
# binGroup = matrix(c(i, i0), nrow = 31) # group gene presence/absence for binning
# d0.tag = cbind(d0[,1:which(colnames(d0)=="SC")], matrix(nrow = nrow(d0), ncol = ncol(binGroup)))
# for(i in 1:ncol(binGroup)){
#   d0.tag[,i+min(binGroup)-1] = tbh(apply(d0[,binGroup[1,i]:binGroup[nrow(binGroup),i]], 1, paste, collapse = ""), sTart = 2)
# };rm(i,i0)

##### Barcoding unique genotypes & construct difference dictionary #####
gNam = colnames(d0)[-(1:mEnd)]
d0$data = apply(d0[-(1:mEnd)],1,paste, collapse = "")
d0.u = cbind(bcod(df = unique(d0[,-(1:mEnd)])), unique(d0[,-(1:mEnd)]))
colnames(d0.u)[1] = "tag"
d0$tag = d0.u[match(d0$data, d0.u$data),1]
d0$data = d0.u$data = NULL

nObs = table(d0$Time) # number of observations in each time points

#g0 = t(t(d0.u[,-c(1,ncol(d0.u))]) - as.numeric(d0.u[1,-c(1,ncol(d0.u))])) # 1/0 difference contrast with reference (1st row)
#d0.u$diff = apply(g0, 1, function(x){paste(which(x!=0), collapse = ";")})
#rm(g0)

d0.u$VT = d0$VT[match(d0.u$tag, d0$tag)]
d0.u$SC = d0$SC[match(d0.u$tag, d0$tag)]
#d0.u$d.pi = d0.u$d.omega = NA

##### Equilibrium calculations #####
eQm = as.data.frame(matrix(nrow = length(unique(d0$Time)), ncol = ncol(d0)-mEnd))
colnames(eQm) = c("Month", colnames(d0)[-c(1:mEnd,ncol(d0))])
eQm$Month = unique(d0$Time)[order(unique(d0$Time))]
eQm[,-1] = t(sapply(split(d0[,-c(1:mEnd,ncol(d0))], d0$Time), colMeans))
eqm.pre = as.numeric(colMeans(eQm[which(eQm[,1] <= 0),-1]))
eQm.date = eQm[which(eQm[,1] > 0), 1]

##### Vaccine duration #####
vCurve = as.numeric(seq_len(eQm[nrow(eQm), 1]) >= as.numeric(g("vaccine start month")))

##### Migration probability (VT & SC balanced, 2017 paper) #####
#mIg0 = table(d0[, c("VT", "SC")])
#mIg0 = as.data.frame(mIg0 / sum(mIg0))
#mIg0$paste = paste(as.character(mIg0[,1]),as.character(mIg0[,2]), sep = ";")
d0.u$paste = paste(as.character(d0.u$VT),as.character(d0.u$SC), sep = ";")
d0.u$migProb = as.numeric(table(factor(d0$tag, levels = d0.u$tag))) / (length(unique(d0$SC)) * as.numeric(table(d0$SC)[as.character(d0.u$SC)])) #mIg0$Freq[match(d0.u$paste, mIg0$paste)]
migIdx = seq_len(nrow(d0.u))

##### Data VT*SC distribution #####
mIg0 = c(); for(i in 1:length(eQm.date)){ mIg0 = c(mIg0, vtsc(d0$tag[which(d0$Time == eQm.date[i])])) };rm(i)
mIg0 = matrix(mIg0, ncol = length(eQm.date))

##### Group genes that face strong / weak selection (constant dataframe) #####
selMode = (colMeans(eQm[eQm[,1] > 0,-1]) - eqm.pre)^2/(1 - eqm.pre * (1 - eqm.pre))
selMode = data.frame(gene = names(selMode), strength = as.numeric(selMode), category = "weak")

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
  distance_threshold_min = ncol(mIg0) * log(2),
  acceptance_rate_min    = 0.01,
  experiment_folderpath  = "../data",
  max_concurrent_jobs    = 1,
  use_lhs_for_first_iter = TRUE,
  verbose                = FALSE,
  progressbar            = TRUE
)
particles  = res$particles
thresholds = res$thresholds
write.csv(d0.u[,c("tag", "diff", "VT", "SC", "migProb")], paste0("../data/nfds_", gsub(" ", "-", date()), "_", sEed, ".csv"), row.names = F, quote = F)


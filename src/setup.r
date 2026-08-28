#!/bin/env Rscript
# author: ph-u, nickjcroucher
# script: setup.r
# desc: NFDS ABCSMC model
# in: Rscript setup.r [../raw/input.csv] [seed entry]
# out: NA
# arg: 1
# date: 20260819

##### env #####
source("src.r");source("nfds.r")
if(!exists("argv")){
  argv = commandArgs(T)
  if(length(argv) != 2){ argv = c("../raw/input.csv", 1) }
}
f.in = read.csv(argv[1], header = T)

k = as.numeric(g("population"))
popRunaway = 10 * k
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
nGen = max(eQm[nrow(eQm), 1],1)

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


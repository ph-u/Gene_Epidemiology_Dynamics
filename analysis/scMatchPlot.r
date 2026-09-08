#!/bin/env Rscript
# author: ph-u
# script: scMatchPlot.r
# desc: Visualizing sequence cluster distribution from NFDS model
# in: Rscript scMatchPlot.r
# out: res/scMatchPlot--*
# arg: 0
# date: 20260903

##### colour #####
cBp = c(); set.seed(123); for(i in c("Okabe-Ito", "alphabet", "polychrome 36", "dark 2", "set 1", "classic tableau")){
  if(i == "Okabe-Ito"){
    cBp = c(cBp, rev(palette.colors(palette = i, alpha=1, recycle = F)))
  }else{
    cBp.t = palette.colors(palette = i, alpha=1, recycle = F)
    cBp = c(cBp, sample(cBp.t, length(cBp.t)))
  }};rm(i, cBp.t);cBp = unique(cBp)

##### env #####
sEed = read.csv("../raw/seed.csv", header = F)
#rAw = read.table("../data/mass.input", header = T, sep = "\t")
parIn = read.csv("../data/all_accepted_particles.csv", header = T)
pOri = getwd(); setwd("../src/"); source("setup.r"); setwd(pOri); rm(pOri)

parIn$gen = sEed$V1[match(parIn$gen, row.names(sEed))]

##### Simulation reconstruction #####
cat(date(),": start model\n")
for(i in seq_len(nrow(parIn))){
  cat(date(),":",i,"/",nrow(parIn),"(",round(i/nrow(parIn)*100,2),"% )       \r")
  set.seed(parIn$gen[i])
  r.nfds = cbind(as.numeric(i), vtsc.lev, read.table(text = vtsc.lev, sep = ";"), m.nfds(parIn$propStrong[i], parIn$fSelected[i], parIn$wSelected[i], parIn$vSelected[i], parIn$migration[i], meanStandardize = F, keepGenotypes = T)$ss)
  colnames(r.nfds) = c("rep", "vtsc", "vt", "sc", paste0("t", eQm.date))
  r.nfds[,(-1:0)+ncol(r.nfds)] = r.nfds[,(-1:0)+ncol(r.nfds)]/colSums(r.nfds[,(-1:0)+ncol(r.nfds)])
  if(i>1){
    rEs = rbind(rEs, r.nfds)
  }else{
    rEs = r.nfds
  }
};rm(i, r.nfds);cat("\n")

r.Data = as.data.frame(table(d0[which(d0$Time>0),c("Time", "VT", "SC")]))
r.Data$vtsc = paste0(r.Data$VT, ";", r.Data$SC)
r.Data = merge(r.Data[which(r.Data$Time==36),-which(colnames(r.Data)=="Time")], r.Data[which(r.Data$Time!=36),-which(colnames(r.Data)=="Time")], by = c("VT", "SC", "vtsc"), all = T)
r.Data[,(-1:0)+ncol(r.Data)] = r.Data[,(-1:0)+ncol(r.Data)]/colSums(r.Data[,(-1:0)+ncol(r.Data)])

##### Compare data & simulation #####
rEs0 = aggregate(cbind(t36, t72) ~ vt + sc + vtsc, data = rEs, quantile, probs = c(.05,.5,.95), simplify = T)
rEs0$d36 = r.Data$Freq.x[match(rEs0$vtsc, r.Data$vtsc)]
rEs0$d72 = r.Data$Freq.y[match(rEs0$vtsc, r.Data$vtsc)]
rEs0[,(-1:0)+ncol(rEs0)][is.na(rEs0[,(-1:0)+ncol(rEs0)])] = 0

##### Plot #####
pLt = cbind(rEs0[,1:3],rep(c(36,72), each = nrow(rEs0)), c(rEs0$t36[,1], rEs0$t72[,1]), c(rEs0$t36[,2], rEs0$t72[,2]), c(rEs0$t36[,3], rEs0$t72[,3]), c(rEs0$d36, rEs0$d72))
colnames(pLt)[-c(1:3)] = c("month", "p.05", "p.50", "p.95", "real")

jpeg("../res/scMatchPlot--compare.jpeg", width = 6000, height = 2100, res = 300)
par(mfrow = c(2,1), mar = c(5,4,3,0)+.1)
i0 = unique(pLt$month); for(i in seq_len(length(i0))){
  pLt.0 = pLt[pLt$month==i0[i],]; pLt.0 = pLt.0[order(as.numeric(pLt.0$sc)),]
  plot(x = seq_len(nrow(pLt.0)), y = pLt.0$real, col = "#000000ff", pch = as.numeric(as.factor(pLt.0$vt))+2, xaxt = "n", xlab = "Sequence Cluster", ylab = "Proportion", ylim = c(0, max(unlist(pLt[,-c(1:4)]))), main = "+ = NVT, x = VT; Black = data, Pink = simulation", cex = 2)
  segments(x0 = seq_len(nrow(pLt.0)), y0 = pLt.0$p.05, y1 = pLt.0$p.95, col = cBp[2])
  points(x = seq_len(nrow(pLt.0)), y = pLt.0$p.50, col = cBp[2], pch = 20)
  axis(1, at = seq_len(nrow(pLt.0)), labels = pLt.0$sc)
  text(x = 10, y = max(unlist(pLt[,-c(1:4)]))*.9, labels = paste0("Month = ",i0[i]))
};rm(i,i0, pLt.0)
invisible(dev.off())

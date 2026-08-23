#!/bin/env Rscript
# author: ph-u
# script: src.r
# desc: technical functions of model.r
# in: source("src.r")
# out: NA
# arg: 0
# date: 20260819

##### get values from input.csv #####
g = function(x, f = f.in){return(f$Value[f$Type==x])}

##### transform between binary and hexadecimals #####
# tbh = function(x = 1, sTart = 2){
#   if(!(sTart %in% c(2,16))){stop("Only bin and hex textstrings allowed.")}
#   if(sTart > 2){
#     return(paste(rev(as.integer(intToBits(strtoi(as.character(x), base = 16)))), collapse = ""))
#   }else{
#     return(as.hexmode(strtoi(as.character(x), base = 2)))
#   }
# }

##### extract row and column from the gene matrix #####
# rcG = function(x, gM = binGroup, colAdj = which(colnames(d0)=="SC")){
#   if(x > max(gM) || x < min(gM)){stop(paste("Invalid gene notation:",min(gM),"-",max(gM),"."))}
#   r0 = (which(gM==x) %% nrow(gM))
#   return(c(ifelse(r0>0,r0,nrow(gM)), ceiling(which(gM==x) / nrow(gM)) + colAdj))
# }

##### barcoding #####
bcod = function(df = unique(d0[,-(1:2)])){
  nC = ceiling(log(nrow(df)) / log(length(LETTERS)))
  b = data.frame(a = LETTERS, b = rep(LETTERS, each = length(LETTERS)))
  if(nC > 2){
    for(i in 1:(nC-2)){
      b = cbind(b, rep(LETTERS, each = nrow(b)))
    }
  }
  return(apply(b,1,paste, collapse = "")[1:nrow(df)])
}

##### reconstruct genome #####
#g.recon = function(x, u = d0.u){
#  i0 = charToRaw(u$data[which(u$diff == "")])
#  i1 = as.numeric(strsplit(x, ";")[[1]])
#  stopifnot(all(i1 >= 1 & i1 <= length(i0)))
#  i0[i1] = as.raw(bitwXor(as.integer(i0[i1]), 1))
#  return(rawToChar(i0))
#}

##### reconstruct genotype differences #####
#gD = function(dIff, sTart, u = d0.u){ # dIff = "1;14;2;3;57"; sTart = "AB"
#  repeat{
#    tAg = paste0(sample(letters, min(nchar(u$tag))+3, replace = T), collapse = "")
#    if(!(any(u$tag %in% tAg))){break}
#  }
#  if(u$diff[u$tag == sTart] != ""){
#    u0 = strsplit(c(dIff, u$diff[u$tag == sTart]), ";")
#    u0 = as.numeric(c(setdiff(u0[[1]], u0[[2]]), setdiff(u0[[2]], u0[[1]])))
#    x = paste(u0[order(u0)], collapse = ";")
#  }else{ x = dIff }
#  if(length(grep(x, u$diff))>0){return(NULL)}else{
#    x0 = g.recon(x)
#    x = c(tAg, as.numeric(strsplit(x0, "")[[1]]), x0, x, u$VT[u$tag == sTart], u$SC[u$tag == sTart])
#    return(c(x, rep(NA, ncol(u)-length(x))))
#  }
#}

##### proportion of genes from tags to ratio #####
#t2r = function(idx, G = Gm){ #function(tags, u = d0.u){
#  return( drop(tabulate(idx, nbins = nrow(G)) %*% G) / length(idx) )
#  return(as.numeric(colMeans(as.data.frame(do.call(rbind, lapply(strsplit(u$data[match(tags, u$tag)], ""), as.integer))))))
#}

##### Calculate pi and omega from reference #####
#piOmega = function(g, sEl, fSel){
#  g = as.numeric(strsplit(g, "")[[1]]) * sEl
#  if(length(fSel)==length(sEl)){g1 = 0}else{g1 = sum(g[-fSel])}
#  if(length(fSel)==0){g0 = 0}else{g0 = sum(g[fSel])}
#  return(c(g0, g1))
#}

##### VT*SC distribution vector #####
#vtsc = function(tags, u = d0.u){
#  vc0 = data.frame(src = unique(u$paste)[order(unique(u$paste))], f = 0)
#  vc1 = table(d0.u$paste[match(tags, d0.u$tag)])
#  vc0$f = vc1[match(vc0$src, names(vc1))]
#  vc0$f[is.na(vc0$f)] = 0
#  return(vc0$f)
#}


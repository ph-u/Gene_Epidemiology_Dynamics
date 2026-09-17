#!/bin/env Rscript
# author: ph-u
# script: other_data_import.r
# desc: import sample data ideas
# in: source("other_data_import.r")
# out: NA
# arg: 0
# date: 20260916

##### Invasiveness dataset #####
a = list.files("../data/invasiveness_clonal_ukhsa/", pattern = ".rda$", full.names=T)
for(i in seq_len(length(a))){
  load(a[i])
};rm(i)
pneu.Adult = merge(S_pneumoniae_adult_serotype, read.csv(sub("rda", "csv", a[grep("adult", a)]), header = T, row.names = 1), all = T)
pneu.Mixed = merge(S_pneumoniae_mixed_strain, read.csv(sub("rda", "csv", a[grep("mixed", a)]), header = T), all = T)
pneu.ifSer = merge(S_pneumoniae_infant_serotype, read.csv(sub("rda", "csv", a[grep("infant_ser", a)]), header = T), all = T)
pneu.ifStr = merge(S_pneumoniae_infant_strain, read.csv(sub("rda", "csv", a[grep("infant_str", a)]), header = T), all = T)
rm(a, S_pneumoniae_adult_serotype, S_pneumoniae_mixed_strain, S_pneumoniae_infant_serotype, S_pneumoniae_infant_strain)

##### Israel dataset prelim #####
qqGpsDedu = read.table("../data/israel_prelim/avai-microreact-table-qq-gps-dedu.csv", header=T, quote = "", comment.char="", sep = "\t")
cog.Freq = read.csv("../data/israel_prelim/Israel_prelim.cog_frequencies.csv", header=T)
genePa = read.table("../data/israel_prelim/Israel_prelim.gene_pa.csv", header=T, sep = "\t", comment.char="", quote="")

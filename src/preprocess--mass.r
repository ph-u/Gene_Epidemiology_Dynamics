#!/bin/env Rscript
# author: ph-u
# script: preprocess--mass.r
# desc: preprocess data/mass.input
# in: Rscript preprocess--mass.r
# out: data/pp-*
# arg: 0
# date: 20261029

set.seed(20261029)
a0 = read.csv("../data/input.csv", header = T)
a1 = read.table("../data/mass.input", header = T, sep = "\t")

##### data.tsv #####
i0 = "../data/pp-mass.tsv"
cat(date(), ": Reorganizing", i0, "\n")
a0$Value[a0$Type=="data"] = i0
gCols = colnames(a1)[(which(colnames(a1)=="SC")+1):ncol(a1)]
a1$Age = sample(3:84, nrow(a1), replace = TRUE) + a1$Time
a1$Age = ifelse(a1$Age > 84, 84, a1$Age) # original dataset at (Huang et al. 2009)
a2 = data.frame(Time=a1$Time, Deme = "d1", Age = a1$Age, VT1 = a1$VT, SC = a1$SC, a1[,gCols, drop = FALSE])

write.table(a2, file = i0, sep = "\t", row.names = F, quote = F)

##### rollout.tsv #####
i0 = "../data/pp-roll.tsv"
cat(date(), ": Formulating", i0, "\n")
a0$Value[a0$Type=="rollout"] = i0
a3 = data.frame(Deme = unique(a2$Deme), Vaccine = colnames(a2)[grep("VT", colnames(a2))], StartMonth = 0, AgeWindow = max(a2$Age), Uptake = 1)

write.table(a3, file = i0, sep = "\t", row.names = F, quote = F)

write.csv(a0, "../data/input.csv", row.names = F, quote = F)
cat(date(), ": Done\n")

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
library(BRREWABC) # https://github.com/GaelBn/BRREWABC
nCPU = as.integer(Sys.getenv("LSB_DJOB_NUMPROC", unset = "1"))
sEed = read.csv("../raw/seed.csv", header = F)[[1]][as.numeric(argv[2])]

##### per-run directory: concurrent seeds must not share ../data/tmp #####
oUtDir = file.path("..", "data", paste0("run_", sEed))
dir.create(file.path(oUtDir, "tmp"), recursive = T, showWarnings = F)
oUtDir = normalizePath(oUtDir, mustWork = T)
Sys.setenv(NFDS_STATE = file.path(oUtDir, "setupState.RData"))

source("setup.r")

save(list = ls(envir = globalenv()), file = Sys.getenv("NFDS_STATE"), envir = globalenv(), compress = T)
set.seed(sEed)

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
  experiment_folderpath  = oUtDir,
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
  saveRDS(list(G = r$G, meta = gT), file.path(oUtDir, paste0("nfdsG_", sP, ".rds")))
};rm(p)
write.csv(do.call(rbind, oUt), file.path(oUtDir, paste0("nfdsGenotypes_", gsub(" ", "-", date()), "_", sEed, ".csv")), row.names = F, quote = F)


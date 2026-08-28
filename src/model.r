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

save(list = ls(envir = globalenv()), file = "../data/setupState.RData", envir = globalenv(), compress = T)

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


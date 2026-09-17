#!/usr/bin/env Rscript
# author: ph-u, nickjcroucher
# script: model.r
# desc: NFDS ABCSMC -- metapopulation, multi-vaccine, graded roll-out, co-infection
# in: Rscript model.r [../raw/input.csv] [seed entry]
# out: ../data/run_<seed>/{setupState.RData, nfdsG_*.rds, nfdsGenotypes_*.csv, res/}
# arg: 1
# date: 20260918 (metapopulation, multi-vaccine), 20260819

argv = (commandArgs(T))
if(length(argv) != 2){ argv = c("../raw/input.csv", 1) }

##### env #####
library(BRREWABC) # https://github.com/GaelBn/BRREWABC
nCPU = as.integer(Sys.getenv("LSB_DJOB_NUMPROC", unset = "1"))
sEed = read.csv("../raw/seed.csv", header = F)[[1]][as.numeric(argv[2])]
stopifnot(length(sEed) == 1, !is.na(sEed))

##### per-run directory: concurrent seeds must not share tmp/ #####
oUtDir = file.path("..", "data", paste0("run_", sEed))
dir.create(file.path(oUtDir, "tmp"), recursive = T, showWarnings = F)
oUtDir = normalizePath(oUtDir, mustWork = T)
Sys.setenv(NFDS_STATE = file.path(oUtDir, "setupState.RData"))

source("setup.r")

## workers are fresh R processes; they restore this rather than re-running setup.r.
## Saved BEFORE set.seed(), so every worker does not inherit one shared RNG state.
save(list = ls(envir = globalenv()), file = Sys.getenv("NFDS_STATE"),
     envir = globalenv(), compress = T)
set.seed(sEed)

##### Parameters: 7 shared + one efficacy per vaccine #####
## all priors are unif(0,1); m.nfds rescales internally
## (vSelected x .5 per month, migration x .2, mWithin x .1)
prior_dist <- list(nfds = c(
  list(c("propStrong", "unif", 0, 1),      # fraction of loci under strong NFDS
       c("fSelected",  "unif", 1e-6, .22), # strong NFDS coefficient
       c("wSelected",  "unif", 1e-6, .15), # weak NFDS coefficient
       c("migration",  "unif", 0, 1),      # external immigration into the metapopulation
       c("coInf",      "unif", 0, 1),      # co-infection: 0 = Poisson, >0 = overdispersed
       c("vSelMod",    "unif", 0, 1),      # cohort roll-out ramp width; 0 = step
       c("mWithin",    "unif", 0, 1)),     # between-deme movement
  lapply(seq_len(nVac), function(v) c(paste0("vSelected", v), "unif", 0, 1))))

##### Pilot: where can the distance actually reach? #####
## Set PILOT=1 for a short diagnostic run; read the quantiles, then set
## distance_threshold_min just under them and run for real with PILOT unset.
if(nzchar(Sys.getenv("PILOT"))){
  pR = prior_dist$nfds
  d  = replicate(as.integer(Sys.getenv("PILOT")), nfds_jsd(
         setNames(lapply(pR, function(z) runif(1, as.numeric(z[3]), as.numeric(z[4]))),
                  vapply(pR, `[`, character(1), 1)), mIg0))
  print(summary(d)); print(quantile(d, c(.01, .05, .10)))
  cat("rejected simulations:", sum(d == ncol(mIg0) * log(2)), "of", length(d), "\n")
  quit(save = "no")
}

##### ABCSMC-NFDS #####
res <- abcsmc(
  model_list             = list(nfds = nfds_jsd),
  prior_dist             = prior_dist,
  ss_obs                 = mIg0,
  nb_threshold           = 1,
  nb_acc_prtcl_per_gen   = 500,
  max_number_of_gen      = 500,   # runs were still descending at 10; 500 is unreachable
  new_threshold_quantile = .9,
  distance_threshold_min = .05,  # set from the PILOT quantiles above
  acceptance_rate_min    = .001,
  experiment_folderpath  = oUtDir,
  max_concurrent_jobs    = nCPU,
  use_lhs_for_first_iter = TRUE,
  verbose                = TRUE,
  progressbar            = FALSE
)

##### Export: re-run each posterior particle, keeping the genotype state #####
pOst = res$particles[res$particles$gen == max(res$particles$gen), ]
oUt  = vector("list", nrow(pOst))
for(p in seq_len(nrow(pOst))){
  sP = sEed + p; set.seed(sP) # recorded, so the re-run is reproducible
  vS = vapply(seq_len(nVac), function(v) pOst[[paste0("vSelected", v)]][p], 0)
  r  = m.nfds(pOst$propStrong[p], pOst$fSelected[p], pOst$wSelected[p], vS,
              pOst$migration[p], pOst$coInf[p], pOst$vSelMod[p], pOst$mWithin[p],
              keepGenotypes = T)
  if(is.null(r)){ next }   # an accepted particle can still fail under a different seed
  gT = r$genotypes; gT$particle = p; gT$seed = sP
  colnames(r$demeCount) = paste0("n_", dLev)
  gT = cbind(gT, as.data.frame(r$demeCount))
  oUt[[p]] = gT
  saveRDS(list(G = r$G, meta = gT, demeCount = r$demeCount,
               par = pOst[p, , drop = FALSE]),
          file.path(oUtDir, paste0("nfdsG_", sP, ".rds")))
};rm(p)
write.csv(do.call(rbind, oUt),
          file.path(oUtDir, paste0("nfdsGenotypes_", gsub(" ", "-", date()), "_", sEed, ".csv")),
          row.names = F, quote = F)


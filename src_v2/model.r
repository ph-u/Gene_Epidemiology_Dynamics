#!/usr/bin/env Rscript
# author: ph-u, nickjcroucher
# script: model.r
# desc: NFDS ABCSMC v2 -- ABC fitting and posterior export
# Features: metapopulation (heterogeneous mWithin matrix) | multi-vaccine |
#           host age structure | co-infection | infection/colonisation history
# in: Rscript model.r [../raw/input.csv] [seed entry]
# out: ../data/run_<seed>/{setupState.RData, tmp/, res/,
#                          nfdsG_*.rds, nfdsHostState_*.rds,
#                          nfdsHHist_*.rds, nfdsGenotypes_*.csv}
# arg: 1
# date: 20260919

argv = (commandArgs(TRUE))
if(length(argv) != 2){ argv = c("../raw/input.csv", 1) }

##### Environment #####
library(BRREWABC)
nCPU = as.integer(Sys.getenv("LSB_DJOB_NUMPROC", unset = "1"))
sEed = read.csv("../raw/seed.csv", header = FALSE)[[1]][as.numeric(argv[2])]
stopifnot(length(sEed) == 1, !is.na(sEed))

##### Per-run directory (concurrent seeds must NOT share tmp/) #####
oUtDir = file.path("..", "data", paste0("run_", sEed))
dir.create(file.path(oUtDir, "tmp"), recursive = TRUE, showWarnings = FALSE)
oUtDir = normalizePath(oUtDir, mustWork = TRUE)
Sys.setenv(NFDS_STATE = file.path(oUtDir, "setupState.RData"))

source("setup.r")

## Save BEFORE set.seed() -- workers must not share one RNG state
save(list = ls(envir = globalenv()), file = Sys.getenv("NFDS_STATE"),
     envir = globalenv(), compress = TRUE)
set.seed(sEed)

##### Priors (all unif(0,1); m.nfds rescales internally) #####
## Parameters:
##   propStrong       -- fraction of loci under strong NFDS
##   fSelected        -- strong NFDS coefficient   (rescaled x0.22)
##   wSelected        -- weak NFDS coefficient     (rescaled x0.15)
##   migration        -- external immigration       (rescaled x0.2)
##   coInf            -- co-infection overdispersion
##   mWithin_d1_d2    -- directed migration d1->d2  (rescaled x0.1)
##   vSelected_v      -- vaccine v efficacy         (rescaled x0.5)
##
## Total parameters: 5 + nD*(nD-1) + nVac
prior_dist <- list(nfds = c(
  list(c("propStrong", "unif", 0, 1),
       c("fSelected",  "unif", 1e-6, .22),
       c("wSelected",  "unif", 1e-6, .15),
       c("migration",  "unif", 0, 1),
       c("coInf",      "unif", 0, 1)),
  ## Heterogeneous migration: one parameter per DIRECTED deme pair
  do.call(c, lapply(seq_len(nD), function(d1)
    lapply(seq_len(nD)[seq_len(nD) != d1], function(d2)
      list(c(paste0("mWithin_", d1, "_", d2), "unif", 0, 1))))),
  ## Per-vaccine efficacy
  lapply(seq_len(nVac), function(v) c(paste0("vSelected", v), "unif", 0, 1))))

##### Pilot mode: estimate achievable distance range before the full fit #####
## Usage: PILOT=200 Rscript model.r ../raw/input.csv 1
## Read the printed quantiles; set distance_threshold_min just below them.
if(nzchar(Sys.getenv("PILOT"))){
  pR = prior_dist$nfds
  d  = replicate(as.integer(Sys.getenv("PILOT")), nfds_jsd(
         setNames(lapply(pR, function(z) runif(1, as.numeric(z[3]), as.numeric(z[4]))),
                  vapply(pR, `[`, character(1), 1)), mIg0))
  print(summary(d)); print(quantile(d, c(.01, .05, .10)))
  cat("rejected simulations:", sum(d == ncol(mIg0) * log(2)), "of", length(d), "\n")
  quit(save = "no")
}

##### ABC-SMC fitting #####
res <- abcsmc(
  model_list             = list(nfds = nfds_jsd),
  prior_dist             = prior_dist,
  ss_obs                 = mIg0,
  nb_threshold           = 1,
  nb_acc_prtcl_per_gen   = 500,
  max_number_of_gen      = 30,
  new_threshold_quantile = .9,
  distance_threshold_min = .05,       # set from PILOT quantiles above
  acceptance_rate_min    = .001,
  experiment_folderpath  = oUtDir,
  max_concurrent_jobs    = nCPU,
  use_lhs_for_first_iter = TRUE,
  verbose                = TRUE,
  progressbar            = FALSE
)

##### Export: re-run each posterior particle to obtain full genotype and host state #####
pOst = res$particles[res$particles$gen == max(res$particles$gen), ]
oUt  = vector("list", nrow(pOst))

for(p in seq_len(nrow(pOst))){
  sP = sEed + p; set.seed(sP)

  ## Assemble per-vaccine vector
  vS = vapply(seq_len(nVac), function(v) pOst[[paste0("vSelected", v)]][p], 0)

  ## Assemble heterogeneous migration matrix
  mW = matrix(0, nD, nD)
  for(d1 in seq_len(nD)) for(d2 in seq_len(nD))
    if(d1 != d2) mW[d1, d2] = pOst[[paste0("mWithin_", d1, "_", d2)]][p]

  r = m.nfds(pOst$propStrong[p], pOst$fSelected[p], pOst$wSelected[p], vS,
             pOst$migration[p], pOst$coInf[p], mW,
             keepGenotypes = TRUE)
  if(is.null(r)){ next }   # accepted particle may fail under a different seed

  gT = r$genotypes; gT$particle = p; gT$seed = sP
  colnames(r$demeCount) = paste0("n_", dLev)
  gT = cbind(gT, as.data.frame(r$demeCount))
  oUt[[p]] = gT

  ## (i) Genotype table + deme counts + history summary + parameters
  saveRDS(list(G          = r$G,
               genotypes  = gT,
               demeCount  = r$demeCount,
               hHistStats = r$hHistStats,
               par        = pOst[p, , drop = FALSE]),
          file.path(oUtDir, paste0("nfdsG_", sP, ".rds")))

  ## (ii) Host state snapshot (one row per bacterium at the final generation)
  saveRDS(r$hostState,
          file.path(oUtDir, paste0("nfdsHostState_", sP, ".rds")))

  ## (iii) Infection/colonisation history (cleared episodes during the simulation)
  saveRDS(r$hHist,
          file.path(oUtDir, paste0("nfdsHHist_", sP, ".rds")))
};rm(p)

write.csv(do.call(rbind, oUt),
          file.path(oUtDir, paste0("nfdsGenotypes_", gsub(" ", "-", date()),
                                   "_", sEed, ".csv")),
          row.names = FALSE, quote = FALSE)


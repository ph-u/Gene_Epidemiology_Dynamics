##### shared fixture for test-src.R / test-nfds.R / test-model.R #####
## testthat sources every helper-*.R before any test-*.R, into the same
## environment, so everything below is visible to the tests.

SRC_DIR <- normalizePath(Sys.getenv("NFDS_SRC", unset = file.path("..", "src")),
                         mustWork = TRUE)

##### small but structurally faithful: 3 demes, 2 nested vaccines, staggered roll-out #####
make_fixture <- function(seed = 42L) {
  dir <- tempfile("nfdsfix"); set.seed(seed)
  for (s in c("raw", "src", "data")) dir.create(file.path(dir, s), recursive = TRUE)

  nLin <- 4L; nVar <- 6L; nDm <- 3L
  vt1  <- c(1L, 0L, 0L, 0L)          # vaccine 1 covers lineage 1
  vt2  <- c(1L, 1L, 0L, 0L)          # vaccine 2 nests vaccine 1, adds lineage 2
  tPt  <- c(-6, 0, 6, 12)
  mix  <- list(c(2,2,2,2), c(2,2,2,2), c(1,2,3,2), c(1,1,3,3))   # VT lineages decline

  rows <- do.call(rbind, lapply(seq_along(tPt), function(ti)
            do.call(rbind, lapply(seq_len(nDm), function(dd)
              data.frame(Time = tPt[ti], Dm = dd,
                         lin  = rep(seq_len(nLin), mix[[ti]]))))))

  ## loci 1-4 are one-hot lineage markers, so genotype determines SC and VT profile,
  ## as ~1090 intermediate-frequency loci do for real
  gene <- cbind(diag(nLin)[rows$lin, , drop = FALSE],
                matrix(rbinom(nrow(rows) * nVar, 1, .5), ncol = nVar),
                1L, 0L)                          # fixed + absent: d0.keep drops both
  colnames(gene) <- sprintf("CLS%05d", seq_len(ncol(gene)))

  meta <- data.frame(Time = rows$Time, Deme = paste0("D", rows$Dm),
                     VT1 = vt1[rows$lin], VT2 = vt2[rows$lin],
                     SC  = rows$lin)             # SC must be the LAST metadata column
  write.table(cbind(meta, as.data.frame(gene)), file.path(dir, "raw", "data.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  ## staggered roll-out; D3 never receives VT2, exercising the vStart = Inf path
  write.table(data.frame(Deme       = c("D1","D1","D2","D2","D3"),
                         Vaccine    = c("VT1","VT2","VT1","VT2","VT1"),
                         StartMonth = c(0, 6, 3, 9, 6)),
              file.path(dir, "raw", "rollout.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  write.csv(data.frame(
    Type  = c("data", "rollout", "population", "percentage initial infected"),
    Value = c(file.path(dir, "raw", "data.tsv"), file.path(dir, "raw", "rollout.tsv"),
              "1200", "100")),
    file.path(dir, "raw", "input.csv"), row.names = FALSE, quote = FALSE)

  write.table(c(11, 22, 33), file.path(dir, "raw", "seed.csv"),
              row.names = FALSE, col.names = FALSE, sep = ",")

  ## link, not copy -- src/ holds the only copy of each script.
  ## Linking the directory would not work: setwd() canonicalises, so ".." would
  ## resolve back to the real repo root.
  for (f in c("src.r", "nfds.r", "setup.r", "model.r")) {
    tgt <- normalizePath(file.path(SRC_DIR, f), mustWork = TRUE)
    ok  <- suppressWarnings(file.symlink(tgt, file.path(dir, "src", f)))
    if (!ok) {
      warning("symlink unavailable; falling back to a copy of ", f, call. = FALSE)
      file.copy(tgt, file.path(dir, "src", f), overwrite = TRUE)
    }
  }
  dir
}

##### stub standing in for BRREWABC::abcsmc (columns must match the fixture's nVac = 2) #####
stub_abcsmc <- function(model_list, prior_dist, ss_obs, ...) {
  assign("ABCSMC_CALL",
         c(list(model_list = model_list, prior_dist = prior_dist, ss_obs = ss_obs),
           list(...)), envir = globalenv())
  list(particles = data.frame(
         gen        = c(1, 1),
         propStrong = c(.20, .30),  fSelected  = c(.10, .12),
         wSelected  = c(.001, .002), migration = c(.02, .03),
         coInf      = c(0, .05),    vSelMod    = c(.1, .2), mWithin = c(.05, .10),
         vSelected1 = c(.10, .15),  vSelected2 = c(.08, .12)),
       thresholds = c(1, .5))
}

##### keeps tests readable and immune to argument-order changes #####
mn <- function(propStrong = .25, fSelected = .1, wSelected = .001, vSel = c(.2, .15),
               migration = .02, coInf = 0, vSelMod = .1, mWithin = .05, ...)
  m.nfds(propStrong, fSelected, wSelected, vSel, migration, coInf, vSelMod, mWithin, ...)

p0 <- function(...) modifyList(
  list(propStrong = .25, fSelected = .1, wSelected = .001, migration = .02,
       coInf = 0, vSelMod = .1, mWithin = .05, vSelected1 = .2, vSelected2 = .15),
  list(...))

##### run setup.r (data prep only) or model.r (full pipeline); both cached #####
load_model <- function(force = FALSE) {
  if (!force && isTRUE(getOption("nfds.loaded"))) return(invisible(TRUE))
  dir <- make_fixture()
  assign("abcsmc", stub_abcsmc, envir = globalenv())
  assign("FIXTURE_DIR", dir, envir = globalenv())
  old <- setwd(file.path(dir, "src")); on.exit(setwd(old))
  stopifnot(identical(normalizePath(".."), normalizePath(dir)))
  Sys.setenv(NFDS_SRC = getwd())
  suppressMessages(source("setup.r"))
  options(nfds.loaded = TRUE)
  invisible(TRUE)
}

load_full <- function() {
  if (isTRUE(getOption("nfds.full"))) return(invisible(TRUE))
  dir <- make_fixture()
  assign("abcsmc", stub_abcsmc, envir = globalenv())
  assign("FIXTURE_DIR", dir, envir = globalenv())
  old <- setwd(file.path(dir, "src")); on.exit(setwd(old))
  Sys.setenv(NFDS_SRC = getwd())
  suppressMessages(source("model.r"))
  options(nfds.full = TRUE, nfds.loaded = TRUE)
  invisible(TRUE)
}


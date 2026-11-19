##### shared fixture for test-src.R / test-nfds.R / test-model.R #####
## testthat sources every helper-*.R before any test-*.R, into the same environment.

SRC_DIR <- normalizePath(Sys.getenv("NFDS_SRC", unset = file.path("..", "src")),
                         mustWork = TRUE)

##### small but structurally faithful: 3 demes, 2 nested vaccines, 4 host ages #####
make_fixture <- function(seed = 42L) {
  dir <- tempfile("nfdsfix"); set.seed(seed)
  for (s in c("raw", "src", "data")) dir.create(file.path(dir, s), recursive = TRUE)

  nLin <- 4L; nVar <- 6L; nDm <- 3L
  vt1  <- c(1L, 0L, 0L, 0L)          # vaccine 1 covers lineage 1
  vt2  <- c(1L, 1L, 0L, 0L)          # vaccine 2 nests vaccine 1, adds lineage 2
  aLv  <- c(6L, 18L, 36L, 60L)       # host age in MONTHS, same unit as AgeWindow
  tPt  <- c(-6, 0, 6, 12)
  mix  <- list(c(2,2,2,2), c(2,2,2,2), c(1,2,3,2), c(1,1,3,3))   # VT lineages decline

  rows <- do.call(rbind, lapply(seq_along(tPt), function(ti)
            do.call(rbind, lapply(seq_len(nDm), function(dd)
              data.frame(Time = tPt[ti], Dm = dd,
                         Age  = rep(aLv, length.out = sum(mix[[ti]])),
                         lin  = rep(seq_len(nLin), mix[[ti]]))))))

  ## loci 1-4 are one-hot lineage markers, so genotype determines SC and VT profile,
  ## as ~1090 intermediate-frequency loci do for real
  gene <- cbind(diag(nLin)[rows$lin, , drop = FALSE],
                matrix(rbinom(nrow(rows) * nVar, 1, .5), ncol = nVar),
                1L, 0L)                          # fixed + absent: d0.keep drops both
  colnames(gene) <- sprintf("CLS%05d", seq_len(ncol(gene)))

  ## column order: Time | Deme | Age | VT1 | VT2 | SC | genes   (SC LAST metadata)
  meta <- data.frame(Time = rows$Time, Deme = paste0("D", rows$Dm), Age = rows$Age,
                     VT1 = vt1[rows$lin], VT2 = vt2[rows$lin], SC = rows$lin)
  write.table(cbind(meta, as.data.frame(gene)), file.path(dir, "data", "data.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  ## VT1 timing varies across demes with uptake held at 1, so D1-vs-D3 isolates
  ## the roll-out date. D3 has no VT2 row, exercising the vStart = Inf path.
  write.table(data.frame(Deme       = c("D1","D1","D2","D2","D3"),
                         Vaccine    = c("VT1","VT2","VT1","VT2","VT1"),
                         StartMonth = c(0, 6, 3, 9, 9),
                         AgeWindow  = c(24, 24, 24, 24, 24),
                         Uptake     = c(1, .8, 1, .7, 1)),
              file.path(dir, "raw", "rollout.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  write.csv(data.frame(
    Type  = c("data", "rollout", "population", "percentage initial infected"),
    Value = c(file.path(dir, "data", "data.tsv"), file.path(dir, "raw", "rollout.tsv"),
              "1200", "100")),
    file.path(dir, "data", "input.csv"), row.names = FALSE, quote = FALSE)

  write.table(c(11, 22, 33), file.path(dir, "raw", "seed.csv"),
              row.names = FALSE, col.names = FALSE, sep = ",")

  ## link, not copy -- src/ holds the only copy of each script. Linking the
  ## directory would fail: setwd() canonicalises, so ".." would resolve to the repo root.
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

##### stub standing in for BRREWABC::abcsmc #####
## FIX B1: mWithin is now a heterogeneous nD x nD matrix; the prior has one
## named parameter mWithin_d1_d2 per directed deme pair (nD*(nD-1) of them).
## nD is available in globalenv by the time abcsmc() is called (model.r sources
## setup.r first), so we generate the column names dynamically.
stub_abcsmc <- function(model_list, prior_dist, ss_obs, ...) {
  assign("ABCSMC_CALL",
         c(list(model_list = model_list, prior_dist = prior_dist, ss_obs = ss_obs),
           list(...)), envir = globalenv())
  mW_nms  <- unlist(lapply(seq_len(nD), function(d1)
    paste0("mWithin_", d1, "_", seq_len(nD)[seq_len(nD) != d1])))
  mW_vals <- setNames(lapply(mW_nms, function(.) c(.05, .08)), mW_nms)
  list(particles = do.call(data.frame, c(
         list(gen        = c(1, 1),
              propStrong = c(.20, .30),  fSelected  = c(.10, .12),
              wSelected  = c(.001,.002), migration  = c(.02, .03),
              coInf      = c(0,   .05)),
         mW_vals,
         list(vSelected1 = c(.10, .15),  vSelected2 = c(.08, .12)))),
       thresholds = c(1, .5))
}

##### keeps tests readable and immune to argument-order changes #####
## FIX B2: m.nfds() expects mWithin as an nD x nD matrix. mn() accepts a scalar
## for convenience and converts it: scalar s -> s * off-diagonal matrix. An actual
## matrix passes through unchanged, so mn(mWithin = matrix(...)) also works.
mn <- function(propStrong = .25, fSelected = .1, wSelected = .001, vSel = c(.2, .15),
               migration = .02, coInf = 0, mWithin = .05, ...)
  m.nfds(propStrong, fSelected, wSelected, vSel, migration, coInf,
         if (is.matrix(mWithin)) mWithin
         else matrix(mWithin, nD, nD) * (1L - diag(nD)), ...)

## FIX B3: nfds_jsd() assembles the mWithin matrix from named parameters
## mWithin_d1_d2. p0() now emits one entry per directed deme pair.
## migration and all mWithin pairs default to 0: mixing noise would swamp
## the weak vaccine signals tested by nfds-suite tests (L6 fix).
p0 <- function(...) modifyList(
  c(list(propStrong = .25, fSelected = .1, wSelected = .001, migration = 0, coInf = 0),
    setNames(as.list(rep(0, nD * (nD - 1))),
             unlist(lapply(seq_len(nD), function(d1)
               paste0("mWithin_", d1, "_", seq_len(nD)[seq_len(nD) != d1])))),
    list(vSelected1 = .2, vSelected2 = .15)),
  list(...))

##### temporarily override a global, restoring it when the test block exits #####
swap <- function(nm, value, env = parent.frame()) {
  old <- get(nm, envir = globalenv())
  do.call(on.exit, list(bquote(assign(.(nm), .(old), envir = globalenv())),
                        add = TRUE), envir = env)
  assign(nm, value, envir = globalenv())
}

##### run setup.r (data prep only) or model.r (full pipeline); both cached #####
load_model <- function(force = FALSE) {
  if (!force && isTRUE(getOption("nfds.loaded"))) return(invisible(TRUE))
  dir <- make_fixture()
  assign("abcsmc", stub_abcsmc, envir = globalenv())
  assign("FIXTURE_DIR", dir, envir = globalenv())
  old <- setwd(file.path(dir, "src")); on.exit(setwd(old))
  stopifnot(identical(normalizePath(".."), normalizePath(dir)))
  Sys.setenv(NFDS_SRC = getwd())
  ## FIX B4: setup.r defaults to ../raw/input.csv but the fixture writes input.csv
  ## to ../data/input.csv. Assign argv before sourcing so the guard is satisfied.
  assign("argv", c("../data/input.csv", "1"), envir = globalenv())
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
  stopifnot(identical(normalizePath(".."), normalizePath(dir)))
  Sys.setenv(NFDS_SRC = getwd())
  ## model.r sets argv itself from commandArgs(); its default is ../data/input.csv
  ## (fixed in src/model.r B5). No argv assignment needed here.
  suppressMessages(source("model.r"))
  options(nfds.full = TRUE, nfds.loaded = TRUE)
  invisible(TRUE)
}

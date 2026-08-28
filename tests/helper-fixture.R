##### shared fixture for test-src.R / test-nfds.R / test-model.R #####
## testthat sources every helper-*.R before any test-*.R, into the same
## environment, so all three definitions below are visible to the tests.

SRC_DIR <- normalizePath(Sys.getenv("NFDS_SRC", unset = file.path("..", "src")),
                         mustWork = TRUE)

##### a tiny but structurally faithful dataset #####
make_fixture <- function(seed = 42L) {
  dir <- tempfile("nfdsfix"); set.seed(seed)
  for (s in c("raw", "src", "data")) dir.create(file.path(dir, s), recursive = TRUE)

  ## Four lineages. Loci 1-4 are one-hot lineage markers, so a genotype always
  ## determines its lineage -- as ~1090 intermediate-frequency loci do for real.
  nLin  <- 4L; nVar <- 6L
  linVT <- c(1L, 0L, 1L, 0L)
  tPoint <- c(-6, 0, 6, 12)
  mix    <- list(c(4,4,4,4), c(4,4,4,4), c(2,4,2,4), c(1,3,1,3))  # VT lineages decline
  lin    <- unlist(lapply(mix, function(m) rep(seq_len(nLin), m)))

  base <- matrix(0L, nLin, nLin); diag(base) <- 1L
  gene <- cbind(base[lin, , drop = FALSE],
                matrix(rbinom(length(lin) * nVar, 1, 0.5), ncol = nVar),
                1L, 0L)                                  # fixed + absent: d0.keep drops both
  colnames(gene) <- sprintf("CLS%05d", seq_len(ncol(gene)))

  meta <- data.frame(Time = rep(tPoint, vapply(mix, sum, 0)),
                     VT   = linVT[lin],
                     SC   = lin)                          # SC must be the last meta column
  write.table(cbind(meta, as.data.frame(gene)), file.path(dir, "raw", "data.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  write.csv(data.frame(
    Type  = c("data", "population", "percentage initial infected", "vaccine start month"),
    Value = c(file.path(dir, "raw", "data.tsv"), "1000", "100", "0")),
    file.path(dir, "raw", "input.csv"), row.names = FALSE, quote = FALSE)

  write.table(c(11, 22, 33), file.path(dir, "raw", "seed.csv"),
              row.names = FALSE, col.names = FALSE, sep = ",")

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
stub_abcsmc <- function(model_list, prior_dist, ss_obs, ...) {
  assign("ABCSMC_CALL",
         c(list(model_list = model_list, prior_dist = prior_dist, ss_obs = ss_obs),
           list(...)),
         envir = globalenv())
  list(particles = data.frame(gen        = c(1, 1),
                              propStrong = c(.20, .30), fSelected = c(.10, .12),
                              wSelected  = c(.001, .002), vSelected = c(.10, .15),
                              migration  = c(.02, .03)),
       thresholds = c(1, .5))
}

##### run model.r once against the fixture; cached across test files #####
# model.r calls source("src.r") without `local=`, so everything must land in
# globalenv() for the functions' lexical scope to reach the derived objects.
load_model <- function(force = FALSE) {
  if (!force && isTRUE(getOption("nfds.loaded"))) return(invisible(TRUE))
  dir <- make_fixture()
  assign("abcsmc", stub_abcsmc, envir = globalenv())
  assign("FIXTURE_DIR", dir, envir = globalenv())
  old <- setwd(file.path(dir, "src")); on.exit(setwd(old))
  stopifnot(identical(normalizePath(".."), normalizePath(dir)))
  Sys.setenv(NFDS_SRC = getwd())          # what a worker would read
  source("setup.r")                        # data prep only; no abcsmc, no set.seed
  options(nfds.loaded = TRUE)
  invisible(TRUE)
}

load_full <- function() {
  dir <- make_fixture()
  assign("abcsmc", stub_abcsmc, envir = globalenv())
  assign("FIXTURE_DIR", dir, envir = globalenv())
  old <- setwd(file.path(dir, "src")); on.exit(setwd(old))
  Sys.setenv(NFDS_SRC = getwd())
  source("model.r")
  invisible(TRUE)
}

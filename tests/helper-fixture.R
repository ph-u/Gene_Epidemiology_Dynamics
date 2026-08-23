##### shared fixture for test-src.R / test-nfds.R / test-model.R #####
## testthat sources every helper-*.R before any test-*.R, into the same
## environment, so all three definitions below are visible to the tests.

SRC_DIR <- normalizePath(Sys.getenv("NFDS_SRC", unset = file.path("..", "src")),
                         mustWork = TRUE)

##### a tiny but structurally faithful dataset #####
make_fixture <- function(seed = 42L) {
  dir <- tempfile("nfdsfix"); set.seed(seed)
  for (s in c("raw", "src", "data")) dir.create(file.path(dir, s), recursive = TRUE)

  nLoci  <- 12L
  tPoint <- c(-6, 0, 6, 12)          # negative month exercises the pre-vaccine path
  nPer   <- c(15, 15, 12, 10)

  meta <- data.frame(Time = rep(tPoint, nPer),
                     VT   = rbinom(sum(nPer), 1, 0.4),
                     SC   = sample(1:3, sum(nPer), replace = TRUE))   # SC must be last
  gene <- matrix(rbinom(nrow(meta) * nLoci, 1, 0.5), ncol = nLoci)
  colnames(gene) <- sprintf("CLS%05d", seq_len(nLoci))
  gene[, 1] <- 1L; gene[, 2] <- 0L   # fixed + absent: d0.keep must drop both
  write.table(cbind(meta, as.data.frame(gene)), file.path(dir, "raw", "data.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  write.csv(data.frame(
    Type  = c("data", "population", "percentage initial infected", "vaccine start month"),
    Value = c(file.path(dir, "raw", "data.tsv"), "1000", "100", "0")),
    file.path(dir, "raw", "input.csv"), row.names = FALSE, quote = FALSE)

  write.table(c(11, 22, 33), file.path(dir, "raw", "seed.csv"),
              row.names = FALSE, col.names = FALSE, sep = ",")

  ## link, do not copy -- src/ keeps the only copy of each script.
  ## NB: linking the *directory* would not work; setwd() canonicalises the
  ## path, so ".." would resolve back to the real repo root.
  for (f in c("src.r", "nfds.r", "model.r")) {
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
  stopifnot(identical(normalizePath(".."), normalizePath(dir)))  # redirection worked
  source("model.r")                      # commandArgs(T) is empty -> defaults are used
  options(nfds.loaded = TRUE)
  invisible(TRUE)
}

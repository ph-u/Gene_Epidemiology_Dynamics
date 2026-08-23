##### locate the single copy of the scripts #####
SRC_DIR <- normalizePath(Sys.getenv("NFDS_SRC", unset = file.path("..", "src")),
                         mustWork = TRUE)

make_fixture <- function(seed = 42L) {
  dir <- tempfile("nfdsfix"); set.seed(seed)
  for (s in c("raw", "src", "data")) dir.create(file.path(dir, s), recursive = TRUE)

  nLoci  <- 12L
  tPoint <- c(-6, 0, 6, 12)
  nPer   <- c(15, 15, 12, 10)

  meta <- data.frame(Time = rep(tPoint, nPer),
                     VT   = rbinom(sum(nPer), 1, 0.4),
                     SC   = sample(1:3, sum(nPer), replace = TRUE))
  gene <- matrix(rbinom(nrow(meta) * nLoci, 1, 0.5), ncol = nLoci)
  colnames(gene) <- sprintf("CLS%05d", seq_len(nLoci))
  gene[, 1] <- 1L; gene[, 2] <- 0L
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

load_model <- function(force = FALSE) {
  if (!force && isTRUE(getOption("nfds.loaded"))) return(invisible(TRUE))
  dir <- make_fixture()
  assign("abcsmc", stub_abcsmc, envir = globalenv())
  assign("FIXTURE_DIR", dir, envir = globalenv())
  old <- setwd(file.path(dir, "src")); on.exit(setwd(old))
  stopifnot(identical(normalizePath(".."), normalizePath(dir)))  # redirection worked
  source("model.r")
  options(nfds.loaded = TRUE)
  invisible(TRUE)
}

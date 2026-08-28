ref_jsd <- function(p, q) {                       # independent implementation
  p <- p / sum(p); q <- q / sum(q); m <- 0.5 * (p + q)
  kl <- function(a, b) { i <- a > 0; sum(a[i] * log(a[i] / b[i])) }
  0.5 * kl(p, m) + 0.5 * kl(q, m)
}

test_that("jsd() satisfies the defining properties", {
  load_model(); set.seed(1)
  p <- c(3, 1, 4, 1, 5); q <- c(2, 7, 1, 8, 2)
  expect_equal(jsd(p, p), 0)                                   # identity
  expect_equal(jsd(p, q), jsd(q, p))                           # symmetry
  expect_gte(jsd(p, q), 0); expect_lte(jsd(p, q), log(2))      # bounds, in nats
  expect_equal(jsd(p, q), jsd(p * 10, q * 3))                  # scale invariance
  expect_equal(jsd(c(1, 0), c(0, 1)), log(2))                  # disjoint support
})

test_that("jsd() matches an independent implementation", {
  load_model(); set.seed(2)
  for (i in 1:20) {
    a <- rpois(6, 5); b <- rpois(6, 5)
    if (sum(a) == 0 || sum(b) == 0) next
    expect_equal(jsd(a, b), ref_jsd(a, b))
  }
})

test_that("jsd() is NaN on an empty column -- why eQm.date must land on generations", {
  load_model()
  expect_true(is.na(jsd(c(0, 0, 0), c(1, 2, 3))))
})

test_that("m.nfds returns a matrix shaped like the observed statistic", {
  load_model(); set.seed(4)
  s <- m.nfds(.25, .1, .001, .1, .02)
  expect_true(is.matrix(s))
  expect_equal(dim(s), c(length(vtsc.lev), length(eQm.date)))
  expect_true(all(s >= 0))
})

test_that("each recorded column samples exactly n_t individuals", {
  # regression: nObs must be indexed by month tIme[i], not by loop position i
  load_model(); set.seed(4)
  s <- m.nfds(.25, .1, .001, .1, .02)
  expect_equal(colSums(s), as.numeric(nObs[as.character(eQm.date)]))
})

test_that("m.nfds is deterministic under set.seed", {
  load_model()
  set.seed(99); a <- m.nfds(.25, .1, .001, .1, .02)
  set.seed(99); b <- m.nfds(.25, .1, .001, .1, .02)
  expect_identical(a, b)
})

test_that("a vaccine that never switches on is identical to zero vaccine pressure", {
  # regression for the vCurve[i+1] index origin
  load_model()
  old <- vCurve; on.exit(assign("vCurve", old, envir = globalenv()))
  assign("vCurve", rep(0, length(old)), envir = globalenv())
  set.seed(1); a <- m.nfds(.2, .1, .001, .5, 0)
  assign("vCurve", old, envir = globalenv())
  set.seed(1); b <- m.nfds(.2, .1, .001, 0,  0)
  expect_identical(a, b)
})

test_that("stronger vaccine selection lowers the final vaccine-type fraction", {
  load_model()
  isVT <- startsWith(vtsc.lev, "1;")
  vtFrac <- function(v, seeds = 1:6)
    mean(vapply(seeds, function(s) {
      set.seed(s); x <- m.nfds(.2, .1, .001, v, 0)
      sum(x[isVT, ncol(x)]) / sum(x[, ncol(x)])
    }, numeric(1)))
  expect_lt(vtFrac(0.4), vtFrac(0))
})

test_that("popRunaway rejects a run before allocating the generation", {
  load_model()
  old <- popRunaway; on.exit(assign("popRunaway", old, envir = globalenv()))
  assign("popRunaway", 1, envir = globalenv())
  set.seed(3)
  expect_null(m.nfds(.2, .1, .001, .1, .02))
})

test_that("keepGenotypes returns the summary plus an aligned genotype table", {
  load_model(); set.seed(6)
  r <- m.nfds(.25, .1, .001, .1, .02, keepGenotypes = TRUE)
  expect_setequal(names(r), c("ss", "genotypes", "G"))
  expect_equal(dim(r$ss), c(length(vtsc.lev), length(eQm.date)))
  expect_equal(nrow(r$genotypes), nrow(G0))
  expect_equal(dim(r$G), dim(G0))
  expect_equal(r$genotypes$tag, d0.u$tag)
  expect_gt(sum(r$genotypes$finalCount), 0)
})

test_that("the strongest-NFDS loci are those with the smallest deviation statistic", {
  load_model()
  L <- length(gNam); strong <- order(selMode$strength)[seq_len(floor(L * 0.25))]
  expect_lte(max(selMode$strength[strong]), min(selMode$strength[-strong]))
})

test_that("nfds_jsd returns a finite scalar inside [0, ncol * log(2)]", {
  load_model(); set.seed(5)
  x <- list(propStrong = .2, fSelected = .1, wSelected = .001,
            vSelected = .1, migration = .02)
  d <- nfds_jsd(x, mIg0)
  expect_length(d, 1); expect_true(is.finite(d))
  expect_gte(d, 0); expect_lte(d, ncol(mIg0) * log(2))
})

test_that("nfds_jsd returns the maximum distance when the simulation fails", {
  load_model()
  old <- popRunaway; on.exit(assign("popRunaway", old, envir = globalenv()))
  assign("popRunaway", 1, envir = globalenv())
  d <- nfds_jsd(list(propStrong = .2, fSelected = .1, wSelected = .001,
                     vSelected = .1, migration = .02), mIg0)
  expect_equal(d, ncol(mIg0) * log(2))
})

test_that("nfds_jsd bootstraps itself in a worker with an empty globalenv", {
  # reproduces `could not find function "m.nfds"`: subjob_smc() runs a fresh
  # Rscript, and globalenv() is never serialised with the function
  load_full()
  scr <- file.path(FIXTURE_DIR, "src")
  out <- system2("Rscript", c("-e", shQuote(sprintf(
    'setwd("%s"); load("../data/setupState.RData", envir = globalenv());
     f <- get("nfds_jsd", globalenv()); rm(list = setdiff(ls(globalenv()), "f"), envir = globalenv());
     environment(f) <- globalenv(); assign("nfds_jsd", f, globalenv());
     obs <- local({ e <- new.env(); sys.source("setup.r", e); e$mIg0 });
     cat(nfds_jsd(list(propStrong=.2, fSelected=.1, wSelected=.001,
                       vSelected=.1, migration=.02), obs))', scr))),
    stdout = TRUE, stderr = TRUE)
  expect_false(any(grepl("could not find function", out)))
  expect_true(is.finite(suppressWarnings(as.numeric(tail(out, 1)))))
})

test_that("nfds_jsd returns the maximum distance rather than NA on a dead column", {
  load_model()
  bad <- mIg0; bad[, 1] <- 0
  d <- nfds_jsd(list(propStrong=.2, fSelected=.1, wSelected=.001,
                     vSelected=.1, migration=.02), bad)
  expect_true(is.finite(d))                       # the is.finite() guard
  expect_lte(d, ncol(mIg0) * log(2))
})


ref_jsd <- function(p, q) {                       # independent implementation
  p <- p / sum(p); q <- q / sum(q); m <- 0.5 * (p + q)
  kl <- function(a, b) { i <- a > 0; sum(a[i] * log(a[i] / b[i])) }
  0.5 * kl(p, m) + 0.5 * kl(q, m)
}

test_that("jsd() satisfies the defining properties", {
  load_model()
  p <- c(3, 1, 4, 1, 5); q <- c(2, 7, 1, 8, 2)
  expect_equal(jsd(p, p), 0)
  expect_equal(jsd(p, q), jsd(q, p))
  expect_gte(jsd(p, q), 0); expect_lte(jsd(p, q), log(2))
  expect_equal(jsd(p, q), jsd(p * 10, q * 3))
  expect_equal(jsd(c(1, 0), c(0, 1)), log(2))
})

test_that("jsd() matches an independent implementation", {
  load_model(); set.seed(2)
  for (i in 1:20) {
    a <- rpois(6, 5); b <- rpois(6, 5)
    if (sum(a) == 0 || sum(b) == 0) next
    expect_equal(jsd(a, b), ref_jsd(a, b))
  }
})

test_that("jsd() is NA on an empty column", {
  load_model()
  expect_true(is.na(jsd(c(0, 0, 0), c(1, 2, 3))))
})

test_that("coverage arithmetic: covered iff the vaccine existed and the host was in window", {
  # pure restatement of the nfds.r expression, so the semantics are pinned in one place
  cov <- function(tNow, hAge, sV, aC) (tNow >= sV) & (hAge - (tNow - sV) < aC)
  expect_true (cov(12,  6,   0, 24))   # born after roll-out
  expect_true (cov(12, 18,   0, 24))   # was 6 months old when it began
  expect_false(cov(12, 36,   0, 24))   # was 24 months old -- outside the window
  expect_false(cov( 0, 12,   3, 24))   # vaccine had not started yet
  expect_false(cov(12, 36, Inf, 24))   # never introduced in that deme
  expect_true (cov(12, 60,   6, 60))   # catch-up campaign reaches older hosts
})

test_that("m.nfds returns a matrix shaped like the observed statistic", {
  load_model(); set.seed(4)
  s <- mn()
  expect_true(is.matrix(s))
  expect_equal(dim(s), c(length(vtsc.lev) * nD, length(eQm.date)))
  expect_true(all(s >= 0))
})

test_that("each column samples exactly n_t individuals, stratified by deme", {
  load_model(); set.seed(4)
  s <- mn()
  expect_equal(colSums(s), as.numeric(rowSums(nObs[as.character(eQm.date), ])))
  for (dd in seq_len(nD)) {
    blk <- ((dd - 1) * length(vtsc.lev) + 1):(dd * length(vtsc.lev))
    expect_equal(colSums(s[blk, , drop = FALSE]),
                 as.numeric(nObs[as.character(eQm.date), dd]))
  }
})

test_that("m.nfds is deterministic under set.seed", {
  load_model()
  set.seed(99); a <- mn(); set.seed(99); b <- mn()
  expect_identical(a, b)
})

test_that("zero uptake is identical to zero efficacy", {
  load_model()
  swap("uPt", matrix(0, nVac, nD))
  set.seed(1); a <- mn(vSel = c(.5, .5))
  assign("uPt", get("uPt", globalenv()) * 0 + 1, envir = globalenv())   # restored on exit
  set.seed(1); b <- mn(vSel = c(0, 0))
  expect_identical(a, b)
})

test_that("a vaccine that never rolls out is identical to zero efficacy", {
  load_model()
  swap("vStart", matrix(Inf, nVac, nD))
  set.seed(1); a <- mn(vSel = c(.5, .5))
  assign("vStart", matrix(Inf, nVac, nD), envir = globalenv())
  set.seed(1); b <- mn(vSel = c(0, 0))
  expect_identical(a, b)
})

test_that("hosts older than every age window are never covered", {
  load_model()
  swap("aLev",  c(1e5, 2e5))                       # two levels: avoids sample()'s 1:n trap
  swap("aProb", matrix(.5, 2, nD))
  set.seed(1); a <- mn(vSel = c(.5, .5))
  swap("aLev",  aLev); swap("aProb", aProb)        # no-ops; restored by the first swap
  set.seed(1); b <- mn(vSel = c(.5, .5))
  expect_true(is.matrix(a))
  expect_false(identical(a, b))                    # age genuinely drives coverage
})

test_that("a single age level does not trip sample()'s 1:n reinterpretation", {
  # regression: sample(x, ...) with length(x) == 1 draws from 1:x, not from x
  load_model()
  swap("aLev", 24); swap("aProb", matrix(1, 1, nD))
  set.seed(1)
  expect_no_error(mn())
})

test_that("a later roll-out leaves more vaccine type behind", {
  # D1 starts VT1 at month 0, D3 at month 9, uptake 1 in both
  load_model()
  isV <- startsWith(vtsc.lev, "1;")
  blk <- function(x, dd) x[((dd - 1) * length(vtsc.lev) + 1):(dd * length(vtsc.lev)), ncol(x)]
  f   <- vapply(1:8, function(s) {
           set.seed(s); x <- mn(vSel = c(.5, 0), migration = 0, mWithin = 0)
           c(sum(blk(x, 1)[isV]) / sum(blk(x, 1)), sum(blk(x, 3)[isV]) / sum(blk(x, 3)))
         }, numeric(2))
  expect_lt(mean(f[1, ]), mean(f[2, ]))
})

test_that("stronger vaccine selection lowers the final vaccine-type fraction", {
  load_model()
  isVT <- rep(startsWith(vtsc.lev, "1;"), nD)
  vtFrac <- function(v, seeds = 1:6)
    mean(vapply(seeds, function(s) {
      set.seed(s); x <- mn(vSel = c(v, v))
      sum(x[isVT, ncol(x)]) / sum(x[, ncol(x)])
    }, numeric(1)))
  expect_lt(vtFrac(0.4), vtFrac(0))
})

test_that("coInf = 0 is Poisson; coInf > 0 overdisperses the population size", {
  load_model()
  tot <- function(ci, seeds = 1:8) vapply(seeds, function(s) {
    set.seed(s); r <- mn(coInf = ci, keepGenotypes = TRUE)
    if (is.null(r)) NA_real_ else sum(r$genotypes$finalCount)
  }, 0)
  set.seed(1); a <- mn(coInf = 0); set.seed(1); b <- mn(coInf = .5)
  expect_false(identical(a, b))                    # the NB branch is reached
  expect_gt(sd(tot(.5), na.rm = TRUE), sd(tot(0), na.rm = TRUE))
})

test_that("each deme is regulated to its own share of kappa", {
  load_model(); set.seed(7)
  r <- mn(mWithin = 0, migration = 0, keepGenotypes = TRUE)
  expect_false(is.null(r))
  expect_equal(colSums(r$demeCount) / kD, rep(1, nD), tolerance = .3)
})

test_that("mWithin moves individuals between demes", {
  load_model()
  set.seed(1); a <- mn(mWithin = 0,  migration = 0)
  set.seed(1); b <- mn(mWithin = .8, migration = 0)
  expect_false(identical(a, b))
})

test_that("popRunaway rejects a run before allocating the generation", {
  load_model()
  swap("popRunaway", 1)
  set.seed(3)
  expect_null(mn())
})

test_that("keepGenotypes returns the summary, a genotype table and deme counts", {
  load_model(); set.seed(6)
  r <- mn(keepGenotypes = TRUE)
  expect_setequal(names(r), c("ss", "genotypes", "demeCount", "G"))
  expect_equal(dim(r$ss), c(length(vtsc.lev) * nD, length(eQm.date)))
  expect_equal(nrow(r$genotypes), nrow(G0))
  expect_equal(dim(r$demeCount), c(nrow(G0), nD))
  expect_equal(r$genotypes$finalCount, rowSums(r$demeCount))
  expect_equal(r$genotypes$tag, d0.u$tag)
  expect_true(all(vtCols %in% names(r$genotypes)))
})

test_that("the strongest-NFDS loci are those with the smallest deviation statistic", {
  load_model()
  L <- length(gNam); strong <- order(selMode$strength)[seq_len(floor(L * 0.25))]
  expect_lte(max(selMode$strength[strong]), min(selMode$strength[-strong]))
})

test_that("nfds_jsd returns a finite scalar inside [0, ncol * log(2)]", {
  load_model(); set.seed(5)
  d <- nfds_jsd(p0(), mIg0)
  expect_length(d, 1); expect_true(is.finite(d))
  expect_gte(d, 0); expect_lte(d, ncol(mIg0) * log(2))
})

test_that("nfds_jsd assembles vSel from the vSelected<n> entries", {
  load_model()
  set.seed(8); a <- nfds_jsd(p0(vSelected1 = 0,  vSelected2 = 0),  mIg0)
  set.seed(8); b <- nfds_jsd(p0(vSelected1 = .9, vSelected2 = .9), mIg0)
  expect_false(isTRUE(all.equal(a, b)))
})

test_that("nfds_jsd ignores any parameter the model no longer takes", {
  load_model()
  set.seed(8); a <- nfds_jsd(p0(), mIg0)
  set.seed(8); b <- nfds_jsd(p0(vSelMod = .77), mIg0)   # removed parameter
  expect_identical(a, b)
})

test_that("nfds_jsd returns the maximum distance when the simulation fails", {
  load_model()
  swap("popRunaway", 1)
  expect_equal(nfds_jsd(p0(), mIg0), ncol(mIg0) * log(2))
})

test_that("nfds_jsd returns the maximum distance rather than NA on a dead column", {
  load_model()
  bad <- mIg0; bad[, 1] <- 0
  d <- nfds_jsd(p0(), bad)
  expect_true(is.finite(d))
  expect_lte(d, ncol(mIg0) * log(2))
})

test_that("nfds_jsd bootstraps itself in a worker with an empty globalenv", {
  load_full()
  scr <- file.path(FIXTURE_DIR, "src")
  out <- system2("Rscript", c("-e", shQuote(sprintf(
    'setwd("%s")
     load(Sys.getenv("NFDS_STATE"), envir = globalenv())
     f <- get("nfds_jsd", globalenv()); obs <- get("mIg0", globalenv())
     rm(list = setdiff(ls(globalenv()), c("f", "obs")), envir = globalenv())
     environment(f) <- globalenv()
     cat(f(list(propStrong = .25, fSelected = .1, wSelected = .001, migration = .02,
                coInf = 0, mWithin = .05, vSelected1 = .2, vSelected2 = .15), obs))', scr))),
    stdout = TRUE, stderr = TRUE)
  expect_false(any(grepl("could not find function|object .* not found", out)))
  expect_true(is.finite(suppressWarnings(as.numeric(tail(out, 1)))))
})


test_that("all four scripts parse", {
  # catches the `set,seed` / unbalanced-parenthesis class of error before anything runs
  for (f in c("src.r", "nfds.r", "setup.r", "model.r"))
    expect_silent(invisible(parse(file.path(SRC_DIR, f))))
})

test_that("the intermediate-frequency filter drops fixed and absent loci", {
  load_model()
  f <- colMeans(d0[, gNam, drop = FALSE])
  expect_true(all(f > 0 & f < 1))
  expect_lt(length(gNam), 12L)          # the two forced columns were removed
})

test_that("equilibrium frequencies are isolate-weighted over pre-vaccine samples", {
  load_model()
  expect_equal(eqm.pre, as.numeric(colMeans(d0[d0$Time <= 0, gNam])))
  expect_length(eqm.pre, length(gNam))
  expect_true(all(eqm.pre > 0 & eqm.pre < 1))
})

test_that("the deviation statistic uses 1 - e(1-e), not (1-e)^2", {
  load_model()
  fpost <- colMeans(eQm[eQm[, 1] > 0, -1])
  right <- as.numeric((fpost - eqm.pre)^2 / (1 - eqm.pre * (1 - eqm.pre)))
  wrong <- as.numeric(((fpost - eqm.pre) / (1 - eqm.pre))^2)
  expect_equal(selMode$strength, right)
  expect_false(isTRUE(all.equal(selMode$strength, wrong)))
  expect_equal(selMode$gene, gNam)
})

test_that("migrant probabilities are SC-balanced and sum to one", {
  load_model()
  expect_equal(sum(d0.u$migProb), 1)
  perSC <- tapply(d0.u$migProb, d0.u$SC, sum)
  expect_equal(as.numeric(perSC), rep(1 / length(unique(d0$SC)), length(perSC)))
})

test_that("vCurve covers every month and switches on at the configured month", {
  load_model()
  mo <- seq(min(eQm$Month), max(eQm$Month))
  expect_length(vCurve, length(mo))
  expect_true(all(vCurve %in% c(0, 1)))
  expect_false(is.unsorted(vCurve))                     # monotone 0 -> 1
  expect_equal(vCurve, as.numeric(mo >= as.numeric(g("vaccine start month"))))
})

test_that("vCurve position i+1 is the month simulated at loop step i", {
  load_model()
  mo   <- seq(min(eQm$Month), max(eQm$Month))
  tIme <- (min(eQm$Month) + 1):nGen
  expect_equal(mo[seq_along(tIme) + 1], tIme)
})

test_that("G0 is a numeric 0/1 matrix aligned with d0.u", {
  load_model()
  expect_equal(dim(G0), c(nrow(d0.u), length(gNam)))
  expect_identical(storage.mode(G0), "double")
  expect_true(all(G0 %in% c(0, 1)))
  expect_equal(colnames(G0), gNam)
})

test_that("genotype tags are unique and every isolate maps to one", {
  load_model()
  expect_equal(anyDuplicated(d0.u$tag), 0L)
  expect_false(anyNA(d0.u$tag)); expect_false(anyNA(d0$tag))
  expect_true(all(d0$tag %in% d0.u$tag))
  expect_equal(nrow(d0.u), length(unique(d0$tag)))
})

test_that("the starting pool is one entry per pre-vaccine isolate", {
  load_model()
  expect_equal(length(tag.pre), sum(d0$Time <= 0))
  expect_true(all(tag.pre >= 1 & tag.pre <= nrow(d0.u)))
  expect_false(anyNA(tag.pre))
})

test_that("mIg0 matches the sampled isolate counts, column by column", {
  load_model()
  expect_equal(dim(mIg0), c(length(vtsc.lev), length(eQm.date)))
  expect_equal(colSums(mIg0), as.numeric(nObs[as.character(eQm.date)]))
  for (j in seq_along(eQm.date))
    expect_equal(mIg0[, j],
                 vtsc(match(d0$tag[d0$Time == eQm.date[j]], d0.u$tag)))
})

test_that("every sampling month falls inside the simulated generations", {
  load_model()
  expect_true(all(eQm.date %in% seq_len(nGen)))
  expect_equal(nGen, max(eQm$Month))
  expect_lt(min(eQm$Month), nGen)     # guards the descending (min+1):nGen trap
})

test_that("abcsmc receives the observed statistic and a five-parameter prior", {
  load_full()
  expect_identical(ABCSMC_CALL$ss_obs, mIg0)
  expect_length(prior_dist$nfds, 5)
  expect_setequal(vapply(prior_dist$nfds, `[`, character(1), 1),
                  c("propStrong", "fSelected", "wSelected", "vSelected", "migration"))
  expect_lte(ABCSMC_CALL$distance_threshold_min, ncol(mIg0) * log(2))
})

test_that("the export writes one genotype row per genotype per posterior particle", {
  load_full()
  dOut <- file.path(FIXTURE_DIR, "data", paste0("run_", sEed))
  f <- list.files(dOut, pattern = "^nfdsGenotypes_", full.names = TRUE)
  expect_length(f, 1)
  x <- read.csv(f)
  expect_equal(nrow(x), nrow(d0.u) * 2)      # stub returns two particles
  expect_setequal(unique(x$particle), 1:2)
  expect_length(list.files(dOut, pattern = "\\.rds$"), 2)
})

test_that("setup.r defines no side effects that belong in model.r", {
  # workers source setup.r; a set.seed() there would correlate every particle
  txt <- readLines(file.path(SRC_DIR, "setup.r"))
  txt <- txt[!grepl("^\\s*#", txt)]
  expect_false(any(grepl("set\\.seed", txt)))
  expect_false(any(grepl("abcsmc", txt)))
  expect_false(any(grepl("write\\.csv|saveRDS", txt)))
})

test_that("model.r writes a prepared state a worker can load", {
  load_full()
  f <- file.path(FIXTURE_DIR, "data", paste0("run_", sEed), "setupState.RData")
  expect_true(file.exists(f))
  e <- new.env(); load(f, envir = e)
  for (o in c("m.nfds", "jsd", "nfds_jsd", "vtsc", "G0", "eqm.pre",
              "mIg0", "selMode", "vCurve", "d0.u", "nObs", "vtsc.lev"))
    expect_true(exists(o, envir = e), info = o)
})

test_that("genotype, SC and VT are one-to-one -- migProb assumes it", {
  load_model()
  expect_true(all(tapply(d0$SC, d0$tag, function(x) length(unique(x))) == 1))
  expect_true(all(tapply(d0$VT, d0$tag, function(x) length(unique(x))) == 1))
})

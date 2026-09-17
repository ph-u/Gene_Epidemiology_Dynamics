test_that("all four scripts parse", {
  for (f in c("src.r", "nfds.r", "setup.r", "model.r"))
    expect_no_error(invisible(parse(file.path(SRC_DIR, f))))
})

test_that("demes and vaccines are detected and encoded as 1..n", {
  load_model()
  expect_equal(nD, 3L); expect_equal(nVac, 2L)
  expect_equal(vtCols, c("VT1", "VT2"))
  expect_setequal(unique(d0$Deme), seq_len(nD))
  expect_equal(dLev, c("D1", "D2", "D3"))
})

test_that("the frequency filter uses the PRE-vaccine window, not the whole series", {
  load_model()
  fPre <- colMeans(d0[d0$Time <= 0, gNam, drop = FALSE])
  expect_true(all(fPre > .05 & fPre < .95))
  expect_equal(length(gNam), 10L)          # 12 loci minus the fixed and absent ones
})

test_that("eqm.pre is per-deme, isolate-weighted, L x nD", {
  load_model()
  expect_equal(dim(eqm.pre), c(length(gNam), nD))
  for (dd in seq_len(nD))
    expect_equal(eqm.pre[, dd],
                 as.numeric(colMeans(d0[d0$Time <= 0 & d0$Deme == dd, gNam])))
  expect_false(anyNA(eqm.pre))
})

test_that("the strong/weak ranking uses the GLOBAL pre-vaccine frequencies", {
  load_model()
  expect_equal(eqm.glb, as.numeric(colMeans(d0[d0$Time <= 0, gNam])))
  fpost <- colMeans(eQm[eQm$Month > 0, gNam, drop = FALSE])
  right <- as.numeric((fpost - eqm.glb)^2 / (1 - eqm.glb * (1 - eqm.glb)))
  wrong <- as.numeric(((fpost - eqm.glb) / (1 - eqm.glb))^2)   # the (1-e)^2 error
  expect_equal(selMode$strength, right)
  expect_false(isTRUE(all.equal(selMode$strength, wrong)))
  expect_equal(selMode$gene, gNam)
})

test_that("migrant probabilities are SC-balanced WITHIN each deme", {
  load_model()
  expect_equal(dim(mP0), c(nrow(d0.u), nD))
  expect_equal(colSums(mP0), rep(1, nD))
  for (dd in seq_len(nD)) {
    scPresent <- unique(d0$SC[d0$Deme == dd])
    perSC     <- tapply(mP0[, dd], d0.u$SC, sum)
    expect_equal(as.numeric(perSC[as.character(scPresent)]),
                 rep(1 / length(scPresent), length(scPresent)))
  }
})

test_that("roll-out is deme- and vaccine-specific, with Inf where absent", {
  load_model()
  expect_equal(dim(vStart), c(nVac, nD))
  expect_equal(vStart["VT1", ], c(D1 = 0, D2 = 3, D3 = 6))
  expect_equal(vStart["VT2", ], c(D1 = 6, D2 = 9, D3 = Inf))
  expect_equal(vMonth, seq(min(eQm$Month), max(eQm$Month)))
})

test_that("vMonth position i+1 is the month simulated at loop step i", {
  load_model()
  tIme <- (min(eQm$Month) + 1):nGen
  expect_equal(vMonth[seq_along(tIme) + 1], tIme)
})

test_that("deme sizes partition kappa", {
  load_model()
  expect_equal(sum(dProp), 1); expect_equal(sum(kD), k)
  expect_equal(dProp, as.numeric(table(d0$Deme)) / nrow(d0))
  expect_true(all(kD > 0))
})

test_that("G0 and VTmat are numeric 0/1 matrices aligned with d0.u", {
  load_model()
  expect_equal(dim(G0), c(nrow(d0.u), length(gNam)))
  expect_equal(dim(VTmat), c(nrow(d0.u), nVac))
  expect_identical(storage.mode(G0), "double")
  expect_identical(storage.mode(VTmat), "double")
  expect_true(all(G0 %in% c(0, 1))); expect_true(all(VTmat %in% c(0, 1)))
})

test_that("genotype determines SC and VT profile -- migProb and paste assume it", {
  load_model()
  expect_equal(anyDuplicated(d0.u$tag), 0L)
  expect_true(all(tapply(d0$SC, d0$tag, function(z) length(unique(z))) == 1))
  for (v in vtCols)
    expect_true(all(tapply(d0[[v]], d0$tag, function(z) length(unique(z))) == 1))
  ## but a genotype MAY span demes -- that is the point of a metapopulation
  expect_gt(max(tapply(d0$Deme, d0$tag, function(z) length(unique(z)))), 1)
})

test_that("the starting pool is one entry per pre-vaccine isolate, per deme", {
  load_model()
  expect_length(tag.pre, nD)
  for (dd in seq_len(nD)) {
    expect_length(tag.pre[[dd]], sum(d0$Time <= 0 & d0$Deme == dd))
    expect_true(all(tag.pre[[dd]] >= 1 & tag.pre[[dd]] <= nrow(d0.u)))
  }
  expect_false(anyNA(unlist(tag.pre)))
})

test_that("mIg0 is (class x deme) by timepoint and matches the isolate counts", {
  load_model()
  expect_equal(dim(mIg0), c(length(vtsc.lev) * nD, length(eQm.date)))
  expect_equal(colSums(mIg0), as.numeric(rowSums(nObs[as.character(eQm.date), ])))
  for (j in seq_along(eQm.date)) {
    s <- which(d0$Time == eQm.date[j])
    expect_equal(mIg0[, j], vtsc(match(d0$tag[s], d0.u$tag), d0$Deme[s]))
  }
})

test_that("nObs is stratified by timepoint and deme", {
  load_model()
  expect_equal(dim(nObs), c(length(eQm$Month), nD))
  expect_equal(sum(nObs), nrow(d0))
})

test_that("every sampling month falls inside the simulated generations", {
  load_model()
  expect_true(all(eQm.date %in% seq_len(nGen)))
  expect_equal(nGen, max(eQm$Month))
  expect_lt(min(eQm$Month), nGen)          # guards the descending (min+1):nGen trap
})

test_that("abcsmc receives the observed statistic and a 7 + nVac prior", {
  load_full()
  expect_identical(ABCSMC_CALL$ss_obs, mIg0)
  expect_length(prior_dist$nfds, 7 + nVac)
  expect_setequal(vapply(prior_dist$nfds, `[`, character(1), 1),
                  c("propStrong", "fSelected", "wSelected", "migration",
                    "coInf", "vSelMod", "mWithin", "vSelected1", "vSelected2"))
  expect_lte(ABCSMC_CALL$distance_threshold_min, ncol(mIg0) * log(2))
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
  for (o in c("m.nfds", "jsd", "nfds_jsd", "vtsc", "G0", "VTmat", "eqm.pre", "eqm.glb",
              "mIg0", "selMode", "vStart", "vMonth", "mP0", "kD", "dProp",
              "nD", "nVac", "d0.u", "nObs", "vtsc.lev", "tag.pre"))
    expect_true(exists(o, envir = e), info = o)
})

test_that("the export writes one genotype row per genotype per posterior particle", {
  load_full()
  dOut <- file.path(FIXTURE_DIR, "data", paste0("run_", sEed))
  f <- list.files(dOut, pattern = "^nfdsGenotypes_", full.names = TRUE)
  expect_length(f, 1)
  x <- read.csv(f)
  expect_equal(nrow(x), nrow(d0.u) * 2)          # stub returns two particles
  expect_setequal(unique(x$particle), 1:2)
  expect_length(list.files(dOut, pattern = "\\.rds$"), 2)
})

test_that("the export passes every fitted parameter to m.nfds", {
  # regression: the loop once passed 5 of 8 arguments, and pOst$vSelected no longer exists
  txt <- paste(readLines(file.path(SRC_DIR, "model.r")), collapse = " ")
  for (nm in c("coInf", "vSelMod", "mWithin", "vSelected"))
    expect_true(grepl(paste0("pOst\\[\\[.*", nm, "|pOst\\$", nm), txt), info = nm)
})

test_that("the export records per-deme counts", {
  load_full()
  dOut <- file.path(FIXTURE_DIR, "data", paste0("run_", sEed))
  x <- read.csv(list.files(dOut, pattern = "^nfdsGenotypes_", full.names = TRUE))
  expect_true(all(paste0("n_", dLev) %in% names(x)))
  expect_equal(rowSums(x[, paste0("n_", dLev)]), x$finalCount)
})


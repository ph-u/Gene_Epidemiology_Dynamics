test_that("g() retrieves values and yields character(0) for a missing key", {
  load_model()
  expect_equal(as.numeric(g("population")), 1200)
  expect_true(file.exists(g("rollout")))
  expect_length(g("no such key"), 0)
})

test_that("bcod() gives one unique, non-NA tag per row", {
  load_model()
  for (n in c(1L, 5L, 26L, 27L, 676L)) {
    tg <- bcod(data.frame(x = seq_len(n)))
    expect_length(tg, n); expect_false(anyNA(tg))
    expect_equal(length(unique(tg)), n)
  }
})

test_that("vtsc() tabulates into (VT profile x SC) x deme, class varying fastest", {
  load_model()
  idx <- c(1L, 1L, 2L, 3L); dm <- c(1L, 2L, 1L, 3L)
  v   <- vtsc(idx, dm)
  expect_length(v, length(vtsc.lev) * nD)
  expect_equal(sum(v), length(idx))
  expect_equal(v, as.vector(table(factor(vtsc.idx[idx], levels = seq_along(vtsc.lev)),
                                  factor(dm, levels = seq_len(nD)))))
  ## position of class c in deme d is (d-1)*nLev + c
  expect_equal(which(v > 0), sort(unique((dm - 1) * length(vtsc.lev) + vtsc.idx[idx])))
})

test_that("vtsc() returns all zeros for an empty sample", {
  load_model()
  expect_equal(vtsc(integer(0), integer(0)), rep(0L, length(vtsc.lev) * nD))
})

test_that("vtsc.idx and vtsc.lev stay consistent with d0.u", {
  load_model()
  expect_equal(vtsc.lev, sort(unique(d0.u$paste)))
  expect_equal(vtsc.idx, match(d0.u$paste, vtsc.lev))
  expect_false(anyNA(vtsc.idx))
})

test_that("each seed gets its own output directory", {
  load_full()
  dOut <- file.path(FIXTURE_DIR, "data", paste0("run_", sEed))
  expect_true(dir.exists(file.path(dOut, "tmp")))
  expect_identical(normalizePath(Sys.getenv("NFDS_STATE")),
                   normalizePath(file.path(dOut, "setupState.RData")))
})


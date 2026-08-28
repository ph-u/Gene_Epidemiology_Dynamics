test_that("g() retrieves values and yields character(0) for a missing key", {
  load_model()
  expect_equal(as.numeric(g("population")), 1000)
  expect_equal(as.numeric(g("vaccine start month")), 0)
  expect_length(g("no such key"), 0)   # the failure mode model.r's stopifnot guards
})

test_that("bcod() gives one unique, non-NA tag per row", {
  load_model()
  for (n in c(1L, 5L, 26L, 27L, 676L)) {
    tg <- bcod(data.frame(x = seq_len(n)))
    expect_length(tg, n)
    expect_false(anyNA(tg))
    expect_equal(length(unique(tg)), n)
  }
})

test_that("vtsc() tabulates individuals into VT;SC classes", {
  load_model()
  idx <- c(1L, 1L, 2L, 3L)
  expect_equal(vtsc(idx),
               as.vector(table(factor(d0.u$paste[idx], levels = vtsc.lev))))
  expect_equal(sum(vtsc(idx)), length(idx))
  expect_length(vtsc(idx), length(vtsc.lev))
})

test_that("vtsc() returns all zeros for an empty sample", {
  load_model()
  expect_equal(vtsc(integer(0)), rep(0L, length(vtsc.lev)))
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

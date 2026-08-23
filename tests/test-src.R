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

test_that("bcod() fails above 676 genotypes -- KNOWN LIMITATION", {
  # nC > 2 does cbind(b, rep(LETTERS, each = nrow(b))), where the vector is 26x
  # longer than b has rows. Delete this test once bcod() is fixed; the real
  # dataset has ~616 unique genotypes, uncomfortably close to the ceiling.
  load_model()
  expect_error(bcod(data.frame(x = seq_len(700))))
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

#The exported `col_w_*()` functions, and the shared input processing behind them
#(`process_mat1()`, `process_mat2()`, `.process_bin_vars()`).

m2 <- function() as.matrix(lalonde[c("age", "educ")])
m3 <- function() as.matrix(lalonde[c("age", "educ", "re74")])

test_that("col_w_mean() and col_w_sd() compute weighted moments", {
  m <- m2()

  #Unweighted values match the base equivalents.
  expect_equal(unname(col_w_mean(m)), unname(colMeans(m)))
  expect_equal(unname(col_w_sd(m)), unname(apply(m, 2L, sd)))

  #Weights shift the mean.
  expect_false(isTRUE(all.equal(col_w_mean(m), col_w_mean(m, weights = w_fixed))))

  #`s.weights` and `weights` multiply together.
  expect_equal(col_w_mean(m, weights = w_fixed, s.weights = sw_fixed),
               col_w_mean(m, weights = w_fixed * sw_fixed))

  #A single vector is treated as one column.
  expect_length(col_w_mean(lalonde$age), 1L)

  #`subset` restricts the units used.
  keep <- lalonde$age < 30
  expect_equal(unname(col_w_mean(m, subset = keep)),
               unname(colMeans(m[keep, , drop = FALSE])))

  #Names are carried through.
  expect_named(col_w_mean(m), c("age", "educ"))
})

test_that("col_w_smd() standardizes by the requested denominator", {
  m <- m2()
  tb <- lalonde$treat

  for (d in c("pooled", "treated", "control", "all", "hedges")) {
    v <- col_w_smd(m, treat = tb, s.d.denom = d)
    expect_length(v, 2L)
    expect_true(all(is.finite(v)))
  }

  #`std = FALSE` gives the raw difference in means.
  raw <- col_w_smd(m, treat = tb, std = FALSE)
  expect_equal(unname(raw),
               unname(colMeans(m[tb == 1, , drop = FALSE]) -
                        colMeans(m[tb == 0, , drop = FALSE])))

  #`std` may be given per column.
  mixed <- col_w_smd(m, treat = tb, std = c(TRUE, FALSE), s.d.denom = "pooled")
  expect_equal(mixed[[2L]], raw[[2L]])

  #`abs` folds to absolute values.
  expect_true(all(col_w_smd(m, treat = tb, abs = TRUE, s.d.denom = "pooled") >= 0))

  #A numeric `s.d.denom` is used directly.
  expect_equal(col_w_smd(m, treat = tb, s.d.denom = c(1, 1)), raw)
})

test_that("col_w_vr(), col_w_ks(), and col_w_ovl() return sensible ranges", {
  m <- m2()
  tb <- lalonde$treat

  expect_true(all(col_w_vr(m, treat = tb) > 0))

  ks <- col_w_ks(m, treat = tb)
  expect_true(all(ks >= 0 & ks <= 1))

  ovl <- col_w_ovl(m, treat = tb)
  expect_true(all(ovl >= 0 & ovl <= 1))

  #Weighted versions still lie in range.
  expect_true(all(col_w_ks(m, treat = tb, weights = w_fixed) <= 1))
  expect_true(all(col_w_ovl(m, treat = tb, weights = w_fixed) <= 1))

  #A perfectly separating covariate has an overlap complement of 1.
  sep <- matrix(ifelse(tb == 1, 100, -100), ncol = 1L,
                dimnames = list(NULL, "sep"))
  expect_equal(unname(col_w_ovl(sep, treat = tb)), 1)
})

test_that("col_w_ovl() honours `integrate`, `steps`, and `bw`", {
  m <- m2()
  tb <- lalonde$treat

  #Both integration methods approximate the same quantity.
  v_int <- col_w_ovl(m, treat = tb, integrate = TRUE)
  v_riem <- col_w_ovl(m, treat = tb, integrate = FALSE)
  expect_equal(v_int, v_riem, tolerance = 1e-2)

  #`steps` was once overwritten internally and so had no effect.
  v5 <- col_w_ovl(m, treat = tb, integrate = FALSE, steps = 5)
  v1001 <- col_w_ovl(m, treat = tb, integrate = FALSE, steps = 1001)
  expect_false(isTRUE(all.equal(v5, v1001)))

  #A finer grid must agree more closely with adaptive quadrature.
  expect_lt(max(abs(v1001 - v_int)), max(abs(v5 - v_int)))

  #`steps` is validated even when `integrate = TRUE`, because it is used if
  #`integrate()` fails and the Riemann sum is the fallback.
  expect_err(col_w_ovl(m, treat = tb, integrate = FALSE, steps = 2),
             "must be greater than or equal to 5")
  expect_err(col_w_ovl(m, treat = tb, integrate = TRUE, steps = 2),
             "must be greater than or equal to 5")

  #Alternative bandwidth rules are accepted.
  for (bw in c("nrd", "nrd0", "SJ")) {
    expect_true(all(is.finite(col_w_ovl(m, treat = tb, bw = bw))))
  }
  expect_err(col_w_ovl(m, treat = tb, bw = "bogus"),
             "is not an acceptable entry to `bw`")
})

test_that("col_w_cov() and col_w_corr() handle continuous treatments", {
  m <- m2()
  tc <- lalonde$re75

  cv <- col_w_cov(m, treat = tc)
  expect_length(cv, 2L)

  #`col_w_corr()` is `col_w_cov(std = TRUE)`.
  expect_equal(col_w_corr(m, treat = tc), col_w_cov(m, treat = tc, std = TRUE))

  #Correlations lie in [-1, 1].
  cr <- col_w_corr(m, treat = tc)
  expect_true(all(abs(cr) <= 1 + 1e-8))

  #Spearman ranks both sides.
  sp <- col_w_corr(m, treat = tc, type = "spearman")
  expect_true(all(abs(sp) <= 1 + 1e-8))
  expect_false(isTRUE(all.equal(cr, sp)))

  expect_true(all(col_w_cov(m, treat = tc, abs = TRUE) >= 0))
})

test_that("col_w_dcov() and col_w_dcorr() work for any number of columns", {
  #Regression test: `std` was validated but never recycled, so any matrix with
  #more than one column failed with "missing value where TRUE/FALSE needed".
  tc <- lalonde$re75

  for (m in list(m2()[, 1L, drop = FALSE], m2(), m3())) {
    v <- col_w_dcov(m, treat = tc)
    expect_length(v, ncol(m))
    expect_true(all(is.finite(v)))
    expect_true(all(v >= 0))

    r <- col_w_dcorr(m, treat = tc)
    expect_length(r, ncol(m))
    expect_true(all(is.finite(r)))
  }

  #`col_w_dcorr()` is `col_w_dcov(std = TRUE)`.
  expect_equal(col_w_dcorr(m2(), treat = tc),
               col_w_dcov(m2(), treat = tc, std = TRUE))

  #`std` may be given per column.
  expect_length(col_w_dcov(m2(), treat = tc, std = c(TRUE, FALSE)), 2L)

  #A data frame with a factor is split first, and `std` is recycled to the
  #post-split width.
  expect_length(col_w_dcov(lalonde[c("age", "race")], treat = tc), 4L)

  expect_err(col_w_dcov(m2(), treat = tc, std = c(TRUE, FALSE, TRUE)),
             "must have length equal to 1 or the number of columns")
})

test_that("data frames with factors are split automatically", {
  tb <- lalonde$treat
  d <- lalonde[c("age", "race")]

  v <- col_w_smd(d, treat = tb, s.d.denom = "pooled")
  expect_named(v, c("age", "race_black", "race_hispan", "race_white"))

  #The result matches the pre-split matrix.
  pre <- as.matrix(splitfactor(d, drop.first = "if2"))
  expect_equal(unname(v), unname(col_w_smd(pre, treat = tb, s.d.denom = "pooled")))
})

test_that("`bin.vars` accepts logical, numeric, and character forms", {
  m <- as.matrix(lalonde[c("age", "married")])

  expect_length(col_w_sd(m, bin.vars = c(FALSE, TRUE)), 2L)
  expect_length(col_w_sd(m, bin.vars = 2L), 2L)
  expect_length(col_w_sd(m, bin.vars = "married"), 2L)
  expect_length(col_w_sd(m, bin.vars = -1L), 2L)

  expect_err(col_w_sd(m, bin.vars = TRUE),
             "it must have length equal to the number of columns")
  expect_err(col_w_sd(m, bin.vars = c(-1L, 2L)),
             "Positive and negative indices cannot be mixed")
  expect_err(col_w_sd(m, bin.vars = 99L),
             "none of its values can exceed the number of columns")
  expect_err(col_w_sd(unname(m), bin.vars = "married"),
             "`mat` must have column names")
  expect_err(col_w_sd(m, bin.vars = "zzz"),
             "all its values must be column names")
})

test_that("`mat` must be a data frame or numeric matrix", {
  expect_err(col_w_mean(matrix("a", 2L, 2L)),
             "must be a data frame or numeric matrix")
  expect_err(col_w_mean(list(1, 2)),
             "must be a data frame or numeric matrix")

  #A 3-dimensional array is rejected rather than silently flattened.
  a <- array(1:8, c(2L, 2L, 2L))
  expect_err(col_w_mean(a), "must be a data frame or numeric matrix")
  expect_err(col_w_sd(a), "must be a data frame or numeric matrix")
  expect_err(col_w_dcov(a, treat = 1:2), "must be a data frame or numeric matrix")
})

test_that("degenerate weights are reported informatively", {
  m <- m2()
  tb <- lalonde$treat
  n <- nrow(m)
  zero <- rep(0, n)
  one_only <- c(1, rep(0, n - 1L))

  expect_err(col_w_mean(m, weights = zero),
             "at least one unit must have a nonzero weight")
  expect_err(col_w_sd(m, weights = one_only),
             "at least two units must have nonzero weights")
  expect_err(col_w_cov(m, treat = lalonde$re75, weights = one_only),
             "at least two units must have nonzero weights")
  expect_err(col_w_smd(m, treat = tb, weights = ifelse(tb == 1, 0, 1)),
             "at least one unit in each level of `treat` must have a nonzero weight")
  expect_err(col_w_vr(m, treat = tb, weights = ifelse(tb == 1, 0, 1)),
             "at least two units in each level of `treat` must have nonzero weights")
  expect_err(col_w_ks(m, treat = tb, weights = ifelse(tb == 1, 0, 1)),
             "at least one unit in each level of `treat` must have a nonzero weight")
  expect_err(col_w_ovl(m, treat = tb, weights = ifelse(tb == 1, 0, 1)),
             "at least one unit in each level of `treat` must have a nonzero weight")
})

test_that("statistics requiring a binary treatment reject other types", {
  m <- m2()

  for (f in list(col_w_smd, col_w_vr, col_w_ks, col_w_ovl)) {
    expect_err(f(m, treat = lalonde$race), "must be a binary variable")
  }
})

test_that("mismatched lengths and invalid `subset` are rejected", {
  m <- m2()
  tb <- lalonde$treat

  expect_err(col_w_smd(m, treat = tb[-1L]), "must have the same number of units")
  expect_err(col_w_mean(m, weights = w_fixed[-1L]),
             "must have the same number of units")
  expect_err(col_w_smd(m, treat = tb, s.d.denom = c(1, 2, 3)),
             "must have length")
  expect_err(col_w_smd(m, treat = tb, std = c(TRUE, FALSE, TRUE)),
             "must have length equal to 1 or the number of columns")
})

test_that("`na.rm` controls how missing values are treated", {
  m <- as.matrix(lalonde_mis[c("age", "re74")])

  #With `na.rm = TRUE` (the default) a value is still produced.
  expect_true(all(is.finite(col_w_mean(m, na.rm = TRUE))))

  #With `na.rm = FALSE` the affected column is NA.
  v <- col_w_mean(m, na.rm = FALSE)
  expect_true(anyNA(v))
})

#The exported `col_w_*()` functions, and the shared input processing behind them
#(`process_mat1()`, `process_mat2()`, `.process_bin_vars()`); plus
#`.bal_tab_col_spec()`, which lays out the columns of every balance table.

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

test_that("col_w_ovl() honors `integrate`, `steps`, and `bw`", {
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

# .bal_tab_col_spec() is the single source for balance-table columns. Three table
# builders create their columns from it, and `print()` and `love.plot()` select
# columns by name, so the layout it produces has to be recoverable from a
# `bal.tab` object's own `print.options` -- not just from the builders' locals.

spec_from_p.ops <- function(p, table = "Balance") {
  no.adj <- !isTRUE(p[["disp.adj"]])
  wn <- if (no.adj) "Adj" else p[["weight.names"]]
  thr.s <- if (no.adj) "Un" else wn

  if (table == "Balance") {
    return(cobalt:::.bal_tab_col_spec(p[["type"]], p[["compute"]], p[["thresholds"]],
                                      samples = c("Un", wn), threshold.samples = thr.s))
  }

  if (table == "Subclass.Balance") {
    return(cobalt:::.bal_tab_col_spec(p[["type"]], p[["compute"]], p[["thresholds"]],
                                      samples = "Adj"))
  }

  agg <- switch(table,
                "Balance.Across.Clusters" = p[["cluster.fun"]],
                "Balance.Across.Imputations" = p[["imp.fun"]],
                "max")

  cobalt:::.bal_tab_col_spec(p[["type"]], p[["compute"]], p[["thresholds"]],
                            samples = c("Un", wn), quantities = NULL,
                            agg.funs = if (isTRUE(p[["quick"]])) agg else c("min", "mean", "max"),
                            threshold.samples = if (length(agg) != 1L) character(0L) else thr.s,
                            threshold.agg.fun = agg,
                            include.times = table == "Balance.Across.Times")
}

expect_spec_matches <- function(b, table = "Balance", element = NULL) {
  p <- attr(b, "print.options")
  tbl <- if (is_null(element)) b[[table]] else b[[table]][[element]]

  expect_identical(names(tbl), spec_from_p.ops(p, table)[["name"]],
                   info = paste(table, element))
}

test_that(".bal_tab_col_spec() reproduces the columns of a leaf balance table", {
  covs <- lalonde[c("age", "educ", "married")]
  t <- lalonde$treat

  #No weights: the "Adj" block is still laid out, unfilled.
  expect_spec_matches(bal.tab(covs, treat = t, s.d.denom = "pooled"))
  expect_spec_matches(bal.tab(covs, treat = t, s.d.denom = "pooled",
                              thresholds = c(m = .1, ks = .1)))

  #One weight set: the threshold column loses its suffix.
  b <- bal.tab(covs, treat = t, s.d.denom = "pooled", weights = w_fixed,
               un = TRUE, thresholds = c(m = .1))
  expect_spec_matches(b)
  expect_true("M.Threshold" %in% names(b$Balance))

  #Two weight sets: it keeps one suffix per set, interleaved in each block.
  b2 <- bal.tab(covs, treat = t, s.d.denom = "pooled",
                weights = data.frame(W1 = w_fixed, W2 = rev(w_fixed)),
                un = TRUE, thresholds = c(m = .1))
  expect_spec_matches(b2)
  expect_identical(names(b2$Balance),
                   c("Type", "Diff.Un", "Diff.W1", "M.Threshold.W1",
                     "Diff.W2", "M.Threshold.W2"))

  #Moments, every statistic, and `quick = FALSE`.
  expect_spec_matches(bal.tab(covs, treat = t, s.d.denom = "pooled", weights = w_fixed,
                              disp = c("means", "sds"), un = TRUE))
  expect_spec_matches(bal.tab(covs, treat = t, s.d.denom = "pooled", weights = w_fixed,
                              stats = cobalt:::all_STATS("bin"), un = TRUE))
  expect_spec_matches(bal.tab(covs, treat = t, s.d.denom = "pooled", weights = w_fixed,
                              quick = FALSE, thresholds = c(m = .1)))

  #Continuous treatments have no group part, and the `.target` statistics are
  #adjusted-only, so they appear in no "Un" block.
  b_c <- bal.tab(covs, treat = lalonde$re75, weights = w_fixed, un = TRUE,
                 stats = cobalt:::all_STATS("cont"), disp = c("means", "sds"))
  expect_spec_matches(b_c)
  expect_false(any(grepl("Target.Un", names(b_c$Balance), fixed = TRUE)))
})

test_that(".bal_tab_col_spec() reproduces the columns of every summary table", {
  covs <- lalonde[c("age", "educ", "married")]
  t <- lalonde$treat

  #A threshold labels one aggregate, so it is dropped when several are displayed.
  b <- bal.tab(covs, treat = t, s.d.denom = "pooled", weights = w_fixed,
               cluster = cl_idx, un = TRUE, thresholds = c(m = .1),
               cluster.summary = TRUE, cluster.fun = "mean")
  expect_spec_matches(b, "Balance.Across.Clusters")
  expect_true("M.Threshold" %in% names(b$Balance.Across.Clusters))

  b3 <- bal.tab(covs, treat = t, s.d.denom = "pooled", weights = w_fixed,
                cluster = cl_idx, un = TRUE, thresholds = c(m = .1),
                cluster.summary = TRUE, cluster.fun = c("min", "mean", "max"))
  expect_spec_matches(b3, "Balance.Across.Clusters")
  expect_false(any(grepl("Threshold", names(b3$Balance.Across.Clusters), fixed = TRUE)))

  #`quick = FALSE` widens the aggregates shown but not the ones asked for, so the
  #threshold survives alongside all three.
  bq <- bal.tab(covs, treat = t, s.d.denom = "pooled", weights = w_fixed,
                cluster = cl_idx, thresholds = c(m = .1), cluster.summary = TRUE,
                cluster.fun = "mean", quick = FALSE)
  expect_spec_matches(bq, "Balance.Across.Clusters")
  expect_true(all(c("Min.Diff.Adj", "Mean.Diff.Adj", "Max.Diff.Adj", "M.Threshold") %in%
                    names(bq$Balance.Across.Clusters)))

  expect_spec_matches(bal.tab(covs, treat = t, s.d.denom = "pooled", weights = w_fixed,
                              imp = imp_idx, un = TRUE, thresholds = c(m = .1)),
                      "Balance.Across.Imputations")

  expect_spec_matches(bal.tab(covs, treat = lalonde$race, s.d.denom = "pooled",
                              weights = w_fixed, thresholds = c(m = .1)),
                      "Balance.Across.Pairs")

  #The msm summary carries a leading `Times` column.
  b_msm <- bal.tab(list(treat ~ age, treat ~ age + educ),
                   data = cbind(lalonde, list(w = w_fixed)), weights = "w",
                   s.d.denom = "pooled", un = TRUE, thresholds = c(m = .1))
  expect_spec_matches(b_msm, "Balance.Across.Times")
  expect_identical(names(b_msm$Balance.Across.Times)[1L], "Times")
})

test_that(".bal_tab_col_spec() reproduces the columns of subclass tables", {
  covs <- lalonde[c("age", "educ", "married")]
  sub <- factor(rep(1:4, length.out = nrow(lalonde)))

  for (q in c(TRUE, FALSE)) {
    b <- bal.tab(covs, treat = lalonde$treat, s.d.denom = "pooled", subclass = sub,
                 which.subclass = .all, subclass.summary = TRUE,
                 thresholds = c(m = .1, ks = .1), quick = q)

    for (i in names(b$Subclass.Balance)) {
      expect_spec_matches(b, "Subclass.Balance", i)
    }

    #The across-subclass table is an ordinary balance table, so its `compute` is
    #intersected with `stats` where the per-subclass blocks' is not (H3).
    expect_identical(names(b$Balance.Across.Subclass),
                     cobalt:::.bal_tab_col_spec(
                       "bin",
                       if (q) c(attr(b, "print.options")$disp, "mean.diffs", "ks.statistics")
                       else c("means", "sds", "mean.diffs", "ks.statistics"),
                       attr(b, "print.options")$thresholds,
                       samples = c("Un", "Adj"), threshold.samples = "Adj")[["name"]])
  }
})

test_that(".bal_tab_col_spec() labels each column it generates", {
  spec <- cobalt:::.bal_tab_col_spec(
    "bin", c("means", "sds", "mean.diffs", "variance.ratios"),
    thresholds = list(mean.diffs = .1),
    samples = c("Un", "W1", "W2"), threshold.samples = c("W1", "W2"))

  expect_named(spec, c("name", "quantity", "stat", "sample", "agg.fun", "group"))
  expect_false(anyDuplicated(spec$name) > 0L)

  #Every row is one of the six kinds of column, and the metadata identifies it
  #without anyone having to parse the name.
  expect_setequal(spec$quantity, c("type", "means", "sds", "stat", "threshold"))
  expect_identical(spec$name[spec$quantity == "means" & spec$sample == "W1"],
                   c("M.0.W1", "M.1.W1"))
  expect_identical(spec$group[spec$quantity == "means" & spec$sample == "W1"],
                   c("0", "1"))
  expect_identical(spec$name[spec$quantity == "threshold"],
                   c("M.Threshold.W1", "M.Threshold.W2"))
  expect_identical(unique(spec$stat[spec$quantity == "stat"]),
                   c("mean.diffs", "variance.ratios"))

  #The `.target` statistics are the adjusted-only ones, so they get no unadjusted
  #column while their non-target counterparts do.
  cont <- cobalt:::.bal_tab_col_spec(
    "cont", c("correlations", "mean.diffs.target"), samples = c("Un", "Adj"))

  expect_identical(cont$name,
                   c("Type", "Corr.Un", "Corr.Adj", "Diff.Target.Adj"))

  #Aggregated summaries prefix each statistic with each function, thresholds last.
  agg <- cobalt:::.bal_tab_col_spec(
    "bin", "mean.diffs", thresholds = list(mean.diffs = .1), samples = c("Un", "Adj"),
    quantities = NULL, agg.funs = c("min", "mean", "max"),
    threshold.samples = "Adj", threshold.agg.fun = "mean")

  expect_identical(agg$name,
                   c("Type", "Min.Diff.Un", "Mean.Diff.Un", "Max.Diff.Un",
                     "Min.Diff.Adj", "Mean.Diff.Adj", "Max.Diff.Adj", "M.Threshold"))
  expect_identical(agg$agg.fun[agg$name == "M.Threshold"], "mean")
  expect_identical(agg$agg.fun[agg$name == "Min.Diff.Un"], "min")
})

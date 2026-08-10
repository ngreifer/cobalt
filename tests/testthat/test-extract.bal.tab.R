#`as.data.frame()` and `format()` for `bal.tab` objects.

covs3 <- function() lalonde[c("age", "educ", "married")]

#The columns every tidy result has. The two threshold columns join them only when a
#threshold is on display; see the test that pins that.
tidy_cols <- c("variable", "type", "sample", "stat", "group", "estimate")

tidy_cols_thr <- c(tidy_cols, "threshold", "threshold.value")

test_that("as.data.frame() returns one row per covariate, sample, and statistic", {
  b <- bal.tab(covs3(), treat = lalonde$treat, s.d.denom = "pooled", weights = w_fixed,
               un = TRUE, stats = c("m", "ks"), thresholds = c(m = .1))

  d <- as.data.frame(b)

  expect_s3_class(d, "data.frame")
  expect_named(d, tidy_cols_thr)

  #3 covariates x 2 samples x 2 statistics.
  expect_identical(nrow(d), 12L)
  expect_setequal(d$variable, rownames(b$Balance))
  expect_setequal(d$sample, c("Unadjusted", "Adj"))
  expect_setequal(d$stat, c("mean.diffs", "ks.statistics"))

  #The values are the balance table's own, not recomputed.
  for (i in seq_len(nrow(d))) {
    col <- paste.(if (d$stat[i] == "mean.diffs") "Diff" else "KS",
                  if (d$sample[i] == "Unadjusted") "Un" else d$sample[i])
    expect_equal(d$estimate[i], b$Balance[[col]][match(d$variable[i], rownames(b$Balance))],
                 info = paste(d$variable[i], col))
  }

  #Only the statistic with a threshold carries a verdict, and only for the sample the
  #threshold column describes -- with weights present that is the adjusted one, which
  #is where `print()` shows it.
  expect_true(all(is.na(d$threshold[d$stat == "ks.statistics"])))
  expect_true(all(!is.na(d$threshold[d$stat == "mean.diffs" & d$sample == "Adj"])))
  expect_true(all(is.na(d$threshold[d$stat == "mean.diffs" & d$sample == "Unadjusted"])))

  #The numeric threshold belongs to the statistic, so it is reported on every row.
  expect_setequal(d$threshold.value[d$stat == "mean.diffs"], .1)
  expect_true(all(is.na(d$threshold.value[d$stat == "ks.statistics"])))

  expect_setequal(d$threshold[d$stat == "mean.diffs" & d$sample == "Adj"],
                  b$Balance$M.Threshold)

  #`type` comes from the table, and rows are grouped by covariate.
  expect_setequal(d$type, c("Contin.", "Binary"))
  expect_identical(unique(d$variable), rownames(b$Balance))
})

test_that("as.data.frame() carries the threshold columns only when there is a threshold", {
  # They say nothing at all when no threshold is on display, so they are not carried.
  # Whether they are is decided by the arguments, not by whether the values happen to come
  # out `NA`: a threshold on one statistic keeps them for every row, `NA` on the rows that
  # statistic does not cover.
  covs <- covs3()
  t <- lalonde$treat

  no.thr <- as.data.frame(bal.tab(covs, treat = t, s.d.denom = "pooled"))

  expect_named(no.thr, tidy_cols)

  b <- bal.tab(covs, treat = t, s.d.denom = "pooled", stats = c("m", "ks"),
               thresholds = c(m = .1), disp = "means")

  thr <- as.data.frame(b)

  expect_named(thr, tidy_cols_thr)

  #Kept for every row, filled in only where the threshold applies.
  expect_false(all(is.na(thr$threshold)))
  expect_true(all(is.na(thr$threshold[thr$stat != "mean.diffs"])))
  expect_true(all(thr$threshold.value[thr$stat == "mean.diffs"] == .1))

  #Suppressing the display of the only threshold drops them again, since `print()` would
  #not show it either.
  expect_named(as.data.frame(b, disp.thresholds = c(m = FALSE)), tidy_cols)
})

test_that("as.data.frame() reports means and SDs as statistics with a group", {
  b <- bal.tab(covs3(), treat = lalonde$treat, s.d.denom = "pooled", weights = w_fixed,
               un = TRUE, disp = c("means", "sds"))

  d <- as.data.frame(b)

  expect_setequal(d$stat, c("mean", "sd", "mean.diffs"))

  #A binary treatment has a mean per group, named the way `bal.tab()` names the group
  #rather than by the position its column takes; a mean difference belongs to neither.
  expect_setequal(d$group[d$stat == "mean"], c("Control", "Treated"))
  expect_true(all(is.na(d$group[d$stat == "mean.diffs"])))

  #Moments carry no threshold, whatever thresholds were requested.
  expect_true(all(is.na(d$threshold[d$stat %in% c("mean", "sd")])))
  expect_true(all(is.na(d$threshold.value[d$stat %in% c("mean", "sd")])))

  m <- d[d$stat == "mean" & d$sample == "Unadjusted" & d$group == "Treated", ]
  expect_equal(m$estimate[match(rownames(b$Balance), m$variable)], b$Balance$M.1.Un)

  #A continuous treatment has no groups, so its moments describe the whole sample.
  b_c <- bal.tab(covs3(), treat = lalonde$re75, weights = w_fixed, un = TRUE,
                 disp = c("means", "sds"))
  d_c <- as.data.frame(b_c)
  expect_true(all(d_c$group[d_c$stat %in% c("mean", "sd")] == "All"))
  expect_true(all(is.na(d_c$group[d_c$stat == "correlations"])))
  expect_setequal(d_c$stat, c("mean", "sd", "correlations"))
})

test_that("as.data.frame() reports a multi-category group's own moments once", {
  # Balance is computed one pair of groups at a time, but a group's own mean is the same
  # in every pair it appears in. Reported per pair it would say the same thing several
  # times over and imply it depended on the comparison, so it is reported once, named for
  # the group and belonging to no pair.
  covs <- lalonde["educ"]
  race <- lalonde$race

  for (pairwise in c(TRUE, FALSE)) {
    b <- bal.tab(covs, treat = race, s.d.denom = "pooled", weights = w_fixed, un = TRUE,
                 disp = c("means", "sds"), pairwise = pairwise)

    d <- as.data.frame(b, which.treat = .all)

    moments <- d[d$stat %in% c("mean", "sd"), , drop = FALSE]
    contrasts <- d[d$stat == "mean.diffs", , drop = FALSE]

    #Every group is named, and named once per variable, sample, and statistic.
    groups <- if (pairwise) levels(race) else c("All", levels(race))

    expect_setequal(moments$group, groups)
    expect_identical(anyDuplicated(moments[c("variable", "sample", "stat", "group")]), 0L)

    #A moment belongs to a group, a mean difference to a pair, and neither to both.
    expect_true(all(is.na(moments$pair)), info = pairwise)
    expect_false(anyNA(contrasts$pair), info = pairwise)
    expect_true(all(is.na(contrasts$group)), info = pairwise)

    #The values are the ones `bal.tab()` computed, and are right for the named group --
    #the thing the positional labels got wrong, since a pair's name and its `0`/`1`
    #columns run in opposite order.
    un <- moments[moments$stat == "mean" & moments$sample == "Unadjusted", ]

    expect_equal(un$estimate[match(levels(race), un$group)],
                 as.vector(tapply(lalonde$educ, race, mean)))

    if (!pairwise) {
      expect_equal(un$estimate[un$group == "All"], mean(lalonde$educ))
    }
  }
})

test_that("as.data.frame() puts each level of segmentation in its own column", {
  # Segmentation is always a column, never a nested list, so the result is one
  # rectangle whatever the shape of the input.
  covs <- covs3()
  t <- lalonde$treat
  sub <- factor(rep(1:4, length.out = nrow(lalonde)))

  shapes <- list(
    cluster = list(b = bal.tab(covs, treat = t, s.d.denom = "pooled", weights = w_fixed,
                               cluster = cl_idx, un = TRUE),
                   key = "cluster", levels = levels(cl_idx)),
    imp = list(b = bal.tab(covs, treat = t, s.d.denom = "pooled", weights = w_fixed,
                           imp = imp_idx, un = TRUE),
               key = "imp", levels = c("1", "2")),
    pair = list(b = bal.tab(covs, treat = lalonde$race, s.d.denom = "pooled",
                            weights = w_fixed),
                key = "pair", levels = NULL),
    subclass = list(b = bal.tab(covs, treat = t, s.d.denom = "pooled", subclass = sub),
                    key = "subclass", levels = levels(sub))
  )

  for (nm in names(shapes)) {
    d <- as.data.frame(shapes[[nm]]$b)

    expect_named(d, c(shapes[[nm]]$key, tidy_cols), info = nm)
    expect_false(anyNA(d[[shapes[[nm]]$key]]), info = nm)

    if (is_not_null(shapes[[nm]]$levels)) {
      expect_setequal(d[[shapes[[nm]]$key]], shapes[[nm]]$levels)
    }
  }

  #Longitudinal treatments segment by time point.
  b_msm <- bal.tab(list(treat ~ age, treat ~ age + educ),
                   data = cbind(lalonde, list(w = w_fixed)), weights = "w",
                   s.d.denom = "pooled", un = TRUE)
  expect_named(as.data.frame(b_msm), c("time", tidy_cols))
  expect_setequal(as.data.frame(b_msm)$time, c("1", "2"))

  #Two levels of segmentation give two columns.
  b_n <- bal.tab(covs, treat = t, s.d.denom = "pooled", subclass = sub, cluster = cl_idx)
  d_n <- as.data.frame(b_n)
  expect_named(d_n, c("cluster", "subclass", tidy_cols))
  expect_setequal(d_n$cluster, levels(cl_idx))
  expect_setequal(d_n$subclass, levels(sub))
  expect_identical(nrow(d_n), nlevels(cl_idx) * nlevels(sub) * 3L)
})

test_that("as.data.frame() honours print()'s display arguments", {
  b <- bal.tab(covs3(), treat = lalonde$treat, s.d.denom = "pooled", weights = w_fixed,
               un = TRUE, stats = c("m", "ks"), thresholds = c(m = .1), quick = FALSE)

  #`stats` narrows which statistics appear.
  expect_setequal(as.data.frame(b, stats = "ks.statistics")$stat, "ks.statistics")

  #`un` drops the unadjusted sample.
  expect_setequal(as.data.frame(b, un = FALSE)$sample, "Adj")

  #`disp` adds the moments.
  expect_true("sd" %in% as.data.frame(b, disp = "sds")$stat)

  #`disp.thresholds` drops the verdict but keeps the estimate.
  d <- as.data.frame(b, disp.thresholds = c(m = FALSE))
  expect_true(all(is.na(d$threshold)))
  expect_true("mean.diffs" %in% d$stat)

  #`imbalanced.only` drops the balanced covariates.
  d_i <- as.data.frame(b, imbalanced.only = TRUE)
  expect_lt(length(unique(d_i$variable)), 3L)
  expect_true(all(unique(d_i$variable) %in% rownames(b$Balance)))

  #A statistic that was never computed cannot be requested, as in print().
  b_q <- bal.tab(covs3(), treat = lalonde$treat, s.d.denom = "pooled", weights = w_fixed)
  expect_err(as.data.frame(b_q, stats = "ks.statistics"), "cannot contain")

  #`.all`/`.none` are accepted as they are by print().
  b_s <- bal.tab(covs3(), treat = lalonde$treat, s.d.denom = "pooled",
                 subclass = factor(rep(1:4, length.out = nrow(lalonde))))
  expect_no_error(as.data.frame(b_s, which.subclass = .all))
})

test_that("as.data.frame(wide = TRUE) gives print()'s layout", {
  b <- bal.tab(covs3(), treat = lalonde$treat, s.d.denom = "pooled", weights = w_fixed,
               un = TRUE, thresholds = c(m = .1))

  d <- as.data.frame(b, wide = TRUE)

  #One row per covariate, with the names promoted out of the row names.
  expect_identical(nrow(d), nrow(b$Balance))
  expect_identical(d$variable, rownames(b$Balance))
  expect_identical(names(d)[1:2], c("variable", "Type"))
  expect_true(all(c("Diff.Un", "Diff.Adj", "M.Threshold") %in% names(d)))
  expect_equal(d$Diff.Adj, b$Balance$Diff.Adj)

  #Segmentation columns are kept.
  b_c <- bal.tab(covs3(), treat = lalonde$treat, s.d.denom = "pooled",
                 weights = w_fixed, cluster = cl_idx)
  expect_identical(names(as.data.frame(b_c, wide = TRUE))[1L], "cluster")

  expect_err(as.data.frame(b, wide = "yes"), "must be")
})

test_that("format() returns exactly the table print() displays", {
  # The promise of format() is that it is print()'s own table, so it is checked by
  # rendering it and looking for that block in print()'s output.
  #Wide enough that print() does not wrap the table into several blocks.
  rlang::local_options(width = 400)

  covs <- covs3()
  t <- lalonde$treat
  mk <- function(...) bal.tab(covs, treat = t, s.d.denom = "pooled", ...)

  cases <- list(
    leaf = mk(weights = w_fixed, un = TRUE, thresholds = c(m = .1),
              disp = c("means", "sds")),
    two_weights = mk(weights = data.frame(W1 = w_fixed, W2 = rev(w_fixed)), un = TRUE,
                     thresholds = c(m = .1)),
    no_weights = mk(thresholds = c(m = .1)),
    all_stats = mk(weights = w_fixed, un = TRUE, stats = cobalt:::all_STATS("bin")),
    cont = bal.tab(covs, treat = lalonde$re75, weights = w_fixed, un = TRUE),
    cluster = mk(weights = w_fixed, cluster = cl_idx, un = TRUE, cluster.summary = TRUE,
                 cluster.fun = "mean", thresholds = c(m = .1)),
    multi = bal.tab(covs, treat = lalonde$race, s.d.denom = "pooled", weights = w_fixed,
                    thresholds = c(m = .1)),
    subclass = mk(subclass = factor(rep(1:4, length.out = nrow(lalonde))),
                  subclass.summary = TRUE, thresholds = c(m = .1))
  )

  for (nm in names(cases)) {
    f <- format(cases[[nm]])

    expect_s3_class(f, "data.frame")
    expect_true(all(vapply(f, is.character, logical(1L))), info = nm)

    rendered <- trimws(capture.output(print.data.frame(f)))
    printed <- trimws(capture.output(print(cases[[nm]])))

    hit <- which(printed == rendered[1L])
    expect_length(hit, 1L)
    expect_identical(printed[hit + seq_along(rendered) - 1L], rendered, info = nm)
  }
})

test_that("format() renders NA as a period and pads the decimals", {
  b <- bal.tab(covs3(), treat = lalonde$treat, s.d.denom = "pooled", weights = w_fixed,
               un = TRUE, disp = "sds")

  f <- format(b)

  #A binary covariate has no SD to report under `binary = "raw"`.
  expect_identical(f["married", "SD.0.Un"], ".")

  #Everything in a column has the same number of decimal places.
  dec <- function(x) {
    x <- trimws(x[x != "."])
    nchar(sub("^[^.]*\\.?", "", x))
  }
  expect_length(unique(dec(f[["Diff.Un"]])), 1L)

  #`digits` is honoured.
  expect_false(identical(format(b, digits = 2L)[["Diff.Un"]], f[["Diff.Un"]]))
})

test_that("format(component = 'observations') gives the sample size table", {
  b <- bal.tab(covs3(), treat = lalonde$treat, s.d.denom = "pooled", weights = w_fixed)

  f <- format(b, component = "observations")

  expect_s3_class(f, "data.frame")
  expect_named(f, c("Control", "Treated"))
  expect_true(all(c("Unadjusted", "Adjusted") %in% rownames(f)))
  expect_true(all(vapply(f, is.character, logical(1L))))

  #A longitudinal object has one table per time point.
  b_msm <- bal.tab(list(treat ~ age, treat ~ age + educ),
                   data = cbind(lalonde, list(w = w_fixed)), weights = "w",
                   s.d.denom = "pooled")
  f_msm <- format(b_msm, component = "observations")
  expect_type(f_msm, "list")
  expect_length(f_msm, 2L)
  expect_s3_class(f_msm[[1L]], "data.frame")

  expect_err(format(b, component = "bogus"), "should be one of")
})

test_that("format() falls back to the segmented tables when there is no summary", {
  #A clustered, subclassified object has one table per cluster and subclass and no
  #summary across them, so there is no single table for print() to show at the top.
  b <- bal.tab(covs3(), treat = lalonde$treat, s.d.denom = "pooled",
               subclass = factor(rep(1:4, length.out = nrow(lalonde))),
               cluster = cl_idx)

  f <- format(b)

  expect_s3_class(f, "data.frame")
  expect_true(all(c("cluster", "subclass", "variable") %in% names(f)))
  expect_identical(nrow(f), nlevels(cl_idx) * 4L * 3L)
  expect_true(all(vapply(f, is.character, logical(1L))))
})

test_that("the extracted table is usable downstream", {
  b <- bal.tab(covs3(), treat = lalonde$treat, s.d.denom = "pooled", weights = w_fixed,
               un = TRUE, stats = c("m", "ks"), thresholds = c(m = .1))

  #The documented quarto route: a character data frame any renderer accepts.
  expect_no_error(knitr::kable(format(b)))

  #The tidy frame supports the obvious aggregations without reshaping.
  d <- as.data.frame(b)
  agg <- tapply(abs(d$estimate[d$stat == "mean.diffs"]),
                d$sample[d$stat == "mean.diffs"], max)
  expect_length(agg, 2L)
  expect_true(all(is.finite(agg)))

  #No list columns, no factors: it round-trips through a csv.
  expect_false(any(vapply(d, is.list, logical(1L))))
  expect_false(any(vapply(d, is.factor, logical(1L))))
})

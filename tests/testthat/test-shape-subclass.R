#`bal.tab()` with subclassification: the arguments handled by
#`base.bal.tab.subclass()` rather than by `print()`.

test_that("subclasses reject multi-category treatments", {
  covs <- lalonde[c("age", "educ")]

  expect_err(bal.tab(covs, treat = lalonde$race, subclass = sub_idx),
             "not currently compatible with subclasses")
})

test_that("subclassification works with a continuous treatment", {
  covs <- lalonde[c("age", "educ", "married", "race")]
  sub_levels <- as.character(sort(unique(sub_idx)))

  b <- bal.tab(covs, treat = lalonde$re75, subclass = sub_idx, un = TRUE,
               subclass.summary = TRUE, thresholds = c(cor = .1),
               disp = c("means", "sds"))

  expect_s3_class(b, "bal.tab.subclass")
  expect_named(b$Subclass.Balance, sub_levels)
  expect_no_error(capture.output(print(b)))

  #Only adjusted values exist within a subclass, because the subclassification is the
  #adjustment; the summary carries both samples. The threshold column belongs to the
  #adjusted sample alone, as it does for a binary treatment -- `print()` predicts the
  #layout from the print options, so an unexpected column shifts every column after it.
  expect_identical(names(b$Subclass.Balance[[1L]]),
                   c("Type", "M.Adj", "SD.Adj", "Corr.Adj", "R.Threshold"))
  expect_identical(names(b$Balance.Across.Subclass),
                   c("Type", "M.Un", "SD.Un", "Corr.Un",
                     "M.Adj", "SD.Adj", "Corr.Adj", "R.Threshold"))

  #The summary weights each subclass by its share of the subclassified units.
  p <- as.numeric(prop.table(table(sub_idx)))
  from_subclasses <- function(col) {
    do.call("cbind", lapply(b$Subclass.Balance, `[[`, col))
  }

  expect_equal(b$Balance.Across.Subclass$Corr.Adj,
               drop(from_subclasses("Corr.Adj") %*% p))

  #Standard deviations are combined in quadrature, so the summary value is the pooled
  #within-subclass standard deviation rather than an average of standard deviations.
  expect_equal(b$Balance.Across.Subclass$SD.Adj,
               sqrt(drop(from_subclasses("SD.Adj")^2 %*% p)))
  expect_true(all(b$Balance.Across.Subclass$SD.Adj <=
                    b$Balance.Across.Subclass$SD.Un))

  #Subclassification does not change the covariate distribution, so the summary means
  #are the means of the original sample.
  expect_equal(b$Balance.Across.Subclass$M.Adj,
               b$Balance.Across.Subclass$M.Un)

  #The tally and max-imbalance tables have one column or row per subclass.
  expect_named(b$Balanced.correlations.Subclass, paste("Subclass", sub_levels))
  expect_identical(rownames(b$Max.Imbalance.correlations.Subclass),
                   paste("Subclass", sub_levels))

  #A continuous treatment's sample-size table is a single `Total` row.
  expect_identical(rownames(b$Observations), "Total")
  expect_identical(colnames(b$Observations), c(sub_levels, "All"))
  expect_equal(b$Observations[["All"]], nrow(lalonde))
  expect_equal(sum(b$Observations[sub_levels]), nrow(lalonde))
})

test_that("continuous subclassification honors the display arguments", {
  covs <- lalonde[c("age", "educ", "married", "race")]

  #`s.weights` make a subclass's share its share of the population, which is what keeps
  #the summary means equal to the unadjusted ones.
  b_sw <- bal.tab(covs, treat = lalonde$re75, subclass = sub_idx, s.weights = sw_fixed,
                  un = TRUE, subclass.summary = TRUE, disp = "means")

  expect_equal(b_sw$Balance.Across.Subclass$M.Adj,
               b_sw$Balance.Across.Subclass$M.Un)

  #`abs` folds the summary statistic, and some correlation here is negative, so this
  #is not vacuous.
  b <- bal.tab(covs, treat = lalonde$re75, subclass = sub_idx,
               subclass.summary = TRUE, un = TRUE)
  b_abs <- bal.tab(covs, treat = lalonde$re75, subclass = sub_idx,
                   subclass.summary = TRUE, un = TRUE, abs = TRUE)

  expect_true(any(b$Balance.Across.Subclass$Corr.Adj < 0))
  expect_equal(b_abs$Balance.Across.Subclass$Corr.Adj,
               abs(b$Balance.Across.Subclass$Corr.Adj))

  #`quick = FALSE` computes the columns otherwise skipped, and the per-subclass tables
  #then carry statistics the summary does not -- `print()` must cope with both.
  b_slow <- bal.tab(covs, treat = lalonde$re75, subclass = sub_idx, quick = FALSE,
                    stats = c("correlations", "spearman.correlations"),
                    subclass.summary = TRUE, un = TRUE)
  expect_no_error(capture.output(print(b_slow)))
  expect_no_error(capture.output(print(b_slow, which.subclass = .all)))
  expect_no_error(capture.output(print(b_slow, stats = "spearman.correlations")))

  #And the other consumers of a subclassified object.
  expect_s3_class(love.plot(b), "love.plot")
  expect_s3_class(as.data.frame(b), "data.frame")
  expect_s3_class(format(b), "data.frame")

  #Nested inside a cluster, as a binary treatment can be.
  b_cl <- bal.tab(covs, treat = lalonde$re75, subclass = sub_idx, cluster = cl_idx)
  expect_s3_class(b_cl, "bal.tab.cluster")
  expect_s3_class(b_cl$Cluster.Balance[[1L]], "bal.tab.subclass")
  expect_no_error(capture.output(print(b_cl)))
})

test_that("which.subclass and disp.subclass are honored at bal.tab() time", {
  covs <- lalonde[c("age", "educ")]

  #The deprecated `disp.subclass` selects every subclass.
  b <- bal.tab(covs, treat = lalonde$treat, subclass = sub_idx,
               s.d.denom = "pooled", disp.subclass = TRUE)
  expect_match(squish(capture.output(print(b))), "Subclass 1")

  #`which.subclass` given to bal.tab() rather than print().
  b2 <- bal.tab(covs, treat = lalonde$treat, subclass = sub_idx,
                s.d.denom = "pooled", which.subclass = 1:2)
  out <- squish(capture.output(print(b2)))
  expect_match(out, "Subclass 1")
  expect_false(grepl("Subclass 4", out, fixed = TRUE))

  #`.all` selects every subclass and turns the summary off.
  b3 <- bal.tab(covs, treat = lalonde$treat, subclass = sub_idx,
                s.d.denom = "pooled", which.subclass = .all)
  expect_match(squish(capture.output(print(b3))), "Subclass 4")

  #Indices that match nothing leave the summary on.
  b4 <- bal.tab(covs, treat = lalonde$treat, subclass = sub_idx,
                s.d.denom = "pooled", which.subclass = c(90L, 91L))
  expect_true(is_not_null(b4[["Balance.Across.Subclass"]]))
})

test_that("thresholds are computed per subclass", {
  covs <- lalonde[c("age", "educ", "married")]

  b <- bal.tab(covs, treat = lalonde$treat, subclass = sub_idx,
               s.d.denom = "pooled", subclass.summary = TRUE,
               thresholds = c(m = .1))

  #A tally and a max-imbalance table across subclasses.
  expect_true(is_not_null(b[["Balanced.mean.diffs.Subclass"]]) ||
                is_not_null(b[["Balanced.mean.diffs"]]))
  expect_match(squish(capture.output(print(b))), "Balance tally for")

  #The threshold column reaches the per-subclass tables.
  expect_true("M.Threshold" %in% names(b[["Subclass.Balance"]][[1L]]))
})

test_that("per-subclass thresholds use each statistic's own absolute value", {
  #Variance ratios are compared on `pmax(x, 1/x)`, not `abs(x)`, so a ratio below
  #1 can still exceed the threshold. The per-subclass tables previously used
  #`abs()` for every statistic and so labeled those as balanced.
  covs <- lalonde[c("age", "educ", "re74")]

  b <- bal.tab(covs, treat = lalonde$treat, subclass = sub_idx,
               s.d.denom = "pooled", stats = c("mean.diffs", "variance.ratios"),
               thresholds = c(v = 2), disp.subclass = TRUE)

  for (i in names(b[["Subclass.Balance"]])) {
    sb <- b[["Subclass.Balance"]][[i]]
    vr <- sb[["V.Ratio.Adj"]]

    expect_identical(sb[["V.Threshold"]],
                     ifelse(pmax(vr, 1 / vr) < 2, "Balanced, <2", "Not Balanced, >2"),
                     info = i)
  }

  #At least one subclass must actually exercise the ratio: a variance ratio below
  #1 whose reciprocal exceeds the threshold. Otherwise the test proves nothing.
  vrs <- unlist(lapply(b[["Subclass.Balance"]], `[[`, "V.Ratio.Adj"))
  expect_true(any(vrs < 1 & 1 / vrs > 2))

  #Mean differences are unaffected -- their registry entry is plain `abs()`.
  b2 <- bal.tab(covs, treat = lalonde$treat, subclass = sub_idx,
                s.d.denom = "pooled", thresholds = c(m = .1), disp.subclass = TRUE)

  sb2 <- b2[["Subclass.Balance"]][[1L]]
  expect_identical(sb2[["M.Threshold"]],
                   ifelse(abs(sb2[["Diff.Adj"]]) < .1,
                          "Balanced, <0.1", "Not Balanced, >0.1"))
})

test_that("quick, s.weights, and raw statistics work with subclasses", {
  covs <- lalonde[c("age", "educ", "married")]

  #`quick = FALSE` computes the columns that are otherwise skipped, and the
  #result must still print.
  b <- bal.tab(covs, treat = lalonde$treat, subclass = sub_idx,
               s.d.denom = "pooled", quick = FALSE, subclass.summary = TRUE)
  expect_no_error(capture.output(print(b)))
  expect_true(is_not_null(b[["Balance.Across.Subclass"]]))

  #Sampling weights change the reported balance.
  b_sw <- bal.tab(covs, treat = lalonde$treat, subclass = sub_idx,
                  s.d.denom = "pooled", s.weights = sw_fixed, un = TRUE)
  b_no <- bal.tab(covs, treat = lalonde$treat, subclass = sub_idx,
                  s.d.denom = "pooled", un = TRUE)
  expect_false(isTRUE(all.equal(b_sw$Balance.Across.Subclass$Diff.Un,
                                b_no$Balance.Across.Subclass$Diff.Un)))

  #Raw statistics skip the standardization factor entirely.
  expect_no_message(bal.tab(covs, treat = lalonde$treat, subclass = sub_idx,
                            binary = "raw", continuous = "raw"))

  #Interactions and polynomials under subclassification.
  b_int <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat,
                   subclass = sub_idx, s.d.denom = "pooled", int = TRUE)
  expect_true(any(startsWith(rownames(b_int$Balance.Across.Subclass), "age * educ")))

  b_poly <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat,
                    subclass = sub_idx, s.d.denom = "pooled", poly = 2)
  expect_gt(nrow(b_poly$Balance.Across.Subclass), 2L)
})

test_that("subclassification from a matchit object is accepted", {
  skip_if_not_installed("MatchIt")

  m <- fx("matchit_sub")

  b <- bal.tab(m, method = "subclassification", subclass.summary = TRUE)
  expect_s3_class(b, "bal.tab.subclass")
  expect_no_error(capture.output(print(b)))
})

test_that("the subclass sample-size table counts discarded units correctly", {
  skip_if_not_installed("MatchIt")

  m <- MatchIt::matchit(treat ~ age + educ + re74, data = lalonde,
                        method = "subclass", subclass = 4, discard = "both")

  nn <- bal.tab(m)$Observations

  expect_true("Discarded" %in% colnames(nn))

  #Every column's Total is the sum of its treatment-group rows. The `Discarded`
  #column previously reported `length(treat)` here instead of `sum(discarded)`.
  groups <- setdiff(rownames(nn), "Total")

  for (col in colnames(nn)) {
    expect_equal(nn["Total", col], sum(nn[groups, col]), info = col)
  }

  expect_equal(nn["Total", "Discarded"], sum(m$discarded))
  expect_equal(nn["Total", "All"], nrow(lalonde))

  #With nothing discarded the column is dropped entirely.
  m2 <- MatchIt::matchit(treat ~ age + educ + re74, data = lalonde,
                         method = "subclass", subclass = 4)

  expect_false("Discarded" %in% colnames(bal.tab(m2)$Observations))
})

test_that("longitudinal treatments accept data frames and addl/distance lists", {
  covs <- lalonde[c("age", "educ")]
  ps <- fitted(glm(treat ~ age + educ, data = lalonde, family = binomial))

  #The data.frame.list method, with a treatment per time point.
  b <- bal.tab(list(covs, covs),
               treat.list = list(lalonde$treat, lalonde$nodegree))
  expect_s3_class(b, "bal.tab.msm")
  expect_length(b$Time.Balance, 2L)

  #`addl.list` and `distance.list` are applied per time point.
  b2 <- bal.tab(list(covs, covs),
                treat.list = list(lalonde$treat, lalonde$nodegree),
                addl.list = list(lalonde["re74"], lalonde["re75"]),
                distance.list = list(ps, ps))
  expect_true("re74" %in% rownames(b2$Time.Balance[[1L]]$Balance))
  expect_true("re75" %in% rownames(b2$Time.Balance[[2L]]$Balance))
  expect_identical(rownames(b2$Time.Balance[[1L]]$Balance)[1L], "distance")
})

test_that("cluster.summary is omitted rather than erroring with subclassification", {
  # Regression test: a subclassified cluster carries no `nweights` in its print
  # options, and `balance_summary()` derived `no.adj` as `nweights == 0`, giving
  # `logical(0)` and failing with "missing value where TRUE/FALSE needed".
  covs <- lalonde[c("age", "educ", "married")]
  sub <- factor(rep(1:4, length.out = nrow(lalonde)))

  for (q in c(TRUE, FALSE)) {
    b <- bal.tab(covs, treat = lalonde$treat, s.d.denom = "pooled", subclass = sub,
                 cluster = cl_idx, cluster.summary = TRUE, cluster.fun = "mean",
                 thresholds = c(m = .1), quick = q)

    #Omitted, as it already is for multiply imputed data.
    expect_false("Balance.Across.Clusters" %in% names(b))
    expect_named(b$Cluster.Balance, levels(cl_idx))
    expect_s3_class(b$Cluster.Balance[[1L]], "bal.tab.subclass")

    #Printing does not error, and still shows the per-cluster subclass tables.
    expect_no_error(out <- capture.output(print(b, cluster.summary = TRUE)))
    expect_true(any(grepl("Cluster: a", out, fixed = TRUE)))
  }

  #Without subclassification the summary is still produced.
  b2 <- bal.tab(covs, treat = lalonde$treat, s.d.denom = "pooled", weights = w_fixed,
                cluster = cl_idx, cluster.summary = TRUE, cluster.fun = "mean")
  expect_s3_class(b2$Balance.Across.Clusters, "data.frame")
})

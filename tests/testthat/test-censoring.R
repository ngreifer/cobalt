#`bal.tab()` with a censoring indicator: the balance of the units still under
#observation against the full at-risk sample.
#
#The oracle for every number here is the stacking documented under "Assessing balance"
#in `WeightIt::.cens()`, computed with cobalt's ordinary binary machinery. Asserting
#against that rather than against values this implementation produced is the point: it
#is an independent statement of what the comparison is supposed to be.

#The censoring weights, as a censoring model would produce them: the inverse
#probability of remaining under observation, and exactly 0 once censored.
cens_w <- function(C = cens_idx) {
  p <- fitted(glm(C ~ age + educ + race, data = lalonde, family = binomial))

  ifelse(C == 0, 1 / (1 - p), 0)
}

#The comparison spelled out by hand, with no censoring code involved.
cens_manual <- function(covs, C = cens_idx, weights = NULL, s.weights = NULL, ...) {
  u <- which(!is.na(C) & C == 0)
  a <- which(!is.na(C))

  bal.tab(rbind(covs[u, , drop = FALSE], covs[a, , drop = FALSE]),
          treat = rep(0:1, times = c(length(u), length(a))),
          weights = if (is_not_null(weights)) c(weights[u], rep.int(1, length(a))),
          s.weights = if (is_not_null(s.weights)) s.weights[c(u, a)],
          estimand = "ATT", s.d.denom = "treated", ...)
}

test_that(".cens() tags an indicator without changing its values", {
  C <- .cens(cens_idx)

  expect_identical(as.numeric(C), as.numeric(cens_idx))
  expect_identical(get.treat.type(C), "censoring")
  expect_s3_class(C, "treat")
  expect_identical(attr(C, "treat.name"), "cens_idx")

  #Logicals and 0/1 factors are accepted and coerced.
  expect_identical(as.numeric(.cens(cens_idx == 1L)), as.numeric(cens_idx))
  expect_identical(as.numeric(.cens(factor(cens_idx))), as.numeric(cens_idx))

  #Missing values are preserved: a unit censored earlier has no later indicator.
  C_na <- cens_idx
  C_na[1:10] <- NA
  expect_identical(which(is.na(.cens(C_na))), 1:10)

  #Anything else is an error, including a coercion that would silently produce NAs.
  expect_err(.cens(lalonde$age), "must contain only the values")
  expect_err(.cens(factor(c("Yes", "No"))), "must contain only the values")
  expect_err(.cens(rep(NA, 10L)), "no non-missing values")
})

test_that("censoring balance is the uncensored sample against the full one", {
  covs <- lalonde[c("age", "educ", "race")]
  w <- cens_w()

  args <- list(un = TRUE, stats = c("mean.diffs", "variance.ratios",
                                    "ks.statistics", "ovl.coefficients"),
               disp = c("means", "sds"))

  b <- do.call(bal.tab, c(list(covs, treat = .cens(cens_idx), weights = w), args))
  manual <- do.call(cens_manual, c(list(covs, weights = w), args))

  expect_s3_class(b, "bal.tab.cens")
  expect_s3_class(b, "bal.tab.bin")

  #Every statistic and every moment, not just the mean differences.
  expect_equal(b$Balance, manual$Balance)

  #Group 0 is the uncensored sample and group 1 the full one, so the unadjusted means
  #of group 1 are the means of every at-risk unit.
  expect_equal(b$Balance$M.1.Un,
               unname(col_w_mean(splitfactor(covs, drop.first = "if2"))))

  #`un = TRUE` is the unweighted uncensored sample against the full one, which is not
  #the same comparison as the weighted one.
  expect_false(isTRUE(all.equal(b$Balance$Diff.Un, b$Balance$Diff.Adj)))
  expect_equal(b$Balance$Diff.Un,
               do.call(cens_manual, c(list(covs), args))$Balance$Diff.Un)
})

test_that("the censoring sample-size table counts the two samples", {
  covs <- lalonde[c("age", "educ")]
  w <- cens_w()

  nn <- bal.tab(covs, treat = .cens(cens_idx), weights = w)$Observations

  expect_identical(colnames(nn), "Total")
  expect_identical(rownames(nn), c("Full", "Uncensored", "Adjusted", "Censored"))

  expect_equal(nn["Full", "Total"], sum(!is.na(cens_idx)))
  expect_equal(nn["Uncensored", "Total"], sum(cens_idx == 0))
  expect_equal(nn["Censored", "Total"], sum(cens_idx == 1))

  #The two samples partition the at-risk units.
  expect_equal(nn["Uncensored", "Total"] + nn["Censored", "Total"],
               nn["Full", "Total"])

  #Weighting can only cost effective sample size.
  expect_lt(nn["Adjusted", "Total"], nn["Uncensored", "Total"])

  #With no weights there is nothing to adjust, and the tag says sizes not ESSs.
  nn_un <- bal.tab(covs, treat = .cens(cens_idx))$Observations
  expect_identical(rownames(nn_un), c("Full", "Uncensored", "Censored"))
  expect_identical(attr(nn_un, "tag"), "Sample sizes")

  #One row per set of weights.
  nn_2 <- bal.tab(covs, treat = .cens(cens_idx),
                  weights = list(a = w, b = sqrt(w)))$Observations
  expect_identical(rownames(nn_2), c("Full", "Uncensored", "a", "b", "Censored"))

  #A `Censored` row nobody reached is dropped, as `Discarded` is.
  nn_0 <- bal.tab(covs, treat = .cens(rep(0L, n_lalonde)))$Observations
  expect_false("Censored" %in% rownames(nn_0))
})

test_that("units with a missing indicator are in neither sample", {
  covs <- lalonde[c("age", "educ")]
  C <- cens_idx
  C[1:40] <- NA
  w <- cens_w(ifelse(is.na(C), 0L, C))

  b <- bal.tab(covs, treat = .cens(C), weights = w, un = TRUE)

  #The at-risk sample is the units with a non-missing indicator, so the whole table is
  #what it would be with those units removed.
  expect_equal(b$Balance, cens_manual(covs, C = C, weights = w, un = TRUE)$Balance)

  nn <- b$Observations
  expect_equal(nn["Full", "Total"], sum(!is.na(C)))
  expect_equal(nn["Uncensored", "Total"] + nn["Censored", "Total"],
               nn["Full", "Total"])
})

test_that("every input interface reaches the same table", {
  covs <- lalonde[c("age", "educ", "race")]
  w <- cens_w()
  d <- transform(lalonde, cens_idx = cens_idx)

  from_df <- bal.tab(covs, treat = .cens(cens_idx), weights = w)
  from_f <- bal.tab(.cens(cens_idx) ~ age + educ + race, data = d, weights = w)

  expect_equal(from_f$Balance, from_df$Balance)

  #A formula written for WeightIt works unchanged, marker and all.
  skip_if_not_installed("WeightIt")

  from_wf <- bal.tab(WeightIt::.cens(cens_idx) ~ age + educ + race,
                     data = d, weights = w)
  expect_equal(from_wf$Balance, from_df$Balance)

  #The treatment keeps the indicator's own name rather than `.cens(cens_idx)`.
  expect_true("age" %in% rownames(from_f$Balance))
  expect_false(any(grepl(".cens", rownames(from_f$Balance), fixed = TRUE)))
})

test_that("a weightit censoring object is accepted", {
  W <- fx("weightit_cens")

  b <- bal.tab(W, un = TRUE, stats = c("mean.diffs", "ks.statistics"))

  expect_s3_class(b, "bal.tab.cens")

  #The censoring propensity score is P(C = 1 | X); balance on it is reported as it is
  #for a binary treatment.
  expect_identical(rownames(b$Balance)[1L], "prop.score")

  #And the covariates match the documented manual stacking on the same object.
  manual <- cens_manual(splitfactor(W$covs, drop.first = "if2"),
                        weights = W$weights, un = TRUE,
                        stats = c("mean.diffs", "ks.statistics"))

  expect_equal(b$Balance[rownames(manual$Balance), ], manual$Balance)

  expect_equal(b$Observations["Uncensored", "Total"], sum(W$treat == 0))
})

test_that("s.d.denom names the two samples and defaults to the full one", {
  covs <- lalonde[c("age", "educ")]
  w <- cens_w()

  b_default <- bal.tab(covs, treat = .cens(cens_idx), weights = w)
  b_full <- bal.tab(covs, treat = .cens(cens_idx), weights = w, s.d.denom = "full")

  expect_equal(b_default$Balance, b_full$Balance)

  #The target is fixed by the design, so nothing is inferred and nothing is announced.
  expect_no_message(bal.tab(covs, treat = .cens(cens_idx), weights = w))

  #"full" is the full sample's SD, which is the stacked treated group's.
  expect_equal(b_full$Balance$Diff.Adj,
               cens_manual(covs, weights = w)$Balance$Diff.Adj)

  #The other values are accepted and genuinely different.
  b_unc <- bal.tab(covs, treat = .cens(cens_idx), weights = w,
                   s.d.denom = "uncensored")
  expect_false(isTRUE(all.equal(b_unc$Balance$Diff.Adj, b_full$Balance$Diff.Adj)))

  for (sdd in c("pooled", "all", "weighted", "hedges")) {
    expect_true(all(is.finite(bal.tab(covs, treat = .cens(cens_idx), weights = w,
                                      s.d.denom = sdd)$Balance$Diff.Adj)),
                label = sdd)
  }

  #Treated/control are not the vocabulary here.
  expect_err(bal.tab(covs, treat = .cens(cens_idx), weights = w,
                     s.d.denom = "treated"),
             '`s.d.denom` should be one of "full", "uncensored"')
})

test_that("censoring composes with clusters and imputations", {
  covs <- lalonde[c("age", "educ", "married")]
  w <- cens_w()

  b_cl <- bal.tab(covs, treat = .cens(cens_idx), weights = w, cluster = cl_idx,
                  cluster.summary = TRUE, cluster.fun = "mean")

  expect_s3_class(b_cl, "bal.tab.cluster")
  expect_named(b_cl$Cluster.Balance, levels(cl_idx))
  expect_s3_class(b_cl$Cluster.Balance[[1L]], "bal.tab.cens")
  expect_s3_class(b_cl$Balance.Across.Clusters, "data.frame")
  expect_no_error(capture.output(print(b_cl)))

  #Each cluster's table is the censoring comparison within that cluster.
  in.a <- cl_idx == "a"
  expect_equal(b_cl$Cluster.Balance[["a"]]$Balance,
               cens_manual(covs[in.a, ], C = cens_idx[in.a],
                           weights = w[in.a])$Balance)

  #The at-risk counts add up across clusters.
  expect_equal(b_cl$Observations["Full", "Total"], n_lalonde)

  b_imp <- bal.tab(covs, treat = .cens(cens_idx), weights = w, imp = imp_idx,
                   imp.summary = TRUE, imp.fun = "mean")

  expect_s3_class(b_imp, "bal.tab.imp")
  expect_s3_class(b_imp$Imputation.Balance[[1L]], "bal.tab.cens")
  expect_no_error(capture.output(print(b_imp)))

  in.1 <- imp_idx == 1L
  expect_equal(b_imp$Imputation.Balance[["1"]]$Balance,
               cens_manual(covs[in.1, ], C = cens_idx[in.1],
                           weights = w[in.1])$Balance)

  #A cluster with no uncensored units has no sample to compare against its own.
  cl_bad <- as.character(cl_idx)
  cl_bad[cens_idx == 1L] <- "d"
  expect_err(bal.tab(covs, treat = .cens(cens_idx), weights = w,
                     cluster = factor(cl_bad)),
             "every unit is censored in at least one cluster")
})

test_that("s.weights, thresholds, and the other display options apply", {
  covs <- lalonde[c("age", "educ", "race")]
  w <- cens_w()

  #Sampling weights apply to both samples, since they describe the population.
  b_sw <- bal.tab(covs, treat = .cens(cens_idx), weights = w, s.weights = sw_fixed,
                  un = TRUE)
  expect_equal(b_sw$Balance,
               cens_manual(covs, weights = w, s.weights = sw_fixed,
                           un = TRUE)$Balance)
  expect_false(isTRUE(all.equal(
    b_sw$Balance$Diff.Adj,
    bal.tab(covs, treat = .cens(cens_idx), weights = w)$Balance$Diff.Adj)))

  #Thresholds, tallies, and the max-imbalance table.
  b_thr <- bal.tab(covs, treat = .cens(cens_idx), weights = w,
                   thresholds = c(m = .05))
  expect_true("M.Threshold" %in% names(b_thr$Balance))
  expect_match(squish(capture.output(print(b_thr))), "Balance tally for")

  #Interactions, polynomials, additional covariates, and a distance measure.
  b_int <- bal.tab(covs[c("age", "educ")], treat = .cens(cens_idx), weights = w,
                   int = TRUE, poly = 2, addl = covs["race"],
                   distance = lalonde$re74)
  expect_true(any(startsWith(rownames(b_int$Balance), "age * educ")))
  expect_true("race_black" %in% rownames(b_int$Balance))

  #`quick = FALSE` and `abs`.
  expect_no_error(capture.output(print(
    bal.tab(covs, treat = .cens(cens_idx), weights = w, quick = FALSE))))
  expect_equal(bal.tab(covs, treat = .cens(cens_idx), weights = w,
                       abs = TRUE)$Balance$Diff.Adj,
               abs(bal.tab(covs, treat = .cens(cens_idx),
                           weights = w)$Balance$Diff.Adj))
})

test_that("the other consumers of a bal.tab accept a censoring one", {
  covs <- lalonde[c("age", "educ", "race")]
  w <- cens_w()

  b <- bal.tab(covs, treat = .cens(cens_idx), weights = w, un = TRUE,
               thresholds = c(m = .05))

  local_null_device()

  expect_s3_class(suppressWarnings(love.plot(b)), "love.plot")

  d <- as.data.frame(b)
  expect_s3_class(d, "data.frame")
  expect_setequal(unique(d$sample), c("Unadjusted", "Adj"))

  expect_s3_class(format(b), "data.frame")
  expect_no_error(capture.output(print(b)))
})

test_that("what does not apply to a censoring indicator says so", {
  covs <- lalonde[c("age", "educ")]
  w <- cens_w()

  #Subclassification is not a way of estimating censoring weights.
  expect_err(bal.tab(covs, treat = .cens(cens_idx), subclass = sub_idx),
             "subclasses are not allowed with a censoring indicator")
  expect_err(bal.tab(covs, treat = .cens(cens_idx), match.strata = sub_idx),
             "matching strata are not allowed with a censoring indicator")

  #`estimand` and `focal` name a treatment group; there is none to name.
  expect_wrn(bal.tab(covs, treat = .cens(cens_idx), weights = w, estimand = "ATT"),
             "`estimand` does not apply to a censoring indicator")
  expect_wrn(bal.tab(covs, treat = .cens(cens_idx), weights = w, focal = "1"),
             "`focal` does not apply to a censoring indicator")

  #With everyone censored there is no sample left.
  expect_err(bal.tab(covs, treat = .cens(rep(1L, n_lalonde))),
             "every unit is censored")

  #`bal.plot()` does not handle the two samples yet, and says so rather than failing
  #somewhere further in.
  expect_err(bal.plot(covs, treat = .cens(cens_idx), weights = w, var.name = "age"),
             "does not yet support censoring indicators")
})

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
#
#The moment columns come back named for the pseudo-treatment's groups -- `M.0` and
#`M.1` -- and are renamed to the two samples, so that the correspondence (group 0 is the
#uncensored sample, group 1 the full one) is part of what the comparison asserts rather
#than something the test works around.
cens_manual <- function(covs, C = cens_idx, weights = NULL, s.weights = NULL,
                        subclass = NULL, ...) {
  u <- which(!is.na(C) & C == 0)
  a <- which(!is.na(C))

  b <- bal.tab(rbind(covs[u, , drop = FALSE], covs[a, , drop = FALSE]),
               treat = rep(0:1, times = c(length(u), length(a))),
               weights = if (is_not_null(weights)) c(weights[u], rep.int(1, length(a))),
               s.weights = if (is_not_null(s.weights)) s.weights[c(u, a)],
               subclass = if (is_not_null(subclass)) subclass[c(u, a)],
               estimand = "ATT", s.d.denom = "treated", ...)

  for (tab in c("Balance", "Balance.Across.Subclass")) {
    if (is.data.frame(b[[tab]])) {
      names(b[[tab]]) <- names(b[[tab]]) |>
        sub("^(M|SD)\\.0\\.", "\\1.Uncensored.", x = _) |>
        sub("^(M|SD)\\.1\\.", "\\1.Full.", x = _)
    }
  }

  b
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
  is.na(C_na[1:10]) <- TRUE
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

  #The moment columns are named for the two samples rather than for a control and a
  #treated group, and the full sample's unadjusted means are the means of every at-risk
  #unit.
  expect_true(all(c("M.Uncensored.Un", "M.Full.Un", "SD.Uncensored.Adj", "SD.Full.Adj")
                  %in% names(b$Balance)))
  expect_false(any(grepl("M.0", names(b$Balance), fixed = TRUE)))

  expect_equal(b$Balance$M.Full.Un,
               unname(col_w_mean(splitfactor(covs, drop.first = "if2"))))

  #Weighting the uncensored units cannot change the target, so the full sample's
  #moments are the same before and after.
  expect_equal(b$Balance$M.Full.Adj, b$Balance$M.Full.Un)

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

  #The heading says the table holds effective sample sizes, so the adjusted row is not
  #starred as well and no footnote explains a star. The only row that could be an
  #effective sample size is that one; the rest are counts, as their names say.
  out <- capture.output(print(bal.tab(covs, treat = .cens(cens_idx), weights = w)))

  expect_false(any(grepl("*", out, fixed = TRUE)))
  expect_false(any(grepl("indicates effective sample size", out, fixed = TRUE)))
  expect_true(any(grepl("Effective sample sizes", out, fixed = TRUE)))
})

test_that("units with a missing indicator are in neither sample", {
  covs <- lalonde[c("age", "educ")]
  C <- cens_idx
  is.na(C[1:40]) <- TRUE
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

  #Matching strata become weights before the two samples exist, and say the same thing
  #given to `subclass`, which is supported.
  expect_err(bal.tab(covs, treat = .cens(cens_idx), match.strata = sub_idx),
             "matching strata are not allowed with a censoring indicator; supply them to `subclass`")

  #`estimand` and `focal` name a treatment group; there is none to name.
  expect_wrn(bal.tab(covs, treat = .cens(cens_idx), weights = w, estimand = "ATT"),
             "`estimand` does not apply to a censoring indicator")
  expect_wrn(bal.tab(covs, treat = .cens(cens_idx), weights = w, focal = "1"),
             "`focal` does not apply to a censoring indicator")

  #With everyone censored there is no sample left.
  expect_err(bal.tab(covs, treat = .cens(rep(1L, n_lalonde))),
             "every unit is censored")
})

test_that("subclassification is itself a censoring solution", {
  covs <- lalonde[c("age", "educ", "race")]
  p <- fitted(glm(cens_idx ~ age + educ + race, data = lalonde, family = binomial))
  sub <- findInterval(p, quantile(p, c(.25, .5, .75))) + 1L

  b <- bal.tab(covs, treat = .cens(cens_idx), subclass = sub, un = TRUE,
               thresholds = c(m = .1))

  expect_s3_class(b, "bal.tab.cens")
  expect_s3_class(b, "bal.tab.subclass")
  expect_named(b$Subclass.Balance, as.character(sort(unique(sub))))
  expect_no_error(capture.output(print(b, which.subclass = .all)))

  #Within each subclass, the uncensored units are compared against every at-risk unit in
  #it -- nothing is weighted, because the subclassification is the adjustment. The oracle
  #is the same documented stacking with the subclasses stacked along with it.
  manual <- cens_manual(covs, subclass = sub, un = TRUE, thresholds = c(m = .1))

  expect_equal(b$Subclass.Balance, manual$Subclass.Balance)

  #The summary across subclasses is subclassification expressed as censoring weights:
  #`n_k / n_{k,uncensored}` for the uncensored units and 1 for the full sample.
  n.k <- table(sub)
  n.k.u <- table(sub[cens_idx == 0])
  w <- ifelse(cens_idx == 0, (n.k / n.k.u)[as.character(sub)], 0)

  expect_equal(b$Balance.Across.Subclass,
               cens_manual(covs, weights = w, un = TRUE,
                           thresholds = c(m = .1))$Balance)

  #Which is also what the same stacking produces when it is subclassified.
  expect_equal(b$Balance.Across.Subclass, manual$Balance.Across.Subclass)

  #The sample sizes divide the two samples up by subclass.
  nn <- b$Observations
  expect_identical(rownames(nn), c("Full", "Uncensored", "Censored"))
  expect_identical(colnames(nn), c(names(b$Subclass.Balance), "All"))
  expect_equal(unlist(nn["Uncensored", ] + nn["Censored", ]), unlist(nn["Full", ]))
  expect_equal(nn[["Full", "All"]], n_lalonde)

  #And it composes with clusters, as the other subclassified shapes do.
  b_cl <- bal.tab(covs, treat = .cens(cens_idx), subclass = sub, cluster = cl_idx)
  expect_s3_class(b_cl, "bal.tab.cluster")
  expect_s3_class(b_cl$Cluster.Balance[[1L]], "bal.tab.cens")
})

test_that("bal.plot() shows the weighted sample against the full one", {
  covs <- lalonde[c("age", "educ", "race")]
  w <- cens_w()
  p <- fitted(glm(cens_idx ~ age + educ + race, data = lalonde, family = binomial))
  sub <- findInterval(p, quantile(p, c(.25, .5, .75))) + 1L

  local_null_device()

  #`bal.plot()` gives each group its own layer, so what it plots is their union. The
  #column naming the panel is `which` when the panels are samples and `subclass` when
  #they are subclasses.
  plot_data <- function(p) {
    do.call("rbind", lapply(p$layers, function(l) {
      panel <- intersect(c("which", "subclass"), names(l$data))[1L]

      setNames(l$data[c("treat", panel, "weights")], c("treat", "panel", "weights"))
    }))
  }

  #The two groups are the two samples, and they are labelled as such rather than as a
  #control and a treated group.
  pl <- bal.plot(covs, treat = .cens(cens_idx), weights = w, var.name = "age",
                 which = "both")
  d <- plot_data(pl)

  expect_setequal(as.character(unique(d$treat)), c("Uncensored", "Full"))
  expect_setequal(as.character(unique(d$panel)),
                  c("Unadjusted Sample", "Adjusted Sample"))

  #The full sample is every at-risk unit; the uncensored sample is the units still under
  #observation.
  is.full <- d$treat == "Full"
  is.adj <- d$panel == "Adjusted Sample"

  expect_equal(sum(is.full & is.adj), sum(!is.na(cens_idx)))
  expect_equal(sum(!is.full & is.adj), sum(cens_idx == 0))

  #The full sample is the target and so is never reweighted, in either panel; the
  #uncensored sample carries the weights in the adjusted panel only. (`bal.plot()`
  #rescales weights for plotting, so what is checked is that they vary, not their values.)
  expect_true(all_the_same(d$weights[is.full]))
  expect_true(all_the_same(d$weights[!is.full & !is.adj]))
  expect_false(all_the_same(d$weights[!is.full & is.adj]))

  #The other plot shapes work too.
  for (a in list(list(var.name = "race"),
                 list(var.name = "age", type = "histogram", mirror = TRUE),
                 list(var.name = "age", type = "ecdf"),
                 list(var.name = "age", disp.means = TRUE))) {
    pl <- do.call(bal.plot, c(list(covs, treat = .cens(cens_idx), weights = w), a))
    expect_no_error(ggplot2::ggplot_build(pl))
  }

  #With subclasses, one panel per subclass.
  pl <- bal.plot(covs, treat = .cens(cens_idx), subclass = sub, var.name = "age")
  expect_setequal(as.character(unique(plot_data(pl)$panel)),
                  paste("Subclass", sort(unique(sub))))

  #And on a weightit object.
  skip_if_not_installed("WeightIt")
  expect_no_error(ggplot2::ggplot_build(bal.plot(fx("weightit_cens"),
                                                var.name = "age", which = "both")))
})

#A longitudinal data set of the shape `weightitMSM()` takes: a treatment, then a
#censoring indicator, then a treatment only the units still under observation have. The
#risk set entering the third model is the units `C` did not censor.
msm_lalonde <- function() {
  d <- transform(lalonde, C = cens_idx)
  is.na(d$nodegree[cens_idx == 1L]) <- TRUE
  d
}

msm_forms <- function() {
  list(treat ~ age + educ,
       .cens(C) ~ age + educ,
       nodegree ~ age + educ)
}

test_that("a longitudinal list gives a table per model, of that model's kind", {
  d <- msm_lalonde()
  covs <- d[c("age", "educ")]
  w <- w_fixed
  risk3 <- which(cens_idx == 0L)

  b <- bal.tab(msm_forms(), data = d, weights = w, un = TRUE)

  expect_length(b$Time.Balance, 3L)
  expect_s3_class(b$Time.Balance[[1L]], "bal.tab.bin")
  expect_s3_class(b$Time.Balance[[2L]], "bal.tab.cens")
  expect_s3_class(b$Time.Balance[[3L]], "bal.tab.bin")

  #Each table is the table that model would have produced on its own, computed on the
  #units still under observation entering it. The censoring one is the manual stacking
  #that every other number in this file is pinned to.
  expect_equal(b$Time.Balance[[1L]]$Balance,
               bal.tab(covs, treat = d$treat, weights = w, s.d.denom = "pooled",
                       un = TRUE)$Balance)

  expect_equal(b$Time.Balance[[2L]]$Balance,
               cens_manual(covs, C = d$C, weights = w, un = TRUE)$Balance)

  expect_equal(b$Time.Balance[[3L]]$Balance,
               bal.tab(covs[risk3, ], treat = d$nodegree[risk3], weights = w[risk3],
                       s.d.denom = "pooled", un = TRUE)$Balance)

  #A list mixing kinds of model has nothing to aggregate across, exactly as a list
  #mixing continuous and binary treatments does not.
  expect_null(b$Balance.Across.Times)

  #Sample sizes are produced anyway, one table per time point, each of its own shape.
  expect_length(b$Observations, 3L)
  expect_identical(names(b$Observations[[2L]]), "Total")
  expect_true(all(c("Full", "Uncensored", "Censored") %in% rownames(b$Observations[[2L]])))
  expect_equal(unname(unlist(b$Observations[[2L]]["Full", ])), sum(!is.na(d$C)))
  expect_equal(sum(unlist(b$Observations[[3L]]["Unadjusted", ])), length(risk3))

  #They are not printed as a table of their own, though: that table belongs to the
  #summary, and each time point's own table has already reported them. So each time point
  #is named once, over its own block, rather than twice.
  out <- capture.output(print(b))

  expect_identical(sum(grepl("Uncensored", out, fixed = TRUE)), 1L)
  expect_identical(sum(grepl("2. Censoring: C", out, fixed = TRUE)), 1L)

  #Asking for the summary of a list that has one does print the collected table, so there
  #each time point is named twice. (No censoring indicator here, so `nodegree` is the
  #recorded one rather than the blanked one `msm_lalonde()` carries.)
  b_sum <- bal.tab(list(treat ~ age + educ, nodegree ~ age + educ), data = lalonde,
                   weights = w, s.d.denom = "pooled", which.time = .all,
                   msm.summary = TRUE)

  expect_identical(sum(grepl("2. Treatment: nodegree",
                             capture.output(print(b_sum)), fixed = TRUE)), 2L)

  #Every entry a censoring model is not a mixture, so that list is summarized.
  C2 <- rep(c(0L, 0L, 0L, 1L, 0L), length.out = n_lalonde)
  is.na(C2[cens_idx == 1L]) <- TRUE

  b_all <- bal.tab(list(.cens(C) ~ age + educ, .cens(C2) ~ age + educ),
                   data = transform(d, C2 = C2), weights = w, un = TRUE,
                   msm.summary = TRUE)

  expect_s3_class(b_all$Balance.Across.Times, "data.frame")
  expect_true("Max.Diff.Adj" %in% names(b_all$Balance.Across.Times))

  #The second indicator's full sample is the risk set entering it, not the cohort.
  expect_equal(unname(unlist(b_all$Observations[[2L]]["Full", ])), sum(cens_idx == 0L))
})

test_that("the risk set is accumulated from the censoring indicators", {
  skip_if_not_installed("WeightIt")

  W <- fx("weightitmsm_cens")
  X <- x2base(process_obj(W))

  #Two independent derivations of the same thing: WeightIt builds `at.risk` while
  #fitting the models, and this rebuilds it from the indicators alone.
  at.risk <- .msm_at_risk(X$treat.list)

  expect_length(at.risk, ncol(W$at.risk))
  for (ti in seq_along(at.risk)) {
    expect_identical(unname(at.risk[[ti]]), unname(W$at.risk[, ti]))
  }

  #A censoring model's own risk set still holds the units it is about to remove -- they
  #are half of what it compares -- and the model after it does not.
  expect_true(all(at.risk[[2L]]))
  expect_identical(base::which(at.risk[[3L]]), base::which(cens_idx == 0L))

  b <- bal.tab(W, un = TRUE)

  expect_length(b$Time.Balance, 3L)
  expect_s3_class(b$Time.Balance[[2L]], "bal.tab.cens")
  expect_identical(names(b$Time.Balance), c("treat", "cens_idx", "nodegree"))
})

test_that("both shapes of a weightitMSM censoring object are read the same way", {
  skip_if_not_installed("WeightIt")

  #WeightIt 2.1.0 treats censoring as a treatment type, so `treat.list`/`covs.list` hold
  #every model in the order it was fit. Development versions between 2.0.0 and 2.1.0
  #segregated the censoring models into `cens.list`/`cens.covs.list`, with `cens.time`
  #recording where each belonged. Both have to give the same tables, so whichever shape
  #the fixture happens to be built with, this converts it to the other and compares.
  W <- fx("weightitmsm_cens")

  segregated <- is_not_null(W[["cens.list"]])

  other <- W

  if (segregated) {
    models <- .weightitMSM_models(W)

    other[["treat.list"]] <- models[["treat.list"]]
    other[["covs.list"]] <- models[["covs.list"]]
    other[c("cens.list", "cens.covs.list", "cens.formula.list", "cens.time")] <- NULL
  }
  else {
    is.cens <- vapply(W[["treat.list"]], .is_cens, logical(1L))

    other[["treat.list"]] <- W[["treat.list"]][!is.cens]
    other[["covs.list"]] <- W[["covs.list"]][!is.cens]
    other[["cens.list"]] <- unname(W[["treat.list"]][is.cens])
    other[["cens.covs.list"]] <- unname(W[["covs.list"]][is.cens])
    other[["cens.time"]] <- setNames(base::which(is.cens),
                                     names(W[["treat.list"]])[is.cens])
  }

  #The shapes really are different, or the comparison below says nothing.
  expect_false(identical(other[["treat.list"]], W[["treat.list"]]))
  expect_identical(is_not_null(other[["cens.list"]]), !segregated)

  expect_equal(bal.tab(other, un = TRUE), bal.tab(W, un = TRUE))
})

test_that("a `weightitMSM` object holds its censoring models among the treatments", {
  #\pkg{WeightIt} 2.1.0 keeps every model in `treat.list`/`covs.list` in the order it was
  #fit, censoring included, rather than apart in `cens.list`/`cens.covs.list` with
  #`cens.time` recording where each sat. Built by hand rather than by fitting so that the
  #shape under test is that one whatever version of WeightIt happens to be installed:
  #`.weightitMSM_models()` still accepts the older shape, and `fx("weightitmsm_cens")`
  #produces whichever the installed version makes.
  d <- msm_lalonde()
  covs <- d[c("age", "educ")]

  W <- list(weights = w_fixed,
            treat.list = list(treat = d$treat,
                              C = .cens(d$C),
                              nodegree = d$nodegree),
            covs.list = list(covs, covs, covs),
            estimand = "ATE",
            call = quote(weightitMSM(f, data = d)))
  class(W) <- c("weightitMSM", "weightit")

  b <- bal.tab(W, un = TRUE)

  #The same three tables, of the same three kinds, as the formula list it stands for
  expect_identical(names(b$Time.Balance), c("treat", "C", "nodegree"))
  expect_s3_class(b$Time.Balance[[2L]], "bal.tab.cens")

  b_f <- bal.tab(msm_forms(), data = d, weights = w_fixed, un = TRUE)

  expect_equal(grab(b$Time.Balance, "Balance"),
               grab(b_f$Time.Balance, "Balance"))
  expect_equal(b$Observations, b_f$Observations)
})

test_that("which units are at risk does not depend on how the data records the rest", {
  #A unit censored at time 2 has no treatment at time 3. Reading the risk set from the
  #indicators rather than from where treatments happen to be missing means it makes no
  #difference whether the data blanks those treatments or records them.
  d_na <- msm_lalonde()
  d_recorded <- transform(lalonde, C = cens_idx)

  b_na <- bal.tab(msm_forms(), data = d_na, weights = w_fixed, un = TRUE)
  b_recorded <- bal.tab(msm_forms(), data = d_recorded, weights = w_fixed, un = TRUE)

  expect_equal(grab(b_na$Time.Balance, "Balance"),
               grab(b_recorded$Time.Balance, "Balance"))
  expect_equal(b_na$Observations, b_recorded$Observations)

  #A unit still under observation whose treatment is unknown is a different thing, and
  #the error says which model it turned up in.
  d_bad <- d_na
  is.na(d_bad$nodegree[base::which(cens_idx == 0L)[1L]]) <- TRUE

  expect_err(bal.tab(msm_forms(), data = d_bad, weights = w_fixed),
             'in time point "nodegree": missing values must not exist in `treat`')

  #Without a censoring indicator in the list a missing treatment stays an error, as it
  #is for a point treatment.
  expect_err(bal.tab(list(treat ~ age + educ, nodegree ~ age + educ),
                     data = d_na, weights = w_fixed),
             "missing values must not exist in `treat`")
})

test_that("a longitudinal censoring model composes with imputations", {
  d <- msm_lalonde()
  risk3 <- which(cens_idx == 0L)

  b <- bal.tab(msm_forms(), data = d, weights = w_fixed, imp = imp_idx, un = TRUE)

  #`msm` ranks above `imp`, so each time point holds the imputations and the censoring
  #model is a leaf inside them.
  expect_s3_class(b$Time.Balance[[2L]], "bal.tab.imp")
  expect_s3_class(b$Time.Balance[[2L]]$Imputation.Balance[[1L]], "bal.tab.cens")

  in.imp <- imp_idx == 1L

  expect_equal(b$Time.Balance[[2L]]$Imputation.Balance[[1L]]$Balance,
               cens_manual(d[in.imp, c("age", "educ")], C = d$C[in.imp],
                           weights = w_fixed[in.imp], un = TRUE)$Balance)

  #Each imputation's last time point sees only the units still under observation in it.
  expect_equal(sum(unlist(b$Time.Balance[[3L]]$Imputation.Balance[[1L]]$Observations["Unadjusted", ])),
               sum(in.imp[risk3]))

  #A mixture still has no summary across time points, and still has sample sizes.
  expect_null(b$Balance.Across.Times)
  expect_length(b$Observations, 3L)
})

test_that("a longitudinal censoring model composes with clusters", {
  d <- msm_lalonde()
  risk3 <- which(cens_idx == 0L)

  b <- bal.tab(msm_forms(), data = d, weights = w_fixed, cluster = cl_idx, un = TRUE)

  #`cluster` ranks above `msm`, so each cluster gets its own longitudinal balance table
  #and the censoring model is a leaf inside it.
  expect_s3_class(b, "bal.tab.cluster")

  in.cluster <- cl_idx == levels(cl_idx)[1L]
  b_cl <- b$Cluster.Balance[[1L]]

  expect_s3_class(b_cl, "bal.tab.msm")
  expect_s3_class(b_cl$Time.Balance[[2L]], "bal.tab.cens")

  #Each cluster's censoring table is that cluster's two samples, and the last time point
  #sees only the units still under observation in it.
  expect_equal(b_cl$Time.Balance[[2L]]$Balance,
               cens_manual(d[in.cluster, c("age", "educ")], C = d$C[in.cluster],
                           weights = w_fixed[in.cluster], un = TRUE)$Balance)

  expect_equal(sum(unlist(b_cl$Observations[[3L]]["Unadjusted", ])),
               sum(in.cluster[risk3]))
})

test_that("a time point is named for its position, its kind, and its variable", {
  d <- msm_lalonde()

  #The form `WeightIt::summary.weightitMSM()` uses, so that the same model goes by the
  #same name in both packages. Numbered by position in the list, so a censoring model has
  #a number too and it is the one `which.time` takes.
  b <- bal.tab(msm_forms(), data = d, weights = w_fixed)

  labels <- c("1. Treatment: treat", "2. Censoring: C", "3. Treatment: nodegree")

  expect_identical(attr(b, "print.options")[["time.labels"]], labels)

  out <- capture.output(print(b))
  for (l in labels) {
    expect_true(any(grepl(l, out, fixed = TRUE)), info = l)
  }

  #The collected sample sizes, which print with the summary, are labelled the same way.
  b_sum <- bal.tab(list(treat ~ age + educ, nodegree ~ age + educ), data = lalonde,
                   weights = w_fixed, s.d.denom = "pooled", msm.summary = TRUE)

  expect_true(any(grepl("2. Treatment: nodegree",
                        capture.output(print(b_sum)), fixed = TRUE)))

  #And so are the facets of a plot. `which.time` still selects on the bare variable name.
  local_null_device()

  p <- bal.plot(msm_forms(), data = d, var.name = "age", which = "unadjusted")
  panels <- do.call("rbind", lapply(p$layers, function(l) l$data))

  expect_identical(levels(panels$time), labels)

  p2 <- bal.plot(msm_forms(), data = d, var.name = "age", which.time = "C",
                 which = "unadjusted")
  panels2 <- do.call("rbind", lapply(p2$layers, function(l) l$data))

  expect_identical(levels(panels2$time), "2. Censoring: C")
})

test_that("bal.plot() draws a longitudinal censoring model at each kind of time point", {
  d <- msm_lalonde()
  w <- cens_w()
  risk3 <- which(cens_idx == 0L)

  local_null_device()

  plot_data <- function(p) do.call("rbind", lapply(p$layers, function(l) l$data))

  #A treatment time point shows only the units still under observation entering it; a
  #censoring one shows the two samples it compares, stacked.
  p1 <- bal.plot(msm_forms(), data = d, weights = w, var.name = "age", which.time = 1,
                 which = "unadjusted")
  expect_identical(nrow(plot_data(p1)), n_lalonde)

  p2 <- bal.plot(msm_forms(), data = d, weights = w, var.name = "age", which.time = 2,
                 which = "both")
  d2 <- plot_data(p2)

  expect_setequal(as.character(unique(d2$treat)), c("Uncensored", "Full"))
  expect_identical(nrow(d2), 2L * (n_lalonde + length(risk3)))

  #The full sample is the target, so it is never reweighted; the uncensored sample
  #carries the weights in the adjusted panel. (`bal.plot()` rescales weights within a
  #panel, so what is checked is that they vary, not their values.)
  is.adj <- d2$which == "Adjusted Sample"

  expect_length(unique(d2$weights[is.adj & d2$treat == "Full"]), 1L)
  expect_gt(length(unique(d2$weights[is.adj & d2$treat == "Uncensored"])), 1L)

  p3 <- bal.plot(msm_forms(), data = d, weights = w, var.name = "age", which.time = 3,
                 which = "unadjusted")
  expect_identical(nrow(plot_data(p3)), length(risk3))

  #Every time point at once, and time points picked out by the name of their model.
  all.times <- plot_data(bal.plot(msm_forms(), data = d, weights = w, var.name = "age",
                                  which = "adjusted"))
  expect_identical(nrow(all.times),
                   as.integer(n_lalonde + (n_lalonde + length(risk3)) + length(risk3)))

  by.name <- plot_data(bal.plot(msm_forms(), data = d, weights = w, var.name = "age",
                                which.time = "C", which = "unadjusted"))
  expect_identical(nrow(by.name), n_lalonde + length(risk3))
  expect_setequal(as.character(unique(by.name$treat)), c("Uncensored", "Full"))
})

test_that("a treat carries its type and group names through subsetting", {
  #`[` on a `treat` keeps every attribute the class exists to hold; `[.factor` alone
  #would drop them, which would silently turn a censoring indicator back into a binary
  #treatment.
  C <- process_treat(.cens(cens_idx))

  expect_identical(get.treat.type(C), "censoring")
  expect_identical(treat_names(C),
                   setNames(c("Uncensored", "Censored"), c("control", "treated")))
  expect_identical(group_labels(C), c("Uncensored", "Censored"))

  for (obj in list(C[1:20], subset_treat(C, 1:20))) {
    expect_identical(get.treat.type(obj), "censoring")
    expect_identical(treat_names(obj), treat_names(C))
    expect_identical(group_labels(obj), group_labels(C))
  }

  #A subset in which nobody is censored is still a censoring indicator, with both
  #levels named.
  none <- subset_treat(C, which(cens_idx == 0))
  expect_identical(get.treat.type(none), "censoring")
  expect_identical(treat_names(none), treat_names(C))

  #An ordinary binary treatment keeps `c("0", "1")` as its column labels -- `M.0` has
  #always meant the control group -- and narrows its levels when subset.
  tb <- process_treat(lalonde$treat)

  expect_identical(group_labels(tb), c("0", "1"))
  expect_identical(treat_names(tb[1:20]), treat_names(tb))

  #Narrowing a binary treatment to a single group leaves nothing to compare, which
  #`subset_treat()` reports rather than silently returning a one-level treatment.
  expect_err(subset_treat(tb, which(lalonde$treat == 1)),
             "the treatment must have at least two unique values")

  #A binary treatment's balance table still uses the positional labels.
  b <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat, s.d.denom = "pooled",
               disp = "means")
  expect_true("M.0.Un" %in% names(b$Balance))
})

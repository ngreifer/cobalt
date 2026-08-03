#The cluster and imputation wrappers: the continuous-treatment standardization
#branch and the threshold summary, neither of which the tests above reach.

test_that("clusters and imputations handle continuous treatments", {
  covs <- lalonde[c("age", "educ", "married")]

  #A continuous treatment with correlation statistics takes a different
  #standardization path from the binary case.
  b_cl <- bal.tab(covs, treat = lalonde$re75, cluster = cl_idx, weights = w_fixed,
                  un = TRUE, binary = "std",
                  stats = c("correlations", "spearman.correlations"))
  expect_s3_class(b_cl, "bal.tab.cluster")
  expect_no_error(capture.output(print(b_cl)))

  b_imp <- bal.tab(covs, treat = lalonde$re75, imp = imp_idx, weights = w_fixed,
                   un = TRUE, binary = "std", stats = "correlations")
  expect_s3_class(b_imp, "bal.tab.imp")
  expect_no_error(capture.output(print(b_imp)))
})

test_that("a single aggregation function produces a threshold summary", {
  covs <- lalonde[c("age", "educ", "married")]

  #`cluster.fun`/`imp.fun` of length 1 adds the balance tally and max-imbalance
  #tables across the grouping dimension.
  b_cl <- bal.tab(covs, treat = lalonde$treat, cluster = cl_idx, weights = w_fixed,
                  s.d.denom = "pooled", cluster.summary = TRUE,
                  cluster.fun = "mean", thresholds = c(m = .1))
  expect_true(any(startsWith(names(b_cl), "Balanced")))
  out <- squish(capture.output(print(b_cl)))
  expect_match(out, "Balance tally for")

  b_imp <- bal.tab(covs, treat = lalonde$treat, imp = imp_idx, weights = w_fixed,
                   s.d.denom = "pooled", imp.summary = TRUE,
                   imp.fun = "mean", thresholds = c(m = .1))
  expect_true(any(startsWith(names(b_imp), "Balanced")))
  expect_match(squish(capture.output(print(b_imp))), "Balance tally for")
})

test_that("errors inside one group are labelled with that group", {
  covs <- lalonde[c("age", "educ")]

  #All-zero weights in one arm of a single cluster or imputation.
  w_cl <- rep(1, n_lalonde)
  w_cl[cl_idx == "b" & lalonde$treat == 1] <- 0
  expect_err(bal.tab(covs, treat = lalonde$treat, cluster = cl_idx,
                     weights = w_cl, s.d.denom = "pooled"),
             "in cluster")

  w_im <- rep(1, n_lalonde)
  w_im[imp_idx == 1 & lalonde$treat == 1] <- 0
  expect_err(bal.tab(covs, treat = lalonde$treat, imp = imp_idx,
                     weights = w_im, s.d.denom = "pooled"),
             "in imputation")
})

test_that("the wrapper classes nest in the documented precedence", {
  covs <- lalonde[c("age", "educ")]

  expect_s3_class(bal.tab(covs, treat = lalonde$treat, cluster = cl_idx,
                          imp = imp_idx, s.d.denom = "pooled"), "bal.tab.cluster")
  expect_s3_class(bal.tab(covs, treat = lalonde$race, cluster = cl_idx),
                  "bal.tab.cluster")
  expect_s3_class(bal.tab(covs, treat = lalonde$treat, subclass = sub_idx,
                          cluster = cl_idx, s.d.denom = "pooled"),
                  "bal.tab.cluster")
  expect_s3_class(bal.tab(list(treat ~ age, treat ~ age + educ), data = lalonde,
                          cluster = cl_idx), "bal.tab.cluster")
  expect_s3_class(bal.tab(list(treat ~ age, treat ~ age + educ), data = lalonde,
                          imp = imp_idx), "bal.tab.msm")
})

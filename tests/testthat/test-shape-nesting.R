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

test_that("arguments supplied once are expanded across imputations", {
  # `length_imp_process()` stacks anything given for a single imputation up to the
  # full stacked data, keeping each imputation's block in the order `imp` gives. It
  # used to read and write its caller's frame by name; it is now pure, so the
  # expansion is worth pinning directly.
  set.seed(4)
  n <- 100L
  d1 <- data.frame(age = rnorm(n), educ = rnorm(n))
  t1 <- rbinom(n, 1L, .5)

  #Two imputations' worth of covariates, one imputation's worth of treatment.
  d2 <- rbind(d1, transform(d1, age = age + 1))

  for (sorted in c(TRUE, FALSE)) {
    imp <- if (sorted) rep(1:2, each = n) else rep(1:2, times = n)

    b <- bal.tab(d2, treat = t1, imp = imp, s.d.denom = "pooled")

    expect_length(b$Imputation.Balance, 2L)

    #Each imputation's table must equal the one computed from that imputation's rows
    #with the treatment as given -- which is only true if `treat` was replicated into
    #the right positions.
    for (i in c("1", "2")) {
      ref <- bal.tab(d2[imp == as.integer(i), , drop = FALSE], treat = t1,
                     s.d.denom = "pooled")

      expect_equal(b$Imputation.Balance[[i]]$Balance, ref$Balance, info = i)
    }
  }

  #Data frames, factors, and weights expand the same way.
  imp <- rep(1:2, times = n)
  b_w <- bal.tab(d2, treat = t1, imp = imp, s.d.denom = "pooled",
                 weights = rep(c(1, 2), n), addl = d1["age"],
                 cluster = factor(rep(c("a", "b"), length.out = n)))
  expect_length(b_w$Cluster.Balance, 2L)

  #A length that matches neither `imp` nor one imputation is an error.
  expect_err(bal.tab(d2[1:50, ], treat = t1, imp = rep(1:2, each = n)),
             "must have the same number of observations as")
})

test_that("covariates expanded across imputations keep their names", {
  # The expansion is done with `[`, which keeps a matrix's `dim` and `dimnames` and
  # drops everything else -- including the `co.names` that records what each column is.
  # A covariate matrix without them reads as having no covariates at all, so supplying
  # the covariates once and `imp` for several imputations produced a balance table with
  # no rows in it.
  set.seed(9)
  n <- 100L
  covs <- data.frame(age = rnorm(n),
                     race = factor(sample(c("a", "b", "c"), n, replace = TRUE)))

  for (sorted in c(TRUE, FALSE)) {
    imp <- if (sorted) rep(1:2, each = n) else rep(1:2, times = n)
    treat <- rbinom(2L * n, 1L, .5)

    b <- bal.tab(covs, treat = treat, imp = imp, s.d.denom = "pooled")

    expect_length(b$Imputation.Balance, 2L)

    # Each imputation is the same covariates against that imputation's treatment, which
    # holds only if they were replicated into the right positions and kept their names
    # on the way.
    for (i in c("1", "2")) {
      ref <- bal.tab(covs, treat = treat[imp == as.integer(i)], s.d.denom = "pooled")

      expect_equal(b$Imputation.Balance[[i]]$Balance, ref$Balance, info = i)
    }

    # A split factor's dummies are named from `co.names`, so they are the rows that go
    # missing first.
    expect_identical(rownames(b$Balance.Across.Imputations),
                     c("age", "race_a", "race_b", "race_c"))
  }
})

test_that("a covariate set left with no columns is subset with everything else", {
  # `.get_C2()` returns a matrix with a row per unit and no columns when nothing is left
  # to assess balance on -- here a factor that takes a single value within every cluster,
  # whose dummies say nothing about any of them. That matrix has length zero, which used
  # to read as "not supplied" and leave it at full length while every other slot was cut
  # down to the cluster.
  b <- expect_no_warning(
    bal.tab(data.frame(g = cl_idx), treat = lalonde$treat, cluster = cl_idx,
            s.d.denom = "pooled")
  )

  expect_identical(nrow(b$Cluster.Balance[[1L]]$Balance), 0L)

  # A covariate that survives alongside it is unaffected.
  b2 <- expect_no_warning(
    bal.tab(data.frame(g = cl_idx, age = lalonde$age), treat = lalonde$treat,
            cluster = cl_idx, s.d.denom = "pooled")
  )

  expect_identical(rownames(b2$Cluster.Balance[[1L]]$Balance), "age")
})

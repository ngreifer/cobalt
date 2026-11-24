skip_on_cran()
test_that("bal.compute() and bal.init() work", {
    set.seed(100)
    data("lalonde")
    cov_names <- c("age", "educ", "race", "married", "nodegree", "re74", 
                   "re75")
    sw <- runif(nrow(lalonde))
    
    tb <- lalonde$treat
    tm <- factor(sample(LETTERS[1:4], nrow(lalonde), TRUE))
    tc <- rnorm(nrow(lalonde))
    
    # Binary treatment
    init <- bal.init(lalonde[cov_names], tb, "smd.mean")
    expect_s3_class(init, "bal.init")
    expect_equal(bal.compute(init), 0.60909674)
    
    
    # Multi-category treatment
    #(The treatment here is `tm`, not `tb`; a multi-category `bal.init()` builds a
    #different object and must be exercised with one.)
    init <- bal.init(lalonde[cov_names], tm, "smd.mean")
    expect_s3_class(init, "bal.init")
    expect_equal(bal.compute(init), 0.10969168)

    # Continuous treatment
    init <- bal.init(lalonde[cov_names], tc, "p.mean")
    expect_s3_class(init, "bal.init")
    expect_equal(bal.compute(init), 0.033749304)


})

# ---------------------------------------------------------------------------
# The blocks above pin three values and cross-check the `smd`/`ks` families against
# `bal.tab()` for `pairwise`. The blocks below extend that cross-check to every
# statistic family and aggregator, and add the algebraic properties each statistic
# is supposed to have.
#
# Cross-checks are preferred over hardcoded values because they state the *defini-
# tion* of the statistic: `smd.rms` is the root-mean-square of the per-covariate
# SMDs that `bal.tab()` reports, and no reorganization of the internals may change
# that.

bc_fixture <- function() {
  data("lalonde", package = "cobalt")

  cov_names <- c("age", "educ", "race", "married", "nodegree", "re74", "re75")

  set.seed(505)

  list(covs = lalonde[cov_names],
       num_covs = lalonde[c("age", "educ", "re74", "re75")],
       tb = lalonde$treat,
       tm = lalonde$race,
       tc = lalonde$re75,
       w = runif(nrow(lalonde)),
       sw = runif(nrow(lalonde)))
}

test_that("the .mean/.max/.rms aggregators summarize the bal.tab() columns", {
  skip_if_not_installed("patrick")

  eps <- if (capabilities("long.double")) 1e-8 else 1e-1

  f <- bc_fixture()

  #`binary = "std"` so every row is on the same scale as the `smd.*` statistics,
  #and `estimand = "ATE"` so `bal.init()` and `bal.tab()` use the same denominator.
  bt <- bal.tab(f$covs, treat = f$tb, binary = "std", estimand = "ATE",
                weights = f$w, un = TRUE,
                stats = c("mean.diffs", "ks.statistics", "ovl.coefficients"))$Balance

  columns <- c(smd = "Diff.Adj", ks = "KS.Adj", ovl = "OVL.Adj")

  patrick::with_parameters_test_that(
    "{family}.{agg}",
    {
      v <- bt[[columns[[family]]]]

      #`mean` and `max` are taken over absolute values; `rms` squares first, so the
      #sign drops out on its own.
      expected <- c(mean = mean(abs(v)), max = max(abs(v)), rms = sqrt(mean(v^2)))

      expect_equal(bal.compute(bal.init(f$covs, f$tb, paste0(family, ".", agg),
                                        estimand = "ATE"),
                               weights = f$w),
                   expected[[agg]], tolerance = eps)
    },
    .cases = expand.grid(agg = c("mean", "max", "rms"), family = names(columns),
                         stringsAsFactors = FALSE)
  )
})

test_that("the continuous aggregators summarize the correlation columns", {
  eps <- if (capabilities("long.double")) 1e-8 else 1e-1

  f <- bc_fixture()

  covs <- f$covs[c("age", "educ", "married", "nodegree", "re74")]

  bt <- bal.tab(covs, treat = f$tc, weights = f$w, un = TRUE,
                stats = c("correlations", "spearman.correlations"))$Balance

  #`p.*` are Pearson correlations, `s.*` Spearman. Mixing the two up is invisible
  #in the class and the shape of the output, so pin both.
  columns <- c(p = "Corr.Adj", s = "S.Corr.Adj")

  #The Spearman statistics must actually differ from the Pearson ones -- they were
  #once computed identically because of a bad loop index, which no shape or class
  #check would have caught.
  expect_false(isTRUE(all.equal(bal.compute(bal.init(covs, f$tc, "s.mean"), weights = f$w),
                                bal.compute(bal.init(covs, f$tc, "p.mean"), weights = f$w))))

  skip_if_not_installed("patrick")

  patrick::with_parameters_test_that(
    "{family}.{agg}",
    {
      v <- bt[[columns[[family]]]]
      expected <- c(mean = mean(abs(v)), max = max(abs(v)), rms = sqrt(mean(v^2)))

      expect_equal(bal.compute(bal.init(covs, f$tc, paste0(family, ".", agg)),
                               weights = f$w),
                   expected[[agg]], tolerance = eps)
    },
    .cases = expand.grid(agg = c("mean", "max", "rms"), family = names(columns),
                         stringsAsFactors = FALSE)
  )
})

test_that("the aggregators summarize the bal.tab() columns when covariates have missing values", {
  eps <- if (capabilities("long.double")) 1e-8 else 1e-1

  #Two continuous covariates and one binary covariate with missing values. Each
  #statistic must use only the units with its covariate observed, in `bal.init()`
  #as in `bal.tab()`. `bal.tab()` adds a row flagging the missingness of each such
  #covariate, which `bal.init()` has no counterpart for, so only the covariates'
  #own rows are compared.
  covs <- lalonde_mis[c("age", "re74", "re75", "married")]

  expect_wrn({
    bt <- bal.tab(covs, treat = lalonde_mis$treat, binary = "std", estimand = "ATE",
                  weights = w_fixed,
                  stats = c("mean.diffs", "ks.statistics", "ovl.coefficients"))
  }, "Missing values exist in the covariates")

  bt <- bt$Balance[names(covs), ]

  columns <- c(smd = "Diff.Adj", ks = "KS.Adj", ovl = "OVL.Adj")

  for (family in names(columns)) {
    stat <- paste0(family, ".max")

    expect_equal(bal.compute(bal.init(covs, lalonde_mis$treat, stat, estimand = "ATE"),
                             weights = w_fixed),
                 max(abs(bt[[columns[[family]]]])), tolerance = eps, info = stat)
  }

  expect_wrn({
    bt <- bal.tab(covs, treat = lalonde_mis$re78, weights = w_fixed,
                  stats = c("correlations", "spearman.correlations"))
  }, "Missing values exist in the covariates")

  bt <- bt$Balance[names(covs), ]

  columns <- c(p = "Corr.Adj", s = "S.Corr.Adj")

  for (family in names(columns)) {
    stat <- paste0(family, ".max")

    expect_equal(bal.compute(bal.init(covs, lalonde_mis$re78, stat), weights = w_fixed),
                 max(abs(bt[[columns[[family]]]])), tolerance = eps, info = stat)
  }
})

test_that("multi-category aggregators summarize the pairwise bal.tab() columns", {
  skip_if_not_installed("patrick")

  eps <- if (capabilities("long.double")) 1e-8 else 1e-1

  f <- bc_fixture()

  covs <- f$covs[c("age", "educ", "married", "nodegree", "re74", "re75")]

  bt <- bal.tab(covs, treat = f$tm, binary = "std", estimand = "ATE",
                weights = f$w, un = TRUE,
                stats = c("mean.diffs", "ks.statistics"))$Pair.Balance

  columns <- c(smd = "Diff.Adj", ks = "KS.Adj")

  patrick::with_parameters_test_that(
    "{family}.{agg}",
    {
      #Pairwise statistics pool every pair's covariates into one vector before
      #aggregating, so `max` is over all pairs at once, not per pair.
      v <- unlist(lapply(bt, function(i) i$Balance[[columns[[family]]]]))
      expected <- c(mean = mean(abs(v)), max = max(abs(v)), rms = sqrt(mean(v^2)))

      expect_equal(bal.compute(bal.init(covs, f$tm, paste0(family, ".", agg),
                                        estimand = "ATE"),
                               weights = f$w),
                   expected[[agg]], tolerance = eps)
    },
    .cases = expand.grid(agg = c("mean", "max", "rms"), family = names(columns),
                         stringsAsFactors = FALSE)
  )
})

test_that("bal.compute() is unchanged by transformations the statistics ignore", {
  f <- bc_fixture()

  #The Mahalanobis distance is invariant to any rescaling and shifting of the
  #covariates, because it standardizes by the covariance itself. This is a genuine
  #mathematical property of the statistic, not just of this implementation.
  rescaled <- f$num_covs
  rescaled$age <- rescaled$age * 100 + 7
  rescaled$re74 <- rescaled$re74 / 1000

  expect_equal(bal.compute(bal.init(rescaled, f$tb, "mahalanobis"), weights = f$w),
               bal.compute(bal.init(f$num_covs, f$tb, "mahalanobis"), weights = f$w))

  #SMDs are scale-free too, being differences divided by a standard deviation.
  expect_equal(bal.compute(bal.init(rescaled, f$tb, "smd.mean"), weights = f$w),
               bal.compute(bal.init(f$num_covs, f$tb, "smd.mean"), weights = f$w))

  #Rank-based and distributional statistics survive any monotone transformation.
  logged <- f$num_covs
  logged$re74 <- log1p(logged$re74)
  logged$age <- exp(logged$age / 50)

  expect_equal(bal.compute(bal.init(logged, f$tb, "ks.max"), weights = f$w),
               bal.compute(bal.init(f$num_covs, f$tb, "ks.max"), weights = f$w))
  expect_equal(bal.compute(bal.init(logged, f$tc, "s.mean"), weights = f$w),
               bal.compute(bal.init(f$num_covs, f$tc, "s.mean"), weights = f$w))

  skip_if_not_installed("patrick")

  #Weights enter every statistic as a ratio, so a constant multiple cannot matter.
  #A refactor that starts normalizing by `sum(w)` in one place but not another
  #shows up here.
  patrick::with_parameters_test_that(
    "weights times a constant: {stat}",
    {
      init <- bal.init(f$covs, f$tb, stat)

      expect_equal(bal.compute(init, weights = f$w),
                   bal.compute(init, weights = f$w * 7))
    },
    stat = c("smd.mean", "ks.max", "ovl.mean", "energy.dist", "mahalanobis")
  )
})

test_that("bal.compute() reaches its optimum at perfect balance", {
  f <- bc_fixture()

  #Every statistic is a non-negative discrepancy measure, minimized when the
  #weighted groups match. Duplicating the sample so each unit appears in both
  #groups makes the two groups identical by construction.
  doubled <- rbind(f$num_covs, f$num_covs)
  paired <- rep(0:1, each = nrow(f$num_covs))

  for (stat in c("smd.mean", "smd.max", "smd.rms", "ks.mean", "ks.max",
                 "mahalanobis", "energy.dist")) {
    expect_equal(bal.compute(bal.init(doubled, paired, stat)), 0,
                 tolerance = 1e-8, info = stat)
  }

  #The `ovl.*` statistics integrate a kernel density estimate, so identical groups
  #give a small numerical residual rather than an exact zero. Bounding it keeps the
  #test meaningful without pretending the integration is exact -- this is why the
  #`ovl.max` case in the "unweighted target" block above is commented out.
  for (stat in c("ovl.mean", "ovl.max", "ovl.rms")) {
    v <- bal.compute(bal.init(doubled, paired, stat))

    expect_gte(v, 0)
    expect_lt(v, 1e-3)
  }

  #All the aggregators are non-negative on real data, and ordered mean <= rms <= max
  #by the power-mean inequality.
  for (family in c("smd", "ks", "ovl")) {
    v <- vapply(c("mean", "rms", "max"),
                function(agg) {
                  bal.compute(bal.init(f$covs, f$tb, paste0(family, ".", agg)),
                              weights = f$w)
                },
                numeric(1L))

    expect_gte(v[["mean"]], 0)
    expect_lte(v[["mean"]], v[["rms"]] + 1e-12)
    expect_lte(v[["rms"]], v[["max"]] + 1e-12)
  }

  #`ks.*` and `ovl.*` are bounded above by 1, being distances between CDFs.
  for (stat in c("ks.max", "ovl.max")) {
    expect_lte(bal.compute(bal.init(f$covs, f$tb, stat), weights = f$w), 1)
  }
})

test_that("bal.compute() uses `s.weights` fixed at bal.init() time", {
  f <- bc_fixture()

  #`s.weights` belong to the sample, so they are baked into the `bal.init()` object
  #and multiply whatever weights are supplied later. For `ks.*`, whose value depends
  #on the weights only through their product, that makes the two routes identical.
  expect_equal(bal.compute(bal.init(f$covs, f$tb, "ks.max", s.weights = f$sw),
                           weights = f$w),
               bal.compute(bal.init(f$covs, f$tb, "ks.max"),
                           weights = f$w * f$sw))

  #For `smd.*` they are *not* interchangeable: the standard deviation in the
  #denominator is computed once from `s.weights` alone and is not reweighted by the
  #balancing weights. Pinning the difference documents which of the two the
  #denominator follows.
  with_sw <- bal.compute(bal.init(f$covs, f$tb, "smd.mean", s.weights = f$sw),
                         weights = f$w)
  folded <- bal.compute(bal.init(f$covs, f$tb, "smd.mean"), weights = f$w * f$sw)

  expect_false(isTRUE(all.equal(with_sw, folded)))

  expect_equal(with_sw,
               mean(abs(col_w_smd(splitfactor(f$covs, drop.first = FALSE),
                                  treat = f$tb, weights = f$w, s.weights = f$sw,
                                  s.d.denom = "pooled", std = TRUE))))

  #Constant `s.weights` are the same as none at all.
  expect_equal(bal.compute(bal.init(f$covs, f$tb, "smd.mean",
                                    s.weights = rep(2, length(f$tb))),
                           weights = f$w),
               bal.compute(bal.init(f$covs, f$tb, "smd.mean"), weights = f$w))
})

test_that("bal.compute() reports the same value however the init is reached", {
  f <- bc_fixture()

  #The `default` method builds an init and computes in one step; doing it in two
  #must give the same answer.
  expect_equal(bal.compute(f$covs, treat = f$tb, stat = "smd.mean", weights = f$w),
               bal.compute(bal.init(f$covs, f$tb, "smd.mean"), weights = f$w))

  #Omitting `weights` is the same as supplying uniform ones.
  init <- bal.init(f$covs, f$tb, "smd.mean")
  expect_equal(bal.compute(init),
               bal.compute(init, weights = rep(1, length(f$tb))))

  #A matrix of pre-split dummies must behave like the data frame it came from.
  expect_equal(bal.compute(bal.init(as.matrix(splitfactor(f$covs, drop.first = FALSE)),
                                    f$tb, "smd.mean"), weights = f$w),
               bal.compute(init, weights = f$w))

  #`l1.med` bins at random, so it is reproducible only under a fixed seed. Two runs
  #from the same seed must agree; that is what makes it usable in a test at all.
  set.seed(31)
  a <- bal.compute(bal.init(f$covs, f$tb, "l1.med"), weights = f$w)
  set.seed(31)
  b <- bal.compute(bal.init(f$covs, f$tb, "l1.med"), weights = f$w)

  expect_equal(a, b)
  expect_gte(a, 0)
  expect_lte(a, 1)
})

test_that("bal.init() records what it was given", {
  f <- bc_fixture()

  init <- bal.init(f$covs, f$tb, "smd.mean", estimand = "ATT", s.weights = f$sw)

  expect_s3_class(init, "bal.init")

  skip_if_not_installed("patrick")

  #`print()` must name the statistic; the phrase table behind it has been missing
  #entries before, which produced an empty line rather than an error.
  treats <- list(binary = f$tb, continuous = f$tc)

  patrick::with_parameters_test_that(
    "print() names {treat.type} {stat}",
    {
      out <- capture.output(print(bal.init(f$num_covs, treats[[treat.type]], stat)))

      expect_true(any(nzchar(trimws(out))))
      expect_false(any(grepl("NA", out, fixed = TRUE)))
    },
    .cases = stat_cases(c("binary", "continuous"))
  )
})

test_that("bal.compute() returns 0 for unweighted target balance statistics", {
  set.seed(1004)
  data("lalonde")
  cov_names <- c("age", "educ", "race", "married", "nodegree", "re74", 
                 "re75")
  sw <- runif(nrow(lalonde))
  
  expect_equal(bal.compute(lalonde[cov_names], stat = "smd.rms"),
               0)
  expect_equal(bal.compute(lalonde[cov_names], stat = "ks.max"),
               0)
  # expect_equal(bal.compute(lalonde[cov_names], stat = "ovl.max"),
  #              0)
  expect_equal(bal.compute(lalonde[cov_names], stat = "energy.dist"),
               0)
  expect_equal(bal.compute(lalonde[cov_names], stat = "mahalanobis"),
               0)
  
  expect_equal(bal.compute(lalonde[cov_names], stat = "smd.rms", s.weights = sw),
               0)
  expect_equal(bal.compute(lalonde[cov_names], stat = "ks.max", s.weights = sw),
               0)
  # expect_equal(bal.compute(lalonde[cov_names], stat = "ovl.max", s.weights = sw),
  #              0)
  expect_equal(bal.compute(lalonde[cov_names], stat = "energy.dist", s.weights = sw),
               0)
  expect_equal(bal.compute(lalonde[cov_names], stat = "mahalanobis", s.weights = sw),
               0)
})

test_that("bal.compute() works with pairwise = FALSE", {
  eps <- if (capabilities("long.double")) 1e-8 else 1e-1
  
  set.seed(100)
  data("lalonde")
  cov_names <- c("age", "educ", "race", "married", "nodegree", "re74", 
                 "re75")
  sw <- runif(nrow(lalonde))
  
  tb <- lalonde$treat
  tm <- factor(sample(LETTERS[1:4], nrow(lalonde), TRUE))

  w <- runif(nrow(lalonde))
  
  # Binary treatment
  ## SMD
  init <- bal.init(lalonde[cov_names], tb, "smd.mean")
  expect_s3_class(init, "bal.init")
  
  baltab <- bal.tab(lalonde[cov_names], treat = tb, binary = "std",
                    stats = "m", estimand = "ATE", weights = w,
                    un = TRUE)$Balance
  
  expect_equal(bal.compute(init),
               mean(abs(baltab$Diff.Un)),
               tolerance = eps)
  
  expect_equal(bal.compute(init, weights = w),
               mean(abs(baltab$Diff.Adj)),
               tolerance = eps)
  
  init <- bal.init(lalonde[cov_names], treat = tb,
                   stat = "smd.mean",
                   pairwise = FALSE)
  expect_s3_class(init, "bal.init")
  
  baltab <- bal.tab(lalonde[cov_names], treat = tb, binary = "std",
                    stats = "m", estimand = "ATE", weights = w,
                    un = TRUE, pairwise = FALSE)$Pair.Balance
  
  expect_equal(bal.compute(init),
               mean(abs(unlist(lapply(baltab, function(i) i$Balance$Diff.Un)))),
               tolerance = eps)
  
  expect_equal(bal.compute(init, weights = w),
               mean(abs(unlist(lapply(baltab, function(i) i$Balance$Diff.Adj)))),
               tolerance = eps)
  
  ## KS
  init <- bal.init(lalonde[cov_names], tb, "ks.mean")
  expect_s3_class(init, "bal.init")
  
  baltab <- bal.tab(lalonde[cov_names], treat = tb, binary = "std",
                    stats = "ks", estimand = "ATE", weights = w,
                    un = TRUE)$Balance
  
  expect_equal(bal.compute(init),
               mean(abs(baltab$KS.Un)),
               tolerance = eps)
  
  expect_equal(bal.compute(init, weights = w),
               mean(abs(baltab$KS.Adj)),
               tolerance = eps)
  
  init <- bal.init(lalonde[cov_names], treat = tb,
                   stat = "ks.mean",
                   pairwise = FALSE)
  expect_s3_class(init, "bal.init")
  
  baltab <- bal.tab(lalonde[cov_names], treat = tb, binary = "std",
                    stats = "ks", estimand = "ATE", weights = w,
                    un = TRUE, pairwise = FALSE)$Pair.Balance
  
  expect_equal(bal.compute(init),
               mean(abs(unlist(lapply(baltab, function(i) i$Balance$KS.Un)))),
               tolerance = eps)
  
  expect_equal(bal.compute(init, weights = w),
               mean(abs(unlist(lapply(baltab, function(i) i$Balance$KS.Adj)))),
               tolerance = eps)
  
  # Multi-category treatment
  ## SMD
  init <- bal.init(lalonde[cov_names], tm, "smd.mean")
  expect_s3_class(init, "bal.init")
  
  baltab <- bal.tab(lalonde[cov_names], treat = tm, binary = "std",
                    stats = "m", estimand = "ATE", weights = w,
                    un = TRUE)$Pair.Balance
  
  expect_equal(bal.compute(init),
               mean(abs(unlist(lapply(baltab, function(i) i$Balance$Diff.Un)))),
               tolerance = eps)
  
  expect_equal(bal.compute(init, weights = w),
               mean(abs(unlist(lapply(baltab, function(i) i$Balance$Diff.Adj)))),
               tolerance = eps)
  
  init <- bal.init(lalonde[cov_names], treat = tm,
                   stat = "smd.mean",
                   pairwise = FALSE)
  expect_s3_class(init, "bal.init")
  
  baltab <- bal.tab(lalonde[cov_names], treat = tm, binary = "std",
                    stats = "m", estimand = "ATE", weights = w,
                    un = TRUE, pairwise = FALSE)$Pair.Balance
  
  expect_equal(bal.compute(init),
               mean(abs(unlist(lapply(baltab, function(i) i$Balance$Diff.Un)))),
               tolerance = eps)
  
  expect_equal(bal.compute(init, weights = w),
               mean(abs(unlist(lapply(baltab, function(i) i$Balance$Diff.Adj)))),
               tolerance = eps)
  
  ## KS
  init <- bal.init(lalonde[cov_names], tm, "ks.mean")
  expect_s3_class(init, "bal.init")
  
  baltab <- bal.tab(lalonde[cov_names], treat = tm, binary = "std",
                    stats = "ks", estimand = "ATE", weights = w,
                    un = TRUE)$Pair.Balance
  
  expect_equal(bal.compute(init),
               mean(abs(unlist(lapply(baltab, function(i) i$Balance$KS.Un)))),
               tolerance = eps)
  
  expect_equal(bal.compute(init, weights = w),
               mean(abs(unlist(lapply(baltab, function(i) i$Balance$KS.Adj)))),
               tolerance = eps)
  
  init <- bal.init(lalonde[cov_names], treat = tm,
                   stat = "ks.mean",
                   pairwise = FALSE)
  expect_s3_class(init, "bal.init")
  
  baltab <- bal.tab(lalonde[cov_names], treat = tm, binary = "std",
                    stats = "ks", estimand = "ATE", weights = w,
                    un = TRUE, pairwise = FALSE)$Pair.Balance
  
  expect_equal(bal.compute(init),
               mean(abs(unlist(lapply(baltab, function(i) i$Balance$KS.Un)))),
               tolerance = eps)
  
  expect_equal(bal.compute(init, weights = w),
               mean(abs(unlist(lapply(baltab, function(i) i$Balance$KS.Adj)))),
               tolerance = eps)
})
test_that("bal.compute() with a matrix of weights matches computing each column", {
  skip_if_not_installed("patrick")

  f <- bc_fixture()

  W <- cbind(f$w, runif(length(f$w)), 1)
  W[1:20, 2L] <- 0

  treats <- list(binary = f$tb, multinomial = f$tm, continuous = f$tc, target = NULL)

  #`energy.dist`, `kernel.dist`, and `distance.cov`/`distance.cor` compute all columns
  #at once and the rest loop over them; either way each value must be the one the
  #column gives on its own.
  patrick::with_parameters_test_that(
    "{treat.type}: {stat}",
    {
      #`tm` is race, so race is dropped from the covariates it is balanced on
      covs <- {
        if (treat.type == "multinomial") f$covs[names(f$covs) != "race"]
        else f$covs
      }

      init <- bal.init(covs, treats[[treat.type]], stat, s.weights = f$sw)

      expect_equal(bal.compute(init, weights = W),
                   apply(W, 2L, function(w) bal.compute(init, weights = w)),
                   tolerance = 1e-10)
    },
    .cases = stat_cases(c("binary", "multinomial", "continuous", "target"))
  )

  patrick::with_parameters_test_that(
    "multinomial energy.dist: {config}",
    {
      init <- do.call("bal.init", c(list(f$num_covs, f$tm, "energy.dist"), args))

      expect_equal(bal.compute(init, weights = W),
                   apply(W, 2L, function(w) bal.compute(init, weights = w)),
                   tolerance = 1e-10)
    },
    patrick::cases(
      list(config = "estimand = NULL", args = list(estimand = NULL)),
      list(config = "ATE, improved = FALSE", args = list(estimand = "ATE", improved = FALSE)),
      list(config = "ATT, focal = black", args = list(estimand = "ATT", focal = "black"))
    )
  )
})

test_that("bal.compute() names matrix results by column and accepts a single column", {
  f <- bc_fixture()

  init <- bal.init(f$covs, f$tb, "energy.dist")

  W <- cbind(a = f$w, b = 1)

  out <- bal.compute(init, weights = W)
  expect_named(out, c("a", "b"))
  expect_length(out, 2L)

  expect_equal(unname(bal.compute(init, weights = W[, "a", drop = FALSE])),
               bal.compute(init, weights = f$w))

  expect_null(names(bal.compute(init, weights = unname(W))))
})

test_that("bal.compute() rejects invalid weights", {
  f <- bc_fixture()
  n <- length(f$w)

  init <- bal.init(f$covs, f$tb, "smd.mean")

  expect_err(bal.compute(init, weights = f$w[-1]), "one value for each unit")
  expect_err(bal.compute(init, weights = cbind(f$w, f$w)[-1, ]), "one row for each unit")
  expect_err(bal.compute(init, weights = matrix(numeric(0), nrow = n)), "at least one column")
  expect_err(bal.compute(init, weights = replace(f$w, 3, NA)), "NA")
  expect_err(bal.compute(init, weights = replace(f$w, 3, Inf)), "infinite")
  expect_err(bal.compute(init, weights = as.character(f$w)), "must be numeric")
  expect_err(bal.compute(init, weights = data.frame(f$w)), "vector or a matrix")
  expect_err(bal.compute(init, weights = array(1, c(n, 2, 2))), "vector or a matrix")

  #The target statistics have the same number of units as the covariates
  init <- bal.init(f$covs, stat = "energy.dist")
  expect_err(bal.compute(init, weights = f$w[-1]), "one value for each unit")
})

test_that("energy.dist, kernel.dist, and distance.cor match their definitions", {
  f <- bc_fixture()

  #All-continuous covariates with no sampling weights, so the scaling `bal.init()`
  #applies is the usual standardization
  X <- scale(as.matrix(f$num_covs))
  D <- as.matrix(dist(X))
  n <- nrow(D)

  w <- f$w
  a1 <- w * (f$tb == 1) / sum(w[f$tb == 1])
  a0 <- w * (f$tb == 0) / sum(w[f$tb == 0])
  s <- rep(1 / n, n)

  edist <- function(u, v) -drop(t(u - v) %*% D %*% (u - v))

  expect_equal(bal.compute(bal.init(f$num_covs, f$tb, "energy.dist"), weights = w),
               edist(a1, a0))

  expect_equal(bal.compute(bal.init(f$num_covs, f$tb, "energy.dist", estimand = "ATE"),
                           weights = w),
               edist(a1, s) + edist(a0, s) + edist(a1, a0))

  expect_equal(bal.compute(bal.init(f$num_covs, f$tb, "energy.dist", estimand = "ATT",
                                    improved = FALSE),
                           weights = w),
               edist(a1, s * (f$tb == 1) * n / sum(f$tb == 1)) +
                 edist(a0, s * (f$tb == 1) * n / sum(f$tb == 1)))

  expect_equal(bal.compute(bal.init(f$num_covs, stat = "energy.dist"), weights = w),
               edist(w / sum(w), s))

  #Kernel distance: Gaussian kernel with the median distance as the bandwidth
  K <- exp(-D^2 / median(D))
  tw <- a1 - a0
  expect_equal(bal.compute(bal.init(f$num_covs, f$tb, "kernel.dist"), weights = w),
               sqrt(drop(t(tw) %*% K %*% tw)))

  #Distance correlation: double-centered distance matrices of the covariates and of
  #the treatment, with the weights in the outer product
  H <- diag(n) - 1 / n
  XX <- H %*% D %*% H
  AA <- H %*% as.matrix(dist(f$tc)) %*% H
  wn <- w / sum(w)
  dcor2 <- drop(t(wn) %*% (XX * AA) %*% wn) /
    sqrt(drop(t(s) %*% (XX^2) %*% s) * drop(t(s) %*% (AA^2) %*% s))

  expect_equal(bal.compute(bal.init(f$num_covs, f$tc, "distance.cor"), weights = w),
               sqrt(dcor2))
})

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
  eps <- if (capabilities("long.double")) 1e-8 else 1e-1

  f <- bc_fixture()

  #`binary = "std"` so every row is on the same scale as the `smd.*` statistics,
  #and `estimand = "ATE"` so `bal.init()` and `bal.tab()` use the same denominator.
  bt <- bal.tab(f$covs, treat = f$tb, binary = "std", estimand = "ATE",
                weights = f$w, un = TRUE,
                stats = c("mean.diffs", "ks.statistics", "ovl.coefficients"))$Balance

  columns <- c(smd = "Diff.Adj", ks = "KS.Adj", ovl = "OVL.Adj")

  for (family in names(columns)) {
    v <- bt[[columns[[family]]]]

    #`mean` and `max` are taken over absolute values; `rms` squares first, so the
    #sign drops out on its own.
    expected <- c(mean = mean(abs(v)), max = max(abs(v)), rms = sqrt(mean(v^2)))

    for (agg in names(expected)) {
      stat <- paste0(family, ".", agg)

      expect_equal(bal.compute(bal.init(f$covs, f$tb, stat, estimand = "ATE"),
                               weights = f$w),
                   expected[[agg]], tolerance = eps, info = stat)
    }
  }
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

  for (family in names(columns)) {
    v <- bt[[columns[[family]]]]
    expected <- c(mean = mean(abs(v)), max = max(abs(v)), rms = sqrt(mean(v^2)))

    for (agg in names(expected)) {
      stat <- paste0(family, ".", agg)

      expect_equal(bal.compute(bal.init(covs, f$tc, stat), weights = f$w),
                   expected[[agg]], tolerance = eps, info = stat)
    }
  }

  #The Spearman statistics must actually differ from the Pearson ones -- they were
  #once computed identically because of a bad loop index, which no shape or class
  #check would have caught.
  expect_false(isTRUE(all.equal(bal.compute(bal.init(covs, f$tc, "s.mean"), weights = f$w),
                                bal.compute(bal.init(covs, f$tc, "p.mean"), weights = f$w))))
})

test_that("multi-category aggregators summarize the pairwise bal.tab() columns", {
  eps <- if (capabilities("long.double")) 1e-8 else 1e-1

  f <- bc_fixture()

  covs <- f$covs[c("age", "educ", "married", "nodegree", "re74", "re75")]

  bt <- bal.tab(covs, treat = f$tm, binary = "std", estimand = "ATE",
                weights = f$w, un = TRUE,
                stats = c("mean.diffs", "ks.statistics"))$Pair.Balance

  columns <- c(smd = "Diff.Adj", ks = "KS.Adj")

  for (family in names(columns)) {
    #Pairwise statistics pool every pair's covariates into one vector before
    #aggregating, so `max` is over all pairs at once, not per pair.
    v <- unlist(lapply(bt, function(i) i$Balance[[columns[[family]]]]))
    expected <- c(mean = mean(abs(v)), max = max(abs(v)), rms = sqrt(mean(v^2)))

    for (agg in names(expected)) {
      stat <- paste0(family, ".", agg)

      expect_equal(bal.compute(bal.init(covs, f$tm, stat, estimand = "ATE"),
                               weights = f$w),
                   expected[[agg]], tolerance = eps, info = stat)
    }
  }
})

test_that("bal.compute() is unchanged by transformations the statistics ignore", {
  f <- bc_fixture()

  #Weights enter every statistic as a ratio, so a constant multiple cannot matter.
  #A refactor that starts normalizing by `sum(w)` in one place but not another
  #shows up here.
  for (stat in c("smd.mean", "ks.max", "ovl.mean", "energy.dist", "mahalanobis")) {
    init <- bal.init(f$covs, f$tb, stat)

    expect_equal(bal.compute(init, weights = f$w),
                 bal.compute(init, weights = f$w * 7),
                 info = stat)
  }

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

  #`print()` must name the statistic; the phrase table behind it has been missing
  #entries before, which produced an empty line rather than an error.
  for (stat in available.stats("binary")) {
    out <- capture.output(print(bal.init(f$num_covs, f$tb, stat)))

    expect_true(any(nzchar(trimws(out))), info = stat)
    expect_false(any(grepl("NA", out, fixed = TRUE)), info = stat)
  }

  for (stat in available.stats("continuous")) {
    out <- capture.output(print(bal.init(f$num_covs, f$tc, stat)))

    expect_true(any(nzchar(trimws(out))), info = stat)
    expect_false(any(grepl("NA", out, fixed = TRUE)), info = stat)
  }
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
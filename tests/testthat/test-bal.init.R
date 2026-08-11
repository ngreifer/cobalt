#Systematic coverage of `bal.init()`/`bal.compute()` across every statistic in
#`available.stats()` and every treatment type, plus `print.bal.init()` and the
#per-statistic `...` options. Complements `test-bal.compute.R`, which pins a few
#specific numeric values.
#
#All of these are pure functions of a covariate matrix and a treatment vector, so
#no package from Suggests is involved.

t_bin <- function() lalonde$treat
t_multi <- function() lalonde$race
t_cont <- function() lalonde$re75

covs_num <- function() lalonde[c("age", "educ", "re74")]

test_that("available.stats() reports the documented statistics per treatment type", {
  expect_type(available.stats("binary"), "character")

  for (tt in c("binary", "multinomial", "continuous", "target")) {
    expect_gt(length(available.stats(tt)), 0L)
  }

  #"multi-category" is accepted as an alias for "multinomial".
  expect_identical(available.stats("multi-category"), available.stats("multinomial"))

  expect_err(available.stats("bogus"), "`treat.type` should be one of")
})

test_that("every binary statistic initialises and computes a finite value", {
  x <- covs_num()
  treat <- t_bin()

  for (s in available.stats("binary")) {
    init <- bal.init(x, treat = treat, stat = s)
    expect_s3_class(init, "bal.init")
    expect_identical(attr(init, "stat"), s)
    expect_identical(attr(init, "treat.type"), "binary")

    val <- bal.compute(init)
    expect_true(is.numeric(val) && length(val) == 1L && is.finite(val),
                label = sprintf("bal.compute(%s)", s))

    #Supplying weights must also work and stay finite.
    val_w <- bal.compute(init, weights = w_fixed)
    expect_true(is.finite(val_w), label = sprintf("bal.compute(%s, weights=)", s))
  }
})

test_that("every multi-category statistic initialises and computes a finite value", {
  x <- covs_num()
  treat <- t_multi()

  for (s in available.stats("multinomial")) {
    init <- bal.init(x, treat = treat, stat = s)
    expect_identical(attr(init, "treat.type"), "multinomial")
    expect_true(is.finite(bal.compute(init)), label = s)
    expect_true(is.finite(bal.compute(init, weights = w_fixed)), label = s)
  }
})

test_that("every continuous statistic initialises and computes a finite value", {
  x <- covs_num()
  treat <- t_cont()

  for (s in available.stats("continuous")) {
    init <- bal.init(x, treat = treat, stat = s)
    expect_identical(attr(init, "treat.type"), "continuous")
    expect_true(is.finite(bal.compute(init)), label = s)
    expect_true(is.finite(bal.compute(init, weights = w_fixed)), label = s)
  }
})

test_that("every target statistic initialises and computes a finite value", {
  x <- covs_num()

  for (s in available.stats("target")) {
    #`treat = NULL` selects target balance.
    init <- bal.init(x, treat = NULL, stat = s)
    expect_identical(attr(init, "treat.type"), "target")
    expect_true(is.finite(bal.compute(init)), label = s)
    expect_true(is.finite(bal.compute(init, weights = w_fixed)), label = s)
  }
})

test_that("bal.compute() dispatches on a raw covariate set as well as a bal.init", {
  x <- covs_num()

  direct <- bal.compute(x, treat = t_bin(), stat = "smd.mean")
  viainit <- bal.compute(bal.init(x, treat = t_bin(), stat = "smd.mean"))

  expect_equal(direct, viainit)

  #`s.weights` are honored on both paths.
  d_sw <- bal.compute(x, treat = t_bin(), stat = "smd.mean", s.weights = sw_fixed)
  i_sw <- bal.compute(bal.init(x, treat = t_bin(), stat = "smd.mean",
                               s.weights = sw_fixed))
  expect_equal(d_sw, i_sw)
  expect_false(isTRUE(all.equal(direct, d_sw)))
})

test_that("print.bal.init() reports the treatment type and statistic", {
  x <- covs_num()

  #Every statistic must have a phrase; `r2.2` and `r2.3` were once missing.
  for (s in available.stats("binary")) {
    out <- squish(capture.output(print(bal.init(x, treat = t_bin(), stat = s))))
    expect_match(out, "A `bal.init` object", fixed = TRUE)
    expect_match(out, "treatment type: binary", fixed = TRUE)
    expect_match(out, sprintf("statistic: %s", s), fixed = TRUE)
  }

  for (s in available.stats("continuous")) {
    expect_no_error(capture.output(print(bal.init(x, treat = t_cont(), stat = s))))
  }
})

test_that("the r2 family accepts poly and int options", {
  x <- covs_num()

  #`r2` uses the covariates as given; `r2.2`/`r2.3` add squares and cubes.
  r2 <- bal.compute(x, treat = t_bin(), stat = "r2")
  r22 <- bal.compute(x, treat = t_bin(), stat = "r2.2")
  r23 <- bal.compute(x, treat = t_bin(), stat = "r2.3")

  for (v in list(r2, r22, r23)) {
    expect_true(is.finite(v) && v >= 0 && v <= 1)
  }

  #Adding terms cannot reduce the model R-squared.
  expect_gte(r22, r2)

  #`r2` with explicit `poly`/`int` reproduces the shortcuts.
  expect_equal(bal.compute(bal.init(x, treat = t_bin(), stat = "r2", poly = 2)), r22)
  expect_equal(bal.compute(bal.init(x, treat = t_bin(), stat = "r2", poly = 3)), r23)

  expect_true(is.finite(bal.compute(bal.init(x, treat = t_bin(), stat = "r2",
                                             int = TRUE))))
})

test_that("smd, ks, and ovl accept estimand, focal, and pairwise", {
  x <- covs_num()

  for (s in c("smd.mean", "ks.max", "ovl.mean")) {
    for (e in c("ATE", "ATT", "ATC")) {
      expect_true(is.finite(bal.compute(bal.init(x, treat = t_bin(), stat = s,
                                                 estimand = e))), label = paste(s, e))
    }
  }

  #For multi-category treatments, `focal` selects the reference group and
  #`pairwise = FALSE` compares each group to the target.
  for (s in c("smd.mean", "ks.mean", "ovl.mean")) {
    expect_true(is.finite(bal.compute(bal.init(x, treat = t_multi(), stat = s,
                                               estimand = "ATT", focal = "white"))))
    expect_true(is.finite(bal.compute(bal.init(x, treat = t_multi(), stat = s,
                                               pairwise = FALSE))))
  }
})

test_that("ovl accepts integrate and steps", {
  x <- covs_num()

  v_int <- bal.compute(bal.init(x, treat = t_bin(), stat = "ovl.mean",
                                integrate = TRUE))
  v_riem <- bal.compute(bal.init(x, treat = t_bin(), stat = "ovl.mean",
                                 integrate = FALSE))

  expect_true(is.finite(v_int))
  expect_true(is.finite(v_riem))

  #The two methods approximate the same integral.
  expect_equal(v_int, v_riem, tolerance = 1e-2)
})

test_that("energy.dist accepts improved and estimand; distance.cov accepts std", {
  x <- covs_num()

  for (imp in c(TRUE, FALSE)) {
    expect_true(is.finite(bal.compute(bal.init(x, treat = t_bin(),
                                               stat = "energy.dist",
                                               improved = imp))))
  }

  expect_true(is.finite(bal.compute(bal.init(x, treat = t_cont(),
                                             stat = "distance.cov"))))
  expect_true(is.finite(bal.compute(bal.init(x, treat = t_cont(),
                                             stat = "distance.cor"))))
})

test_that("l1.med accepts its binning options", {
  x <- covs_num()

  expect_true(is.finite(bal.compute(bal.init(x, treat = t_bin(), stat = "l1.med"))))
  expect_true(is.finite(bal.compute(bal.init(x, treat = t_bin(), stat = "l1.med",
                                             l1.min.bin = 2, l1.max.bin = 6,
                                             l1.n = 10))))
})

test_that("factor covariates are split before computation", {
  #A data frame with a factor gives the same answer as the pre-split matrix.
  with_factor <- lalonde[c("age", "race")]
  pre_split <- splitfactor(with_factor, drop.first = "if2")

  expect_equal(bal.compute(with_factor, treat = t_bin(), stat = "smd.mean"),
               bal.compute(as.matrix(pre_split), treat = t_bin(), stat = "smd.mean"))
})

test_that("statistics are rejected for the wrong treatment type", {
  x <- covs_num()

  #`match_arg()` reports the statistics that *are* available.
  expect_err(bal.init(x, treat = t_cont(), stat = "smd.mean"),
             "`stat` should be one of")
  expect_err(bal.init(x, treat = t_bin(), stat = "p.mean"),
             "`stat` should be one of")
  expect_err(bal.init(x, treat = t_multi(), stat = "kernel.dist"),
             "`stat` should be one of")
  expect_err(bal.init(x, treat = t_multi(), stat = "r2"),
             "`stat` should be one of")
  expect_err(bal.init(x, treat = t_bin(), stat = "bogus"),
             "`stat` should be one of")
})

test_that("statistics that cannot handle missing values say so", {
  x_mis <- lalonde_mis[c("age", "re74")]

  for (s in c("mahalanobis", "energy.dist", "kernel.dist", "l1.med", "r2.2")) {
    expect_err(bal.init(x_mis, treat = t_bin(), stat = s),
               "cannot be used when there are missing values in the covariates")
  }
})

test_that("bal.init() validates its covariates and treatment", {
  x <- covs_num()

  expect_err(bal.init(array(1, c(4, 2, 2)), treat = rep(0:1, 2), stat = "smd.mean"),
             "must be a data frame or numeric matrix")
  expect_err(bal.init(x, treat = rep(1, nrow(lalonde)), stat = "smd.mean"),
             "treatment must have at least two unique values")
  expect_err(bal.init(x, treat = t_bin()[-1L], stat = "smd.mean"),
             "must have the same number of units")
  expect_err(bal.init(x, treat = list(1, 2), stat = "smd.mean"),
             "`treat` must be")
})

test_that("bal.compute() validates weights against the init object", {
  init <- bal.init(covs_num(), treat = t_bin(), stat = "smd.mean")

  expect_err(bal.compute(init, weights = w_fixed[-1L]),
             "must have the same number of units")
})

# ---------------------------------------------------------------------------
# The remaining `...` options, and the `init = NULL` fallback that every entry of
# the statistic registry carries.

test_that("each statistic's function builds its own init when not given one", {
  #`bal.init()` always constructs the init and stores the computing function on
  #the object, so `bal.compute()` never exercises the fallback inside each
  #registry entry. Calling it directly must give the same answer.
  x <- covs_num()

  treats <- list(binary = t_bin(), multinomial = t_multi(), continuous = t_cont(),
                 target = NULL)

  for (tt in names(treats)) {
    treat <- treats[[tt]]

    for (s in available.stats(tt)) {
      init <- bal.init(x, treat = treat, stat = s)
      fun <- attr(init, "fun")

      #`l1.med` bins the covariates at random, so both calls need the same seed
      #for their results to be comparable at all.
      set.seed(11L)
      via_init <- fun(init = init, weights = NULL)

      #Without it, so the entry builds its own.
      set.seed(11L)
      built <- fun(covs = x, treat = treat, weights = NULL)

      expect_equal(built, via_init,
                   info = sprintf("%s / %s", tt, s))
    }
  }
})

test_that("s.weights are honored by every statistic, not just smd", {
  x <- covs_num()

  for (s in c("ks.max", "ovl.mean", "energy.dist", "mahalanobis", "r2.2")) {
    plain <- bal.compute(x, treat = t_bin(), stat = s)
    weighted <- bal.compute(x, treat = t_bin(), stat = s, s.weights = sw_fixed)
    expect_true(is.finite(weighted), label = s)
    expect_false(isTRUE(all.equal(plain, weighted)), label = s)
  }

  for (s in c("p.mean", "s.max", "distance.cov")) {
    expect_true(is.finite(bal.compute(x, treat = t_cont(), stat = s,
                                      s.weights = sw_fixed)), label = s)
  }
})

test_that("focal applies to binary treatments and estimand to energy.dist", {
  x <- covs_num()

  #`focal` with a binary treatment names the group to compare against.
  for (f in as.character(unique(lalonde$treat))) {
    expect_true(is.finite(bal.compute(bal.init(x, treat = t_bin(),
                                              stat = "smd.mean",
                                              estimand = "ATT", focal = f))),
                label = f)
  }

  for (e in c("ATE", "ATT", "ATC")) {
    expect_true(is.finite(bal.compute(bal.init(x, treat = t_bin(),
                                              stat = "energy.dist",
                                              estimand = e))), label = e)
  }
})

test_that("a multi-category ATC is an ATT and requires focal", {
  x <- covs_num()

  #With more than two groups there is no single control group, so an ATC names its
  #reference group through `focal` exactly as an ATT does, and both insist on one.
  for (e in c("ATT", "ATC")) {
    expect_err(bal.init(x, treat = t_multi(), stat = "smd.mean", estimand = e),
               sprintf('estimand = "%s"', e))
  }

  att <- bal.init(x, treat = t_multi(), stat = "smd.mean",
                  estimand = "ATT", focal = "white")
  atc <- bal.init(x, treat = t_multi(), stat = "smd.mean",
                  estimand = "ATC", focal = "white")

  expect_identical(atc, att)
})

test_that("distance.cov honors std and matches distance.cor", {
  x <- covs_num()

  raw <- bal.compute(bal.init(x, treat = t_cont(), stat = "distance.cov",
                              std = FALSE))
  std <- bal.compute(bal.init(x, treat = t_cont(), stat = "distance.cov",
                              std = TRUE))
  expect_true(is.finite(raw))
  expect_true(is.finite(std))
  expect_false(isTRUE(all.equal(raw, std)))

  #`distance.cor` is the standardized version.
  expect_equal(bal.compute(x, treat = t_cont(), stat = "distance.cor"), std)
})

test_that("l1.med accepts its remaining options", {
  x <- covs_num()

  expect_true(is.finite(bal.compute(bal.init(x, treat = t_bin(), stat = "l1.med",
                                            l1.n = 5))))
  #`.covs` overrides which covariates the binning uses.
  expect_true(is.finite(bal.compute(bal.init(x, treat = t_bin(), stat = "l1.med",
                                            .covs = x[1L]))))
})

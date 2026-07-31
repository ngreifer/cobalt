#`ps` and `iptw` objects fitted with twang's xgboost engine
#(`version = "xgboost"` in the call to `twang::ps()` or `twang::iptw()`).
#
#These fits store the model differently from the default gbm engine, and cobalt
#branches on that difference in `x2base.ps()` and `x2base.iptw()`:
#
#  1. a gbm fit records `gbm.obj$var.names`, the covariates as the user named them;
#  2. an xgboost fit has no `var.names`, only `gbm.obj$feature_names`, which are the
#     columns of the *model matrix* -- so a factor covariate appears as one one-hot
#     column per level;
#  3. cobalt therefore falls back to `feature_names`, but only when every one of
#     them is still a column of the stored data. When a factor was expanded they
#     are not, and cobalt cannot recover the original covariates, so it asks for
#     `formula`/`covs` (or `formula.list`/`covs.list`) instead.
#
#Two kinds of test below. The first group fits real xgboost models and skips if
#that is not possible on this machine. The second group takes a real gbm fit and
#replaces only `gbm.obj` with an xgboost-shaped one, which reproduces exactly the
#condition cobalt branches on and so covers branches 2 and 3 regardless of whether
#twang and xgboost are currently compatible. That matters: twang 2.6.2 calls
#`xgboost(data=, label=, params=, feval=)`, all of which xgboost 3.x removed or
#renamed, so the real fits skip against a current xgboost.

skip_on_cran()

#Give `gbm.obj` the shape of a fitted xgboost model: no `var.names`, and
#`feature_names` holding model-matrix column names.
#
#This is the layout cobalt reads, and the one xgboost used when twang's xgboost
#support was written. xgboost 3.x no longer stores feature names in the booster
#list -- `m$feature_names` is NULL and they are retrieved with
#`getinfo(m, "feature_name")` -- so `feature_names = NULL` is tested separately
#below as the case a current xgboost would present.
as_xgboost_obj <- function(feature_names = NULL) {
  structure(if (is.null(feature_names)) list()
            else list(feature_names = feature_names),
            class = "xgb.Booster")
}

# ---- real xgboost fits ------------------------------------------------------

test_that("bal.tab() works on a real xgboost ps fit with numeric covariates", {
  ps <- fx("ps_xgboost_num")

  #Every feature name is still a column of the data, so the covariates are
  #recovered without help.
  b <- bal.tab(ps)
  expect_s3_class(b, "bal.tab.bin")
  expect_true(all(c("age", "educ", "re74") %in% rownames(b$Balance)))

  expect_length(get.w(ps), n_lalonde)
  expect_no_error(capture.output(print(b)))
})

test_that("a real xgboost ps fit with a factor covariate needs formula or covs", {
  ps <- fx("ps_xgboost")

  expect_err(bal.tab(ps),
             'when `version = "xgboost"` in the call to `ps()` and any variables are categorical')

  b <- bal.tab(ps, formula = treat ~ age + educ + race, data = lalonde)
  expect_s3_class(b, "bal.tab.bin")
  expect_true("race_black" %in% rownames(b$Balance))
})

test_that("bal.tab() works on a real xgboost iptw fit given covs.list", {
  ip <- fx("iptw_xgboost")

  b <- bal.tab(ip, covs.list = list(lalonde[c("age", "race")],
                                    lalonde[c("age", "race", "treat")]))
  expect_s3_class(b, "bal.tab.msm")
  expect_length(b$Time.Balance, 2L)
})

# ---- xgboost-shaped objects derived from a real gbm fit ---------------------

test_that("ps: feature_names are used when they are all columns of the data", {
  x <- fx("ps")
  x$gbm.obj <- as_xgboost_obj(c("age", "educ"))

  b <- bal.tab(x)
  expect_s3_class(b, "bal.tab.bin")

  #Only the two named features are reported, alongside one propensity score per
  #requested stop method.
  covs_shown <- grep("^prop.score", rownames(b$Balance), invert = TRUE, value = TRUE)
  expect_setequal(covs_shown, c("age", "educ"))

  expect_no_error(capture.output(print(b)))
})

test_that("ps: expanded feature_names require formula or covs", {
  x <- fx("ps")

  #What xgboost would produce for `race`: one column per level, none of which is
  #a column of the stored data.
  x$gbm.obj <- as_xgboost_obj(c("age", "educ", "raceblack", "racehispan",
                               "racewhite"))

  expect_err(bal.tab(x),
             'when `version = "xgboost"` in the call to `ps()` and any variables are categorical')

  #`formula` recovers the original covariates, and cobalt splits the factor itself.
  b_f <- bal.tab(x, formula = treat ~ age + educ + race, data = lalonde)
  expect_true(all(c("age", "educ", "race_black", "race_hispan", "race_white") %in%
                    rownames(b_f$Balance)))

  #`covs` is equivalent.
  b_c <- bal.tab(x, covs = lalonde[c("age", "educ", "race")])
  expect_equal(b_f$Balance, b_c$Balance)

  #The weights do not depend on the covariate names, so `get.w()` is unaffected.
  expect_identical(get.w(x), get.w(fx("ps")))
})

test_that("ps: the other cobalt entry points work once covariates are supplied", {
  x <- fx("ps")
  x$gbm.obj <- as_xgboost_obj(c("age", "educ", "raceblack", "racehispan",
                               "racewhite"))

  f <- treat ~ age + educ + race

  expect_s3_class(love.plot(x, formula = f, data = lalonde, stats = "mean.diffs",
                            binary = "std"),
                  "love.plot")

  p <- bal.plot(x, formula = f, data = lalonde, var.name = "age")
  expect_no_condition(ggplot2::ggplot_build(p))

  p <- bal.plot(x, formula = f, data = lalonde, var.name = "race")
  expect_no_condition(ggplot2::ggplot_build(p))
})

test_that("ps: a real gbm fit still uses var.names", {
  x <- fx("ps")

  #The control for the tests above: with `var.names` present, no formula or covs
  #is needed even though the covariates include a factor.
  expect_true(is_not_null(x$gbm.obj$var.names))
  b <- bal.tab(x)
  expect_true("race_black" %in% rownames(b$Balance))
})

test_that("iptw: feature_names are used when they are all columns of the data", {
  x <- fx("iptw_fac")

  for (i in seq_along(x$psList)) {
    x$psList[[i]]$gbm.obj <- as_xgboost_obj("age")
  }

  b <- bal.tab(x)
  expect_s3_class(b, "bal.tab.msm")
  expect_length(b$Time.Balance, 2L)
  expect_setequal(setdiff(rownames(b$Time.Balance[[1L]]$Balance), "prop.score"),
                  "age")
})

test_that("iptw: expanded feature_names require formula.list or covs.list", {
  x <- fx("iptw_fac")

  for (i in seq_along(x$psList)) {
    x$psList[[i]]$gbm.obj <- as_xgboost_obj(c("age", "raceblack", "racehispan",
                                             "racewhite"))
  }

  expect_err(bal.tab(x),
             'when `version = "xgboost"` in the call to `iptw()` and any variables are categorical')

  fl <- list(treat ~ age + race, nodegree ~ age + race + treat)

  b_f <- bal.tab(x, formula.list = fl, data = lalonde)
  expect_s3_class(b_f, "bal.tab.msm")
  expect_true(all(c("race_black", "race_hispan", "race_white") %in%
                    rownames(b_f$Time.Balance[[1L]]$Balance)))

  b_c <- bal.tab(x, covs.list = list(lalonde[c("age", "race")],
                                     lalonde[c("age", "race", "treat")]))
  expect_equal(b_f$Time.Balance[[1L]]$Balance, b_c$Time.Balance[[1L]]$Balance)
})

test_that("iptw: formula.list is validated", {
  x <- fx("iptw_fac")

  for (i in seq_along(x$psList)) {
    x$psList[[i]]$gbm.obj <- as_xgboost_obj(c("age", "raceblack"))
  }

  expect_err(bal.tab(x, formula.list = list(treat ~ age)),
             "must have as many entries as time points")
  expect_err(bal.tab(x, formula.list = list(1, 2)),
             "must be a list of formulas identifying the covariates")

  #`covs.list` is validated the same way.
  expect_err(bal.tab(x, covs.list = list(lalonde["age"])),
             "must have")
})

test_that("iptw: a real gbm fit still uses var.names", {
  x <- fx("iptw_fac")

  expect_true(is_not_null(x$psList[[1L]]$gbm.obj$var.names))
  b <- bal.tab(x)
  expect_true("race_black" %in% rownames(b$Time.Balance[[1L]]$Balance))
})

test_that("a model object recording no feature names asks for formula or covs", {
  #`all()` of an empty set is TRUE, so a `gbm.obj` with neither `var.names` nor
  #`feature_names` used to enter the feature-names branch and fail inside
  #`reformulate(NULL)` with "'termlabels' must be a character vector". It must
  #reach the informative error instead. This is the shape an xgboost 3.x booster
  #presents, since it no longer stores feature names in the object's list.
  x <- fx("ps")
  x$gbm.obj <- as_xgboost_obj()

  expect_err(bal.tab(x),
             'when `version = "xgboost"` in the call to `ps()` and any variables are categorical')

  #Supplying the covariates still works.
  b <- bal.tab(x, formula = treat ~ age + educ + race, data = lalonde)
  expect_true("race_black" %in% rownames(b$Balance))

  xi <- fx("iptw_fac")
  for (i in seq_along(xi$psList)) {
    xi$psList[[i]]$gbm.obj <- as_xgboost_obj()
  }

  expect_err(bal.tab(xi),
             'when `version = "xgboost"` in the call to `iptw()` and any variables are categorical')

  b <- bal.tab(xi, covs.list = list(lalonde[c("age", "race")],
                                    lalonde[c("age", "race", "treat")]))
  expect_s3_class(b, "bal.tab.msm")
})

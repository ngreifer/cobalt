#Argument handling in the `get.w()` methods. `test-input-guards.R` already calls
#`get.w()` once per supported class with no arguments, which reaches every method
#but none of their options.

skip_on_cran()

test_that("s.weights = TRUE multiplies the sampling weights back in", {
  w <- fx("weightit_sw")

  plain <- get.w(w)
  with_sw <- get.w(w, s.weights = TRUE)

  expect_false(isTRUE(all.equal(plain, with_sw)))
  expect_equal(with_sw, plain * w$s.weights)
})

test_that("get.w.ps() accepts a vector estimand", {
  ps <- fx("ps")

  #One estimand per requested stop method.
  w2 <- get.w(ps, estimand = c("ATE", "ATT"))
  expect_s3_class(w2, "data.frame")
  expect_identical(ncol(w2), 2L)

  #The two columns differ, since the estimands differ.
  expect_false(isTRUE(all.equal(w2[[1L]], w2[[2L]])))

  #More estimands than stop methods keeps the leading ones.
  expect_s3_class(get.w(ps, estimand = c("ATE", "ATT", "ATC")), "data.frame")

  #An estimand vector shorter than the stop methods but longer than 1 is an
  #error; that needs at least three stop methods, so it is not exercised here.

  expect_type(get.w(ps, s.weights = TRUE), "list")
})

test_that("get.w.mnps() handles the ATT estimand and several stop methods", {
  m <- fx("mnps_att")

  w <- get.w(m)
  expect_s3_class(w, "data.frame")
  expect_identical(ncol(w), 2L)
  expect_false(anyNA(w))

  #Selecting one stop method gives a vector.
  expect_type(get.w(m, stop.method = "es.mean"), "double")

  #Sampling weights are honoured.
  expect_no_error(get.w(m, s.weights = TRUE))
})

test_that("get.w.CBPS() can recompute the weights from the propensity scores", {
  cb <- fx("cbps")

  from_obj <- get.w(cb)
  recomputed <- get.w(cb, use.weights = FALSE)

  expect_length(recomputed, n_lalonde)
  expect_false(anyNA(recomputed))
  expect_true(all(recomputed >= 0))

  #The estimand may be given explicitly rather than inferred.
  for (e in c("att", "ate")) {
    w <- get.w(cb, use.weights = FALSE, estimand = e)
    expect_length(w, n_lalonde)
    expect_false(anyNA(w))
  }

  #Both routes describe the same design, so they should agree closely.
  expect_equal(cor(from_obj, recomputed), 1, tolerance = .2)
})

test_that("estimand is honoured by the matching-based methods", {
  for (nm in c("optmatch", "cem_match")) {
    obj <- fx(nm)
    for (e in c("ATT", "ATC", "ATE")) {
      w <- suppressWarnings(get.w(obj, estimand = e))
      expect_length(w, n_lalonde)
      expect_false(anyNA(w))
    }
  }

  dm <- fx("designmatch")
  for (e in c("ATT", "ATE")) {
    expect_length(get.w(dm, treat = lalonde$treat, estimand = e), n_lalonde)
  }
})

test_that("get.w.cem.match() can derive the weights from the matching strata", {
  cm <- fx("cem_match")

  from_w <- get.w(cm)
  from_strata <- get.w(cm, use.match.strata = TRUE)

  expect_length(from_strata, n_lalonde)
  expect_false(anyNA(from_strata))

  #Both identify the same matched units.
  expect_identical(from_w > 0, from_strata > 0)
})

test_that("get.w.ebalance() checks the treatment length", {
  eb <- fx("ebalance")

  expect_length(get.w(eb, treat = lalonde$treat), n_lalonde)

  #A treatment with more controls than the object has weights is rejected.
  expect_err(get.w(eb, treat = c(lalonde$treat, 0)),
             "more control units")

  #A factor treatment is processed like a numeric one.
  expect_length(get.w(eb, treat = factor(lalonde$treat)), n_lalonde)
})

test_that("get.w.mimids()/wimids() read the current MatchThem layout", {
  mi <- fx("mimids")

  w <- get.w(mi)
  expect_length(w, 2L * n_lalonde)
  expect_false(anyNA(w))

  #The pre-0.9 layout had no "approach" element; dropping it exercises the
  #backward-compatible branch.
  old <- mi
  old[["approach"]] <- NULL
  expect_no_error(get.w(old))
})

test_that("get.w() handles a list of cem.match objects", {
  cm <- fx("cem_match")

  cml <- structure(list(cm, cm), class = c("cem.match.list", "list"))

  w <- get.w(cml)
  expect_length(w, 2L * n_lalonde)
  expect_false(anyNA(w))

  #Deriving from the matching strata takes a separate branch for lists.
  w2 <- get.w(cml, use.match.strata = TRUE)
  expect_length(w2, 2L * n_lalonde)
  expect_identical(w > 0, w2 > 0)
})

test_that("get.w.designmatch() rejects a non-1:1 match", {
  #`group_id` must account for every treated and control unit.
  dm <- list(id_1 = NULL, t_id = 1:5, c_id = 6:10, group_id = 1:8,
             obj_total = 0, time = 0)

  expect_err(get.w(structure(dm, class = "designmatch"), treat = rep(0:1, 5)),
             "without 1:1 matching cannot be used")
})

test_that("get.w.CBPS() infers each estimand from the supplied weights", {
  cb <- fx("cbps")

  #The estimand is read off the weights: constant among the treated means ATT,
  #constant among the controls means ATC, neither means ATE.
  t1 <- cb$y == 1
  for (e in c("att", "atc", "ate")) {
    x <- cb
    x$weights <- switch(e,
                        att = ifelse(t1, 1, runif(length(t1), .5, 2)),
                        atc = ifelse(t1, runif(sum(t1), .5, 2), 1),
                        ate = runif(length(t1), .5, 2))
    w <- get.w(x, use.weights = FALSE)
    expect_length(w, length(t1))
    expect_false(anyNA(w))
  }
})

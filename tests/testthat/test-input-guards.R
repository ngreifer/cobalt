#Smoke tests and argument guards for every input class cobalt supports.
#
#`x2base.R` repeats three guards across nearly every method -- `subclass`,
#`match.strata`, and (for most classes) `focal` are "not allowed with <cls>
#objects". Driving them from a table covers roughly fifty `arg::err()` call sites
#for a few dozen lines of test code, and the accompanying `bal.tab()`/`get.w()`
#calls exercise each adapter's main path.
#
#One `test_that()` is generated per class at file top level, NOT a loop inside a
#single `test_that()`: `skip_if_not_installed()` inside `fx()` aborts the enclosing
#block, so a loop would stop at the first missing package.

skip_on_cran()

#Classes whose `x2base` method also rejects `focal`. The rest accept it (either
#because `focal` is meaningful for them or because they have no such guard).
.focal_rejected <- c("weightitmsm", "ps", "iptw", "iptw_fac", "ps_cont", "cbps",
                     "cbps_multi", "cbps_cont", "Match", "optmatch", "ebalance",
                     "sbwcau", "designmatch", "mimids",
                     "ps_xgboost", "ps_xgboost_num", "iptw_xgboost")

#Extra arguments some methods require: the covariates cannot be recovered from the
#object alone, or `data` is mandatory.
.extra_args <- function(nm) {
  if (identical(nm, "cem_match")) {
    return(list(data = lalonde))
  }

  if (nm %in% c("Match", "optmatch", "ebalance", "designmatch")) {
    return(list(formula = treat ~ age + educ + re74, data = lalonde))
  }

  list()
}

#`optweight` has no `x2base` method of its own; it is handled by the default
#method and so has none of these guards.
.guarded_classes <- setdiff(fx_object_names(), "optweight")

for (.nm in .guarded_classes) {
  local({
    nm <- .nm

    test_that(sprintf("bal.tab() works and rejects subclass/match.strata for %s objects", nm), {
      obj <- fx(nm)
      extra <- .extra_args(nm)

      b <- suppressMessages(suppressWarnings(
        do.call(bal.tab, c(list(obj), extra))
      ))
      expect_s3_class(b, "bal.tab")

      #Whatever the shape, at least one balance component must be populated.
      #(The component is `Balance` for point treatments and
      #`{Pair,Time,Imputation,Cluster,Subclass}.Balance` for the wrapper classes.)
      bal <- b[endsWith(names(b), "Balance")]
      expect_gt(length(bal), 0L)
      expect_true(any(vapply(bal, function(z) NROW(z) > 0L, logical(1L))))

      #Printing must work for every supported input.
      expect_no_error(capture.output(print(b)))

      sub <- rep(1:4, length.out = n_lalonde)

      expect_err(do.call(bal.tab, c(list(obj), extra, list(subclass = sub))),
                 "subclasses are not allowed with")
      expect_err(do.call(bal.tab, c(list(obj), extra, list(match.strata = sub))),
                 "matching strata are not allowed with")
    })
  })
}

for (.nm in intersect(.focal_rejected, .guarded_classes)) {
  local({
    nm <- .nm

    test_that(sprintf("`focal` is rejected for %s objects", nm), {
      obj <- fx(nm)

      expect_err(do.call(bal.tab, c(list(obj), .extra_args(nm), list(focal = "1"))),
                 "is not allowed with")
    })
  })
}

#`get.w()` methods. `ebalance` and `designmatch` need the treatment supplied
#because their objects do not record it.
for (.nm in setdiff(.guarded_classes, c("ebalance", "designmatch"))) {
  local({
    nm <- .nm

    test_that(sprintf("get.w() returns weights for %s objects", nm), {
      obj <- fx(nm)

      w <- suppressMessages(suppressWarnings(get.w(obj)))

      #A single set of weights is a vector; several sets come back as a data frame.
      if (is.data.frame(w)) {
        expect_gt(ncol(w), 0L)
        expect_true(all(vapply(w, is.numeric, logical(1L))))
        expect_false(anyNA(w))
      }
      else {
        expect_type(w, "double")
        expect_gt(length(w), 0L)
        expect_false(anyNA(w))
        expect_true(all(w >= 0))
      }
    })
  })
}

test_that("get.w() requires `treat` for ebalance and designmatch objects", {
  #These objects record only the weights, not which units are treated.
  eb <- fx("ebalance")
  expect_err(get.w(eb), "must be supplied")
  expect_type(get.w(eb, treat = lalonde$treat), "double")
})

test_that("get.w() requires `treat` for designmatch objects", {
  dm <- fx("designmatch")
  expect_err(get.w(dm), "must be supplied")
  expect_length(get.w(dm, treat = lalonde$treat), n_lalonde)
})

test_that("get.w.ps() honours stop.method, estimand, and s.weights", {
  ps <- fx("ps")

  #Two stop methods were requested, so a data frame comes back.
  w_all <- get.w(ps)
  expect_s3_class(w_all, "data.frame")
  expect_identical(ncol(w_all), 2L)

  #Selecting one stop method gives a vector.
  w_one <- get.w(ps, stop.method = "es.mean")
  expect_type(w_one, "double")
  expect_length(w_one, n_lalonde)

  #All three estimands are reachable and give different weights. "ATC" was once
  #unreachable because its branch tested for "ATT".
  ws <- lapply(c("ATE", "ATT", "ATC"),
               function(e) get.w(ps, stop.method = "es.mean", estimand = e))
  expect_false(isTRUE(all.equal(ws[[1L]], ws[[2L]])))
  expect_false(isTRUE(all.equal(ws[[2L]], ws[[3L]])))
  expect_false(isTRUE(all.equal(ws[[1L]], ws[[3L]])))

  #The ATT weights are 1 for treated units; the ATC weights are 1 for controls.
  treated <- ps$treat == 1
  expect_true(all(ws[[2L]][treated] == 1))
  expect_true(all(ws[[3L]][!treated] == 1))

  expect_err(get.w(ps, estimand = "ATO"), "`estimand` must be")
})

test_that("bal.tab() accepts the arguments each adapter documents", {
  #`matchit` objects can be re-interpreted as weighting or subclassification.
  m <- fx("matchit_sub")
  expect_s3_class(bal.tab(m, method = "subclassification"), "bal.tab.subclass")
  expect_s3_class(bal.tab(m, method = "weighting"), "bal.tab.bin")

  #`ps` objects accept `stop.method` positionally.
  ps <- fx("ps")
  expect_s3_class(bal.tab(ps, stop.method = "es.mean"), "bal.tab.bin")

  #`CBPS` objects reject `estimand`, which is fixed at fit time.
  expect_err(bal.tab(fx("cbps"), estimand = "ATT"),
             "`estimand` is not allowed with")

  #`CBMSM` objects reject sampling weights.
  #(covered here rather than in a separate file because CBPS is the only source)
  expect_s3_class(bal.tab(fx("cbps_cont")), "bal.tab.cont")
})

test_that("cem.match objects require `data`", {
  cm <- fx("cem_match")

  expect_err(bal.tab(cm), "must be specified with")
  expect_s3_class(bal.tab(cm, data = lalonde), "bal.tab.bin")
})

test_that("multiply imputed objects produce imputation tables", {
  mi <- fx("mimids")

  b <- bal.tab(mi)
  expect_s3_class(b, "bal.tab.imp")
  expect_gt(length(b$Imputation.Balance), 1L)

  wi <- fx("wimids")
  expect_s3_class(bal.tab(wi), "bal.tab.imp")

  #`imp` is not an argument for these methods; the imputations come from the object.
  expect_length(get.w(mi), 2L * n_lalonde)
})

test_that("longitudinal objects produce time-point tables", {
  expect_s3_class(bal.tab(fx("weightitmsm")), "bal.tab.msm")
  expect_s3_class(bal.tab(fx("iptw")), "bal.tab.msm")
})

test_that("optweight objects go through the default method", {
  ow <- fx("optweight")

  #No `x2base.optweight` exists, so the duck-typed default method handles it and
  #none of the class-specific guards apply.
  b <- suppressMessages(bal.tab(ow))
  expect_s3_class(b, "bal.tab")
})

#Fixtures for cobalt's tests.
#
#Fitted objects from the packages in Suggests are built lazily and memoized, then
#reused across every test file in the run. `fx()` must be called from inside
#`test_that()`: it calls `skip_if_not_installed()` before touching the cache, so a
#missing package skips the calling test rather than aborting the run (a
#`skip_if_not_installed()` raised at helper top level is an error, not a skip).
#
#Fixtures are NOT saved to disk. cobalt reads other packages' object internals, so
#serialized fixtures would keep passing while real users break on the next upstream
#release. The one exception is `cem`, which is a frequent install failure; see
#`fixtures/` and `make-fixtures.R`.

.fx_cache <- new.env(parent = emptyenv())

#Deterministic stand-ins for randomly generated inputs. Nothing here uses the RNG,
#so results are reproducible regardless of which test file runs first, and
#snapshots stay stable.
#
#Two traps these avoid:
#  1. `lalonde` is sorted with all 185 treated units first. Building an index with
#     `rep(..., each = )` therefore puts a contiguous run in one group, yielding a
#     cluster or imputation in which the treatment takes only one value. Always
#     interleave with `length.out = `.
#  2. Two indices built with the same period are perfectly collinear (subclass i
#     coincides exactly with cluster i), which empties cells of their crosstab.
#     The periods below (2, 3, 4) are chosen so that does not happen.
n_lalonde <- 614L

w_fixed   <- rep(c(1, 2), length.out = n_lalonde)
sw_fixed  <- rep(c(1, 1.5, 2), length.out = n_lalonde)
cl_idx    <- factor(rep(c("a", "b", "c"), length.out = n_lalonde))
sub_idx   <- rep(1:4, length.out = n_lalonde)
imp_idx   <- rep(1:2, length.out = n_lalonde)

cov_names <- c("age", "educ", "race", "married", "nodegree", "re74", "re75")

#`twang` fits are pinned to a small number of trees. The quality of the boosted
#model is irrelevant here; all that matters is a valid object with a populated
#`$ps` and `$w`. Two stop methods are requested so that `.process_stop_method()`
#and the multiple-weight code paths are reached.
.fx_registry <- function() {
  list(
    ## ---- no optional dependencies ------------------------------------------
    covs = list(pkg = character(), build = function() {
      lalonde[cov_names]
    }),

    ## ---- MatchIt ------------------------------------------------------------
    matchit = list(pkg = "MatchIt", build = function() {
      MatchIt::matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
                       data = lalonde)
    }),
    matchit_full = list(pkg = "MatchIt", build = function() {
      MatchIt::matchit(treat ~ age + educ + race + re74, data = lalonde,
                       method = "full")
    }),
    matchit_sub = list(pkg = "MatchIt", build = function() {
      MatchIt::matchit(treat ~ age + educ + race + re74, data = lalonde,
                       method = "subclass", subclass = 4)
    }),
    matchit_cal = list(pkg = "MatchIt", build = function() {
      MatchIt::matchit(treat ~ age + educ + race + re74, data = lalonde,
                       caliper = .1)
    }),

    ## ---- WeightIt -----------------------------------------------------------
    weightit = list(pkg = "WeightIt", build = function() {
      WeightIt::weightit(treat ~ age + educ + race + married + nodegree + re74 + re75,
                         data = lalonde, method = "glm", estimand = "ATT")
    }),
    weightit_multi = list(pkg = "WeightIt", build = function() {
      WeightIt::weightit(race ~ age + educ + married + re74, data = lalonde,
                         method = "glm", estimand = "ATE")
    }),
    weightit_cont = list(pkg = "WeightIt", build = function() {
      WeightIt::weightit(re75 ~ age + educ + married + re74, data = lalonde,
                         method = "glm")
    }),
    weightitmsm = list(pkg = "WeightIt", build = function() {
      WeightIt::weightitMSM(
        list(treat ~ age + educ,
             nodegree ~ age + educ + treat),
        data = lalonde, method = "glm")
    }),

    ## ---- twang --------------------------------------------------------------
    ps = list(pkg = "twang", build = function() {
      twang::ps(treat ~ age + educ + race + married + nodegree + re74 + re75,
                data = lalonde, estimand = "ATT",
                stop.method = c("es.mean", "ks.max"),
                n.trees = 200, interaction.depth = 2, verbose = FALSE)
    }),
    mnps = list(pkg = "twang", build = function() {
      twang::mnps(race ~ age + educ + married + re74, data = lalonde,
                  estimand = "ATE", stop.method = "es.mean",
                  n.trees = 200, verbose = FALSE)
    }),
    iptw = list(pkg = "twang", build = function() {
      twang::iptw(list(treat ~ age + educ,
                       nodegree ~ age + educ + treat),
                  data = lalonde, stop.method = "es.mean",
                  n.trees = 200, verbose = FALSE)
    }),
    ps_cont = list(pkg = "twangContinuous", build = function() {
      twangContinuous::ps.cont(re75 ~ age + educ + married + re74,
                               data = lalonde, n.trees = 200, verbose = FALSE)
    }),

    ## ---- twang with the xgboost engine --------------------------------------
    #cobalt reads a different component for xgboost fits: `gbm.obj$var.names` is
    #absent, so it falls back to `gbm.obj$feature_names`. These fixtures attempt
    #the real fit; `fx()` turns a failure into an informative skip, which is what
    #happens when the installed xgboost is too new for twang (twang 2.6.2 calls
    #`xgboost(data=, label=, params=, feval=)`, all removed or renamed in
    #xgboost 3.x, and fails with "Passed invalid 'nthreads': NA").
    ps_xgboost = list(pkg = c("twang", "xgboost"), build = function() {
      twang::ps(treat ~ age + educ + race, data = lalonde, estimand = "ATT",
                stop.method = "es.mean", n.trees = 100, verbose = FALSE,
                version = "xgboost")
    }),
    ps_xgboost_num = list(pkg = c("twang", "xgboost"), build = function() {
      #All-numeric covariates, so xgboost creates no one-hot columns and every
      #`feature_names` entry is still a column of the stored data.
      twang::ps(treat ~ age + educ + re74, data = lalonde, estimand = "ATT",
                stop.method = "es.mean", n.trees = 100, verbose = FALSE,
                version = "xgboost")
    }),
    iptw_xgboost = list(pkg = c("twang", "xgboost"), build = function() {
      twang::iptw(list(treat ~ age + race,
                       nodegree ~ age + race + treat),
                  data = lalonde, stop.method = "es.mean", n.trees = 100,
                  verbose = FALSE, version = "xgboost")
    }),

    #A real `iptw` fit with a factor covariate, used as the base for the
    #xgboost-shaped objects in test-input-twang-xgboost.R.
    iptw_fac = list(pkg = "twang", build = function() {
      twang::iptw(list(treat ~ age + race,
                       nodegree ~ age + race + treat),
                  data = lalonde, stop.method = "es.mean",
                  n.trees = 200, verbose = FALSE)
    }),

    ## ---- CBPS ---------------------------------------------------------------
    cbps = list(pkg = "CBPS", build = function() {
      CBPS::CBPS(treat ~ age + educ + re74, data = lalonde, ATT = 1,
                 standardize = FALSE)
    }),
    cbps_multi = list(pkg = "CBPS", build = function() {
      #`race` is already a factor; wrapping it in `as.factor()` in the formula
      #leaves a term cobalt cannot resolve back to a variable.
      CBPS::CBPS(race ~ age + educ + re74, data = lalonde, ATT = 0,
                 standardize = FALSE)
    }),
    cbps_cont = list(pkg = "CBPS", build = function() {
      CBPS::CBPS(re75 ~ age + educ + re74, data = lalonde, ATT = 0,
                 standardize = FALSE)
    }),

    ## ---- Matching / optmatch ------------------------------------------------
    Match = list(pkg = "Matching", build = function() {
      covs <- lalonde[c("age", "educ", "re74")]
      p <- fitted(glm(treat ~ age + educ + re74, data = lalonde,
                      family = binomial))
      Matching::Match(Tr = lalonde$treat, X = p, estimand = "ATT")
    }),
    optmatch = list(pkg = "optmatch", build = function() {
      optmatch::fullmatch(treat ~ age + educ + re74, data = lalonde)
    }),

    ## ---- ebal / sbw / designmatch / optweight -------------------------------
    ebalance = list(pkg = "ebal", build = function() {
      #`drop.first = TRUE` (not "if2") is required: keeping every level of `race`
      #makes the covariate matrix collinear and ebalance() refuses to run.
      covs <- splitfactor(lalonde[c("age", "educ", "race", "re74")],
                          drop.first = TRUE)
      ebal::ebalance(Treatment = lalonde$treat, X = covs)
    }),
    sbwcau = list(pkg = "sbw", build = function() {
      d <- lalonde[c("treat", "age", "educ", "re74")]
      sbw::sbw(d, ind = "treat",
               bal = list(bal_cov = c("age", "educ", "re74"),
                          bal_alg = FALSE, bal_tol = .1),
               par = list(par_est = "att"))
    }),
    designmatch = list(pkg = "designmatch", build = function() {
      #`lalonde` is already sorted with the treated units first, as bmatch()
      #requires. `t_max` is mandatory; without it the solver errors.
      cv <- as.matrix(lalonde[c("age", "educ", "re74", "re75")])
      designmatch::bmatch(
        lalonde$treat,
        total_groups = sum(lalonde$treat == 1),
        mom = list(covs = cv,
                   tols = designmatch::absstddif(cv, lalonde$treat, .05)),
        solver = list(name = "highs", approximate = 0, t_max = 60, trace = 0))
    }),
    optweight = list(pkg = "optweight", build = function() {
      optweight::optweight(treat ~ age + educ + re74, data = lalonde,
                           estimand = "ATT", tols = .01)
    }),

    ## ---- mice / MatchThem ---------------------------------------------------
    mids = list(pkg = "mice", build = function() {
      #m = 2 is the minimum that makes nlevels(X$imp) > 1 in .assign_X_class().
      mice::mice(lalonde_mis[c("treat", "age", "educ", "race", "married", "re74")],
                 m = 2, maxit = 1, printFlag = FALSE)
    }),
    mimids = list(pkg = c("mice", "MatchThem"), build = function() {
      MatchThem::matchthem(treat ~ age + educ + race + married + re74,
                           datasets = fx("mids"), approach = "within",
                           method = "nearest")
    }),
    wimids = list(pkg = c("mice", "MatchThem"), build = function() {
      MatchThem::weightthem(treat ~ age + educ + race + married + re74,
                            datasets = fx("mids"), approach = "within",
                            method = "glm")
    }),

    ## ---- cem (loaded from disk; see fixtures/make-fixtures.R) ---------------
    cem_match = list(pkg = character(), build = function() {
      path <- test_path("fixtures", "cem_match.rds")
      if (!file.exists(path)) {
        stop("fixture file 'fixtures/cem_match.rds' is missing", call. = FALSE)
      }
      readRDS(path)
    })
  )
}

#Retrieve a memoized fixture. Skips the calling test when a required package is
#unavailable or when the fixture cannot be built.
fx <- function(name) {
  reg <- .fx_registry()[[name]]

  if (is.null(reg)) {
    #A mistyped fixture name must be a hard error, not a silent skip.
    stop("no fixture named '", name, "'", call. = FALSE)
  }

  for (p in reg$pkg) {
    skip_if_not_installed(p)
  }

  if (!exists(name, envir = .fx_cache, inherits = FALSE)) {
    built <- tryCatch(suppressMessages(suppressWarnings(reg$build())),
                      error = function(e) e)
    assign(name, built, envir = .fx_cache)
  }

  out <- get(name, envir = .fx_cache, inherits = FALSE)

  #Failures are cached as conditions so a broken upstream costs one attempt and
  #produces informative skips rather than repeated cryptic errors.
  if (inherits(out, "condition")) {
    skip(sprintf("fixture '%s' could not be built: %s", name, conditionMessage(out)))
  }

  out
}

#Names of fixtures that are fitted objects accepted by `bal.tab()`/`get.w()`, used
#to drive the data-driven guard tests. Excludes the plain data fixtures and the
#`mids` object, which is passed via `data`, not as `x`.
fx_object_names <- function() {
  setdiff(names(.fx_registry()), c("covs", "mids"))
}

#Is every package a fixture needs installed?
fx_available <- function(name) {
  pkgs <- .fx_registry()[[name]]$pkg

  length(pkgs) == 0L ||
    all(vapply(pkgs, requireNamespace, logical(1L), quietly = TRUE))
}

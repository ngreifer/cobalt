skip_on_cran()

test_that("`s.d.denom` processes correctly", {
    data("lalonde")
    
    weights <- list(wATE = runif(nrow(lalonde)),
                    wATT = ifelse(lalonde$treat == 1, 1, runif(sum(lalonde$treat != 1))))
    
    # Binary tretament
    cov_names <- c("age", "educ", "race", "married", "nodegree", "re74", "re75")
    
    f <- reformulate(cov_names, "treat")
    
    expect_msg(bal.tab(f, data = lalonde, binary = "std", continuous = "std"), 'Note: `s.d.denom` not specified; assuming "pooled".')
    
    #No message when no variables are to be standardized
    expect_no_message(bal.tab(f, data = lalonde, binary = "raw", continuous = "raw"))
    
    #No message when SMDs not requested
    expect_no_message(bal.tab(f, data = lalonde, binary = "raw", continuous = "raw",
                              stats = "ks"))
    
    expect_msg(bal.tab(f, data = lalonde, binary = "std", continuous = "std",
                           s.d.denom = "weighted"), 'Note: `s.d.denom` specified as "weighted", but no weights supplied; setting to "all".')
    
    expect_equal(bal.tab(f, data = lalonde, binary = "std", continuous = "std", s.d.denom = "pooled")$Balance$Diff.Un,
                 unname(col_w_smd(lalonde[cov_names], treat = lalonde$treat, s.d.denom = "pooled")))
    
    expect_equal(bal.tab(f, data = lalonde, binary = "std", continuous = "std", s.d.denom = "treated")$Balance$Diff.Un,
                 unname(col_w_smd(lalonde[cov_names], treat = lalonde$treat, s.d.denom = "treated")))
    
    expect_equal(bal.tab(f, data = lalonde, binary = "std", continuous = "std", s.d.denom = "treated"),
                 bal.tab(f, data = lalonde, binary = "std", continuous = "std", s.d.denom = "1"))
    
    expect_msg(bal.tab(f, data = lalonde, binary = "std", continuous = "std",
                           weights = weights$wATE), 'Note: `s.d.denom` not specified; assuming "pooled".')
    
    expect_no_message(bal.tab(f, data = lalonde, binary = "std", continuous = "std",
                              weights = weights$wATT))
    
    expect_msg(bal.tab(f, data = lalonde, binary = "std", continuous = "std",
                           weights = weights), 'Note: `s.d.denom` not specified; assuming "pooled" for `wATE` and "treated" for `wATT`.')
    
    #bal.tab() used the first s.d.denom for the unadjusted differences
    expect_equal(suppressMessages(bal.tab(f, data = lalonde, binary = "std", continuous = "std", 
                                          weights = weights, un = TRUE)$Balance$Diff.Un),
                 unname(col_w_smd(lalonde[cov_names], treat = lalonde$treat, s.d.denom = "pooled")))
    
    expect_equal(suppressMessages(bal.tab(f, data = lalonde, binary = "std", continuous = "std", 
                                          weights = weights[2:1], un = TRUE)$Balance$Diff.Un),
                 unname(col_w_smd(lalonde[cov_names], treat = lalonde$treat, s.d.denom = "treated")))
    
    # Continuous treatment
    cov_names <- c("age", "educ", "race", "married", "nodegree", "re74")
    
    f <- reformulate(cov_names, "re75")
    
    expect_msg(bal.tab(f, data = lalonde, s.d.denom = "weighted"), 'Note: `s.d.denom` specified as "weighted", but no weights supplied; setting to "all".')
    
    expect_equal(bal.tab(f, data = lalonde)$Balance$Corr.Un,
                 unname(col_w_corr(lalonde[cov_names], treat = lalonde$re75)))
    
    expect_equal(bal.tab(f, data = lalonde, weights = weights$wATE)$Balance$Corr.Adj,
                 unname(col_w_corr(lalonde[cov_names], treat = lalonde$re75,
                                   weights = weights$wATE)))
    
    # Multicategory treatment
    cov_names <- c("age", "educ", "married", "nodegree", "re74", "re75")
    
    f <- reformulate(cov_names, "race")
    
    expect_msg(bal.tab(f, data = lalonde, binary = "std", continuous = "std"), 'Note: `s.d.denom` not specified; assuming "pooled".')
    
    #No message when no variables are to be standardized
    expect_no_message(bal.tab(f, data = lalonde, binary = "raw", continuous = "raw"))
    
    #No message when SMDs not requested
    expect_no_message(bal.tab(f, data = lalonde, binary = "raw", continuous = "raw",
                              stats = "ks"))
    
    expect_msg(bal.tab(f, data = lalonde, binary = "std", continuous = "std",
                           s.d.denom = "weighted"), 'Note: `s.d.denom` specified as "weighted", but no weights supplied; setting to "all".')
    
    expect_equal(bal.tab(f, data = lalonde, binary = "std", continuous = "std", s.d.denom = "pooled")$Pair.Balance[["hispan vs. black"]]$Balance$Diff.Un,
                 unname(col_w_smd(lalonde[cov_names], treat = lalonde$race, s.d.denom = "pooled",
                                  subset = lalonde$race %in% c("hispan", "black"))))
    expect_equal(suppressMessages(bal.tab(f, data = lalonde, binary = "std", continuous = "std")),
                 bal.tab(f, data = lalonde, binary = "std", continuous = "std", s.d.denom = "pooled"))
    
    # Clustered data
    
    cov_names <- c("age", "educ", "race", "married", "nodegree", "re74", "re75")
    
    f <- reformulate(cov_names, "treat")
    
    expect_msg(bal.tab(f, data = lalonde, binary = "std", continuous = "std",
                           cluster = "race"), 'Note: `s.d.denom` not specified; assuming "pooled".')
    
    expect_no_message(bal.tab(f, data = lalonde, binary = "raw", continuous = "raw",
                              cluster = "race"))
    
    expect_equal(bal.tab(f, data = lalonde, binary = "std", continuous = "std", s.d.denom = "pooled",
                         cluster = "race")$Cluster.Balance$black$Balance$Diff.Un,
                 unname(col_w_smd(lalonde[lalonde$race == "black", setdiff(cov_names, "race")],
                                  treat = lalonde$treat[lalonde$race == "black"], s.d.denom = "pooled")))

})

# ---------------------------------------------------------------------------
# The block above checks the messages, and cross-checks `"pooled"` and `"treated"`
# against `col_w_smd()`. The blocks below do the same for the remaining accepted
# values, check that they are genuinely distinct from one another, and check the
# per-weight-set form.

sdd_fixture <- function() {
    data("lalonde", package = "cobalt")

    cov_names <- c("age", "educ", "race", "married", "nodegree", "re74", "re75")

    set.seed(606)

    list(cov_names = cov_names,
         covs = lalonde[cov_names],
         f = reformulate(cov_names, "treat"),
         treat = lalonde$treat,
         w = runif(nrow(lalonde)),
         sw = runif(nrow(lalonde)))
}

test_that("every accepted `s.d.denom` matches col_w_smd()", {
    fx <- sdd_fixture()

    #`splitfactor()` reproduces the layout `bal.tab()` builds internally. With both
    #`binary` and `continuous` set to "std", every row is standardized, so `std` is
    #TRUE throughout.
    sp <- splitfactor(fx$covs, drop.first = FALSE)

    for (sd in c("pooled", "treated", "control", "all", "weighted", "hedges")) {
        b <- bal.tab(fx$covs, treat = fx$treat, binary = "std", continuous = "std",
                     s.d.denom = sd, weights = fx$w, s.weights = fx$sw, un = TRUE)

        expect_equal(b$Balance$Diff.Adj,
                     unname(col_w_smd(sp, treat = fx$treat, weights = fx$w,
                                      s.weights = fx$sw, s.d.denom = sd, std = TRUE)),
                     info = sd)
        #The unadjusted column uses no balancing weights in the numerator, but
        #under `"weighted"` its denominator still comes from the *weighted* sample.
        #That is what `weighted.weights` exists for; it is ignored otherwise.
        expect_equal(b$Balance$Diff.Un,
                     unname(col_w_smd(sp, treat = fx$treat, s.weights = fx$sw,
                                      s.d.denom = sd, std = TRUE,
                                      weighted.weights = fx$w)),
                     info = sd)
    }

    #`weighted.weights` really is inert for every other denominator.
    expect_equal(col_w_smd(sp, treat = fx$treat, s.weights = fx$sw,
                           s.d.denom = "pooled", std = TRUE, weighted.weights = fx$w),
                 col_w_smd(sp, treat = fx$treat, s.weights = fx$sw,
                           s.d.denom = "pooled", std = TRUE))

    #And `s.weights` must not be applied to the `"weighted"` denominator twice:
    #the denominator is the SD of the full sample under the *product* of the two
    #weight sets. Binary columns use sqrt(p(1-p)), so `bin.vars` has to say which
    #they are.
    bin <- vapply(sp, function(x) all(x %in% c(0, 1)), logical(1L))

    expect_equal(unname(col_w_smd(sp, treat = fx$treat, weights = fx$w,
                                  s.weights = fx$sw, s.d.denom = "weighted",
                                  std = TRUE)),
                 unname((col_w_mean(sp, subset = fx$treat == 1, weights = fx$w * fx$sw) -
                             col_w_mean(sp, subset = fx$treat == 0, weights = fx$w * fx$sw)) /
                            col_w_sd(sp, weights = fx$w * fx$sw, bin.vars = bin)))

    #The default `binary = "raw"` leaves the binary rows unstandardized, so the
    #same call with `std` following `Type` is what matches then.
    b_raw <- bal.tab(fx$covs, treat = fx$treat, s.d.denom = "pooled",
                     weights = fx$w, s.weights = fx$sw, un = TRUE)

    expect_equal(b_raw$Balance$Diff.Adj,
                 unname(col_w_smd(sp, treat = fx$treat, weights = fx$w,
                                  s.weights = fx$sw, s.d.denom = "pooled",
                                  std = b_raw$Balance$Type == "Contin.")))
})

test_that("the `s.d.denom` values are distinct from one another", {
    fx <- sdd_fixture()

    #If a refactor started ignoring `s.d.denom` and always pooled, every
    #cross-check above would still pass for `"pooled"` alone. These assertions are
    #what make the choice observable.
    diffs <- lapply(c("pooled", "treated", "control", "all", "hedges"),
                    function(sd) {
                        bal.tab(fx$f, data = lalonde, binary = "std", continuous = "std",
                                s.d.denom = sd)$Balance$Diff.Un
                    })

    for (i in seq_along(diffs)[-1L]) {
        for (j in seq_len(i - 1L)) {
            expect_false(isTRUE(all.equal(diffs[[i]], diffs[[j]])))
        }
    }

    #The numeric aliases name the treatment *values*, so "0" is the control group
    #and "1" the treated one.
    expect_equal(bal.tab(fx$f, data = lalonde, binary = "std", continuous = "std",
                         s.d.denom = "0")$Balance,
                 bal.tab(fx$f, data = lalonde, binary = "std", continuous = "std",
                         s.d.denom = "control")$Balance)

    #Denominators that condition on one group must equal that group's own SD.
    sp <- splitfactor(fx$covs, drop.first = FALSE)
    b_t <- bal.tab(fx$covs, treat = fx$treat, binary = "std", continuous = "std",
                   s.d.denom = "treated")
    cont <- b_t$Balance$Type == "Contin."

    md <- unname(col_w_mean(sp, subset = fx$treat == 1) -
                     col_w_mean(sp, subset = fx$treat == 0))
    sd_t <- unname(col_w_sd(sp, subset = fx$treat == 1))

    expect_equal(b_t$Balance$Diff.Un[cont], (md / sd_t)[cont])
})

test_that("`s.d.denom` can differ per weight set", {
    fx <- sdd_fixture()

    #Each weight set gets its own denominator, and the resulting column must equal
    #the table produced by that weight set and denominator alone.
    b <- bal.tab(fx$f, data = lalonde, binary = "std", continuous = "std",
                 weights = list(a = fx$w, b = fx$sw),
                 s.d.denom = c("treated", "control"), un = TRUE)

    b_a <- bal.tab(fx$f, data = lalonde, binary = "std", continuous = "std",
                   weights = fx$w, s.d.denom = "treated", un = TRUE)
    b_b <- bal.tab(fx$f, data = lalonde, binary = "std", continuous = "std",
                   weights = fx$sw, s.d.denom = "control", un = TRUE)

    expect_equal(b$Balance$Diff.a, b_a$Balance$Diff.Adj)
    expect_equal(b$Balance$Diff.b, b_b$Balance$Diff.Adj)

    #The shared unadjusted column follows the *first* denominator, which is what
    #the block above asserts for the auto-detected case.
    expect_equal(b$Balance$Diff.Un, b_a$Balance$Diff.Un)
})

test_that("`s.d.denom` is irrelevant to statistics that do not standardize", {
    fx <- sdd_fixture()

    #KS statistics and variance ratios have no standard deviation in them, so
    #changing the denominator must not move them. A refactor that threaded
    #`s.d.denom` too far down would show up here.
    b_p <- bal.tab(fx$covs, treat = fx$treat, s.d.denom = "pooled", weights = fx$w,
                   un = TRUE, stats = c("m", "v", "ks"))
    b_t <- bal.tab(fx$covs, treat = fx$treat, s.d.denom = "treated", weights = fx$w,
                   un = TRUE, stats = c("m", "v", "ks"))

    expect_equal(b_t$Balance$KS.Adj, b_p$Balance$KS.Adj)
    expect_equal(b_t$Balance$V.Ratio.Adj, b_p$Balance$V.Ratio.Adj)

    #And with nothing standardized, the whole table is identical.
    b_p_raw <- bal.tab(fx$covs, treat = fx$treat, binary = "raw", continuous = "raw",
                       s.d.denom = "pooled", weights = fx$w, un = TRUE)
    b_t_raw <- bal.tab(fx$covs, treat = fx$treat, binary = "raw", continuous = "raw",
                       s.d.denom = "treated", weights = fx$w, un = TRUE)

    expect_equal(b_t_raw$Balance, b_p_raw$Balance)

    #`binary`/`continuous` decide which rows are standardized at all.
    b_std <- bal.tab(fx$covs, treat = fx$treat, binary = "std", continuous = "std",
                     s.d.denom = "pooled", weights = fx$w, un = TRUE)

    binary_rows <- b_std$Balance$Type == "Binary"
    expect_true(any(binary_rows))
    expect_false(isTRUE(all.equal(b_std$Balance$Diff.Adj[binary_rows],
                                  b_p_raw$Balance$Diff.Adj[binary_rows])))
    expect_equal(b_std$Balance$Diff.Adj[!binary_rows],
                 b_p$Balance$Diff.Adj[!binary_rows])
})

test_that("`s.d.denom` is honoured at each time point of a longitudinal treatment", {
    fx <- sdd_fixture()

    set.seed(28)
    d <- data.frame(lalonde[c("age", "educ", "re74")],
                    t1 = lalonde$treat,
                    t2 = rbinom(nrow(lalonde), 1L, 0.4),
                    c1 = lalonde$re75,
                    c2 = lalonde$re78,
                    w = fx$w)

    fs <- list(t1 ~ age + educ, t2 ~ age + educ + t1)

    #A longitudinal treatment defaults to the ATE's denominator at every time point,
    #which is what `estimand = "ATE"` gives for a point treatment.
    b_default <- bal.tab(fs, data = d, weights = "w")
    b_ate <- bal.tab(fs, data = d, weights = "w", s.d.denom = "pooled")

    expect_equal(b_default$Balance.Across.Times, b_ate$Balance.Across.Times)

    #But the default is only a default: another value changes the numbers rather than
    #being discarded, which is what used to happen.
    b_treated <- bal.tab(fs, data = d, weights = "w", s.d.denom = "treated")

    expect_false(isTRUE(all.equal(b_treated$Balance.Across.Times,
                                  b_default$Balance.Across.Times)))
    expect_equal(b_treated$Time.Balance[[1L]]$Balance$Diff.Adj,
                 bal.tab(t1 ~ age + educ, data = d, weights = "w",
                         s.d.denom = "treated")$Balance$Diff.Adj)

    #The same for a continuous longitudinal treatment, whose ATE default is "all".
    fs_c <- list(c1 ~ age + educ, c2 ~ age + educ + c1)

    expect_equal(bal.tab(fs_c, data = d, weights = "w")$Balance.Across.Times,
                 bal.tab(fs_c, data = d, weights = "w",
                         s.d.denom = "all")$Balance.Across.Times)

    b_w <- bal.tab(fs_c, data = d, weights = "w", s.d.denom = "weighted")
    expect_false(isTRUE(all.equal(b_w$Balance.Across.Times,
                                  bal.tab(fs_c, data = d,
                                          weights = "w")$Balance.Across.Times)))

    #An unusable value is now rejected rather than ignored.
    expect_err(bal.tab(fs, data = d, weights = "w", s.d.denom = "nope"),
               "`s.d.denom` should be one of")
})

test_that("an unrecognized `estimand` warns instead of silently becoming ATE", {
    fx <- sdd_fixture()

    #`estimand` picks the standardization factor, so a typo previously changed the
    #numbers with no indication: "ATTT" produced the same table as "ATE".
    expect_wrn(bal.tab(fx$covs, treat = fx$treat, estimand = "ATTT"),
               '`estimand` should be "ATT", "ATC", "ATE", "ATO", or "ATM"')
    expect_wrn(bal.tab(fx$covs, treat = fx$treat, estimand = "nope"), "ignoring it")

    #The value is still ignored rather than fatal, and the fallback is unchanged.
    expect_equal(suppressWarnings(suppressMessages(
        bal.tab(fx$covs, treat = fx$treat, estimand = "nope")))$Balance,
        suppressMessages(bal.tab(fx$covs, treat = fx$treat))$Balance)

    #Recognized values, including abbreviations, stay silent.
    for (e in c("ATT", "att", "ATC", "ATE", "ATO", "ATM")) {
        expect_no_warning(suppressMessages(
            bal.tab(fx$covs, treat = fx$treat, estimand = e)))
    }

    #An ATT or ATC with a multi-category treatment is legitimate -- the focal group
    #supplies the denominator -- and must not warn.
    covs_m <- fx$covs[setdiff(fx$cov_names, "race")]

    expect_no_warning(suppressMessages(
        bal.tab(covs_m, treat = lalonde$race, estimand = "ATT", focal = "white")))
    expect_no_warning(suppressMessages(
        bal.tab(covs_m, treat = lalonde$race, estimand = "ATT")))

    #An explicit `s.d.denom` short-circuits the estimand entirely.
    expect_no_warning(
        bal.tab(fx$covs, treat = fx$treat, estimand = "nope", s.d.denom = "pooled"))

    #`col_w_smd()` and friends pass `quietly = TRUE` and must stay silent.
    expect_no_warning(col_w_smd(fx$covs[c("age", "educ")], treat = fx$treat,
                                s.d.denom = "pooled", std = TRUE))
})

# The tests above exercise `s.d.denom` through `bal.tab()`. These call
# `.get_s.d.denom()` directly, one test per branch, because it is the function that
# decides the standardization factor and the branch that produced a given answer is
# not otherwise observable. Written before the function was restructured, so they pin
# the behaviour rather than describe the implementation.

sdd <- function(...) cobalt:::.get_s.d.denom(...)

t_bin <- function() cobalt:::process_treat(lalonde$treat)
t_multi <- function() cobalt:::process_treat(lalonde$race)

#Weights whose constant group identifies the estimand: an ATT design leaves the
#treated at 1, an ATC design leaves the control at 1, an ATE design varies both.
w_att <- function() ifelse(lalonde$treat == 1, 1, seq_len(nrow(lalonde)) / nrow(lalonde))
w_atc <- function() ifelse(lalonde$treat == 0, 1, seq_len(nrow(lalonde)) / nrow(lalonde))
w_ate <- function() seq_len(nrow(lalonde)) / nrow(lalonde)

test_that(".get_s.d.denom() honours an explicit `s.d.denom`", {
  tb <- t_bin()

  #The four denominators that need no treatment information.
  for (d in c("pooled", "all", "weighted", "hedges")) {
    expect_identical(as.vector(sdd(d, treat = tb, quietly = TRUE)), d, info = d)
  }

  #`treated`/`control` are resolved to the treatment's own values.
  expect_identical(as.vector(sdd("treated", treat = tb, quietly = TRUE)), "1")
  expect_identical(as.vector(sdd("control", treat = tb, quietly = TRUE)), "0")

  #A treatment value may be named directly.
  expect_identical(as.vector(sdd("1", treat = tb, quietly = TRUE)), "1")

  #`focal` is only allowable when a focal group was supplied, and resolves to it.
  expect_identical(as.vector(sdd("focal", treat = tb, focal = "1", quietly = TRUE)), "1")
  expect_err(sdd("focal", treat = tb, quietly = TRUE), "`s.d.denom` should be one of")

  expect_err(sdd("bogus", treat = tb, quietly = TRUE), "`s.d.denom` should be one of")

  #The result is marked as checked, and a checked value is returned untouched. This is
  #what stops a wrapper's denominator from being re-derived by each of its children.
  out <- sdd("hedges", treat = tb, quietly = TRUE)
  expect_true(isTRUE(attr(out, "checked")))
  expect_identical(sdd(out, treat = t_multi(), quietly = TRUE), out)
})

test_that(".get_s.d.denom() maps `estimand` to a denominator", {
  tb <- t_bin()

  expect_identical(as.vector(sdd(estimand = "ATT", treat = tb, quietly = TRUE)), "1")
  expect_identical(as.vector(sdd(estimand = "ATC", treat = tb, quietly = TRUE)), "0")
  expect_identical(as.vector(sdd(estimand = "ATE", treat = tb, quietly = TRUE)), "pooled")
  expect_identical(as.vector(sdd(estimand = "ATO", treat = tb, quietly = TRUE)), "weighted")
  expect_identical(as.vector(sdd(estimand = "ATM", treat = tb, quietly = TRUE)), "weighted")

  #An ATT with a multi-category treatment has no "treated" group, so the focal group
  #decides instead -- and without one it falls through to the weights.
  expect_identical(as.vector(sdd(estimand = "ATT", treat = t_multi(), focal = "black",
                                 quietly = TRUE)),
                   "black")
  expect_identical(as.vector(sdd(estimand = "ATT", treat = t_multi(), quietly = TRUE)),
                   "pooled")

  #An unrecognized estimand warns and is ignored rather than silently reading as ATE.
  expect_wrn(sdd(estimand = "ATTT", treat = tb),
             '`estimand` should be "ATT", "ATC", "ATE", "ATO", or "ATM"; ignoring it')
  expect_identical(as.vector(sdd(estimand = "ATTT", treat = tb, quietly = TRUE)), "pooled")

  #`s.d.denom` wins over `estimand`.
  expect_identical(as.vector(sdd("hedges", estimand = "ATT", treat = tb, quietly = TRUE)),
                   "hedges")
})

test_that(".get_s.d.denom() infers a denominator when nothing is specified", {
  tb <- t_bin()

  #No weights and no subclasses: nothing to infer from.
  expect_identical(as.vector(sdd(treat = tb, quietly = TRUE)), "pooled")

  #A focal group is used when there is at most one set of weights.
  expect_identical(as.vector(sdd(treat = tb, focal = "1", quietly = TRUE)), "1")

  #Otherwise the weights are read: the group left unweighted is the denominator.
  expect_identical(as.vector(sdd(treat = tb, weights = data.frame(Adj = w_att()),
                                 quietly = TRUE)),
                   "1")
  expect_identical(as.vector(sdd(treat = tb, weights = data.frame(Adj = w_atc()),
                                 quietly = TRUE)),
                   "0")
  expect_identical(as.vector(sdd(treat = tb, weights = data.frame(Adj = w_ate()),
                                 quietly = TRUE)),
                   "pooled")

  #Each set of weights is read separately.
  expect_identical(as.vector(sdd(treat = tb,
                                 weights = data.frame(A = w_att(), B = w_atc()),
                                 quietly = TRUE)),
                   c("1", "0"))

  #With several sets of weights a focal group is not enough; the weights decide.
  expect_identical(as.vector(sdd(treat = tb, focal = "1",
                                 weights = data.frame(A = w_ate(), B = w_ate()),
                                 quietly = TRUE)),
                   c("pooled", "pooled"))

  #Subclasses of equal size favour no group.
  expect_identical(as.vector(sdd(treat = tb,
                                 subclass = factor(rep(1:4, length.out = nrow(lalonde))),
                                 quietly = TRUE)),
                   "pooled")
})

test_that(".get_s.d.denom() recycles and names per set of weights", {
  tb <- t_bin()
  w2 <- data.frame(W1 = w_ate(), W2 = w_ate())

  #One value is used for every set of weights, and the result takes their names.
  out <- sdd("pooled", treat = tb, weights = w2, quietly = TRUE)
  expect_identical(as.vector(out), c("pooled", "pooled"))
  expect_named(out, c("W1", "W2"))

  #One value per set of weights is kept as given.
  expect_identical(as.vector(sdd(c("pooled", "all"), treat = tb, weights = w2,
                                 quietly = TRUE)),
                   c("pooled", "all"))
  expect_identical(as.vector(sdd(estimand = c("ATT", "ATE"), treat = tb, weights = w2,
                                 quietly = TRUE)),
                   c("1", "pooled"))

  #Any other length is an error, and `s.d.denom` and `estimand` say so separately.
  expect_err(sdd(c("pooled", "all"), treat = tb, weights = data.frame(Adj = w_ate()),
                 quietly = TRUE),
             "`s.d.denom` must have length 1 or equal to the number of valid sets of weights, which is 1")
  expect_err(sdd(estimand = c("ATT", "ATE"), treat = tb,
                 weights = data.frame(Adj = w_ate()), quietly = TRUE),
             "`estimand` must have length 1 or equal to the number of valid sets of weights, which is 1")
})

test_that(".get_s.d.denom() reports what it assumed", {
  tb <- t_bin()

  #An inferred denominator is announced, but only when it is not a treatment group:
  #reading "treated" off ATT weights is unambiguous, guessing "pooled" is not.
  expect_msg(sdd(treat = tb), 'Note: `s.d.denom` not specified; assuming "pooled".')
  expect_no_message(sdd(treat = tb, weights = data.frame(Adj = w_att())))

  #With several sets of weights the note names each one.
  expect_msg(sdd(treat = tb, weights = data.frame(ATT = w_att(), ATE = w_ate())),
             'Note: `s.d.denom` not specified; assuming "treated" for `ATT` and "pooled" for `ATE`.')

  #An explicitly given denominator is never announced.
  expect_no_message(sdd("pooled", treat = tb))
  expect_no_message(sdd(estimand = "ATE", treat = tb))

  #`weighted` needs weights to mean anything; without them it is the same as `all`,
  #and the note says so.
  expect_msg(sdd("weighted", treat = tb),
             'Note: `s.d.denom` specified as "weighted", but no weights supplied; setting to "all".')
  expect_no_message(sdd("weighted", treat = tb, weights = data.frame(Adj = w_ate())))

  #`quietly` suppresses all of it.
  expect_no_message(sdd(treat = tb, quietly = TRUE))
  expect_no_warning(sdd(estimand = "ATTT", treat = tb, quietly = TRUE))
})

test_that('"weighted" without weights computes the same factor as "all"', {
  # This is why the note above can claim it is "setting to \"all\"" without the
  # returned value changing: `.compute_s.d.denom()` reaches the same number either way.
  m <- as.matrix(lalonde[c("age", "educ")])
  tb <- t_bin()

  expect_equal(cobalt:::.compute_s.d.denom(m, tb, s.d.denom = "weighted",
                                           bin.vars = c(FALSE, FALSE)),
               cobalt:::.compute_s.d.denom(m, tb, s.d.denom = "all",
                                           bin.vars = c(FALSE, FALSE)))
})

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

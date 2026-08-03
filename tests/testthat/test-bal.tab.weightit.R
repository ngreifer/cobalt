skip_on_cran()
test_that("bal.tab() works with WeightIt", {
    data("lalonde")
    cov_names <- c("age", "educ", "race", "married", "nodegree", "re74", 
                   "re75")
    sw <- runif(nrow(lalonde))
    
    w <- WeightIt::weightit(treat ~age + educ + race + married + nodegree + re74 + re75,
                            data = lalonde)
    
    expect_s3_class(bal.tab(w), "bal.tab")
    
    w <- WeightIt::weightit(race ~age + educ + married + nodegree + re74 + re75,
                            data = lalonde)
    
    expect_s3_class(bal.tab(w), "bal.tab.multi")
    
    w <- WeightIt::weightit(re75 ~age + educ + married + nodegree + re74,
                            data = lalonde, density = "kernel")
    
    expect_s3_class(bal.tab(w), "bal.tab.cont")
    
    w <- WeightIt::weightit(treat ~ age + splines::ns(educ, 3) + race + married + nodegree + re74 + re75,
                            data = lalonde)

    expect_s3_class(bal.tab(w), "bal.tab")
})

# ---------------------------------------------------------------------------
# The block above checks dispatch. The blocks below check that the object method
# extracts the right pieces out of the fitted object.
#
# `bal.tab.weightit()` is a thin adapter: it pulls `weights`, `s.weights`, `covs`,
# `treat`, and `ps` off the object, derives `s.d.denom` from `estimand`, and hands
# everything to the data.frame method. So the strongest possible test is that the
# two routes agree *exactly* -- any refactor that mis-wires the adapter breaks it,
# and it needs no hardcoded numbers.

skip_if_not_installed("WeightIt")

wi_covs <- c("age", "educ", "race", "married", "nodegree", "re74", "re75")

test_that("bal.tab.weightit() agrees with the data.frame method it delegates to", {
    data("lalonde")

    w <- WeightIt::weightit(reformulate(wi_covs, "treat"), data = lalonde)

    b_obj <- bal.tab(w, un = TRUE, stats = c("m", "v", "ks"), disp = c("m", "sd"))

    #`weightit()` names its propensity score `prop.score`, so name it the same way
    #on the manual route for the tables to line up row for row.
    b_man <- bal.tab(lalonde[wi_covs], treat = lalonde$treat, weights = w$weights,
                     s.d.denom = "pooled", distance = data.frame(prop.score = w$ps),
                     un = TRUE, stats = c("m", "v", "ks"), disp = c("m", "sd"))

    expect_equal(b_obj$Balance, b_man$Balance)
    expect_equal(b_obj$Observations, b_man$Observations)

    expect_identical(rownames(b_obj$Balance)[1L], "prop.score")
    expect_identical(b_obj$Balance[["Type"]][1L], "Distance")
})

test_that("bal.tab.weightit() derives `s.d.denom` from the estimand", {
    data("lalonde")

    #This mapping is the one piece of logic the adapter adds, and getting it wrong
    #changes every standardized difference silently rather than erroring.
    expected_denom <- c(ATE = "pooled", ATT = "treated", ATC = "control")

    for (estimand in names(expected_denom)) {
        w <- WeightIt::weightit(reformulate(wi_covs, "treat"), data = lalonde,
                                estimand = estimand)
        expect_identical(w$estimand, estimand)

        b_obj <- bal.tab(w, un = TRUE)
        b_man <- bal.tab(lalonde[wi_covs], treat = lalonde$treat, weights = w$weights,
                         s.d.denom = expected_denom[[estimand]],
                         distance = data.frame(prop.score = w$ps), un = TRUE)

        expect_equal(b_obj$Balance, b_man$Balance, info = estimand)

        #And the three denominators really are different, so the check above has
        #something to catch.
        if (estimand != "ATE") {
            b_wrong <- bal.tab(lalonde[wi_covs], treat = lalonde$treat,
                               weights = w$weights, s.d.denom = "pooled",
                               distance = data.frame(prop.score = w$ps), un = TRUE)
            expect_false(isTRUE(all.equal(b_obj$Balance$Diff.Adj,
                                          b_wrong$Balance$Diff.Adj)))
        }
    }
})

test_that("bal.tab.weightit() carries the sampling weights through", {
    data("lalonde")

    set.seed(2024)
    sw <- runif(nrow(lalonde))

    w <- WeightIt::weightit(reformulate(wi_covs, "treat"), data = lalonde,
                            s.weights = sw)
    expect_equal(w$s.weights, sw)

    b_obj <- bal.tab(w, un = TRUE, stats = c("m", "ks"))
    b_man <- bal.tab(lalonde[wi_covs], treat = lalonde$treat, weights = w$weights,
                     s.weights = sw, s.d.denom = "pooled",
                     distance = data.frame(prop.score = w$ps), un = TRUE,
                     stats = c("m", "ks"))

    expect_equal(b_obj$Balance, b_man$Balance)

    #Both sample-size rows become effective sample sizes once `s.weights` are in
    #play; dropping them on the floor would leave the unadjusted row at 429/185.
    ess <- function(x) sum(x)^2 / sum(x^2)
    expect_equal(b_obj$Observations[["Control"]],
                 c(ess(sw[lalonde$treat == 0]),
                   ess((w$weights * sw)[lalonde$treat == 0])))
    expect_false(isTRUE(all.equal(b_obj$Observations[["Control"]][1L],
                                  sum(lalonde$treat == 0))))
})

test_that("bal.tab.weightit() works for continuous and multi-category treatments", {
    data("lalonde")

    #Continuous: correlations, no propensity score, no `s.d.denom` to derive.
    cont_covs <- c("age", "educ", "married", "nodegree", "re74")

    wc <- WeightIt::weightit(reformulate(cont_covs, "re75"), data = lalonde,
                             density = "kernel")
    expect_null(wc$ps)

    bc_obj <- bal.tab(wc, un = TRUE)
    bc_man <- bal.tab(lalonde[cont_covs], treat = lalonde$re75,
                      weights = wc$weights, un = TRUE)

    expect_identical(colnames(bc_obj$Balance), c("Type", "Corr.Un", "Corr.Adj"))
    expect_equal(bc_obj$Balance, bc_man$Balance)
    expect_equal(bc_obj$Balance$Corr.Adj,
                 unname(col_w_corr(lalonde[cont_covs], treat = lalonde$re75,
                                   weights = wc$weights)))

    #Multi-category: one table per pair plus the across-pairs summary.
    multi_covs <- c("age", "educ", "married", "nodegree", "re74", "re75")

    wm <- WeightIt::weightit(reformulate(multi_covs, "race"), data = lalonde)

    bm_obj <- bal.tab(wm, un = TRUE, stats = c("m", "ks"))
    bm_man <- bal.tab(lalonde[multi_covs], treat = lalonde$race,
                      weights = wm$weights, s.d.denom = "pooled", un = TRUE,
                      stats = c("m", "ks"))

    expect_identical(names(bm_obj$Pair.Balance), names(bm_man$Pair.Balance))
    expect_equal(bm_obj$Balance.Across.Pairs, bm_man$Balance.Across.Pairs)

    for (p in names(bm_obj$Pair.Balance)) {
        expect_equal(bm_obj$Pair.Balance[[p]]$Balance,
                     bm_man$Pair.Balance[[p]]$Balance, info = p)
    }
})

test_that("bal.tab.weightit() names non-syntactic covariates as supplied", {
    data("lalonde")

    #A spline term arrives as a matrix column; its name must survive intact rather
    #than being mangled or silently dropped.
    w <- WeightIt::weightit(treat ~ age + splines::ns(educ, 3) + race, data = lalonde)

    b <- bal.tab(w, un = TRUE)

    expect_identical(rownames(b$Balance),
                     c("prop.score", "age", "educ", "race_black", "race_hispan",
                       "race_white"))

    #`weightit()` reduces the spline basis back to `educ` in `$covs`, so the row is
    #the balance of the original variable.
    expect_equal(b$Balance["educ", "Diff.Adj"],
                 unname(col_w_smd(lalonde["educ"], treat = lalonde$treat,
                                  weights = w$weights, s.d.denom = "pooled",
                                  std = TRUE)))
})

test_that("bal.tab() accepts the weights out of a weightitMSM object", {
    data("msmdata", package = "WeightIt")

    wm <- WeightIt::weightitMSM(
        list(A_1 ~ X1_0 + X2_0,
             A_2 ~ X1_1 + X2_1 + A_1 + X1_0 + X2_0,
             A_3 ~ X1_2 + X2_2 + A_2 + X1_1 + X2_1 + A_1 + X1_0 + X2_0),
        data = msmdata, method = "glm")

    b <- bal.tab(wm, un = TRUE)

    expect_s3_class(b, "bal.tab.msm")
    expect_length(b$Time.Balance, 3L)

    #Each period's table must equal a standalone `bal.tab()` on that period's
    #covariates with the same MSM weights.
    b1 <- bal.tab(msmdata[c("X1_0", "X2_0")], treat = msmdata$A_1,
                  weights = wm$weights, s.d.denom = "pooled", un = TRUE)

    expect_equal(b$Time.Balance[[1L]]$Balance[c("X1_0", "X2_0"), ], b1$Balance)
})
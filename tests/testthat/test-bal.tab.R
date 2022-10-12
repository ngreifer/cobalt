covs <- subset(lalonde, select = -c(re78, treat))
lalonde_ <- splitfactor(lalonde, "race", drop.first = F)
covs_ <- subset(lalonde_, select = -c(re78, treat))

test_that("No adjustment", {
    data("lalonde", package = "cobalt")
    
    expect_snapshot(
        bal.tab(
            covs,
            treat = lalonde$treat,
            m.threshold = .1,
            v.threshold = 2,
            imbalanced.only = T
        )
    )
    expect_snapshot(bal.tab(
        f.build("treat", covs),
        data = lalonde,
        ks.threshold = .1
    ))
})

test_that("MatchIt with PS", {
    skip_if_not_installed("MatchIt")
    m1 <- MatchIt::matchit(
            treat ~ log(age) * married + educ + race + nodegree + re74 + re75,
            data = lalonde,
            replace = T,
            ratio = 2,
            discard = "both",
            reestimate = TRUE
        )
    expect_snapshot(
        bal.tab(
            m1,
            int = T,
            v.threshold = 2,
            imbalanced.only = F
        )
    )
})

test_that("MatchIt with distance", {
    skip_if_not_installed("MatchIt")
    set.seed(123)
    m1.1 <- MatchIt::matchit(
        f.build("treat", covs),
        data = lalonde,
        distance = runif(nrow(lalonde))
    )
    expect_snapshot(bal.tab(m1.1, data = lalonde))
})

test_that("MatchIt with complicated formula", {
    skip_if_not_installed("MatchIt")
    skip_if_not_installed("rms")
    m1.2 <- MatchIt::matchit(
        treat ~ log(age)*married + poly(educ,2) + rms::rcs(age,3) + race + factor(nodegree) + re74 + re75, 
        data = lalonde, 
        replace = T, 
        ratio = 2
    )
    
    expect_snapshot(bal.tab(m1.2, addl = ~ rms::rcs(educ,3) + poly(age,2)))
})

test_that("MatchIt with Malahanobis", {
    skip_if_not_installed("MatchIt")
    m2 <- MatchIt::matchit(f.build("treat", covs), data = lalonde, distance = "mahalanobis")
    
    expect_snapshot(
        bal.tab(m2, data = lalonde, int = T, v.threshold = 2, addl = "race")
    )
})

test_that("MatchIt with subclassification", {
    skip_if_not_installed("MatchIt")
    m3 <- MatchIt::matchit(f.build("treat", covs), data = lalonde, method = "subclass")
    
    expect_snapshot(
        bal.tab(m3, int = F, v.threshold = 2, disp.subclass = T, ks.threshold = .1)
    )
})

test_that("MatchIt: full matching", {
    skip_if_not_installed("MatchIt")
    m4 <- MatchIt::matchit(f.build("treat", covs), data = lalonde, method = "full", distance = "glm", link = "probit", estimand = "ATE")
    expect_snapshot(
        bal.tab(m4, int = T, ks.threshold = .05)
    )
    expect_snapshot(
        bal.tab(m4, pairwise = FALSE, un = T)
    )
})

test_that("MatchIt: genetic matching, using matching weights", {
    skip_if_not_installed("MatchIt")
    skip_if_not_installed("rgenoud")
    m5 <- MatchIt::matchit(f.build("treat", covs), data = lalonde, method = "genetic", replace = F,
                  ratio = 1, pop.size = 50)
   expect_snapshot(
       bal.tab(m5, method = "m", estimand = "ATT")
   )
})

skip_on_cran()
test_that("bal.tab() works with data.frames", {
    data("lalonde")
    cov_names <- c("age", "educ", "race", "married", "nodegree", "re74", 
                   "re75")
    covs <- lalonde[cov_names]
    w <- runif(nrow(lalonde))
    sw <- runif(nrow(lalonde))
    dist <- runif(nrow(lalonde))
    
    expect_s3_class(bal.tab(covs, treat = lalonde$treat), "bal.tab")
    expect_s3_class(bal.tab(covs, treat = lalonde$treat, distance = dist,
                            weights = w, s.weights = sw), "bal.tab")
    expect_s3_class(bal.tab(covs, treat = lalonde$treat, distance = dist,
                            weights = w, s.weights = sw, cluster = lalonde$race), "bal.tab.cluster")
    expect_s3_class(bal.tab(covs, treat = lalonde$race, distance = dist,
                            weights = w, s.weights = sw), "bal.tab.multi")
})
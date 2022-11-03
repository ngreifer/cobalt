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
                            data = lalonde)
    
    expect_s3_class(bal.tab(w), "bal.tab.cont")
    
    w <- WeightIt::weightit(treat ~ age + splines::ns(educ, 3) + race + married + nodegree + re74 + re75,
                            data = lalonde)
    
    expect_s3_class(bal.tab(w), "bal.tab")
})
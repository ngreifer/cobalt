skip_on_cran()
test_that("bal.compute() and bal.init() work", {
    data("lalonde")
    cov_names <- c("age", "educ", "race", "married", "nodegree", "re74", 
                   "re75")
    sw <- runif(nrow(lalonde))
    
    tb <- lalonde$treat
    tm <- factor(sample(LETTERS[1:4], nrow(lalonde), TRUE))
    tc <- rnorm(nrow(lalonde))
    
    # Binary treatment
    init <- bal.init("smd.mean", tb, lalonde[cov_names])
    expect_s3_class(init, "bal.init")
    
    
    # Multi-category treatment
    
    # Continuous treatment
    

    
})

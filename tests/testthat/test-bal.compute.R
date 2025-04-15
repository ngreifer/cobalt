skip_on_cran()
test_that("bal.compute() and bal.init() work", {
    set.seed(100)
    data("lalonde")
    cov_names <- c("age", "educ", "race", "married", "nodegree", "re74", 
                   "re75")
    sw <- runif(nrow(lalonde))
    
    tb <- lalonde$treat
    tm <- factor(sample(LETTERS[1:4], nrow(lalonde), TRUE))
    tc <- rnorm(nrow(lalonde))
    
    # Binary treatment
    init <- bal.init(lalonde[cov_names], tb, "smd.mean")
    expect_s3_class(init, "bal.init")
    expect_equal(bal.compute(init), 0.60909674)
    
    
    # Multi-category treatment
    init <- bal.init(lalonde[cov_names], tb, "smd.mean")
    expect_s3_class(init, "bal.init")

    # Continuous treatment
    init <- bal.init(lalonde[cov_names], tc, "p.mean")
    expect_s3_class(init, "bal.init")
    expect_equal(bal.compute(init), 0.033749304)

    
})

test_that("bal.compute() returns 0 for unweighted target balance statistics", {
  set.seed(1004)
  data("lalonde")
  cov_names <- c("age", "educ", "race", "married", "nodegree", "re74", 
                 "re75")
  sw <- runif(nrow(lalonde))
  
  expect_equal(bal.compute(lalonde[cov_names], stat = "smd.rms"),
               0)
  expect_equal(bal.compute(lalonde[cov_names], stat = "ks.max"),
               0)
  # expect_equal(bal.compute(lalonde[cov_names], stat = "ovl.max"),
  #              0)
  expect_equal(bal.compute(lalonde[cov_names], stat = "energy.dist"),
               0)
  expect_equal(bal.compute(lalonde[cov_names], stat = "mahalanobis"),
               0)
  
  expect_equal(bal.compute(lalonde[cov_names], stat = "smd.rms", s.weights = sw),
               0)
  expect_equal(bal.compute(lalonde[cov_names], stat = "ks.max", s.weights = sw),
               0)
  # expect_equal(bal.compute(lalonde[cov_names], stat = "ovl.max", s.weights = sw),
  #              0)
  expect_equal(bal.compute(lalonde[cov_names], stat = "energy.dist", s.weights = sw),
               0)
  expect_equal(bal.compute(lalonde[cov_names], stat = "mahalanobis", s.weights = sw),
               0)
})
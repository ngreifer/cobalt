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

test_that("bal.compute() works with pairwise = FALSE", {
  eps <- if (capabilities("long.double")) 1e-8 else 1e-1
  
  set.seed(100)
  data("lalonde")
  cov_names <- c("age", "educ", "race", "married", "nodegree", "re74", 
                 "re75")
  sw <- runif(nrow(lalonde))
  
  tb <- lalonde$treat
  tm <- factor(sample(LETTERS[1:4], nrow(lalonde), TRUE))

  w <- runif(nrow(lalonde))
  
  # Binary treatment
  ## SMD
  init <- bal.init(lalonde[cov_names], tb, "smd.mean")
  expect_s3_class(init, "bal.init")
  
  baltab <- bal.tab(lalonde[cov_names], treat = tb, binary = "std",
                    stats = "m", estimand = "ATE", weights = w,
                    un = TRUE)$Balance
  
  expect_equal(bal.compute(init),
               mean(abs(baltab$Diff.Un)),
               tolerance = eps)
  
  expect_equal(bal.compute(init, weights = w),
               mean(abs(baltab$Diff.Adj)),
               tolerance = eps)
  
  init <- bal.init(lalonde[cov_names], treat = tb,
                   stat = "smd.mean",
                   pairwise = FALSE)
  expect_s3_class(init, "bal.init")
  
  baltab <- bal.tab(lalonde[cov_names], treat = tb, binary = "std",
                    stats = "m", estimand = "ATE", weights = w,
                    un = TRUE, pairwise = FALSE)$Pair.Balance
  
  expect_equal(bal.compute(init),
               mean(abs(unlist(lapply(baltab, function(i) i$Balance$Diff.Un)))),
               tolerance = eps)
  
  expect_equal(bal.compute(init, weights = w),
               mean(abs(unlist(lapply(baltab, function(i) i$Balance$Diff.Adj)))),
               tolerance = eps)
  
  ## KS
  init <- bal.init(lalonde[cov_names], tb, "ks.mean")
  expect_s3_class(init, "bal.init")
  
  baltab <- bal.tab(lalonde[cov_names], treat = tb, binary = "std",
                    stats = "ks", estimand = "ATE", weights = w,
                    un = TRUE)$Balance
  
  expect_equal(bal.compute(init),
               mean(abs(baltab$KS.Un)),
               tolerance = eps)
  
  expect_equal(bal.compute(init, weights = w),
               mean(abs(baltab$KS.Adj)),
               tolerance = eps)
  
  init <- bal.init(lalonde[cov_names], treat = tb,
                   stat = "ks.mean",
                   pairwise = FALSE)
  expect_s3_class(init, "bal.init")
  
  baltab <- bal.tab(lalonde[cov_names], treat = tb, binary = "std",
                    stats = "ks", estimand = "ATE", weights = w,
                    un = TRUE, pairwise = FALSE)$Pair.Balance
  
  expect_equal(bal.compute(init),
               mean(abs(unlist(lapply(baltab, function(i) i$Balance$KS.Un)))),
               tolerance = eps)
  
  expect_equal(bal.compute(init, weights = w),
               mean(abs(unlist(lapply(baltab, function(i) i$Balance$KS.Adj)))),
               tolerance = eps)
  
  # Multi-category treatment
  ## SMD
  init <- bal.init(lalonde[cov_names], tm, "smd.mean")
  expect_s3_class(init, "bal.init")
  
  baltab <- bal.tab(lalonde[cov_names], treat = tm, binary = "std",
                    stats = "m", estimand = "ATE", weights = w,
                    un = TRUE)$Pair.Balance
  
  expect_equal(bal.compute(init),
               mean(abs(unlist(lapply(baltab, function(i) i$Balance$Diff.Un)))),
               tolerance = eps)
  
  expect_equal(bal.compute(init, weights = w),
               mean(abs(unlist(lapply(baltab, function(i) i$Balance$Diff.Adj)))),
               tolerance = eps)
  
  init <- bal.init(lalonde[cov_names], treat = tm,
                   stat = "smd.mean",
                   pairwise = FALSE)
  expect_s3_class(init, "bal.init")
  
  baltab <- bal.tab(lalonde[cov_names], treat = tm, binary = "std",
                    stats = "m", estimand = "ATE", weights = w,
                    un = TRUE, pairwise = FALSE)$Pair.Balance
  
  expect_equal(bal.compute(init),
               mean(abs(unlist(lapply(baltab, function(i) i$Balance$Diff.Un)))),
               tolerance = eps)
  
  expect_equal(bal.compute(init, weights = w),
               mean(abs(unlist(lapply(baltab, function(i) i$Balance$Diff.Adj)))),
               tolerance = eps)
  
  ## KS
  init <- bal.init(lalonde[cov_names], tm, "ks.mean")
  expect_s3_class(init, "bal.init")
  
  baltab <- bal.tab(lalonde[cov_names], treat = tm, binary = "std",
                    stats = "ks", estimand = "ATE", weights = w,
                    un = TRUE)$Pair.Balance
  
  expect_equal(bal.compute(init),
               mean(abs(unlist(lapply(baltab, function(i) i$Balance$KS.Un)))),
               tolerance = eps)
  
  expect_equal(bal.compute(init, weights = w),
               mean(abs(unlist(lapply(baltab, function(i) i$Balance$KS.Adj)))),
               tolerance = eps)
  
  init <- bal.init(lalonde[cov_names], treat = tm,
                   stat = "ks.mean",
                   pairwise = FALSE)
  expect_s3_class(init, "bal.init")
  
  baltab <- bal.tab(lalonde[cov_names], treat = tm, binary = "std",
                    stats = "ks", estimand = "ATE", weights = w,
                    un = TRUE, pairwise = FALSE)$Pair.Balance
  
  expect_equal(bal.compute(init),
               mean(abs(unlist(lapply(baltab, function(i) i$Balance$KS.Un)))),
               tolerance = eps)
  
  expect_equal(bal.compute(init, weights = w),
               mean(abs(unlist(lapply(baltab, function(i) i$Balance$KS.Adj)))),
               tolerance = eps)
})
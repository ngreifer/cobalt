test_that("splictfactor() and unsplitfactor() work", {
  set.seed(123)
  n <- 100
  d <- data.frame(x1 = factor(sample(c("A", "B"), n, TRUE)),
                  x2 = factor(sample(c("A", "B", "C"), n, TRUE)),
                  x3 = rbinom(n, 1, .5),
                  x4a = factor(rep("A", n), levels = c("A")),
                  x4b = factor(rep("A", n), levels = c("A", "B")),
                  x5 = factor(sample(c("A", "B", "C", "D"), n, TRUE)))
  
  
  expect_no_condition({
    d_s0 <- splitfactor(d)
  })
  expect_identical(names(d_s0),
                   c("x1_B", "x2_B", "x2_C", "x3", "x4a_A", "x4b_A", "x5_B", "x5_C", 
                     "x5_D"))
  
  expect_no_condition({
    d_s <- splitfactor(d, drop.first = TRUE)
  })
  expect_identical(d_s, d_s0)
  expect_no_condition({
    d_u <- unsplitfactor(d_s, dropped.level = "A")
  })
  expect_equal(d_u, d, ignore_attr = TRUE)
  
  expect_no_condition({
    d_s <- splitfactor(d, drop.first = FALSE)
  })
  expect_identical(names(d_s),
                   c("x1_A", "x1_B", "x2_A", "x2_B", "x2_C", "x3", "x4a_A", "x4b_A", 
                     "x5_A", "x5_B", "x5_C", "x5_D"))
  expect_no_condition({
    d_u <- unsplitfactor(d_s)
  })
  expect_equal(d_u, d, ignore_attr = TRUE)
  
  expect_no_condition({
    d_s <- splitfactor(d, drop.first = "if2")
  })
  expect_identical(names(d_s),
                   c("x1_B", "x2_A", "x2_B", "x2_C", "x3", "x4a_A", "x4b_A", "x5_A", 
                     "x5_B", "x5_C", "x5_D"))
  
  #
  expect_no_condition({
    d_s <- splitfactor(d, drop.first = TRUE, drop.singleton = TRUE)
  })
  expect_identical(names(d_s),
                   c("x1_B", "x2_B", "x2_C", "x3", "x5_B", "x5_C", "x5_D"))
  
  expect_no_condition({
    d_s <- splitfactor(d, drop.first = FALSE, drop.singleton = TRUE)
  })
  expect_identical(names(d_s),
                   c("x1_A", "x1_B", "x2_A", "x2_B", "x2_C", "x3", "x5_A", "x5_B", 
                     "x5_C", "x5_D"))
  expect_no_condition({
    d_u <- unsplitfactor(d_s)
  })
  expect_equal(d_u, d[-(4:5)], ignore_attr = TRUE)
  
  expect_no_condition({
    d_s <- splitfactor(d, drop.first = "if2", drop.singleton = TRUE)
  })
  expect_identical(names(d_s),
                   c("x1_B", "x2_A", "x2_B", "x2_C", "x3", "x5_A", "x5_B", "x5_C", 
                     "x5_D"))
  
  expect_no_condition({
    d_s <- splitfactor(d, drop.first = "if2", drop.singleton = TRUE,
                       sep = "|")
  })
  expect_identical(names(d_s),
                   c("x1|B", "x2|A", "x2|B", "x2|C", "x3", "x5|A", "x5|B", "x5|C", 
                     "x5|D"))
  
  #
  expect_no_condition({
    d_s <- splitfactor(d, drop.first = "if2", drop.singleton = TRUE,
                       replace = FALSE)
  })
  expect_identical(names(d_s),
                   c("x1", "x2", "x3", "x5", "x1_B", "x2_A", "x2_B", "x2_C", "x5_A", 
                     "x5_B", "x5_C", "x5_D"))
  
  #
  expect_no_condition({
    d_s <- splitfactor(d, "x2")
  })
  expect_identical(names(d_s),
                   c("x1", "x2_B", "x2_C", "x3", "x4a", "x4b", "x5"))
  expect_no_condition({
    d_u <- unsplitfactor(d_s, "x2", dropped.level = levels(d$x2)[1L])
  })
  expect_equal(d_u, d, ignore_attr = TRUE)
  
  expect_no_condition({
    d_s <- splitfactor(d, "x2", drop.level = "C")
  })
  expect_identical(names(d_s),
                   c("x1", "x2_A", "x2_B", "x3", "x4a", "x4b", "x5"))
  expect_no_condition({
    d_u <- unsplitfactor(d_s, "x2", dropped.level = "C")
  })
  expect_equal(d_u |> transform(x2 = as.character(x2)),
               d |> transform(x2 = as.character(x2)),
               ignore_attr = TRUE)
  
  expect_no_condition({
    d_s <- splitfactor(d, "x2", drop.first = "if2")
  })
  expect_identical(names(d_s),
                   c("x1", "x2_A", "x2_B", "x2_C", "x3", "x4a", "x4b", "x5"))
  expect_no_condition({
    d_u <- unsplitfactor(d_s, "x2")
  })
  expect_equal(d_u, d, ignore_attr = TRUE)
  
  expect_no_condition({
    d_s <- splitfactor(d, "x1", drop.first = FALSE)
  })
  expect_identical(names(d_s),
                   c("x1_A", "x1_B", "x2", "x3", "x4a", "x4b", "x5"))
  expect_no_condition({
    d_u <- unsplitfactor(d_s, "x1")
  })
  expect_equal(d_u, d, ignore_attr = TRUE)
  
  expect_no_condition({
    d_s <- splitfactor(d, c("x1", "x2"), drop.first = FALSE)
  })
  expect_identical(names(d_s),
                   c("x1_A", "x1_B", "x2_A", "x2_B", "x2_C", "x3", "x4a", "x4b", 
                     "x5"))
  expect_no_condition({
    d_u <- unsplitfactor(d_s)
  })
  expect_equal(d_u, d, ignore_attr = TRUE)
  
  expect_no_condition({
    d_s <- splitfactor(d, c("x1", "x2"), drop.first = "if2")
  })
  expect_identical(names(d_s),
                   c("x1_B", "x2_A", "x2_B", "x2_C", "x3", "x4a", "x4b", "x5"))
  expect_no_condition({
    d_u <- unsplitfactor(d_s, dropped.level = "A")
  })
  expect_equal(d_u, d, ignore_attr = TRUE)
  
  #Bad inputs
  expect_warning({
    d_s <- splitfactor(d, c("x1", "x2", "bad"), drop.first = "if2")
  }, .w('"bad" is not the name of a factor variable in `data` and will not be split.'))
  expect_identical(names(d_s),
                   c("x1_B", "x2_A", "x2_B", "x2_C", "x3", "x4a", "x4b", "x5"))
  
  expect_warning({
    d_s <- splitfactor(d, c("x1", "x2", "x3"), drop.first = "if2")
  }, .w('"x3" is not the name of a factor variable in `data` and will not be split.'))
  expect_identical(names(d_s),
                   c("x1_B", "x2_A", "x2_B", "x2_C", "x3", "x4a", "x4b", "x5"))
  
  expect_warning({
    d_s <- splitfactor(d, c("x1", "x2", "x3", "bad"), drop.first = "if2")
  }, .w('"x3" and "bad" are not the names of factor variables in `data` and will not be split.'))
  expect_identical(names(d_s),
                   c("x1_B", "x2_A", "x2_B", "x2_C", "x3", "x4a", "x4b", "x5"))
  
  expect_error({
    d_s <- splitfactor(d, 1)
  }, .w("`var.name` must be a character vector of the names of one or more factor variables in `data`."))
  
  expect_error({
    d_s <- splitfactor(d, "bad_x")
  }, .w("No names in `var.name` are names of factor variables in `data`."))
  
  expect_warning({
    d_s <- splitfactor(d, c("x1", "x2"), drop.level = "B")
  }, .w("`drop.level` cannot be used with multiple entries to `var.name`. Ignoring `drop.level`."))
  expect_identical(names(d_s),
                   c("x1_B", "x2_B", "x2_C", "x3", "x4a", "x4b", "x5"))
  
  expect_error({
    d_s <- splitfactor(d, drop.first = 3)
  }, .w('`drop.first` must be `TRUE`, `FALSE`, or "if2".'))
  
  expect_no_condition({
    d_s <- splitfactor(d, "x1", drop.level = "B", drop.first = 3)
  })
  
  m <- matrix(rbinom(5 * n, 1, .5), nrow = n)
  
  expect_warning({
    d_s <- splitfactor(m)
  }, .w('There are no factor variables to split in `data`.'))
  expect_s3_class(d_s, "data.frame")
  
  m[] <- as.character(m)
  
  expect_no_condition({
    d_s <- splitfactor(m, drop.first = FALSE)
  })
  expect_s3_class(d_s, "data.frame")
  expect_identical(names(d_s),
                   c("V1_0", "V1_1", "V2_0", "V2_1", "V3_0", "V3_1", "V4_0", "V4_1", 
                     "V5_0", "V5_1"))
  
})

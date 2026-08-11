#Tests for specific github issues
skip_on_cran()
test_that("(#71) `addl` argument not throwing correct error when variable not available", {
  #Note: actual error was something else, but addl should have thrown a more informative error
  skip_if_not_installed("MatchIt")
  
  data("lalonde")
  
  m <- MatchIt::matchit(treat ~ age + educ, data = lalonde,
                        distance = "scaled_euclidean")
  
  expect_err(bal.tab(m, addl = "race"), "The variable \"race\" cannot be found. Be sure it is entered correctly or supply a dataset that contains this variable to `data`.")

  expect_no_condition(b <- bal.tab(m, addl = "race", data = lalonde))

  #Supplying `data` must actually add the variable, not just stop erroring. A
  #`scaled_euclidean` distance is not a propensity score, so no distance row.
  expect_identical(rownames(b$Balance),
                   c("age", "educ", "race_black", "race_hispan", "race_white"))

  #And the added rows must hold the right numbers: matching weights on the race
  #dummies of the full data set.
  w <- MatchIt::match.data(m, data = lalonde, weights = ".w", drop.unmatched = FALSE)[[".w"]]

  expect_equal(b$Balance[c("race_black", "race_hispan", "race_white"), "Diff.Adj"],
               unname(col_w_smd(splitfactor(lalonde["race"], drop.first = FALSE),
                                treat = lalonde$treat, weights = w,
                                s.d.denom = "treated", std = FALSE)))
})

test_that("(#74) bal.tab() works with mnps objects", {
  skip_if_not_installed("twang")
  
  data("lalonde")
  
  mnps <- twang::mnps(race ~ age + educ + married + re74, 
                      data = lalonde,
                      estimand = "ATT",
                      treatATT = "white",
                      verbose = FALSE,
                      stop.method = c("es.mean"), 
                      n.trees = 1000)
  
  expect_no_condition(
    b <- bal.tab(mnps, un = TRUE, disp = c("m", "sd"))
  )

  #Every pair of treatment levels gets a table, regardless of which one is focal.
  expect_s3_class(b, "bal.tab.multi")
  expect_identical(names(b$Pair.Balance),
                   c("hispan vs. black", "white vs. black", "white vs. hispan"))

  #`disp = c("m", "sd")` must produce mean and SD columns for both groups.
  expect_true(all(c("M.0.Un", "SD.0.Un", "M.1.Un", "SD.1.Un",
                    "M.0.Adj", "SD.0.Adj", "M.1.Adj", "SD.1.Adj") %in%
                      colnames(b$Pair.Balance[[1L]]$Balance)))

  #The weights must be the ones the stop method chose, not uniform ones: the
  #unadjusted and adjusted columns must differ.
  expect_false(isTRUE(all.equal(b$Pair.Balance[[1L]]$Balance$Diff.Un,
                                b$Pair.Balance[[1L]]$Balance$Diff.Adj)))

  #`get.w()` must return those same weights, and `bal.tab()` on them must
  #reproduce the object method's table. `treatATT = "white"` makes white the focal
  #group, which is what the object method infers.
  w <- get.w(mnps)
  expect_length(w, nrow(lalonde))

  b_man <- bal.tab(lalonde[c("age", "educ", "married", "re74")],
                   treat = lalonde$race, weights = w, s.d.denom = "focal",
                   focal = "white", un = TRUE, disp = c("m", "sd"))

  for (p in names(b$Pair.Balance)) {
    expect_equal(b$Pair.Balance[[p]]$Balance, b_man$Pair.Balance[[p]]$Balance,
                 info = p)
  }

  #The focal group's weights are all 1 under an ATT estimand.
  expect_true(all(w[lalonde$race == "white"] == 1))
})

test_that("(#76) bal.tab() doesn't produce an error with missing covariates", {
  
  data("lalonde_mis")
  
  expect_wrn(
    b <- bal.tab(treat ~ age + educ + race + married + nodegree + re74 + re75,
                 data = lalonde_mis, s.d.denom = "pooled"), "Missing values exist in the covariates. Displayed values omit these observations.")
  
  expect_no_error(
    suppressWarnings(
      b <- bal.tab(treat ~ age + educ + race + married + nodegree + re74 + re75,
                   data = lalonde_mis, s.d.denom = "pooled")
    )
  )
  
  expect_identical(rownames(b$Balance),
                   c("age", "educ", "race_black", "race_hispan", "race_white",
                     "married", "married:<NA>", "nodegree",
                     "re74", "re74:<NA>", "re75", "re75:<NA>"))

  #The `:<NA>` rows are missingness indicators, so their group means are the
  #proportions missing. Getting these wrong would silently misreport missingness.
  b2 <- suppressWarnings(
    bal.tab(treat ~ age + educ + race + married + nodegree + re74 + re75,
            data = lalonde_mis, s.d.denom = "pooled", disp = "m", un = TRUE)
  )

  for (v in c("married", "re74", "re75")) {
    row <- paste0(v, ":<NA>")

    expect_identical(b2$Balance[row, "Type"], "Binary")
    expect_equal(b2$Balance[row, "M.0.Un"],
                 mean(is.na(lalonde_mis[[v]][lalonde_mis$treat == 0])), info = row)
    expect_equal(b2$Balance[row, "M.1.Un"],
                 mean(is.na(lalonde_mis[[v]][lalonde_mis$treat == 1])), info = row)
  }

  #The observed-value rows omit the missing observations rather than treating them
  #as zeroes, which is what the warning promises.
  expect_equal(b2$Balance["re74", "M.0.Un"],
               mean(lalonde_mis$re74[lalonde_mis$treat == 0], na.rm = TRUE))

  #A variable with no missingness gets no indicator row.
  expect_false("age:<NA>" %in% rownames(b2$Balance))
  expect_false(anyNA(lalonde_mis$age))
})

test_that("(#77) no conflict with `cluster` argument to bal.tab() when {caret} is loaded", {
  skip_if_not_installed("caret")
  
  suppressPackageStartupMessages(library(caret))
  
  data("lalonde")
  
  expect_no_condition(bal.tab(treat ~ age + educ + race, data = lalonde, cluster = "married",
                              s.d.denom = "pooled"))
})

test_that("(#82) love.plot() doesn't throw any ggplot2-related warnings with multiple stats", {
  skip_if_not_installed("MatchIt")
  
  data("lalonde")
  
  m.out <- MatchIt::matchit(treat ~ age + educ + race + married +
                              nodegree + re74 + re75,
                            data = lalonde)
  
  local_null_device()

  expect_no_condition({
    p <- love.plot(m.out,
                   stats = c("mean.diffs", "variance.ratios"),
                   binary = "std")
  })

  #With two statistics `love.plot()` returns an assembled gtable rather than a
  #single ggplot, and it must contain a panel for each.
  expect_s3_class(p, "love.plot")
  expect_s3_class(p, "gtable")

  #Drawing it is what would surface a deferred ggplot2 warning, so do that too.
  expect_no_condition(print(p))

  #One statistic still gives a plain ggplot, and the values it plots must be the
  #ones `bal.tab()` computed -- a plot that draws the wrong column is the failure
  #mode no "no condition" check can see.
  p1 <- love.plot(m.out, stats = "mean.diffs", binary = "std")

  expect_s3_class(p1, "ggplot")

  pd <- ggplot2::ggplot_build(p1)$plot$data
  b <- bal.tab(m.out, binary = "std", un = TRUE)

  for (samp in c(Unadjusted = "Diff.Un", Adjusted = "Diff.Adj")) {
    nm <- names(which(c(Unadjusted = "Diff.Un", Adjusted = "Diff.Adj") == samp))
    rows <- pd[pd$Sample == nm, ]

    expect_equal(rows$stat[match(rownames(b$Balance), as.character(rows$var))],
                 b$Balance[[samp]], info = nm)
  }
})

test_that("(#89) love.plot() doesn't throw any error when manually removing rows from bal.tab() while supplying var.names", {
  data("lalonde_mis")
  
  v <- data.frame(old = c("age", "educ", "race_black", "race_hispan", 
                          "race_white", "married", "nodegree", "re74", "re75", "distance"),
                  new = c("Age", "Years of Education", "Black", 
                          "Hispanic", "White", "Married", "No Degree Earned", 
                          "Earnings 1974", "Earnings 1975", "Propensity Score"))
  
  covs <- subset(lalonde_mis, select = -c(treat, re78, nodegree, married))
  
  expect_wrn({
    b <- bal.tab(treat ~ covs, data = lalonde_mis, binary = "std", estimand = "ATE")
  }, "Missing values exist in the covariates. Displayed values omit these observations")
  
  b$Balance <- b$Balance[!endsWith(rownames(b$Balance), "<NA>"),]

  local_null_device()

  expect_no_condition({
    p <- love.plot(b, var.names = v)
  })

  #`var.names` lists variables that are no longer in the table (`married`,
  #`nodegree`, `distance`); the rows that remain must be relabeled and the surplus
  #entries quietly ignored.
  labs <- levels(ggplot2::ggplot_build(p)$plot$data$var)

  expect_setequal(labs,
                  c("Age", "Years of Education", "Black", "Hispanic", "White",
                    "Earnings 1974", "Earnings 1975"))
  expect_length(labs, nrow(b$Balance))

  #No original name may survive the relabeling.
  expect_false(any(rownames(b$Balance) %in% labs))

  #A variable absent from `var.names` keeps its original name.
  b2 <- b
  p2 <- love.plot(b2, var.names = v[v$old != "age", ])

  expect_true("age" %in% levels(ggplot2::ggplot_build(p2)$plot$data$var))
})

#Tests for specific github issues
skip_on_cran()
test_that("(#71) `addl` argument not throwing correct error when variable not available", {
    #Note: actual error was something else, but addl should have thrown a more informative error
    skip_if_not_installed("MatchIt")
    
    data("lalonde")
    
    m <- MatchIt::matchit(treat ~ age + educ, data = lalonde, distance = "scaled_euclidean")
    
    expect_error(bal.tab(m, addl = "race"),
                 "The variable \"race\" cannot be found. Be sure it is entered correctly or supply a dataset that contains this varialble to `data`.",
                 fixed = TRUE)
    
    expect_no_condition(bal.tab(m, addl = "race", data = lalonde))
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
        bal.tab(mnps, un = TRUE, disp = c("m", "sd"))
    )
})

test_that("(#76) bal.tab() doesn't produce an error with missing covariates", {
    
    data("lalonde_mis")
    
    expect_warning(
        b <- bal.tab(treat ~ age + educ + race + married + nodegree + re74 + re75,
                     data = lalonde_mis, s.d.denom = "pooled"),
        "Missing values exist in the covariates. Displayed values omit these observations.",
        fixed = TRUE
    )
    
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
})

test_that("(#77) no conflict with `cluster` argument to bal.tab() when {caret} is loaded", {
    skip_if_not_installed("caret")
    
    library(caret)
    
    data("lalonde")
    
    expect_no_condition(bal.tab(treat ~ age + educ + race, data = lalonde, cluster = "married",
                                s.d.denom = "pooled"))
})

test_that("(#82) love.plot() doesn't throw any ggplo2-related warnings with multiple stats", {
    skip_if_not_installed("MatchIt")
    
    data("lalonde")
    
    m.out <- MatchIt::matchit(treat ~ age + educ + race + married +
                                  nodegree + re74 + re75,
                              data = lalonde)
    
    expect_no_condition(love.plot(m.out, stats = c("mean.diffs", "variance.ratios"),
                                  binary = "std"))
})

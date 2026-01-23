skip_on_cran()

test_that("`s.d.denom` processes correctly", {
    data("lalonde")
    
    weights <- list(wATE = runif(nrow(lalonde)),
                    wATT = ifelse(lalonde$treat == 1, 1, runif(sum(lalonde$treat != 1))))
    
    # Binary tretament
    cov_names <- c("age", "educ", "race", "married", "nodegree", "re74", "re75")
    
    f <- reformulate(cov_names, "treat")
    
    expect_message(bal.tab(f, data = lalonde, binary = "std", continuous = "std"),
                   .w('Note: `s.d.denom` not specified; assuming "pooled".'),
                   perl = TRUE)
    
    #No message when no variables are to be standardized
    expect_no_message(bal.tab(f, data = lalonde, binary = "raw", continuous = "raw"))
    
    #No message when SMDs not requested
    expect_no_message(bal.tab(f, data = lalonde, binary = "raw", continuous = "raw",
                              stats = "ks"))
    
    expect_message(bal.tab(f, data = lalonde, binary = "std", continuous = "std",
                           s.d.denom = "weighted"),
                   .w('Note: `s.d.denom` specified as "weighted", but no weights supplied; setting to "all".'),
                   perl = TRUE)
    
    expect_equal(bal.tab(f, data = lalonde, binary = "std", continuous = "std", s.d.denom = "pooled")$Balance$Diff.Un,
                 unname(col_w_smd(lalonde[cov_names], treat = lalonde$treat, s.d.denom = "pooled")))
    
    expect_equal(bal.tab(f, data = lalonde, binary = "std", continuous = "std", s.d.denom = "treated")$Balance$Diff.Un,
                 unname(col_w_smd(lalonde[cov_names], treat = lalonde$treat, s.d.denom = "treated")))
    
    expect_equal(bal.tab(f, data = lalonde, binary = "std", continuous = "std", s.d.denom = "treated"),
                 bal.tab(f, data = lalonde, binary = "std", continuous = "std", s.d.denom = "1"))
    
    expect_message(bal.tab(f, data = lalonde, binary = "std", continuous = "std",
                           weights = weights$wATE),
                   .w('Note: `s.d.denom` not specified; assuming "pooled".'),
                   perl = TRUE)
    
    expect_no_message(bal.tab(f, data = lalonde, binary = "std", continuous = "std",
                              weights = weights$wATT))
    
    expect_message(bal.tab(f, data = lalonde, binary = "std", continuous = "std",
                           weights = weights),
                   .w('Note: `s.d.denom` not specified; assuming "pooled" for `wATE` and "treated" for `wATT`.'),
                   perl = TRUE)
    
    #bal.tab() used the first s.d.denom for the unadjusted differences
    expect_equal(suppressMessages(bal.tab(f, data = lalonde, binary = "std", continuous = "std", 
                                          weights = weights, un = TRUE)$Balance$Diff.Un),
                 unname(col_w_smd(lalonde[cov_names], treat = lalonde$treat, s.d.denom = "pooled")))
    
    expect_equal(suppressMessages(bal.tab(f, data = lalonde, binary = "std", continuous = "std", 
                                          weights = weights[2:1], un = TRUE)$Balance$Diff.Un),
                 unname(col_w_smd(lalonde[cov_names], treat = lalonde$treat, s.d.denom = "treated")))
    
    # Continuous treatment
    cov_names <- c("age", "educ", "race", "married", "nodegree", "re74")
    
    f <- reformulate(cov_names, "re75")
    
    expect_message(bal.tab(f, data = lalonde, s.d.denom = "weighted"),
                   .w('Note: `s.d.denom` specified as "weighted", but no weights supplied; setting to "all".'),
                   perl = TRUE)
    
    expect_equal(bal.tab(f, data = lalonde)$Balance$Corr.Un,
                 unname(col_w_corr(lalonde[cov_names], treat = lalonde$re75)))
    
    expect_equal(bal.tab(f, data = lalonde, weights = weights$wATE)$Balance$Corr.Adj,
                 unname(col_w_corr(lalonde[cov_names], treat = lalonde$re75,
                                   weights = weights$wATE)))
    
    # Multicategory treatment
    cov_names <- c("age", "educ", "married", "nodegree", "re74", "re75")
    
    f <- reformulate(cov_names, "race")
    
    expect_message(bal.tab(f, data = lalonde, binary = "std", continuous = "std"),
                   .w('Note: `s.d.denom` not specified; assuming "pooled".'),
                   perl = TRUE)
    
    #No message when no variables are to be standardized
    expect_no_message(bal.tab(f, data = lalonde, binary = "raw", continuous = "raw"))
    
    #No message when SMDs not requested
    expect_no_message(bal.tab(f, data = lalonde, binary = "raw", continuous = "raw",
                              stats = "ks"))
    
    expect_message(bal.tab(f, data = lalonde, binary = "std", continuous = "std",
                           s.d.denom = "weighted"),
                   .w('Note: `s.d.denom` specified as "weighted", but no weights supplied; setting to "all".'),
                   perl = TRUE)
    
    expect_equal(bal.tab(f, data = lalonde, binary = "std", continuous = "std", s.d.denom = "pooled")$Pair.Balance[["hispan vs. black"]]$Balance$Diff.Un,
                 unname(col_w_smd(lalonde[cov_names], treat = lalonde$race, s.d.denom = "pooled",
                                  subset = lalonde$race %in% c("hispan", "black"))))
    expect_equal(suppressMessages(bal.tab(f, data = lalonde, binary = "std", continuous = "std")),
                 bal.tab(f, data = lalonde, binary = "std", continuous = "std", s.d.denom = "pooled"))
    
    # Clustered data
    
    cov_names <- c("age", "educ", "race", "married", "nodegree", "re74", "re75")
    
    f <- reformulate(cov_names, "treat")
    
    expect_message(bal.tab(f, data = lalonde, binary = "std", continuous = "std",
                           cluster = "race"),
                   .w('Note: `s.d.denom` not specified; assuming "pooled".'),
                   perl = TRUE)
    
    expect_no_message(bal.tab(f, data = lalonde, binary = "raw", continuous = "raw",
                              cluster = "race"))
    
    expect_equal(bal.tab(f, data = lalonde, binary = "std", continuous = "std", s.d.denom = "pooled",
                         cluster = "race")$Cluster.Balance$black$Balance$Diff.Un,
                 unname(col_w_smd(lalonde[lalonde$race == "black", setdiff(cov_names, "race")],
                                  treat = lalonde$treat[lalonde$race == "black"], s.d.denom = "pooled")))
    
})
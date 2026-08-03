skip_on_cran()
test_that("bal.tab() works with data.frames", {
    data("lalonde")
    cov_names <- c("age", "educ", "race", "married", "nodegree", "re74", 
                   "re75")
    covs <- lalonde[cov_names]
    w <- runif(nrow(lalonde))
    sw <- runif(nrow(lalonde))
    dist <- runif(nrow(lalonde))
    
    expect_s3_class(bal.tab(covs, treat = lalonde$treat, s.d.denom = "pooled"), "bal.tab")
    expect_s3_class(bal.tab(covs, treat = lalonde$treat, s.d.denom = "pooled", distance = dist,
                            weights = w, s.weights = sw), "bal.tab")
    expect_s3_class(bal.tab(covs, treat = lalonde$treat, s.d.denom = "pooled", distance = dist,
                            weights = w, s.weights = sw, cluster = lalonde$race), "bal.tab.cluster")
    expect_s3_class(bal.tab(covs[-3], treat = lalonde$race, s.d.denom = "pooled", distance = dist,
                            weights = w, s.weights = sw), "bal.tab.multi")
})

# ---------------------------------------------------------------------------
# The block above checks that each input shape dispatches to the right class. The
# blocks below check that the numbers inside are right, which is what a refactor
# of the internals can break without changing any class.
#
# Every assertion ties a reported value to an independently computed one -- either
# a `col_w_*()` function (the documented definition of the statistic) or a
# `bal.tab()` call by a different route that must agree. None of them hardcodes a
# number, so they survive changes to the fixture but not to the arithmetic.

#Shared fixture. Seeded, so `w` and `sw` are the same on every run and in every
#block; `dist` is deliberately unrelated to `treat` so it behaves like any other
#continuous covariate.
btf <- function() {
  data("lalonde", package = "cobalt")

  cov_names <- c("age", "educ", "race", "married", "nodegree", "re74", "re75")

  set.seed(4321)

  list(covs = lalonde[cov_names],
       cov_names = cov_names,
       treat = lalonde$treat,
       race = lalonde$race,
       w = runif(nrow(lalonde)),
       sw = runif(nrow(lalonde)),
       dist = runif(nrow(lalonde)))
}

test_that("bal.tab() reports the statistics its column names claim", {
  f <- btf()

  #`splitfactor()` gives the same one-column-per-level layout `bal.tab()` uses
  #internally, so the `col_w_*()` results line up row for row with `$Balance`.
  sp <- splitfactor(f$covs, drop.first = FALSE)

  b <- bal.tab(f$covs, treat = f$treat, s.d.denom = "pooled",
               weights = f$w, s.weights = f$sw, un = TRUE,
               disp = c("m", "sd"), stats = c("m", "v", "ks"))
  B <- b$Balance

  expect_identical(rownames(B), names(sp))
  expect_identical(B$Type, ifelse(names(sp) %in% c("age", "educ", "re74", "re75"),
                                  "Contin.", "Binary"))

  cont <- B$Type == "Contin."

  #Sampling weights alone define the "unadjusted" sample; the balancing weights
  #multiply them for the "adjusted" one. Getting this product wrong is the classic
  #weighting refactor bug, and it would show up in all four means below.
  expect_equal(B$M.0.Un, unname(col_w_mean(sp, subset = f$treat == 0, weights = f$sw)))
  expect_equal(B$M.1.Un, unname(col_w_mean(sp, subset = f$treat == 1, weights = f$sw)))
  expect_equal(B$M.0.Adj, unname(col_w_mean(sp, subset = f$treat == 0, weights = f$w * f$sw)))
  expect_equal(B$M.1.Adj, unname(col_w_mean(sp, subset = f$treat == 1, weights = f$w * f$sw)))

  #SDs and variance ratios are reported only for continuous variables; the binary
  #rows are blank by design rather than by accident.
  expect_equal(B$SD.0.Un[cont], unname(col_w_sd(sp, subset = f$treat == 0, weights = f$sw))[cont])
  expect_equal(B$SD.1.Adj[cont], unname(col_w_sd(sp, subset = f$treat == 1, weights = f$w * f$sw))[cont])
  expect_true(all(is.na(B$SD.0.Un[!cont])))
  expect_true(all(is.na(B$SD.1.Adj[!cont])))

  expect_equal(B$V.Ratio.Adj[cont],
               unname(col_w_vr(sp, treat = f$treat, weights = f$w, s.weights = f$sw))[cont])
  expect_true(all(is.na(B$V.Ratio.Adj[!cont])))

  #`continuous = "std"` and `binary = "raw"` are the defaults, so exactly the
  #continuous rows are standardized. `std` is the switch that decides it.
  expect_equal(B$Diff.Un,
               unname(col_w_smd(sp, treat = f$treat, s.weights = f$sw,
                                s.d.denom = "pooled", std = cont)))
  expect_equal(B$Diff.Adj,
               unname(col_w_smd(sp, treat = f$treat, weights = f$w, s.weights = f$sw,
                                s.d.denom = "pooled", std = cont)))

  expect_equal(B$KS.Un, unname(col_w_ks(sp, treat = f$treat, s.weights = f$sw)))
  expect_equal(B$KS.Adj,
               unname(col_w_ks(sp, treat = f$treat, weights = f$w, s.weights = f$sw)))
})

test_that("bal.tab() reports sample sizes as effective sample sizes", {
  f <- btf()

  ess <- function(x) sum(x)^2 / sum(x^2)

  #Without sampling weights the unadjusted row is a plain count.
  b <- bal.tab(f$covs, treat = f$treat, s.d.denom = "pooled", weights = f$w, un = TRUE)

  expect_identical(rownames(b$Observations), c("Unadjusted", "Adjusted"))
  expect_identical(colnames(b$Observations), c("Control", "Treated"))
  expect_identical(attr(b$Observations, "ss.type"), c("ss", "ess"))

  expect_equal(b$Observations[["Control"]],
               c(sum(f$treat == 0), ess(f$w[f$treat == 0])))
  expect_equal(b$Observations[["Treated"]],
               c(sum(f$treat == 1), ess(f$w[f$treat == 1])))

  #With sampling weights *both* rows are effective sample sizes: the unadjusted
  #one for `s.weights`, the adjusted one for the product.
  b2 <- bal.tab(f$covs, treat = f$treat, s.d.denom = "pooled",
                weights = f$w, s.weights = f$sw, un = TRUE)

  expect_equal(b2$Observations[["Control"]],
               c(ess(f$sw[f$treat == 0]), ess((f$w * f$sw)[f$treat == 0])))
  expect_equal(b2$Observations[["Treated"]],
               c(ess(f$sw[f$treat == 1]), ess((f$w * f$sw)[f$treat == 1])))

  #Constant weights leave the sample size alone.
  b3 <- bal.tab(f$covs, treat = f$treat, s.d.denom = "pooled",
                weights = rep(3, length(f$treat)), un = TRUE)

  expect_equal(b3$Observations[["Control"]], rep(sum(f$treat == 0), 2L))
  expect_equal(b3$Observations[["Treated"]], rep(sum(f$treat == 1), 2L))
})

test_that("bal.tab() gives the same answer through every input interface", {
  f <- btf()

  args <- list(treat = f$treat, s.d.denom = "pooled", weights = f$w,
               s.weights = f$sw, un = TRUE, stats = c("m", "v", "ks"),
               disp = c("m", "sd"))

  b_df <- do.call(bal.tab, c(list(f$covs), args))

  #A formula plus `data` must reach the identical table.
  b_fo <- bal.tab(reformulate(f$cov_names, "treat"), data = lalonde,
                  s.d.denom = "pooled", weights = f$w, s.weights = f$sw,
                  un = TRUE, stats = c("m", "v", "ks"), disp = c("m", "sd"))
  expect_equal(b_fo$Balance, b_df$Balance)
  expect_equal(b_fo$Observations, b_df$Observations)

  #So must a pre-split numeric matrix, whose columns are already the dummies
  #`bal.tab()` would have created.
  b_mx <- do.call(bal.tab, c(list(as.matrix(splitfactor(f$covs, drop.first = FALSE))), args))
  expect_identical(rownames(b_mx$Balance), rownames(b_df$Balance))
  expect_equal(b_mx$Balance$Diff.Adj, b_df$Balance$Diff.Adj)
  expect_equal(b_mx$Balance$KS.Adj, b_df$Balance$KS.Adj)

  #Naming the treatment column instead of passing the vector must not matter.
  b_chr <- bal.tab(f$covs, treat = "treat", data = lalonde, s.d.denom = "pooled",
                   weights = f$w, s.weights = f$sw, un = TRUE,
                   stats = c("m", "v", "ks"), disp = c("m", "sd"))
  expect_equal(b_chr$Balance, b_df$Balance)
})

test_that("bal.tab() results do not depend on row order", {
  f <- btf()

  #Every statistic here is a function of the multiset of (covariate, treat, weight)
  #triples, so shuffling the rows must leave the table bit-for-bit identical. This
  #catches any refactor that starts relying on the data arriving pre-sorted -- a
  #real hazard with `lalonde`, which is sorted treated-first.
  set.seed(77)
  p <- sample(length(f$treat))

  args <- list(s.d.denom = "pooled", un = TRUE, stats = c("m", "v", "ks"),
               disp = c("m", "sd"))

  b <- do.call(bal.tab, c(list(f$covs, treat = f$treat, weights = f$w,
                               s.weights = f$sw), args))
  b_p <- do.call(bal.tab, c(list(f$covs[p, ], treat = f$treat[p], weights = f$w[p],
                                 s.weights = f$sw[p]), args))

  expect_equal(b_p$Balance, b$Balance)
  expect_equal(b_p$Observations, b$Observations)
})

test_that("bal.tab() options that only affect presentation do not change values", {
  f <- btf()

  args <- list(f$covs, treat = f$treat, s.d.denom = "pooled", weights = f$w,
               un = TRUE, stats = c("m", "v", "ks"))

  #`quick = FALSE` computes more columns; the ones it shares with `quick = TRUE`
  #must be unchanged.
  b_q <- do.call(bal.tab, c(args, list(quick = TRUE)))
  b_f <- do.call(bal.tab, c(args, list(quick = FALSE)))

  shared <- intersect(colnames(b_q$Balance), colnames(b_f$Balance))
  expect_true(all(c("Diff.Adj", "KS.Adj", "V.Ratio.Adj") %in% shared))
  expect_equal(b_f$Balance[shared], b_q$Balance[shared])

  #`abs = TRUE` is a display transform of the same numbers.
  b_abs <- do.call(bal.tab, c(args, list(abs = TRUE)))
  expect_equal(b_abs$Balance$Diff.Adj, abs(b_q$Balance$Diff.Adj))
  expect_equal(b_abs$Balance$KS.Adj, abs(b_q$Balance$KS.Adj))

  #Thresholds add a verdict column without touching the statistic.
  b_thr <- do.call(bal.tab, c(args, list(thresholds = c(m = .1))))
  expect_equal(b_thr$Balance$Diff.Adj, b_q$Balance$Diff.Adj)
  expect_identical(b_thr$Balance$M.Threshold,
                   ifelse(abs(b_q$Balance$Diff.Adj) < .1, "Balanced, <0.1", "Not Balanced, >0.1"))
})

test_that("bal.tab() puts `distance` first and treats it as a covariate", {
  f <- btf()

  b <- bal.tab(f$covs, treat = f$treat, s.d.denom = "pooled", distance = f$dist,
               weights = f$w, un = TRUE)

  expect_identical(rownames(b$Balance)[1L], "distance")
  expect_identical(b$Balance[["Type"]][1L], "Distance")

  #A distance is always standardized, whatever `continuous` says, so its row must
  #match the standardized difference of the supplied values.
  expect_equal(b$Balance[["Diff.Adj"]][1L],
               unname(col_w_smd(data.frame(distance = f$dist), treat = f$treat,
                                weights = f$w, s.d.denom = "pooled", std = TRUE)))

  #Naming the distance changes the row label and nothing else.
  b2 <- bal.tab(f$covs, treat = f$treat, s.d.denom = "pooled",
                distance = data.frame(prop.score = f$dist), weights = f$w, un = TRUE)
  expect_identical(rownames(b2$Balance)[1L], "prop.score")
  expect_equal(unname(unlist(b2$Balance[-1L])), unname(unlist(b$Balance[-1L])))
})

test_that("bal.tab() with `cluster` computes each cluster on its own subset", {
  f <- btf()

  covs <- f$covs[setdiff(f$cov_names, "race")]

  b <- bal.tab(covs, treat = f$treat, s.d.denom = "pooled", weights = f$w,
               cluster = f$race, un = TRUE, cluster.summary = TRUE,
               stats = c("m", "ks"))

  expect_identical(names(b$Cluster.Balance), levels(f$race))

  #A cluster's table must equal a standalone `bal.tab()` on that cluster's rows --
  #including the standard deviations, which are computed within the cluster.
  for (g in levels(f$race)) {
    in_g <- f$race == g

    b_g <- bal.tab(covs[in_g, ], treat = f$treat[in_g], s.d.denom = "pooled",
                   weights = f$w[in_g], un = TRUE, stats = c("m", "ks"))

    expect_equal(b$Cluster.Balance[[g]]$Balance, b_g$Balance, info = g)

    #The per-cluster tables label their sample-size columns with the treatment
    #values rather than "Control"/"Treated", so compare the numbers only.
    expect_equal(unname(unlist(b$Cluster.Balance[[g]]$Observations)),
                 unname(unlist(b_g$Observations)), info = g)
  }

  #The across-cluster summary is the row-wise min/mean/max over the clusters, on
  #the *signed* statistics -- so `Max` can be smaller than `Min` in absolute value.
  per_cluster <- do.call(cbind, lapply(b$Cluster.Balance,
                                       function(z) z$Balance$Diff.Adj))

  expect_equal(b$Balance.Across.Clusters[["Min.Diff.Adj"]],
               apply(per_cluster, 1L, min), ignore_attr = TRUE)
  expect_equal(b$Balance.Across.Clusters[["Mean.Diff.Adj"]],
               rowMeans(per_cluster), ignore_attr = TRUE)
  expect_equal(b$Balance.Across.Clusters[["Max.Diff.Adj"]],
               apply(per_cluster, 1L, max), ignore_attr = TRUE)

  #`abs = TRUE` is what makes the summary aggregate magnitudes instead.
  b_abs <- bal.tab(covs, treat = f$treat, s.d.denom = "pooled", weights = f$w,
                   cluster = f$race, un = TRUE, cluster.summary = TRUE,
                   stats = c("m", "ks"), abs = TRUE)

  expect_equal(b_abs$Balance.Across.Clusters[["Max.Diff.Adj"]],
               apply(abs(per_cluster), 1L, max), ignore_attr = TRUE)
})

test_that("bal.tab() with a multi-category treatment pools SDs across all groups", {
  f <- btf()

  covs <- f$covs[setdiff(f$cov_names, "race")]

  b <- bal.tab(covs, treat = f$race, s.d.denom = "pooled", weights = f$w,
               un = TRUE, binary = "std", continuous = "std", stats = c("m", "ks"))

  expect_identical(names(b$Pair.Balance),
                   c("hispan vs. black", "white vs. black", "white vs. hispan"))

  #KS statistics need no standard deviation, so a pair's column is exactly what a
  #standalone two-group `bal.tab()` on those rows reports.
  in_pair <- f$race %in% c("hispan", "black")
  b_pair <- bal.tab(covs[in_pair, ], treat = factor(f$race[in_pair]),
                    s.d.denom = "pooled", weights = f$w[in_pair], un = TRUE,
                    binary = "std", continuous = "std", stats = c("m", "ks"))

  expect_equal(b$Pair.Balance[["hispan vs. black"]]$Balance$KS.Adj,
               b_pair$Balance$KS.Adj)

  #The SMDs, by contrast, use a denominator pooled over *all three* groups, not
  #just the two being compared. `col_w_smd()` given the full treatment does the
  #same, so it agrees; the two-group table above deliberately does not.
  expect_equal(b$Pair.Balance[["hispan vs. black"]]$Balance$Diff.Adj,
               unname(col_w_smd(covs, treat = f$race, weights = f$w,
                                s.d.denom = "pooled", subset = in_pair, std = TRUE)))
  expect_false(isTRUE(all.equal(b$Pair.Balance[["hispan vs. black"]]$Balance$Diff.Adj,
                                b_pair$Balance$Diff.Adj)))

  #Each group's sample size is its own, and the summary is the pairwise maximum.
  expect_equal(unlist(b$Observations["Unadjusted", ]),
               vapply(levels(f$race), function(g) sum(f$race == g), numeric(1L)))
  expect_equal(b$Balance.Across.Pairs[["Max.Diff.Adj"]],
               do.call(pmax, lapply(b$Pair.Balance,
                                    function(z) abs(z$Balance$Diff.Adj))),
               ignore_attr = TRUE)
})

test_that("bal.tab() handles several weight sets as independent columns", {
  f <- btf()

  #Each named weight set must be summarized exactly as if it had been passed
  #alone, so a refactor cannot cross-contaminate the columns.
  ws <- list(a = f$w, b = f$sw)

  b <- bal.tab(f$covs, treat = f$treat, s.d.denom = "pooled", weights = ws,
               un = TRUE, stats = c("m", "ks"))

  expect_identical(attr(b, "print.options")$weight.names, c("a", "b"))
  expect_identical(colnames(b$Observations), c("Control", "Treated"))

  #With more than one weight set the shared unadjusted row is labelled "All".
  expect_identical(rownames(b$Observations), c("All", "a", "b"))
  expect_identical(unname(attr(b$Observations, "ss.type")), c("ss", "ess", "ess"))

  for (nm in names(ws)) {
    b_one <- bal.tab(f$covs, treat = f$treat, s.d.denom = "pooled",
                     weights = ws[[nm]], un = TRUE, stats = c("m", "ks"))

    expect_equal(b$Balance[[paste0("Diff.", nm)]], b_one$Balance$Diff.Adj, info = nm)
    expect_equal(b$Balance[[paste0("KS.", nm)]], b_one$Balance$KS.Adj, info = nm)
  }

  #The unadjusted column is shared and must match the no-weights table.
  b_un <- bal.tab(f$covs, treat = f$treat, s.d.denom = "pooled", stats = c("m", "ks"))
  expect_equal(b$Balance$Diff.Un, b_un$Balance$Diff.Un)
})

test_that("bal.tab() adds `addl` and interaction terms without altering the originals", {
  f <- btf()

  covs <- f$covs[c("age", "educ")]

  b <- bal.tab(covs, treat = f$treat, s.d.denom = "pooled", weights = f$w, un = TRUE)

  #`addl` appends rows; the rows already there keep their values.
  b_addl <- bal.tab(covs, treat = f$treat, s.d.denom = "pooled", weights = f$w,
                    un = TRUE, addl = f$covs["race"])
  expect_identical(rownames(b_addl$Balance),
                   c("age", "educ", "race_black", "race_hispan", "race_white"))
  expect_equal(b_addl$Balance[1:2, ], b$Balance)

  #So does `int = TRUE`, which adds the cross term. Squares are `poly`'s job, not
  #`int`'s -- asserting both here pins the division of labour between them.
  b_int <- bal.tab(covs, treat = f$treat, s.d.denom = "pooled", weights = f$w,
                   un = TRUE, int = TRUE)
  expect_equal(b_int$Balance[rownames(b$Balance), ], b$Balance)
  expect_identical(rownames(b_int$Balance), c("age", "educ", "age * educ"))

  b_poly <- bal.tab(covs, treat = f$treat, s.d.denom = "pooled", weights = f$w,
                    un = TRUE, poly = 2)
  expect_identical(rownames(b_poly$Balance), c("age", "educ", "age²", "educ²"))

  b_both <- bal.tab(covs, treat = f$treat, s.d.denom = "pooled", weights = f$w,
                    un = TRUE, int = TRUE, poly = 2)
  expect_identical(rownames(b_both$Balance),
                   c("age", "educ", "age²", "educ²", "age * educ"))

  #The square row is the balance of the squared covariate.
  expect_equal(b_poly$Balance["age²", "Diff.Adj"],
               unname(col_w_smd(data.frame(x = covs$age^2), treat = f$treat,
                                weights = f$w, s.d.denom = "pooled", std = TRUE)))

  #The interaction row is the balance of the product, computed like any covariate.
  expect_equal(b_int$Balance["age * educ", "Diff.Adj"],
               unname(col_w_smd(data.frame(x = covs$age * covs$educ), treat = f$treat,
                                weights = f$w, s.d.denom = "pooled", std = TRUE)))
})
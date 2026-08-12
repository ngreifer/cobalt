#`bal.tab()` for the `formula`, `data.frame`, and `matrix` methods. These are the
#entry points that need no package from Suggests, and they drive most of
#`x2base.data.frame()` and `get_covs_from_formula()`/`get_treat_from_formula()`.

test_that("the formula, data.frame, and matrix methods agree", {
  covs <- fx("covs")

  b_df <- bal.tab(covs, treat = lalonde$treat, s.d.denom = "pooled")
  b_f <- bal.tab(reformulate(cov_names, "treat"), data = lalonde,
                 s.d.denom = "pooled")

  #The stored call differs; everything else must match.
  expect_equal(b_df$Balance, b_f$Balance)
  expect_equal(b_df$Observations, b_f$Observations)

  #A matrix of numeric covariates dispatches to the same method.
  m <- as.matrix(lalonde[c("age", "educ", "re74")])
  b_m <- bal.tab(m, treat = lalonde$treat, s.d.denom = "pooled")
  expect_s3_class(b_m, "bal.tab.bin")
  expect_identical(rownames(b_m$Balance), c("age", "educ", "re74"))

})

test_that("treatment type determines the bal.tab subclass", {
  covs <- fx("covs")

  expect_s3_class(bal.tab(covs, treat = lalonde$treat, s.d.denom = "pooled"), "bal.tab.bin")
  expect_s3_class(bal.tab(covs, treat = lalonde$race, s.d.denom = "pooled"),
                  "bal.tab.multi")
  expect_s3_class(bal.tab(covs, treat = lalonde$re75), "bal.tab.cont")

  #A two-level factor is still binary.
  expect_s3_class(bal.tab(covs, treat = factor(lalonde$treat), s.d.denom = "pooled"),
                  "bal.tab.bin")
})

test_that("`treat` can be named in `data` or supplied as a vector", {
  covs <- fx("covs")

  b_vec <- bal.tab(covs, treat = lalonde$treat, s.d.denom = "pooled")
  b_nm <- bal.tab(covs, treat = "treat", data = lalonde, s.d.denom = "pooled")

  expect_equal(b_vec$Balance, b_nm$Balance)
})

test_that("formulas accept expressions, interactions, and data frames on the RHS", {
  #A function of a variable.
  b <- bal.tab(treat ~ age + log(re74 + 1), data = lalonde, s.d.denom = "pooled")
  expect_identical(rownames(b$Balance), c("age", "log(re74 + 1)"))

  #An explicit interaction. `married` is 0/1, so cobalt treats it as a binary
  #variable and produces one interaction row per level.
  b <- bal.tab(treat ~ age * married, data = lalonde, s.d.denom = "pooled")
  expect_true(all(c("age", "married") %in% rownames(b$Balance)))
  expect_true(any(startsWith(rownames(b$Balance), "age * married")))

  #`:` gives the interaction without the main effects.
  b <- bal.tab(treat ~ age:married, data = lalonde, s.d.denom = "pooled")
  expect_false(any(rownames(b$Balance) == "age"))
  expect_true(all(startsWith(rownames(b$Balance), "age * married")))

  #`I()`.
  b <- bal.tab(treat ~ I(age^2), data = lalonde, s.d.denom = "pooled")
  expect_identical(rownames(b$Balance), "I(age^2)")

  #A data frame as a single RHS term is expanded into its columns.
  covs <- fx("covs")
  b <- bal.tab(treat ~ covs, data = lalonde, s.d.denom = "pooled")
  expect_true(all(c("age", "educ", "married") %in% rownames(b$Balance)))

  #`.` expands to every other column of `data`.
  b <- bal.tab(treat ~ ., data = lalonde[c("treat", "age", "educ")],
               s.d.denom = "pooled")
  expect_identical(rownames(b$Balance), c("age", "educ"))
})

test_that("factors are split into dummies and named with the separator", {
  b <- bal.tab(treat ~ race, data = lalonde, s.d.denom = "pooled")
  expect_identical(rownames(b$Balance), c("race_black", "race_hispan", "race_white"))

  #A two-level factor keeps only its second level.
  d <- transform(lalonde, f2 = factor(ifelse(married == 1, "yes", "no")))
  b <- bal.tab(treat ~ f2, data = d, s.d.denom = "pooled")
  expect_identical(rownames(b$Balance), "f2_yes")
})

test_that("`int` and `poly` add interaction and polynomial rows", {
  covs <- lalonde[c("age", "educ")]

  b1 <- bal.tab(covs, treat = lalonde$treat, s.d.denom = "pooled")
  b_int <- bal.tab(covs, treat = lalonde$treat, int = TRUE, s.d.denom = "pooled")
  b_p2 <- bal.tab(covs, treat = lalonde$treat, poly = 2, s.d.denom = "pooled")
  b_p3 <- bal.tab(covs, treat = lalonde$treat, poly = 3, s.d.denom = "pooled")

  expect_gt(nrow(b_int$Balance), nrow(b1$Balance))
  expect_gt(nrow(b_p2$Balance), nrow(b1$Balance))
  expect_gt(nrow(b_p3$Balance), nrow(b_p2$Balance))

  #`int = TRUE` implies the squared terms as well as the cross-product.
  expect_true("age * educ" %in% rownames(b_int$Balance))

  #A single covariate has no interactions to form, and `poly = 1` adds no
  #polynomials; neither may error (both once returned NULL internally).
  expect_s3_class(bal.tab(lalonde["age"], treat = lalonde$treat, int = TRUE,
                          s.d.denom = "pooled"), "bal.tab")
  expect_s3_class(bal.tab(lalonde["married"], treat = lalonde$treat, int = TRUE,
                          s.d.denom = "pooled"), "bal.tab")

  #Binary covariates are excluded from polynomials, so an all-binary set adds none.
  expect_s3_class(bal.tab(lalonde[c("married", "nodegree")], treat = lalonde$treat,
                          poly = 2, s.d.denom = "pooled"), "bal.tab")
})

test_that("`weights` is accepted as a vector, list, data frame, or variable name", {
  covs <- fx("covs")
  d <- transform(lalonde, w1 = w_fixed, w2 = sw_fixed)

  b_vec <- bal.tab(covs, treat = lalonde$treat, weights = w_fixed, s.d.denom = "pooled")
  b_nm <- bal.tab(covs, treat = lalonde$treat, weights = "w1", data = d,
                  s.d.denom = "pooled")
  expect_equal(b_vec$Balance, b_nm$Balance)

  #Multiple weight sets produce one adjusted column per set.
  b_list <- bal.tab(covs, treat = lalonde$treat,
                    weights = list(a = w_fixed, b = sw_fixed), s.d.denom = "pooled")
  expect_true(all(c("Diff.a", "Diff.b") %in% names(b_list$Balance)))

  b_dfw <- bal.tab(covs, treat = lalonde$treat,
                   weights = data.frame(a = w_fixed, b = sw_fixed), s.d.denom = "pooled")
  expect_equal(b_list$Balance, b_dfw$Balance)

  #Weight sets named in `data` take their column names.
  b_2nm <- bal.tab(covs, treat = lalonde$treat, weights = c("w1", "w2"), data = d,
                   s.d.denom = "pooled")
  expect_true(all(c("Diff.w1", "Diff.w2") %in% names(b_2nm$Balance)))
})

test_that("`addl` and `distance` add rows without entering the main covariates", {
  covs <- lalonde[c("age", "educ")]

  b <- bal.tab(covs, treat = lalonde$treat, addl = lalonde["re74"],
               s.d.denom = "pooled")
  expect_true("re74" %in% rownames(b$Balance))

  b <- bal.tab(covs, treat = lalonde$treat, addl = "race", data = lalonde,
               s.d.denom = "pooled")
  expect_true("race_black" %in% rownames(b$Balance))

  b <- bal.tab(covs, treat = lalonde$treat, addl = ~ re74 + race, data = lalonde,
               s.d.denom = "pooled")
  expect_true(all(c("re74", "race_black") %in% rownames(b$Balance)))

  #`distance` rows sort to the top and are typed "Distance".
  ps <- fitted(glm(treat ~ age + educ + re74, data = lalonde, family = binomial))
  b <- bal.tab(covs, treat = lalonde$treat, distance = ps, s.d.denom = "pooled")
  expect_identical(rownames(b$Balance)[1L], "distance")
  expect_identical(as.character(b$Balance$Type[1L]), "Distance")

  b <- bal.tab(covs, treat = lalonde$treat, distance = data.frame(myps = ps),
               s.d.denom = "pooled")
  expect_identical(rownames(b$Balance)[1L], "myps")

  #A distance measure taking only two values is split into one dummy per value.
  b <- bal.tab(covs, treat = lalonde$treat, distance = w_fixed, s.d.denom = "pooled")
  expect_identical(rownames(b$Balance)[1:2], c("distance_1", "distance_2"))
})

test_that("`s.weights` are applied to both the unadjusted and adjusted samples", {
  covs <- lalonde[c("age", "educ")]

  b_no <- bal.tab(covs, treat = lalonde$treat, un = TRUE, s.d.denom = "pooled")
  b_sw <- bal.tab(covs, treat = lalonde$treat, s.weights = sw_fixed, un = TRUE,
                  s.d.denom = "pooled")

  #Sampling weights change the unadjusted differences.
  expect_false(isTRUE(all.equal(b_no$Balance$Diff.Un, b_sw$Balance$Diff.Un)))

  #Named in `data` gives the same answer as passed as a vector.
  d <- transform(lalonde, sw = sw_fixed)
  b_nm <- bal.tab(covs, treat = lalonde$treat, s.weights = "sw", data = d,
                  un = TRUE, s.d.denom = "pooled")
  expect_equal(b_sw$Balance, b_nm$Balance)
})

test_that("`subset` selects units logically and numerically", {
  covs <- lalonde[c("age", "educ")]
  keep <- lalonde$age < 30

  b_log <- bal.tab(covs, treat = lalonde$treat, subset = keep, s.d.denom = "pooled")
  b_num <- bal.tab(covs, treat = lalonde$treat, subset = which(keep), s.d.denom = "pooled")
  b_neg <- bal.tab(covs, treat = lalonde$treat, subset = -which(!keep),
                   s.d.denom = "pooled")

  expect_equal(b_log$Balance, b_num$Balance)
  expect_equal(b_log$Balance, b_neg$Balance)
  expect_equal(sum(b_log$Observations), sum(keep))
})

test_that("`abs`, `quick`, and `un` control what is computed and reported", {
  covs <- lalonde[c("age", "educ", "married")]

  b <- bal.tab(covs, treat = lalonde$treat, abs = TRUE, un = TRUE,
               s.d.denom = "pooled")
  expect_true(all(b$Balance$Diff.Un >= 0))

  #`un` is a display option: the unadjusted column stays in `$Balance` either way,
  #but it is only printed when `un = TRUE`.
  b_un <- bal.tab(covs, treat = lalonde$treat, weights = w_fixed, un = TRUE,
                  s.d.denom = "pooled")
  b_no <- bal.tab(covs, treat = lalonde$treat, weights = w_fixed, un = FALSE,
                  s.d.denom = "pooled")
  expect_true("Diff.Un" %in% names(b_no$Balance))
  expect_true(attr(b_un, "print.options")$un)
  expect_false(attr(b_no, "print.options")$un)
  expect_output(print(b_un), "Diff.Un")
  expect_false(any(grepl("Diff.Un", capture.output(print(b_no)), fixed = TRUE)))

  #`quick = FALSE` computes the columns that are otherwise skipped.
  b_q <- bal.tab(covs, treat = lalonde$treat, quick = TRUE, s.d.denom = "pooled")
  b_nq <- bal.tab(covs, treat = lalonde$treat, quick = FALSE, s.d.denom = "pooled")
  expect_true(ncol(b_nq$Balance) >= ncol(b_q$Balance))
})

test_that("`stats` and `thresholds` add columns and balance tallies", {
  covs <- lalonde[c("age", "educ", "married")]

  b <- bal.tab(covs, treat = lalonde$treat,
               stats = c("mean.diffs", "variance.ratios", "ks.statistics"),
               s.d.denom = "pooled")
  expect_true(all(c("Diff.Un", "V.Ratio.Un", "KS.Un") %in% names(b$Balance)))

  #Threshold columns follow the sample they describe: with no weights only the
  #unadjusted sample exists, so they are suffixed ".Un".
  b <- bal.tab(covs, treat = lalonde$treat, stats = c("m", "v"),
               thresholds = c(m = .1, v = 2), s.d.denom = "pooled")
  expect_true(all(c("M.Threshold.Un", "V.Threshold.Un") %in% names(b$Balance)))

  #With weights they describe the adjusted sample and carry no suffix.
  bw <- bal.tab(covs, treat = lalonde$treat, weights = w_fixed, stats = c("m", "v"),
                thresholds = c(m = .1, v = 2), s.d.denom = "pooled")
  expect_true(all(c("M.Threshold", "V.Threshold") %in% names(bw$Balance)))
  expect_true(is_not_null(bw$Balanced.mean.diffs))
  expect_true(is_not_null(bw$Max.Imbalance.mean.diffs))

  #Continuous treatments use the correlation statistics.
  b <- bal.tab(covs, treat = lalonde$re75,
               stats = c("correlations", "spearman.correlations"))
  expect_true(all(c("Corr.Un", "S.Corr.Un") %in% names(b$Balance)))
})

test_that("`method` selects how weights are interpreted", {
  covs <- lalonde[c("age", "educ")]

  b_w <- bal.tab(covs, treat = lalonde$treat, weights = w_fixed,
                 method = "weighting", s.d.denom = "pooled")
  b_m <- bal.tab(covs, treat = lalonde$treat, weights = w_fixed,
                 method = "matching", s.d.denom = "pooled")

  #The balance statistics are identical; only the reported sample sizes differ.
  expect_equal(b_w$Balance, b_m$Balance)
  expect_false(identical(rownames(b_w$Observations), rownames(b_m$Observations)))

  #`match.strata` is converted to weights.
  b_ms <- bal.tab(covs, treat = lalonde$treat, match.strata = sub_idx,
                  s.d.denom = "pooled")
  expect_s3_class(b_ms, "bal.tab.bin")

  b_sub <- bal.tab(covs, treat = lalonde$treat, subclass = sub_idx,
                   s.d.denom = "pooled")
  expect_s3_class(b_sub, "bal.tab.subclass")
})

test_that("`estimand` and `focal` select the reference group", {
  covs <- lalonde[c("age", "educ")]

  #For binary treatments, `estimand` determines `s.d.denom` without a message
  #(the note is only emitted when neither is supplied).
  denoms <- c(ATT = "treated", ATC = "control", ATE = "pooled")
  for (e in names(denoms)) {
    expect_equal(bal.tab(covs, treat = lalonde$treat, estimand = e)$Balance$Diff.Un,
                 unname(col_w_smd(covs, treat = lalonde$treat, s.d.denom = denoms[[e]])),
                 info = e)
  }

  #Without either, a note reports the assumed denominator.
  expect_msg(bal.tab(covs, treat = lalonde$treat), "assuming \"pooled\"")

  #`focal` does not drop pairs when `pairwise = TRUE`; every pair is still shown.
  b <- bal.tab(covs, treat = lalonde$race, focal = "white")
  expect_s3_class(b, "bal.tab.multi")
  expect_length(b$Pair.Balance, 3L)

  #With `pairwise = FALSE`, only comparisons against the focal group remain.
  b_np <- bal.tab(covs, treat = lalonde$race, focal = "white", pairwise = FALSE)
  expect_setequal(names(b_np$Pair.Balance), c("white vs. black", "white vs. hispan"))

  #`focal` may be given as an index.
  b_i <- bal.tab(covs, treat = lalonde$race, focal = 3L, pairwise = FALSE)
  expect_identical(names(b_i$Pair.Balance), names(b_np$Pair.Balance))
})

test_that("missing values produce NA rows and a warning", {
  expect_wrn(
    b <- bal.tab(treat ~ age + re74, data = lalonde_mis, s.d.denom = "pooled"),
    "Missing values exist in the covariates")

  #Each variable with missingness gains a `:<NA>` indicator row.
  expect_true("re74:<NA>" %in% rownames(b$Balance))
})

test_that("invalid covariates and treatments are rejected", {
  covs <- lalonde[c("age", "educ")]

  expect_err(bal.tab(1:10), "input object must be an appropriate list, data frame, formula")
  expect_err(bal.tab(list(treat = lalonde$treat)), "`covs` data frame must be specified")
  expect_err(bal.tab(covs, treat = rep(1, nrow(lalonde))),
             "treatment must have at least two unique values")
  expect_err(bal.tab(data.frame(a = c(Inf, rnorm(nrow(lalonde) - 1L))),
                     treat = lalonde$treat, s.d.denom = "pooled"),
             "non-finite values, which are not allowed")

  #A character vector is not a supported `x`; it falls through to the default
  #method. Covariates are named via `addl`/`distance` instead, where the
  #character-vector path is honored.
  expect_err(bal.tab(c("age", "educ"), treat = lalonde$treat, data = lalonde),
             "input object must be an appropriate list, data frame, formula")
  expect_err(bal.tab(covs, treat = lalonde$treat, addl = "nope", data = lalonde),
             "cannot be found")
})

test_that("malformed formulas are rejected informatively", {
  expect_err(bal.tab(nope ~ age, data = lalonde),
             "is not a variable in `data` or the global environment")
  expect_err(bal.tab(treat ~ nosuchvar, data = lalonde),
             "cannot be found")
  expect_err(bal.tab(mean ~ age, data = lalonde),
             "Invalid type (function) for variable")
  expect_err(bal.tab(treat ~ mean, data = lalonde),
             "Invalid type (function) for variable")
  expect_err(bal.tab(treat ~ age, data = 1:10),
             "is not a variable in the global environment")

  #A data frame inside an interaction is not supported.
  covs <- fx("covs")
  expect_err(bal.tab(treat ~ covs * age, data = lalonde),
             "with data frames are not allowed in the input formula")
})

test_that("invalid weights, s.weights, and distance are rejected", {
  covs <- lalonde[c("age", "educ")]
  n <- nrow(lalonde)

  expect_err(bal.tab(covs, treat = lalonde$treat, weights = c(NA, w_fixed[-1L])),
             "NAs are not allowed in the weights")
  expect_err(bal.tab(covs, treat = lalonde$treat, weights = c(Inf, w_fixed[-1L])),
             "Infinite weights are not allowed")
  expect_err(bal.tab(covs, treat = lalonde$treat,
                     weights = ifelse(lalonde$treat == 1, 0, 1)),
             "All weights are zero when the treatment is")
  expect_err(bal.tab(covs, treat = lalonde$treat, s.weights = "nosuchvar"),
             "must be a vector of sampling weights")
  expect_err(bal.tab(covs, treat = lalonde$treat, distance = list(1),
                     s.d.denom = "pooled"),
             "must be a formula or variable containing the distance values")
  expect_err(bal.tab(covs, treat = lalonde$treat,
                     distance = c(NA, w_fixed[-1L]), s.d.denom = "pooled"),
             "Missing values are not allowed in the distance measure")
  expect_err(bal.tab(covs, treat = lalonde$treat, addl = list(a = 1),
                     s.d.denom = "pooled"),
             "must be a formula or variable containing the additional covariates")

  #Length mismatches are caught.
  expect_err(bal.tab(covs, treat = lalonde$treat[-1L]),
             "must have the same number of observations as `treat`")
})

test_that("invalid subset, focal, thresholds, and method are rejected", {
  covs <- lalonde[c("age", "educ")]
  n <- nrow(lalonde)

  expect_err(bal.tab(covs, treat = lalonde$treat, subset = rep(FALSE, n),
                     s.d.denom = "pooled"),
             "All `subset` set to FALSE")
  expect_err(bal.tab(covs, treat = lalonde$treat, subset = "a"),
             "`subset` must be a logical vector or a whole numeric vector")
  expect_err(bal.tab(covs, treat = lalonde$treat, subset = 1e6),
             "cannot be larger than the number of units")
  expect_err(bal.tab(covs, treat = lalonde$treat, subset = c(-1, 2)),
             "Positive and negative indices cannot be mixed")

  expect_err(bal.tab(covs, treat = lalonde$race, focal = "zzz"),
             "is not the name of any treatment group")
  expect_err(bal.tab(covs, treat = lalonde$race, focal = 9),
             "but there are only 3 treatment groups")
  expect_err(bal.tab(covs, treat = lalonde$treat, focal = "zzz"),
             "is not the name of any treatment group")

  expect_err(bal.tab(covs, treat = lalonde$treat, thresholds = c(m = "a"),
                     s.d.denom = "pooled"),
             "`thresholds` must be numeric")

  #Only weights may be combined with more than one method.
  expect_err(bal.tab(covs, treat = lalonde$treat,
                     weights = list(a = w_fixed, b = sw_fixed),
                     method = c("weighting", "subclassification")),
             "Subclassification cannot be specified along with other methods")
  expect_err(bal.tab(covs, treat = lalonde$treat,
                     method = c("weighting", "matching"), match.strata = sub_idx),
             "Only weights can be specified with multiple methods")

  #Subclasses cannot be turned into weights for continuous treatments.
  expect_err(bal.tab(covs, treat = lalonde$re75, match.strata = sub_idx),
             "cannot be turned into weights for continuous treatments")
})

test_that("`s.d.denom` length must match the number of weight sets", {
  covs <- lalonde[c("age", "educ")]

  expect_err(bal.tab(covs, treat = lalonde$treat,
                     weights = list(a = w_fixed, b = sw_fixed),
                     s.d.denom = c("all", "pooled", "treated")),
             "must have length 1 or equal to the number of valid sets of weights")
})

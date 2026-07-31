#`print.bal.tab()` across all seven bal.tab shapes.
#
#`print.bal.tab.R` is almost entirely display logic: `print_process.*` resolves which
#sections and columns to show, and `bal.tab_print.*` emits them. So most of these
#tests assert the *decision* (which columns appear) with `expect_output()` rather
#than snapshotting whole tables, which would be unreviewable and would expose the
#suite to platform differences in the printed digits.
#
#Every object here is built from `lalonde` plus a hard-coded weight vector, so
#nothing depends on an iterative optimizer whose result can vary across platforms.

#Capture printed output as a single whitespace-collapsed string.
printed <- function(x, ...) {
  squish(capture.output(print(x, ...)))
}

b_bin <- function() {
  bal.tab(lalonde[c("age", "educ", "married", "race")], treat = lalonde$treat,
          weights = w_fixed, s.d.denom = "pooled", un = TRUE,
          stats = c("mean.diffs", "variance.ratios"),
          thresholds = c(m = .1, v = 2))
}

test_that("the default print shows the balance table and sample sizes", {
  b <- b_bin()
  out <- printed(b)

  expect_match(out, "Balance Measures")
  expect_match(out, "Effective sample sizes")
  expect_match(out, "Diff.Un")
  expect_match(out, "Diff.Adj")
  expect_match(out, "V.Ratio.Adj")

  #Every covariate appears, with factors split.
  for (v in c("age", "educ", "married", "race_black")) {
    expect_match(out, v, fixed = TRUE)
  }
})

test_that("`un` and `disp.bal.tab` control which sections print", {
  b <- b_bin()

  expect_false(grepl("Diff.Un", printed(b, un = FALSE), fixed = TRUE))
  expect_true(grepl("Diff.Un", printed(b, un = TRUE), fixed = TRUE))

  #`disp.bal.tab = FALSE` drops the covariate table but keeps the tallies, the
  #max-imbalance rows (which still name the statistic column), and sample sizes.
  out <- printed(b, disp.bal.tab = FALSE)
  expect_match(out, "Effective sample sizes")
  expect_false(grepl("age Contin.", out, fixed = TRUE))
})

test_that("`disp.call` prints the stored call when the object has one", {
  #The data.frame and formula methods do not record a call, so `disp.call` warns.
  expect_wrn(print(b_bin(), disp.call = TRUE),
             "cannot be set to TRUE if the input object does not have a call")

  #Objects from supported packages do carry one.
  skip_if_not_installed("MatchIt")
  b <- bal.tab(fx("matchit"))
  expect_match(printed(b, disp.call = TRUE), "Call")
  expect_false(grepl("Call", printed(b, disp.call = FALSE), fixed = TRUE))
})

test_that("`stats` and `disp.thresholds` select columns", {
  b <- b_bin()

  #Only the requested statistic is shown.
  out <- printed(b, stats = "mean.diffs")
  expect_match(out, "Diff.Adj")
  expect_false(grepl("V.Ratio", out, fixed = TRUE))

  #Thresholds can be suppressed per statistic.
  out <- printed(b, disp.thresholds = c(m = FALSE))
  expect_false(grepl("M.Threshold", out, fixed = TRUE))
  expect_match(out, "Diff.Adj")
})

test_that("`imbalanced.only` restricts the table to imbalanced covariates", {
  b <- b_bin()

  out <- printed(b, imbalanced.only = TRUE)
  #The tally line is always shown; it names the threshold being applied.
  expect_match(out, "Balance tally")
})

test_that("`disp` adds means and standard deviations when they were computed", {
  b <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat, weights = w_fixed,
               s.d.denom = "pooled", un = TRUE, disp = c("means", "sds"))

  out <- printed(b)
  expect_match(out, "M.0.Un")
  expect_match(out, "SD.0.Un")

  #A statistic that was not computed cannot be displayed.
  b2 <- b_bin()
  expect_err(print(b2, disp = "bogus"), "not allowed in `disp`")
  expect_err(print(b2, stats = "ks.statistics"),
             "cannot contain")
})

test_that("`digits` changes the printed precision", {
  b <- b_bin()

  expect_false(identical(printed(b, digits = 2L), printed(b, digits = 5L)))
})

test_that("logical print arguments are validated", {
  b <- b_bin()

  expect_err(print(b, un = "yes"), "must be a logical value")
  expect_err(print(b, disp.bal.tab = "yes"), "must be a logical value")
  expect_err(print(b, imbalanced.only = "yes"), "must be a logical value")
  expect_err(print(b, disp.call = "yes"), "must be a logical value")
})

test_that("a threshold cannot be requested for a statistic that lacks one", {
  b <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat,
               s.d.denom = "pooled")

  expect_err(print(b, m.threshold = .1), "must be")
})

test_that("continuous treatments print correlations", {
  b <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$re75, weights = w_fixed,
               un = TRUE, stats = c("correlations", "spearman.correlations"))

  out <- printed(b)
  expect_match(out, "Corr.Un")
  expect_match(out, "S.Corr.Adj")
  expect_match(out, "Effective sample sizes")
})

test_that("multi-category treatments print pairwise tables and a summary", {
  b <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$race,
               weights = w_fixed, un = TRUE, multi.summary = TRUE)

  out <- printed(b)
  expect_match(out, "Balance summary across all treatment pairs")

  #`which.treat` selects individual pairs. The `.all`/`.none` shorthands are
  #rewritten by `bal.tab()` itself, not by `print()`, which takes NULL and NA.
  out_all <- printed(b, which.treat = NULL)
  expect_match(out_all, "black")
  expect_match(out_all, "hispan")

  #By default no individual pair is shown, and NA suppresses them explicitly.
  expect_false(grepl("vs.", printed(b), fixed = TRUE))
  expect_false(grepl("vs.", printed(b, which.treat = NA), fixed = TRUE))

  #A single named group keeps only the pairs involving it.
  expect_match(printed(b, which.treat = "black"), "black")
})

test_that("imputations print per-imputation and across-imputation tables", {
  b <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat, imp = imp_idx,
               weights = w_fixed, s.d.denom = "pooled", un = TRUE,
               imp.summary = TRUE)

  out <- printed(b)
  expect_match(out, "Balance summary across all imputations")
  expect_match(out, "Average effective sample sizes across imputations")

  out <- printed(b, which.imp = 1L)
  expect_match(out, "Imputation 1")

  #`imp.fun` picks a single summary function.
  expect_match(printed(b, imp.fun = "mean"), "Mean")
  expect_err(print(b, imp.fun = "bogus"), "`imp.fun` must be")
})

test_that("clusters print per-cluster and across-cluster tables", {
  b <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat, cluster = cl_idx,
               weights = w_fixed, s.d.denom = "pooled", un = TRUE,
               cluster.summary = TRUE)

  out <- printed(b)
  expect_match(out, "Balance by cluster")
  expect_match(out, "Balance summary across all clusters")
  expect_match(out, "Total effective sample sizes across clusters")

  out <- printed(b, which.cluster = "a")
  expect_match(out, "Cluster")

  expect_match(printed(b, cluster.fun = "mean"), "Mean")
  expect_err(print(b, cluster.fun = "bogus"), "`cluster.fun` must be")
})

test_that("subclasses print per-subclass and across-subclass tables", {
  b <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat, subclass = sub_idx,
               s.d.denom = "pooled", un = TRUE, subclass.summary = TRUE)

  out <- printed(b)
  expect_match(out, "Balance measures across subclasses")
  expect_match(out, "Sample sizes by subclass")

  out <- printed(b, which.subclass = 1L)
  expect_match(out, "Subclass 1")

  #NULL shows every subclass; the default and NA show none.
  expect_true(grepl("Subclass 4", printed(b, which.subclass = NULL), fixed = TRUE))
  expect_false(grepl("Subclass 1", printed(b), fixed = TRUE))
  expect_false(grepl("Subclass 1", printed(b, which.subclass = NA), fixed = TRUE))

  expect_err(print(b, disp = "bogus"), "not allowed in `disp`")
})

test_that("longitudinal treatments print per-time-point tables", {
  b <- bal.tab(list(treat ~ age + educ,
                    nodegree ~ age + educ + treat),
               data = lalonde, un = TRUE, msm.summary = TRUE)

  out <- printed(b)
  expect_match(out, "Balance summary across all time points")

  out <- printed(b, which.time = 1L)
  expect_match(out, "Time")
})

test_that("nested shapes print without error", {
  local_null_device()

  combos <- list(
    `cluster+imp` = bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat,
                            cluster = cl_idx, imp = imp_idx, s.d.denom = "pooled"),
    `cluster+multi` = bal.tab(lalonde[c("age", "educ")], treat = lalonde$race,
                              cluster = cl_idx),
    `cluster+msm` = bal.tab(list(treat ~ age, treat ~ age + educ), data = lalonde,
                            cluster = cl_idx),
    `cluster+subclass` = bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat,
                                 subclass = sub_idx, cluster = cl_idx,
                                 s.d.denom = "pooled"),
    `imp+msm` = bal.tab(list(treat ~ age, treat ~ age + educ), data = lalonde,
                        imp = imp_idx)
  )

  for (nm in names(combos)) {
    expect_no_error(capture.output(print(combos[[nm]])))
  }
})

test_that("error messages from print are stable", {
  b <- b_bin()

  #These are the cli-rendered messages users see; pin them exactly.
  expect_snapshot(error = TRUE, print(b, disp = "bogus"))
  expect_snapshot(error = TRUE, print(b, stats = "ks.statistics"))
  expect_snapshot(error = TRUE, print(b, un = "yes"))
})

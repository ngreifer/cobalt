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

# ---------------------------------------------------------------------------
# Coverage of the deprecated `...` arguments, the `which.*` validation branches in
# each wrapper method, and the display paths that need a purpose-built object.

test_that("deprecated disp.means / disp.sds select the mean and SD columns", {
  b <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat, weights = w_fixed,
               s.d.denom = "pooled", un = TRUE, disp = c("means", "sds"),
               quick = FALSE)

  out <- printed(b, disp.means = TRUE)
  expect_match(out, "M.0.Adj")

  out <- printed(b, disp.sds = TRUE)
  expect_match(out, "SD.0.Adj")

  #Setting them FALSE removes the corresponding columns.
  expect_false(grepl("M.0.Adj", printed(b, disp.means = FALSE, disp.sds = FALSE),
                     fixed = TRUE))

  #Requesting a quantity that `quick = TRUE` skipped warns.
  bq <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat, weights = w_fixed,
                s.d.denom = "pooled")
  expect_wrn(print(bq, disp.means = TRUE), "cannot be set to")
})

test_that("deprecated disp.<stat> toggles a statistic on and off", {
  b <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat, weights = w_fixed,
               s.d.denom = "pooled", quick = FALSE,
               stats = c("mean.diffs", "variance.ratios", "ks.statistics"))

  #Each toggle both adds and removes its own column.
  toggles <- list(disp.v.ratio = "V.Ratio", disp.ks = "KS")
  for (nm in names(toggles)) {
    col <- toggles[[nm]]
    on <- do.call(printed, c(list(b), setNames(list(TRUE), nm)))
    off <- do.call(printed, c(list(b), setNames(list(FALSE), nm)))
    expect_true(grepl(col, on, fixed = TRUE), info = nm)
    expect_false(grepl(col, off, fixed = TRUE), info = nm)
  }

  #A statistic that was not computed cannot be switched on.
  bq <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat, weights = w_fixed,
                s.d.denom = "pooled", stats = "mean.diffs")
  expect_wrn(print(bq, disp.ks = TRUE), "cannot be set to")
})

test_that("a deprecated <stat>.threshold of NULL suppresses the threshold column", {
  b <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat, weights = w_fixed,
               s.d.denom = "pooled", thresholds = c(m = .1))

  expect_match(printed(b), "M.Threshold")
  expect_false(grepl("M.Threshold", printed(b, m.threshold = NULL), fixed = TRUE))
})

test_that("un = TRUE on an object built with un = FALSE warns", {
  b <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat, weights = w_fixed,
               s.d.denom = "pooled", un = FALSE)

  expect_wrn(print(b, un = TRUE), "cannot be set to")
})

test_that("imbalanced.only without a threshold warns", {
  b <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat, weights = w_fixed,
               s.d.denom = "pooled")

  expect_wrn(print(b, imbalanced.only = TRUE), "threshold must be specified")
})

test_that("imbalanced.only reports when everything is balanced", {
  #A very loose threshold leaves nothing imbalanced.
  b <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat, weights = w_fixed,
               s.d.denom = "pooled", thresholds = c(m = 100))

  expect_match(printed(b, imbalanced.only = TRUE), "All covariates are balanced")
})

test_that("disp.thresholds accepts unnamed and partially-named input", {
  b <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat, weights = w_fixed,
               s.d.denom = "pooled", stats = c("mean.diffs", "variance.ratios"),
               thresholds = c(m = .1, v = 2))

  #A single unnamed value applies to every threshold.
  expect_false(grepl("Threshold", printed(b, disp.thresholds = FALSE), fixed = TRUE))

  #More entries than there are thresholds is an error.
  expect_err(print(b, disp.thresholds = c(FALSE, FALSE, FALSE)),
             "more entries were given")

  #A name that is not a threshold warns.
  expect_wrn(print(b, disp.thresholds = c(ks = FALSE)), "threshold")
})

test_that("two weight sets with different methods are reported side by side", {
  b <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat,
               weights = data.frame(m = w_fixed, w = sw_fixed),
               method = c("matching", "weighting"), s.d.denom = "pooled")

  out <- printed(b)
  expect_match(out, "Diff.m")
  expect_match(out, "Diff.w")

  #One sample-size row per weight set, plus the unadjusted row.
  expect_setequal(rownames(b$Observations), c("All", "m", "w"))
})

test_that("pairwise = FALSE prints balance by treatment group", {
  b <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$race, weights = w_fixed,
               pairwise = FALSE, un = TRUE)

  out <- printed(b, which.treat = NULL)
  expect_match(out, "Balance by treatment group")
  expect_match(out, "vs.")
  expect_setequal(names(b$Pair.Balance),
                  c("All vs. black", "All vs. hispan", "All vs. white"))
})

test_that("which.cluster accepts every documented form and warns otherwise", {
  b <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat, cluster = cl_idx,
               weights = w_fixed, s.d.denom = "pooled", cluster.summary = TRUE)

  expect_match(printed(b, which.cluster = 1L), "Cluster")
  expect_match(printed(b, which.cluster = c("a", "b")), "Cluster")

  expect_wrn(print(b, which.cluster = 99L), "are cluster indices")
  expect_wrn(print(b, which.cluster = "zzz"), "are cluster names")
  expect_wrn(print(b, which.cluster = TRUE), "must be")
})

test_that("which.imp accepts every documented form and warns otherwise", {
  b <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat, imp = imp_idx,
               weights = w_fixed, s.d.denom = "pooled", imp.summary = TRUE)

  expect_match(printed(b, which.imp = 1L), "Imputation")
  expect_wrn(print(b, which.imp = 99L), "are imputation numbers")
  expect_wrn(print(b, which.imp = TRUE), "must be")
})

test_that("which.treat accepts every documented form and warns otherwise", {
  b <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$race, weights = w_fixed,
               multi.summary = TRUE)

  expect_match(printed(b, which.treat = 1L), "vs.")
  expect_match(printed(b, which.treat = "black"), "black")
  expect_wrn(print(b, which.treat = 99L), "correspond to treatment values")
  expect_wrn(print(b, which.treat = "zzz"), "correspond to treatment values")
  expect_wrn(print(b, which.treat = TRUE), "must be")
})

test_that("which.time accepts every documented form and warns otherwise", {
  b <- bal.tab(list(treat ~ age + educ, nodegree ~ age + educ + treat),
               data = lalonde, msm.summary = TRUE)

  expect_match(printed(b, which.time = 1L), "Time")
  expect_wrn(print(b, which.time = 99L), "are treatment time points")
  expect_wrn(print(b, which.time = TRUE), "must be")
})

test_that("a subclass bal.tab with thresholds prints tallies and max imbalance", {
  b <- bal.tab(lalonde[c("age", "educ", "married")], treat = lalonde$treat,
               subclass = sub_idx, s.d.denom = "pooled", un = TRUE,
               subclass.summary = TRUE, thresholds = c(m = .1))

  out <- printed(b)
  expect_match(out, "Balance measures across subclasses")
  expect_match(out, "Balance tally for")
  expect_match(out, "Variable with the greatest")

  #The per-subclass tables and their threshold column.
  out <- printed(b, which.subclass = NULL)
  expect_match(out, "Subclass 1")
  expect_match(out, "M.Threshold")
})

test_that("the subclass print method honours stats, disp.call, and disp.thresholds", {
  b <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat, subclass = sub_idx,
               s.d.denom = "pooled", un = TRUE, subclass.summary = TRUE,
               stats = c("mean.diffs", "variance.ratios"),
               thresholds = c(m = .1), quick = FALSE)

  out <- printed(b, stats = "mean.diffs")
  expect_match(out, "Diff.Adj")
  expect_false(grepl("V.Ratio", out, fixed = TRUE))

  b_quick <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat,
                     subclass = sub_idx, s.d.denom = "pooled",
                     subclass.summary = TRUE, stats = "mean.diffs")
  expect_err(print(b_quick, stats = "ks.statistics"), "cannot contain")

  expect_false(grepl("M.Threshold", printed(b, disp.thresholds = c(m = FALSE)),
                     fixed = TRUE))

  #`disp.call` warns when the object carries no call.
  expect_wrn(print(b, disp.call = TRUE), "does not have a call")

  #`disp.thresholds` and `disp.call` together used to abort: the loop over
  #`names(disp.thresholds)` overwrote `x`, and `disp.call` then read `x$call` off
  #a character vector. Either argument alone was fine, which is why the two tests
  #above did not catch it.
  expect_wrn(out2 <- printed(b, disp.thresholds = c(m = FALSE), disp.call = TRUE),
             "does not have a call")
  expect_false(grepl("M.Threshold", out2, fixed = TRUE))
  expect_match(out2, "Diff.Adj")

  #`disp` selects the mean columns on the subclass path too.
  expect_match(printed(b, disp = "means"), "M.0.Adj")

  #`which.subclass` validation.
  expect_wrn(print(b, which.subclass = 99L), "subclass")
  expect_wrn(print(b, which.subclass = "a"), "subclass")
})

test_that("print() rewrites the .all and .none shorthands", {
  b <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$race, weights = w_fixed,
               multi.summary = TRUE)

  #The rewrite inspects print()'s own call, so these must be called directly
  #rather than through a wrapper: passing `.all` down through another frame
  #leaves it to be evaluated as an ordinary object, and it does not exist.
  expect_match(squish(capture.output(print(b, which.treat = .all))), "vs.")
  expect_false(grepl("vs.", squish(capture.output(print(b, which.treat = .none))),
                     fixed = TRUE))

  #NULL and NA are the equivalent literal values.
  expect_match(printed(b, which.treat = NULL), "vs.")
  expect_false(grepl("vs.", printed(b, which.treat = NA), fixed = TRUE))
})

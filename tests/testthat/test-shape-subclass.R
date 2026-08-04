#`bal.tab()` with subclassification: the arguments handled by
#`base.bal.tab.subclass()` rather than by `print()`.

test_that("subclasses require a binary treatment", {
  covs <- lalonde[c("age", "educ")]

  expect_err(bal.tab(covs, treat = lalonde$re75, subclass = sub_idx),
             "not yet compatible with continuous treatments")
  expect_err(bal.tab(covs, treat = lalonde$race, subclass = sub_idx),
             "not currently compatible with subclasses")
})

test_that("which.subclass and disp.subclass are honoured at bal.tab() time", {
  covs <- lalonde[c("age", "educ")]

  #The deprecated `disp.subclass` selects every subclass.
  b <- bal.tab(covs, treat = lalonde$treat, subclass = sub_idx,
               s.d.denom = "pooled", disp.subclass = TRUE)
  expect_match(squish(capture.output(print(b))), "Subclass 1")

  #`which.subclass` given to bal.tab() rather than print().
  b2 <- bal.tab(covs, treat = lalonde$treat, subclass = sub_idx,
                s.d.denom = "pooled", which.subclass = 1:2)
  out <- squish(capture.output(print(b2)))
  expect_match(out, "Subclass 1")
  expect_false(grepl("Subclass 4", out, fixed = TRUE))

  #`.all` selects every subclass and turns the summary off.
  b3 <- bal.tab(covs, treat = lalonde$treat, subclass = sub_idx,
                s.d.denom = "pooled", which.subclass = .all)
  expect_match(squish(capture.output(print(b3))), "Subclass 4")

  #Indices that match nothing leave the summary on.
  b4 <- bal.tab(covs, treat = lalonde$treat, subclass = sub_idx,
                s.d.denom = "pooled", which.subclass = c(90L, 91L))
  expect_true(is_not_null(b4[["Balance.Across.Subclass"]]))
})

test_that("thresholds are computed per subclass", {
  covs <- lalonde[c("age", "educ", "married")]

  b <- bal.tab(covs, treat = lalonde$treat, subclass = sub_idx,
               s.d.denom = "pooled", subclass.summary = TRUE,
               thresholds = c(m = .1))

  #A tally and a max-imbalance table across subclasses.
  expect_true(is_not_null(b[["Balanced.mean.diffs.Subclass"]]) ||
                is_not_null(b[["Balanced.mean.diffs"]]))
  expect_match(squish(capture.output(print(b))), "Balance tally for")

  #The threshold column reaches the per-subclass tables.
  expect_true("M.Threshold" %in% names(b[["Subclass.Balance"]][[1L]]))
})

test_that("quick, s.weights, and raw statistics work with subclasses", {
  covs <- lalonde[c("age", "educ", "married")]

  #`quick = FALSE` computes the columns that are otherwise skipped, and the
  #result must still print.
  b <- bal.tab(covs, treat = lalonde$treat, subclass = sub_idx,
               s.d.denom = "pooled", quick = FALSE, subclass.summary = TRUE)
  expect_no_error(capture.output(print(b)))
  expect_true(is_not_null(b[["Balance.Across.Subclass"]]))

  #Sampling weights change the reported balance.
  b_sw <- bal.tab(covs, treat = lalonde$treat, subclass = sub_idx,
                  s.d.denom = "pooled", s.weights = sw_fixed, un = TRUE)
  b_no <- bal.tab(covs, treat = lalonde$treat, subclass = sub_idx,
                  s.d.denom = "pooled", un = TRUE)
  expect_false(isTRUE(all.equal(b_sw$Balance.Across.Subclass$Diff.Un,
                                b_no$Balance.Across.Subclass$Diff.Un)))

  #Raw statistics skip the standardization factor entirely.
  expect_no_message(bal.tab(covs, treat = lalonde$treat, subclass = sub_idx,
                            binary = "raw", continuous = "raw"))

  #Interactions and polynomials under subclassification.
  b_int <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat,
                   subclass = sub_idx, s.d.denom = "pooled", int = TRUE)
  expect_true(any(startsWith(rownames(b_int$Balance.Across.Subclass), "age * educ")))

  b_poly <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat,
                    subclass = sub_idx, s.d.denom = "pooled", poly = 2)
  expect_gt(nrow(b_poly$Balance.Across.Subclass), 2L)
})

test_that("subclassification from a matchit object is accepted", {
  skip_if_not_installed("MatchIt")

  m <- fx("matchit_sub")

  b <- bal.tab(m, method = "subclassification", subclass.summary = TRUE)
  expect_s3_class(b, "bal.tab.subclass")
  expect_no_error(capture.output(print(b)))
})

test_that("the subclass sample-size table counts discarded units correctly", {
  skip_if_not_installed("MatchIt")

  m <- MatchIt::matchit(treat ~ age + educ + re74, data = lalonde,
                        method = "subclass", subclass = 4, discard = "both")

  nn <- bal.tab(m)$Observations

  expect_true("Discarded" %in% colnames(nn))

  #Every column's Total is the sum of its treatment-group rows. The `Discarded`
  #column previously reported `length(treat)` here instead of `sum(discarded)`.
  groups <- setdiff(rownames(nn), "Total")

  for (col in colnames(nn)) {
    expect_equal(nn["Total", col], sum(nn[groups, col]), info = col)
  }

  expect_equal(nn["Total", "Discarded"], sum(m$discarded))
  expect_equal(nn["Total", "All"], nrow(lalonde))

  #With nothing discarded the column is dropped entirely.
  m2 <- MatchIt::matchit(treat ~ age + educ + re74, data = lalonde,
                         method = "subclass", subclass = 4)

  expect_false("Discarded" %in% colnames(bal.tab(m2)$Observations))
})

test_that("longitudinal treatments accept data frames and addl/distance lists", {
  covs <- lalonde[c("age", "educ")]
  ps <- fitted(glm(treat ~ age + educ, data = lalonde, family = binomial))

  #The data.frame.list method, with a treatment per time point.
  b <- bal.tab(list(covs, covs),
               treat.list = list(lalonde$treat, lalonde$nodegree))
  expect_s3_class(b, "bal.tab.msm")
  expect_length(b$Time.Balance, 2L)

  #`addl.list` and `distance.list` are applied per time point.
  b2 <- bal.tab(list(covs, covs),
                treat.list = list(lalonde$treat, lalonde$nodegree),
                addl.list = list(lalonde["re74"], lalonde["re75"]),
                distance.list = list(ps, ps))
  expect_true("re74" %in% rownames(b2$Time.Balance[[1L]]$Balance))
  expect_true("re75" %in% rownames(b2$Time.Balance[[2L]]$Balance))
  expect_identical(rownames(b2$Time.Balance[[1L]]$Balance)[1L], "distance")
})

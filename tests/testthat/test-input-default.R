#`bal.tab.default()`, which duck-types an arbitrary list by the names and types of
#its components. This is the documented escape hatch for objects from packages
#cobalt does not support directly, and it needs no package from Suggests.
#
#The accepted names per component are listed in the `Q` table in `x2base.default()`.

ps_score <- function() {
  fitted(glm(treat ~ age + educ + re74, data = lalonde, family = binomial))
}

test_that("a list of `treat` and `covs` behaves like the data.frame method", {
  covs <- fx("covs")

  b_list <- bal.tab(list(treat = lalonde$treat, covs = covs), s.d.denom = "pooled")
  b_df <- bal.tab(covs, treat = lalonde$treat, s.d.denom = "pooled")

  expect_s3_class(b_list, "bal.tab.bin")
  expect_equal(b_list$Balance, b_df$Balance)
})

test_that("`treat` and `covs` accept their documented aliases", {
  covs <- lalonde[c("age", "educ")]
  ref <- bal.tab(list(treat = lalonde$treat, covs = covs), s.d.denom = "pooled")$Balance

  for (a in c("treat", "tr")) {
    L <- list(covs = covs)
    L[[a]] <- lalonde$treat
    expect_equal(bal.tab(L, s.d.denom = "pooled")$Balance, ref, info = a)
  }

  for (a in c("covs", "covariates", "x")) {
    L <- list(treat = lalonde$treat)
    L[[a]] <- covs
    expect_equal(bal.tab(L, s.d.denom = "pooled")$Balance, ref, info = a)
  }
})

test_that("`weights` and `s.weights` accept their aliases", {
  covs <- lalonde[c("age", "educ")]

  ref_w <- bal.tab(list(treat = lalonde$treat, covs = covs, weights = w_fixed),
                   s.d.denom = "pooled")$Balance
  for (a in c("weights", "w", "wts")) {
    L <- list(treat = lalonde$treat, covs = covs)
    L[[a]] <- w_fixed
    expect_equal(bal.tab(L, s.d.denom = "pooled")$Balance, ref_w, info = a)
  }

  ref_sw <- bal.tab(list(treat = lalonde$treat, covs = covs, s.weights = sw_fixed),
                    un = TRUE, s.d.denom = "pooled")$Balance
  for (a in c("s.weights", "sw", "sweights", "sampw")) {
    L <- list(treat = lalonde$treat, covs = covs)
    L[[a]] <- sw_fixed
    expect_equal(bal.tab(L, un = TRUE, s.d.denom = "pooled")$Balance, ref_sw, info = a)
  }

  #Sampling weights must actually change the unadjusted differences.
  plain <- bal.tab(list(treat = lalonde$treat, covs = covs), un = TRUE,
                   s.d.denom = "pooled")$Balance
  expect_false(isTRUE(all.equal(ref_sw$Diff.Un, plain$Diff.Un)))
})

test_that("a `distance` in the input list is used and keeps its component name", {
  #Regression test: the distance component used to be dropped silently because
  #the name attribute was lost while coercing it to a data frame, and
  #`setNames(., NULL)` then left the column unnamed.
  covs <- lalonde[c("age", "educ")]
  ps <- ps_score()

  for (a in c("distance", "ps", "pscore", "p.score", "propensity.score")) {
    L <- list(treat = lalonde$treat, covs = covs)
    L[[a]] <- ps
    b <- bal.tab(L, s.d.denom = "pooled")
    expect_identical(rownames(b$Balance), c(a, "age", "educ"), info = a)
    expect_identical(as.character(b$Balance$Type[1L]), "Distance", info = a)
  }

  #Passing it through `...` instead still names the row "distance".
  b <- bal.tab(list(treat = lalonde$treat, covs = covs), distance = ps,
               s.d.denom = "pooled")
  expect_identical(rownames(b$Balance), c("distance", "age", "educ"))
})

test_that("`subclass`, `match.strata`, and `estimand` aliases are honored", {
  covs <- lalonde[c("age", "educ")]

  for (a in c("subclass", "strata")) {
    L <- list(treat = lalonde$treat, covs = covs)
    L[[a]] <- sub_idx
    expect_s3_class(bal.tab(L, s.d.denom = "pooled"), "bal.tab.subclass")
  }

  L <- list(treat = lalonde$treat, covs = covs, match.strata = sub_idx)
  expect_s3_class(bal.tab(L, s.d.denom = "pooled"), "bal.tab.bin")

  #`estimand` selects the default `s.d.denom`.
  for (a in c("estimand", "target")) {
    L <- list(treat = lalonde$treat, covs = covs)
    L[[a]] <- "ATT"
    expect_equal(bal.tab(L)$Balance$Diff.Un,
                 unname(col_w_smd(covs, treat = lalonde$treat, s.d.denom = "treated")),
                 info = a)
  }
})

test_that("`focal` accepts its documented `treatATT` alias", {
  covs <- lalonde[c("age", "educ")]

  #The component names are lowercased before the aliases are matched, so a
  #mixed-case alias could never match and `treatATT` was silently ignored.
  ref <- suppressMessages(
    bal.tab(list(treat = lalonde$race, covs = covs, weights = w_fixed,
                 focal = "white"))
  )
  alias <- suppressMessages(
    bal.tab(list(treat = lalonde$race, covs = covs, weights = w_fixed,
                 treatATT = "white"))
  )
  none <- suppressMessages(
    bal.tab(list(treat = lalonde$race, covs = covs, weights = w_fixed))
  )

  expect_equal(alias, ref)
  expect_false(isTRUE(all.equal(alias, none)))
})

test_that("`focal` restricts non-pairwise comparisons for multi-category treatments", {
  covs <- lalonde[c("age", "educ")]

  L <- list(treat = lalonde$race, covs = covs, focal = "white")
  b <- bal.tab(L, pairwise = FALSE)
  expect_setequal(names(b$Pair.Balance), c("white vs. black", "white vs. hispan"))

  #Without `focal`, `pairwise = FALSE` compares each group to the whole sample.
  b0 <- bal.tab(list(treat = lalonde$race, covs = covs), pairwise = FALSE)
  expect_true(all(startsWith(names(b0$Pair.Balance), "All vs. ")))
})

test_that("a list carrying `formula` and `data` is accepted", {
  ref <- bal.tab(treat ~ age + educ, data = lalonde, s.d.denom = "pooled")$Balance

  for (a in c("formula", "form")) {
    L <- list(data = lalonde)
    L[[a]] <- treat ~ age + educ
    expect_equal(bal.tab(L, s.d.denom = "pooled")$Balance, ref, info = a)
  }
})

test_that("a fitted glm is accepted through the default method", {
  g <- glm(treat ~ age + educ + re74, data = lalonde, family = binomial)

  b <- bal.tab(g, s.d.denom = "pooled")
  expect_s3_class(b, "bal.tab.bin")
  expect_identical(rownames(b$Balance), c("age", "educ", "re74"))
})

test_that("longitudinal treatments are accepted via covs.list and treat.list", {
  covs <- lalonde[c("age", "educ")]

  b <- bal.tab(list(covs.list = list(covs, covs),
                    treat.list = list(lalonde$treat, lalonde$nodegree)))
  expect_s3_class(b, "bal.tab.msm")
  expect_length(b$Time.Balance, 2L)

  #A distance in the object applies to every time point.
  ps <- ps_score()
  b <- bal.tab(list(covs.list = list(covs, covs),
                    treat.list = list(lalonde$treat, lalonde$nodegree),
                    distance = ps))
  expect_identical(rownames(b$Time.Balance[[1L]]$Balance), c("distance", "age", "educ"))
  expect_identical(rownames(b$Time.Balance[[2L]]$Balance), c("distance", "age", "educ"))
})

test_that("arguments in `...` take precedence over components of the object", {
  covs <- lalonde[c("age", "educ")]

  #The object carries one set of weights; `...` supplies another.
  L <- list(treat = lalonde$treat, covs = covs, weights = w_fixed)

  b_obj <- bal.tab(L, s.d.denom = "pooled")
  b_dots <- bal.tab(L, weights = sw_fixed, s.d.denom = "pooled")

  expect_false(isTRUE(all.equal(b_obj$Balance$Diff.Adj, b_dots$Balance$Diff.Adj)))
  expect_equal(b_dots$Balance$Diff.Adj,
               bal.tab(list(treat = lalonde$treat, covs = covs),
                       weights = sw_fixed, s.d.denom = "pooled")$Balance$Diff.Adj)
})

test_that("the default method rejects malformed input", {
  covs <- lalonde[c("age", "educ")]

  expect_err(bal.tab(1:10),
             "input object must be an appropriate list, data frame, formula")
  expect_err(bal.tab("a"),
             "input object must be an appropriate list, data frame, formula")

  #`bal.tab.default()` has 19 named formals, so a few unnamed arguments are
  #matched positionally (here to `stats`) and rejected by that argument's own
  #check. The "must be named" guard is only reached once the formals are
  #exhausted.
  expect_err(bal.tab(list(treat = lalonde$treat, covs = covs), 1),
             "`stats` must be a character vector")
  expect_err(do.call(bal.tab, c(list(list(treat = lalonde$treat, covs = covs)),
                                as.list(1:20))),
             "must be named")

  #A list with no usable covariates.
  expect_err(bal.tab(list(treat = lalonde$treat)),
             "`covs` data frame must be specified")
})

test_that("subclasses, match.strata, and focal are rejected for longitudinal input", {
  covs <- lalonde[c("age", "educ")]
  L <- list(covs.list = list(covs, covs),
            treat.list = list(lalonde$treat, lalonde$nodegree))

  expect_err(bal.tab(L, subclass = sub_idx),
             "subclasses are not allowed with longitudinal treatments")
  expect_err(bal.tab(L, match.strata = sub_idx),
             "matching strata are not allowed with longitudinal treatments")
  expect_err(bal.tab(L, focal = "1"),
             "not allowed with longitudinal treatments")
})

test_that("malformed covs.list and treat.list are rejected", {
  covs <- lalonde[c("age", "educ")]

  #Entries that are not data frames fail the component's type check, which clears
  #`covs.list` entirely and reports it as missing.
  expect_err(bal.tab(list(covs.list = list(1:10),
                          treat.list = list(lalonde$treat))),
             "`covs.list` must be specified")
  expect_err(bal.tab(list(covs.list = list(covs, 1:10),
                          treat.list = list(lalonde$treat, lalonde$treat))),
             "`covs.list` must be specified")

  #A `treat.list` of the wrong length is caught explicitly.
  expect_err(bal.tab(list(covs.list = list(covs, covs),
                          treat.list = list(lalonde$treat))),
             "must be a list of treatment statuses at each time point")

  #Entries of the wrong length are caught downstream.
  expect_err(bal.tab(list(covs.list = list(covs, covs),
                          treat.list = list(lalonde$treat, 1:10))),
             "must have the same number of units")
})

#`get_covs_from_formula()` turns a formula into the covariate matrix balance is
#assessed on, and attaches `co.names`: for each column, its name decomposed into
#components and their types (`base`, `fsep`, `level`, `isep`, `na`). Everything that
#labels a variable reads that structure -- `.get_C2()`, `love.plot()`, `var.names()`,
#`.baltal()` -- so these tests pin the decomposition, not just the column names.

skip_on_cran()

gcf <- function(...) cobalt:::get_covs_from_formula(...)

#The parsed name of one column, as a "component|type" string per part.
parsed <- function(covs, col) {
  x <- attr(covs, "co.names")[[col]]
  paste(x[["component"]], x[["type"]], sep = "|", collapse = " ")
}

types <- function(covs, col) {
  paste(attr(covs, "co.names")[[col]][["type"]], collapse = "/")
}

test_that("get_covs_from_formula() handles the basic forms", {
  d <- lalonde

  #Plain numeric variables pass through as single `base` components.
  covs <- gcf(treat ~ age + educ, d)
  expect_identical(colnames(covs), c("age", "educ"))
  expect_identical(parsed(covs, "age"), "age|base")
  expect_equal(covs[, "age"], d$age, ignore_attr = TRUE)

  #The response is excluded and the intercept dropped.
  expect_false("treat" %in% colnames(covs))
  expect_false("(Intercept)" %in% colnames(covs))

  #A one-sided formula works.
  expect_identical(colnames(gcf(~ age + educ, d)), c("age", "educ"))

  #`.` expands to every other column.
  expect_identical(colnames(gcf(treat ~ ., d[c("treat", "age", "educ")])),
                   c("age", "educ"))

  #And a term can be removed from it.
  expect_identical(colnames(gcf(treat ~ . - age, d[c("treat", "age", "educ")])), "educ")

  #A formula given as a string is accepted.
  expect_identical(colnames(gcf("treat ~ age", d)), "age")

  #An explicit intercept-suppression makes no difference; there was never one.
  expect_identical(colnames(gcf(treat ~ 0 + age, d)), "age")
})

test_that("get_covs_from_formula() splits factors into named levels", {
  d <- lalonde
  d$ch <- as.character(d$race)
  d$lg <- d$married == 1

  #A factor becomes one column per level, named `variable<sep>level`, and the
  #decomposition records which part is which.
  covs <- gcf(treat ~ race, d)
  expect_identical(colnames(covs), c("race_black", "race_hispan", "race_white"))
  expect_identical(parsed(covs, "race_black"), "race|base _|fsep black|level")

  #The dummy really is the indicator.
  expect_equal(covs[, "race_black"], as.numeric(d$race == "black"), ignore_attr = TRUE)

  #Character and logical variables split the same way.
  expect_identical(colnames(gcf(treat ~ ch, d)), c("ch_black", "ch_hispan", "ch_white"))
  expect_identical(colnames(gcf(treat ~ lg, d)), c("lg_FALSE", "lg_TRUE"))
  expect_identical(types(gcf(treat ~ lg, d), "lg_TRUE"), "base/fsep/level")
})

test_that("get_covs_from_formula() records interactions with a separator component", {
  d <- lalonde

  #`:` gives the interaction alone; `*` gives the main effects too.
  covs <- gcf(treat ~ age:educ, d)
  expect_identical(colnames(covs), "age * educ")
  expect_identical(parsed(covs, "age * educ"), "age|base  * |isep educ|base")
  expect_equal(covs[, "age * educ"], d$age * d$educ, ignore_attr = TRUE)

  expect_identical(colnames(gcf(treat ~ age * educ, d)),
                   c("age", "educ", "age * educ"))

  #Interacting with a factor interacts with each of its dummies, and the
  #decomposition carries both separators.
  covs_f <- gcf(treat ~ age:race, d)
  expect_identical(colnames(covs_f),
                   c("age * race_black", "age * race_hispan", "age * race_white"))
  expect_identical(types(covs_f, "age * race_black"), "base/isep/base/fsep/level")
})

test_that("get_covs_from_formula() keeps transformations as opaque names", {
  d <- lalonde

  #A transformation is one column whose name is the call, recorded as a single `base`.
  for (f in list(treat ~ I(educ^2), treat ~ log(age))) {
    covs <- gcf(f, d)
    expect_length(colnames(covs), 1L)
    expect_identical(types(covs, colnames(covs)), "base")
  }

  expect_identical(colnames(gcf(treat ~ I(educ^2), d)), "I(educ^2)")
  expect_equal(gcf(treat ~ I(educ^2), d)[, 1L], d$educ^2, ignore_attr = TRUE)
  expect_equal(gcf(treat ~ log(age), d)[, 1L], log(d$age), ignore_attr = TRUE)

  #A multi-column term is numbered with the factor separator.
  covs_p <- gcf(treat ~ poly(age, 2), d)
  expect_identical(colnames(covs_p), c("poly(age, 2)_1", "poly(age, 2)_2"))
  expect_identical(types(covs_p, "poly(age, 2)_1"), "base")

  #A matrix column contributes its own columns.
  d$M <- as.matrix(d[c("age", "educ")])
  expect_identical(colnames(gcf(treat ~ M, d)), c("age", "educ"))
})

test_that("get_covs_from_formula() adds a missingness indicator", {
  d <- lalonde
  d$na <- replace(d$age, 1:5, NA)

  covs <- gcf(treat ~ na, d)

  #The variable keeps its column, and a companion indicator is added after it.
  expect_identical(colnames(covs), c("na", "na:<NA>"))
  expect_identical(parsed(covs, "na:<NA>"), "na|base :<NA>|na")
  expect_equal(covs[, "na:<NA>"], as.numeric(is.na(d$na)), ignore_attr = TRUE)

  #The variable's own missing values are filled so the column is usable.
  expect_false(anyNA(covs[, "na:<NA>"]))
})

test_that("get_covs_from_formula() honors the separators", {
  d <- lalonde

  covs <- gcf(treat ~ race + age:educ, d, factor_sep = ".", int_sep = " x ")

  expect_identical(colnames(covs),
                   c("race.black", "race.hispan", "race.white", "age x educ"))

  #The separators are recorded, so a later stage can substitute different ones.
  expect_identical(attr(attr(covs, "co.names"), "seps"),
                   c(factor = ".", int = " x "))

  #The separator is its own component, which is what makes substitution possible.
  expect_identical(parsed(covs, "race.black"), "race|base .|fsep black|level")
})

test_that("get_covs_from_formula() rejects what it cannot process", {
  d <- lalonde

  expect_err(gcf(treat ~ nope, d), 'The variable "nope" cannot be found')
  expect_err(gcf(treat ~ age, "nope"), "must be a data frame")
  expect_err(gcf(treat ~ ., NULL), "'.' is not allowed in formulas")

  #Non-finite values would make every balance statistic meaningless. The message
  #agrees with itself about number: `contain{?s}` had the marker on the wrong word, so
  #one variable used to "contain" and several "contains".
  d$z <- replace(d$age, 1L, Inf)
  d$z2 <- replace(d$educ, 1L, -Inf)

  expect_err(gcf(treat ~ z, d),
             "The variable `z` contains non-finite values, which are not allowed")
  expect_err(gcf(treat ~ z + z2, d),
             "The variables `z` and `z2` contain non-finite values, which are not allowed")

  #A function name used as a variable is caught rather than silently evaluated.
  expect_err(gcf(treat ~ mean, d), "invalid type (function) for variable")
})

test_that("co.names round-trips through the whole pipeline", {
  # The point of the decomposition is that downstream stages can rebuild a label from
  # it. This checks the two that do: `.get_C2()` substitutes new separators, and
  # `var.names(minimal = TRUE)` collapses split dummies back to the base variable.
  d <- lalonde

  b <- bal.tab(treat ~ race + age:educ, data = d, s.d.denom = "pooled",
               int_sep = " X ", factor_sep = "|")

  expect_true(all(c("race|black", "age X educ") %in% rownames(b$Balance)))

  b2 <- bal.tab(treat ~ race + age, data = d, s.d.denom = "pooled")
  expect_true("race" %in% var.names(b2, type = "vec", minimal = TRUE))
  expect_false("race_black" %in% var.names(b2, type = "vec", minimal = TRUE))
})

test_that("the balance tally counts against the labels in the table", {
  # `.baltal()` used to recover the numeric threshold by regex from the label
  # `"Balanced, <0.1"` that `.threshold_label()` had just generated. Both now derive the
  # labels from `.threshold_verdicts()`, so the tally's row names and the table's
  # verdicts come from one place.
  b <- bal.tab(lalonde[c("age", "educ", "married")], treat = lalonde$treat,
               s.d.denom = "pooled", weights = w_fixed, thresholds = c(m = .1))

  tally <- b$Balanced.mean.diffs

  expect_identical(rownames(tally), c("Balanced, <0.1", "Not Balanced, >0.1"))

  #The counts agree with the table they summarize.
  verdicts <- b$Balance$M.Threshold
  expect_identical(tally[["count"]],
                   c(sum(verdicts == "Balanced, <0.1"),
                     sum(verdicts == "Not Balanced, >0.1")))
  expect_identical(sum(tally[["count"]]), sum(nzchar(verdicts)))

  #A threshold needing rounding is rounded the same way in both.
  b3 <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat, s.d.denom = "pooled",
                weights = w_fixed, thresholds = c(m = 0.123456))
  expect_identical(rownames(b3$Balanced.mean.diffs),
                   c("Balanced, <0.123", "Not Balanced, >0.123"))
  expect_true(all(b3$Balance$M.Threshold %in%
                    c("", "Balanced, <0.123", "Not Balanced, >0.123")))
})

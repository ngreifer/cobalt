#The small standalone exports: `f.build()`, `var.names()`,
#`set.cobalt.options()`, and `get.cobalt.options()`.

test_that("f.build() assembles a formula from a response and covariates", {
  #From a character vector of names.
  f <- f.build("treat", c("age", "educ"))
  expect_s3_class(f, "formula")
  expect_identical(deparse(f), "treat ~ age + educ")

  #From a data frame, using its column names.
  f <- f.build("treat", lalonde[c("age", "educ")])
  expect_identical(deparse(f), "treat ~ age + educ")

  #Names needing backticks are quoted.
  d <- data.frame(`a b` = 1:2, check.names = FALSE)
  expect_s3_class(f.build("y", d), "formula")

  #With no response, a one-sided formula results.
  f <- f.build(rhs = c("age", "educ"))
  expect_length(f, 2L)

  #The result is usable by bal.tab().
  expect_s3_class(bal.tab(f.build("treat", c("age", "educ")), data = lalonde,
                          s.d.denom = "pooled"),
                  "bal.tab")
})

test_that("f.build() rejects invalid arguments", {
  expect_err(f.build("y", 1:3),
             "must be a vector of variable names or a data set with named variables")
  expect_err(f.build(list("y"), "x"),
             "must be a string containing the response variable")
})

test_that("var.names() extracts names as a vector or data frame", {
  b <- bal.tab(lalonde[c("age", "educ", "race")], treat = lalonde$treat,
               s.d.denom = "pooled")

  v <- var.names(b, type = "vec")
  expect_type(v, "character")
  expect_named(v)
  expect_true(all(c("age", "educ", "race_black") %in% v))

  d <- var.names(b, type = "df")
  expect_s3_class(d, "data.frame")
  expect_named(d, c("old", "new"))
  expect_identical(d$old, d$new)

  #The default is a vector when no file is given.
  expect_type(var.names(b), "character")

  #`minimal` collapses the split dummies back to their base variable.
  vm <- var.names(b, type = "vec", minimal = TRUE)
  expect_true("race" %in% vm)
  expect_false("race_black" %in% vm)
})

test_that("var.names() writes a csv and validates the filename", {
  b <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat,
               s.d.denom = "pooled")

  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)

  #A file name implies `type = "df"`.
  out <- var.names(b, file = f)
  expect_true(file.exists(f))
  expect_s3_class(out, "data.frame")

  back <- utils::read.csv(f)
  expect_named(back, c("old", "new"))

  #`type = "vec"` with a file warns and returns the vector without writing.
  f2 <- tempfile(fileext = ".csv")
  on.exit(unlink(f2), add = TRUE)
  expect_wrn(var.names(b, type = "vec", file = f2),
             "is compatible with a file name")

  expect_err(var.names(b, file = "out.txt"), "must end in")
})

test_that("var.names() rejects objects that are not bal.tabs", {
  expect_err(var.names(1), "no variable names were found in the object")
  expect_err(var.names(lalonde), "no variable names were found in the object")
})

test_that("set.cobalt.options() sets and get.cobalt.options() returns options", {
  #Scope every cobalt_* option to this test so it cannot leak into other files.
  local_cobalt_options()

  set.cobalt.options(binary = "std")
  expect_identical(getOption("cobalt_binary"), "std")
  expect_identical(get.cobalt.options("binary")[["binary"]], "std")

  #The option really changes bal.tab()'s behaviour.
  b_std <- bal.tab(lalonde[c("married")], treat = lalonde$treat,
                   s.d.denom = "pooled")
  set.cobalt.options(binary = "raw")
  b_raw <- bal.tab(lalonde[c("married")], treat = lalonde$treat,
                   s.d.denom = "pooled")
  expect_false(isTRUE(all.equal(b_std$Balance$Diff.Un, b_raw$Balance$Diff.Un)))

  #Several options at once, including a multi-valued one.
  set.cobalt.options(un = TRUE, disp = c("means", "sds"), cluster.fun = "mean")
  expect_true(getOption("cobalt_un"))
  expect_setequal(getOption("cobalt_disp"), c("means", "sds"))
  expect_identical(getOption("cobalt_cluster.fun"), "mean")

  #With no arguments, every set option is returned.
  all_opts <- get.cobalt.options()
  expect_type(all_opts, "list")
  expect_true("un" %in% names(all_opts))

  #`default = TRUE` clears the options that are not named.
  set.cobalt.options(binary = "std", default = TRUE)
  expect_identical(getOption("cobalt_binary"), "std")
  expect_null(getOption("cobalt_un"))

  #Setting an option to NULL removes it.
  set.cobalt.options(binary = NULL)
  expect_null(getOption("cobalt_binary"))
})

test_that("set.cobalt.options() rejects invalid input", {
  local_cobalt_options()

  expect_err(set.cobalt.options("std"), "must be named")
  expect_err(set.cobalt.options(binary = "std", binary = "raw"),
             "is present more than once")
  expect_err(set.cobalt.options(binary = "bogus"), "No options will be set")
  expect_err(set.cobalt.options(un = "yes"), "No options will be set")

  #`stats` only accepts "mean.diffs".
  expect_err(set.cobalt.options(stats = "ks.statistics"), "No options will be set")
  expect_no_error(set.cobalt.options(stats = "mean.diffs"))

  #An unrecognized name warns rather than erroring.
  expect_wrn(set.cobalt.options(nosuchoption = TRUE))

  #The free-text options accept any string.
  expect_no_error(set.cobalt.options(int_sep = " x ", factor_sep = "."))
})

test_that("get.cobalt.options() validates its arguments", {
  expect_err(get.cobalt.options(1),
             "must be strings containing the name of an option")
  expect_err(get.cobalt.options("bogus"), "not an acceptable option")
})

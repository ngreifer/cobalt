#`var.names`: the names the covariates are displayed under, wherever a `bal.tab` object
#is displayed. `love.plot()` has always taken it; `bal.tab()`, `print()`, `format()`, and
#`as.data.frame()` take it too, and all five resolve it the same way.

covs4 <- function() lalonde[c("age", "educ", "race", "re74")]

vn <- c(age = "Age", race = "Race/Eth", re74 = "Earnings 1974")

#What a plot displays, from the outside in the same order a table lists it.
plot_labels <- function(p) {
  rev(levels(factor(ggplot2::ggplot_build(p)$plot$data[["var"]])))
}

test_that("var.names() in bal.tab() renames the balance table without touching co.names", {
  b <- bal.tab(covs4(), treat = lalonde$treat, s.d.denom = "pooled", var.names = vn)

  expect_identical(rownames(b$Balance),
                   c("age", "educ", "race_black", "race_hispan", "race_white", "re74"))

  expect_identical(rownames(format(b)),
                   c("Age", "educ", "Race/Eth_black", "Race/Eth_hispan",
                     "Race/Eth_white", "Earnings 1974"))

  #`co.names` is the record of what each covariate is, and is what the replacements are
  #resolved against, so it keeps the stored names.
  expect_identical(names(attr(b, "print.options")[["co.names"]]),
                   rownames(b$Balance))

  #And so the object itself is unchanged by having been printed.
  expect_output(print(b))
  expect_identical(rownames(b$Balance),
                   c("age", "educ", "race_black", "race_hispan", "race_white", "re74"))
})

test_that("var.names() reaches every name a base variable is a component of", {
  b <- bal.tab(covs4(), treat = lalonde$treat, s.d.denom = "pooled", int = TRUE,
               var.names = vn)

  nms <- rownames(format(b))

  #Factor levels, and both sides of an interaction.
  expect_true("Race/Eth_black" %in% nms)
  expect_true("Age * Race/Eth_black" %in% nms)
  expect_true("educ * Earnings 1974" %in% nms)

  #A name given for a single expansion renames only that one.
  b2 <- bal.tab(covs4(), treat = lalonde$treat, s.d.denom = "pooled",
                var.names = c(race_black = "Black"))

  expect_identical(rownames(format(b2)),
                   c("age", "educ", "Black", "race_hispan", "race_white", "re74"))
})

test_that("var.names() accepts the same structures everywhere it is taken", {
  b <- bal.tab(covs4(), treat = lalonde$treat, s.d.denom = "pooled")

  expected <- c("Age", "educ", "race_black", "race_hispan", "race_white", "E74")

  #Named vector, named list, one-column data frame with row names, two-column data frame.
  expect_identical(rownames(format(b, var.names = c(age = "Age", re74 = "E74"))),
                   expected)

  expect_identical(rownames(format(b, var.names = list(age = "Age", re74 = "E74"))),
                   expected)

  expect_identical(rownames(format(b, var.names = data.frame(new = c("Age", "E74"),
                                                             row.names = c("age", "re74")))),
                   expected)

  expect_identical(rownames(format(b, var.names = data.frame(old = c("age", "re74"),
                                                             new = c("Age", "E74")))),
                   expected)

  #An entry naming nothing in the table is ignored, and an `NA` replacement is a no-op.
  expect_identical(rownames(format(b, var.names = c(nope = "X", age = "Age", re74 = "E74"))),
                   expected)

  expect_identical(rownames(format(b, var.names = c(educ = NA, age = "Age", re74 = "E74"))),
                   expected)
})

test_that("var.names() is rejected the same way in bal.tab() and love.plot()", {
  b <- bal.tab(covs4(), treat = lalonde$treat, s.d.denom = "pooled")

  for (bad in list(c("Age", "Educ"), 1:3)) {
    expect_err(bal.tab(covs4(), treat = lalonde$treat, s.d.denom = "pooled",
                       var.names = bad),
               "must be named")
    expect_err(love.plot(b, var.names = bad), "must be named")
  }

  expect_err(love.plot(b, var.names = mean),
             "not one of the accepted structures")

  #Two covariates under one name would be indistinguishable in a table and would sit on
  #top of each other in a plot.
  expect_err(print(b, var.names = c(age = "Same", educ = "Same")),
             'gives the name "Same" to more than one variable')

  expect_err(love.plot(b, var.names = c(age = "Same", educ = "Same")),
             'gives the name "Same" to more than one variable')
})

test_that("var.names() given to a display function adds to what bal.tab() was given", {
  b <- bal.tab(covs4(), treat = lalonde$treat, s.d.denom = "pooled", binary = "std",
               var.names = vn)

  #`educ` is added; `age` replaces the stored entry; the rest of the stored ones stand.
  expected <- c("AGE", "Education", "Race/Eth_black", "Race/Eth_hispan",
                "Race/Eth_white", "Earnings 1974")

  added <- c(educ = "Education", age = "AGE")

  expect_identical(rownames(format(b, var.names = added)), expected)

  expect_identical(unique(as.data.frame(b, var.names = added)[["variable"]]), expected)

  expect_identical(plot_labels(love.plot(b, var.names = added)), expected)

  expect_output(print(b, var.names = added), "AGE", fixed = TRUE)
})

test_that("var.names() carries to the balance tally and the greatest imbalance", {
  b <- bal.tab(covs4(), treat = lalonde$treat, s.d.denom = "pooled", var.names = vn,
               thresholds = c(m = .1))

  out <- capture.output(print(b))

  #A tally counts variables rather than naming them, so it is unaffected; the greatest
  #imbalance names one, and must name it the way the table above it does.
  expect_true(any(grepl("Race/Eth_black", out, fixed = TRUE)))
  expect_false(any(grepl("race_black", out, fixed = TRUE)))

  imbal <- b[["Max.Imbalance.mean.diffs"]]
  expect_identical(imbal[["Variable"]], "race_black")
})

test_that("var.names() applies to every segment of a segmented object", {
  b <- bal.tab(covs4(), treat = lalonde$treat, s.d.denom = "pooled", binary = "std",
               weights = w_fixed, un = TRUE, cluster = cl_idx, which.cluster = .all,
               cluster.summary = TRUE, var.names = vn)

  expected <- c("Age", "educ", "Race/Eth_black", "Race/Eth_hispan", "Race/Eth_white",
                "Earnings 1974")

  #The summary across clusters, and each cluster's own table, which `print()` reaches by
  #printing the child object.
  expect_identical(rownames(format(b)), expected)

  expect_setequal(as.data.frame(b)[["variable"]], expected)

  out <- capture.output(print(b))
  expect_false(any(grepl("race_black", out, fixed = TRUE)))
  expect_true(any(grepl("Race/Eth_black", out, fixed = TRUE)))

  #Each cluster's table is reached by printing the child, which carries the names of
  #its own accord.
  for (cl in names(b$Cluster.Balance)) {
    expect_identical(rownames(format(b$Cluster.Balance[[cl]])), expected, info = cl)
  }

  expect_identical(plot_labels(love.plot(b, which.cluster = .none, agg.fun = "max")),
                   expected)
})

test_that("var.names() reaches a covariate that appears only at a later time point", {
  #A longitudinal object's time points need not share covariates, and only the first
  #one's names are recorded at the top level.
  b <- bal.tab(list(treat ~ age + educ, treat ~ re74 + re75),
               data = lalonde, s.d.denom = "pooled", which.time = .all,
               var.names = c(age = "Age", re75 = "Earnings 1975"))

  expect_setequal(as.data.frame(b)[["variable"]],
                  c("Age", "educ", "re74", "Earnings 1975"))

  out <- capture.output(print(b))
  expect_true(any(grepl("Earnings 1975", out, fixed = TRUE)))
  expect_false(any(grepl("re75", out, fixed = TRUE)))
})

test_that("var.names() applies to a subclassified object", {
  b <- bal.tab(covs4(), treat = lalonde$treat, s.d.denom = "pooled", binary = "std",
               un = TRUE, subclass = sub_idx, which.subclass = .all,
               subclass.summary = TRUE, var.names = vn)

  expected <- c("Age", "educ", "Race/Eth_black", "Race/Eth_hispan", "Race/Eth_white",
                "Earnings 1974")

  expect_identical(rownames(format(b)), expected)
  expect_setequal(as.data.frame(b)[["variable"]], expected)
  expect_identical(plot_labels(love.plot(b)), expected)
})

test_that("var.names() reports the names already applied, ready to be edited", {
  b <- bal.tab(covs4(), treat = lalonde$treat, s.d.denom = "pooled", int = TRUE,
               var.names = vn)

  #Keyed by the stored names, which is what the replacements are resolved against, and
  #valued at the names on display.
  v <- var.names(b)

  expect_named(v, rownames(b$Balance))
  expect_identical(unname(v), rownames(format(b)))

  d <- var.names(b, type = "df")
  expect_identical(d[["old"]], names(v))
  expect_identical(d[["new"]], unname(v))

  #`minimal` reports the base variables under the names given for them.
  vm <- var.names(b, minimal = TRUE)
  expect_identical(vm, c(age = "Age", educ = "educ", race = "Race/Eth",
                         re74 = "Earnings 1974"))

  #Passed back unedited, either form reproduces the names it came from.
  expect_identical(rownames(format(b, var.names = v)), rownames(format(b)))
  expect_identical(rownames(format(b, var.names = vm)), rownames(format(b)))

  #Edited, each changes what its form names: a base variable everywhere it appears...
  vm["educ"] <- "Education"

  expect_identical(rownames(format(b, var.names = vm)),
                   sub("educ", "Education", rownames(format(b)), fixed = TRUE))

  #...and a whole name only where it is the whole name, since the full form gives every
  #displayed name an entry of its own.
  v["educ"] <- "Education"

  expect_identical(rownames(format(b, var.names = v)),
                   replace(rownames(format(b)), rownames(b$Balance) == "educ",
                           "Education"))
})

test_that("var.names() reports the stored names when none have been applied", {
  b <- bal.tab(covs4(), treat = lalonde$treat, s.d.denom = "pooled", int = TRUE)

  v <- var.names(b)

  expect_identical(unname(v), names(v))
  expect_identical(names(v), rownames(b$Balance))

  d <- var.names(b, type = "df")
  expect_identical(d[["old"]], d[["new"]])
})

test_that("var.names() reports a covariate that appears only at a later time point", {
  b <- bal.tab(list(treat ~ age + educ, treat ~ re74 + re75),
               data = lalonde, s.d.denom = "pooled",
               var.names = c(re75 = "Earnings 1975"))

  expect_identical(var.names(b),
                   c(age = "age", educ = "educ", re74 = "re74",
                     re75 = "Earnings 1975"))
})

test_that("var.names() gives bal.tab() and love.plot() the same names", {
  b <- bal.tab(covs4(), treat = lalonde$treat, s.d.denom = "pooled", binary = "std",
               weights = w_fixed, un = TRUE, int = TRUE, var.names = vn)

  expect_identical(plot_labels(love.plot(b)), rownames(format(b)))

  #And the names travel with the object into a `love.plot()` written as a `bal.tab()`
  #call, which re-evaluates the call rather than reading the object.
  p <- love.plot(bal.tab(covs4(), treat = lalonde$treat, s.d.denom = "pooled",
                         binary = "std", weights = w_fixed, int = TRUE, var.names = vn))

  expect_identical(plot_labels(p), rownames(format(b)))
})

# `factor_sep` and `int_sep`: the separators the names are displayed with. They are taken
# in the same places `var.names` is, and each being one value rather than a set of
# entries, one given later replaces the one in force.

test_that("the separators can be chosen at display time", {
  b <- bal.tab(covs4(), treat = lalonde$treat, s.d.denom = "pooled", binary = "std",
               int = TRUE)

  stored <- rownames(b$Balance)

  expect_identical(rownames(format(b, factor_sep = ": ")),
                   gsub("race_", "race: ", stored, fixed = TRUE))

  expect_identical(rownames(format(b, int_sep = " x ")),
                   gsub(" * ", " x ", stored, fixed = TRUE))

  both <- gsub(" * ", " x ", gsub("race_", "race: ", stored, fixed = TRUE), fixed = TRUE)

  expect_identical(rownames(format(b, factor_sep = ": ", int_sep = " x ")), both)

  expect_identical(unique(as.data.frame(b, factor_sep = ": ", int_sep = " x ")[["variable"]]),
                   both)

  expect_identical(plot_labels(love.plot(b, factor_sep = ": ", int_sep = " x ")), both)

  expect_output(print(b, factor_sep = ": "), "race: black", fixed = TRUE)

  #The object itself keeps the names it was built with.
  expect_identical(rownames(b$Balance), stored)
})

test_that("a separator given at display time replaces the one bal.tab() was given", {
  b <- bal.tab(covs4(), treat = lalonde$treat, s.d.denom = "pooled", binary = "std",
               int = TRUE, factor_sep = ": ", int_sep = " x ")

  #`bal.tab()`'s separators still decide the names the covariates are stored under.
  expect_true(all(c("race: black", "age x educ") %in% rownames(b$Balance)))

  #Displaying with the defaults puts them back.
  expect_true(all(c("race_black", "age * educ") %in% rownames(format(b, factor_sep = "_",
                                                                    int_sep = " * "))))

  #And with none given, the ones it was built with stand.
  expect_identical(rownames(format(b)), rownames(b$Balance))
})

test_that("the separators and var.names apply together", {
  b <- bal.tab(covs4(), treat = lalonde$treat, s.d.denom = "pooled", binary = "std",
               int = TRUE, var.names = vn)

  #A replacement is keyed by the name the covariate is stored under, whichever separators
  #are on display, and applies within a name displayed with the new ones.
  expected <- gsub("Race/Eth_", "Race/Eth: ", rownames(format(b)), fixed = TRUE)

  expect_identical(rownames(format(b, factor_sep = ": ")), expected)

  expect_identical(plot_labels(love.plot(b, factor_sep = ": ")), expected)

  #A name given for a whole variable replaces the separator along with the rest of it.
  expect_identical(rownames(format(b, factor_sep = ": ",
                                   var.names = c(race_black = "Black")))[3L],
                   "Black")
})

test_that("the separators carry to every segment and to the greatest imbalance", {
  b <- bal.tab(covs4(), treat = lalonde$treat, s.d.denom = "pooled", binary = "std",
               un = TRUE, weights = w_fixed, cluster = cl_idx, which.cluster = .all,
               cluster.summary = TRUE, thresholds = c(m = .1))

  out <- capture.output(print(b, factor_sep = ": "))

  expect_true(any(grepl("race: black", out, fixed = TRUE)))
  expect_false(any(grepl("race_black", out, fixed = TRUE)))

  for (cl in names(b$Cluster.Balance)) {
    expect_true("race: black" %in% rownames(format(b$Cluster.Balance[[cl]],
                                                   factor_sep = ": ")),
                info = cl)
  }
})

test_that("the separators must each be a single string", {
  b <- bal.tab(covs4(), treat = lalonde$treat, s.d.denom = "pooled")

  expect_err(format(b, factor_sep = c("a", "b")), "must be a string")
  expect_err(print(b, int_sep = 1), "must be a string")
  expect_err(love.plot(b, factor_sep = NA), "must be a string")
})

test_that("the separators can be set globally", {
  local_cobalt_options()

  set.cobalt.options(factor_sep = ": ")

  b <- bal.tab(covs4(), treat = lalonde$treat, s.d.denom = "pooled")

  expect_true("race: black" %in% rownames(b$Balance))
})

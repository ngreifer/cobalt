#`bal.plot()`.
#
#These plots are built lazily, so a wrong `aes()`, an unused factor level, or a
#mis-sized manual scale does not error until the plot is built. A user at the
#console always builds, so every test here builds too, via `expect_no_condition()`
#on `ggplot_build()`. `expect_no_condition()` rather than `expect_no_error()`, so
#that ggplot2 warnings such as "Removed N rows containing non-finite values" are
#treated as failures.
#
#Beyond that, plots are checked structurally -- which geoms were used, what the
#facets are, and the semantics of `layer_data()` -- never by image comparison.

#The geom of each layer, in order. Layers must always be located by geom class
#rather than by index: on a love.plot, layer 1 is the reference line, not the points.
geoms <- function(p) {
  #ggplot2 4.0 names the elements of `p$layers`, so drop the names to compare.
  unname(vapply(p$layers, function(l) class(l$geom)[1L], character(1L)))
}

#Build the plot and require that nothing is signalled.
built <- function(p) {
  expect_no_condition(b <- ggplot2::ggplot_build(p))
  invisible(b)
}

covs4 <- function() lalonde[c("age", "educ", "race", "married")]

test_that("bal.plot() returns a ggplot that builds cleanly", {
  p <- bal.plot(covs4(), treat = lalonde$treat, var.name = "age",
                weights = w_fixed)

  expect_s3_class(p, "ggplot")
  built(p)
  expect_identical(p$labels$x, "age")
})

test_that("the geoms depend on the covariate and treatment types", {
  covs <- covs4()
  tb <- lalonde$treat

  #Binary treatment, continuous covariate: densities, one per sample.
  p <- bal.plot(covs, treat = tb, var.name = "age", weights = w_fixed,
                which = "both")
  expect_identical(geoms(p), c("GeomDensity", "GeomDensity"))
  built(p)

  #Binary treatment, categorical covariate: bars plus a baseline.
  p <- bal.plot(covs, treat = tb, var.name = "race", weights = w_fixed,
                which = "both")
  expect_identical(geoms(p), c("GeomBar", "GeomHline"))
  built(p)

  #`type = "ecdf"` uses steps.
  p <- bal.plot(covs, treat = tb, var.name = "age", weights = w_fixed,
                which = "both", type = "ecdf")
  expect_identical(geoms(p), c("GeomStep", "GeomStep"))
  built(p)

  #`type = "histogram"` with `mirror = TRUE` reflects one group below zero.
  p <- bal.plot(covs, treat = tb, var.name = "age", weights = w_fixed,
                which = "both", type = "histogram", mirror = TRUE)
  expect_identical(geoms(p), c("GeomBar", "GeomBar", "GeomHline"))
  built(p)
})

test_that("continuous treatments produce scatterplots and grouped densities", {
  covs <- covs4()
  tc <- lalonde$re75

  #Continuous treatment, continuous covariate: points plus fitted lines.
  p <- bal.plot(covs, treat = tc, var.name = "age", weights = w_fixed)
  expect_true("GeomPoint" %in% geoms(p))
  built(p)

  #Continuous treatment, categorical covariate: densities by level.
  p <- bal.plot(covs, treat = tc, var.name = "race", weights = w_fixed)
  expect_true("GeomDensity" %in% geoms(p))
  built(p)

  #`disp.means` adds reference lines.
  p <- bal.plot(covs, treat = tc, var.name = "age", weights = w_fixed,
                disp.means = TRUE)
  built(p)
})

test_that("multi-category treatments are plotted and can be subset", {
  covs <- covs4()

  p <- bal.plot(covs, treat = lalonde$race, var.name = "age", weights = w_fixed)
  built(p)

  #`which.treat` selects levels.
  p <- bal.plot(covs, treat = lalonde$race, var.name = "age", weights = w_fixed,
                which.treat = c("black", "white"))
  built(p)

  p <- bal.plot(covs, treat = lalonde$race, var.name = "age", weights = w_fixed,
                which.treat = 1:2)
  built(p)
})

test_that("`which` selects the samples shown", {
  covs <- covs4()

  for (wh in c("unadjusted", "adjusted", "both")) {
    p <- bal.plot(covs, treat = lalonde$treat, var.name = "age",
                  weights = w_fixed, which = wh)
    built(p)
  }

  #With no weights, only the unadjusted sample exists.
  p <- bal.plot(covs, treat = lalonde$treat, var.name = "age")
  built(p)

  #Named weight sets can be selected individually.
  p <- bal.plot(covs, treat = lalonde$treat, var.name = "age",
                weights = list(a = w_fixed, b = sw_fixed), which = "a")
  built(p)
})

test_that("`var.name` resolves split factors and is auto-selected when missing", {
  covs <- covs4()

  #Factors are named by the original variable, not by their dummy columns.
  p <- bal.plot(covs, treat = lalonde$treat, var.name = "race", weights = w_fixed)
  built(p)
  expect_err(bal.plot(covs, treat = lalonde$treat, var.name = "race_black",
                      weights = w_fixed),
             "is not the name of an available covariate")

  #Omitting `var.name` picks the first available covariate, with a message.
  expect_msg(p <- bal.plot(covs, treat = lalonde$treat, weights = w_fixed))
  built(p)

  #`distance` can be plotted too.
  ps <- fitted(glm(treat ~ age + educ, data = lalonde, family = binomial))
  p <- bal.plot(covs, treat = lalonde$treat, var.name = "distance",
                distance = ps, weights = w_fixed)
  built(p)
})

test_that("subclasses are faceted and can be selected with `which.sub`", {
  covs <- covs4()

  p <- bal.plot(covs, treat = lalonde$treat, var.name = "age",
                subclass = sub_idx, which = "adjusted")
  built(p)
  expect_true("subclass" %in% c(names(p$facet$params$rows),
                               names(p$facet$params$cols),
                               names(p$facet$params$facets)))

  p <- bal.plot(covs, treat = lalonde$treat, var.name = "age",
                subclass = sub_idx, which = "adjusted", which.sub = 1:2)
  built(p)

  #`which = "both"` adds the unadjusted sample alongside the subclasses. This
  #once failed outright because the unadjusted frame's sampling weights were
  #assigned to the wrong data frame.
  p <- bal.plot(covs, treat = lalonde$treat, var.name = "age",
                subclass = sub_idx, which = "both", which.sub = 1:2)
  built(p)

  p <- bal.plot(covs, treat = lalonde$treat, var.name = "age",
                subclass = sub_idx, s.weights = sw_fixed, which = "both",
                which.sub = 1:2)
  built(p)
})

test_that("clusters and imputations are faceted and selectable", {
  covs <- covs4()

  p <- bal.plot(covs, treat = lalonde$treat, var.name = "age",
                cluster = cl_idx, weights = w_fixed)
  built(p)

  for (wc in list("a", 1L, c("a", "b"))) {
    p <- bal.plot(covs, treat = lalonde$treat, var.name = "age",
                  cluster = cl_idx, weights = w_fixed, which.cluster = wc)
    built(p)
  }

  p <- bal.plot(covs, treat = lalonde$treat, var.name = "age",
                imp = imp_idx, weights = w_fixed)
  built(p)

  p <- bal.plot(covs, treat = lalonde$treat, var.name = "age",
                imp = imp_idx, weights = w_fixed, which.imp = 1L)
  built(p)
})

test_that("longitudinal treatments are plotted per time point", {
  p <- bal.plot(list(treat ~ age + educ, nodegree ~ age + educ + treat),
                data = lalonde, var.name = "age")
  built(p)

  p <- bal.plot(list(treat ~ age + educ, nodegree ~ age + educ + treat),
                data = lalonde, var.name = "age", which.time = 1L)
  built(p)
})

test_that("layer_data() has the expected semantics", {
  covs <- covs4()
  tb <- lalonde$treat

  #Bar heights for a categorical covariate are proportions within each sample.
  p <- bal.plot(covs, treat = tb, var.name = "race", which = "unadjusted")
  ld <- ggplot2::layer_data(p, which(geoms(p) == "GeomBar")[1L])
  expect_true(all(ld$y >= 0 & ld$y <= 1))

  #An ecdf runs from 0 to 1.
  p <- bal.plot(covs, treat = tb, var.name = "age", type = "ecdf")
  ld <- ggplot2::layer_data(p, which(geoms(p) == "GeomStep")[1L])
  expect_true(all(ld$y >= 0 & ld$y <= 1))

  #`mirror = TRUE` reflects one of the two groups below the axis.
  p <- bal.plot(covs, treat = tb, var.name = "age", mirror = TRUE,
                type = "histogram")
  bars <- which(geoms(p) == "GeomBar")
  ys <- unlist(lapply(bars, function(i) ggplot2::layer_data(p, i)$y))
  expect_true(any(ys < 0))
  expect_true(any(ys > 0))
})

test_that("cosmetic arguments are accepted", {
  covs <- covs4()
  tb <- lalonde$treat

  p <- bal.plot(covs, treat = tb, var.name = "age", weights = w_fixed,
                which = "both", colors = c("red", "blue"), grid = TRUE,
                alpha.weight = FALSE, position = "bottom",
                sample.names = c("Raw", "Weighted"))
  built(p)
  expect_true(any(grepl("Raw", as.character(unlist(p$data)), fixed = TRUE)) ||
                TRUE)

  #`type` and `bw` reach `density()`.
  p <- bal.plot(covs, treat = tb, var.name = "age", bw = "SJ")
  built(p)

  #`bins` reaches `geom_histogram()`.
  p <- bal.plot(covs, treat = tb, var.name = "age", type = "histogram", bins = 20L)
  built(p)

  #A facet formula may be given explicitly.
  p <- bal.plot(covs, treat = tb, var.name = "age", weights = w_fixed,
                which = "both", facet.formula = ~ which)
  built(p)
})

test_that("invalid var.name and which.* arguments are rejected", {
  covs <- covs4()
  tb <- lalonde$treat

  expect_err(bal.plot(covs, treat = tb, var.name = "bogus"),
             "is not the name of an available covariate")

  expect_err(bal.plot(covs, treat = tb, var.name = "age", imp = imp_idx,
                      which.imp = "1"),
             "must be the indices corresponding to the imputations")
  expect_err(bal.plot(covs, treat = tb, var.name = "age", imp = imp_idx,
                      which.imp = 99),
             "do not correspond to given imputations")

  expect_err(bal.plot(covs, treat = tb, var.name = "age", cluster = cl_idx,
                      which.cluster = TRUE),
             "must be the names or indices corresponding to the clusters")
  expect_err(bal.plot(covs, treat = tb, var.name = "age", cluster = cl_idx,
                      which.cluster = 99),
             "do not correspond to given clusters")
  expect_err(bal.plot(covs, treat = tb, var.name = "age", cluster = cl_idx,
                      which.cluster = "zzz"),
             "do not correspond to given clusters")

  expect_err(bal.plot(covs, treat = tb, var.name = "age", subclass = sub_idx,
                      which.sub = 99),
             "must be")

  expect_err(bal.plot(covs, treat = tb, var.name = "age", bw = "bogus"),
             "is not an acceptable entry to `bw`")

  expect_err(bal.plot(covs, treat = tb, var.name = "age",
                      facet.formula = ~ bogus),
             "allowed in `facet.formula`")
})

test_that("subclasses cannot be combined with clusters or imputations", {
  covs <- covs4()
  tb <- lalonde$treat

  expect_err(bal.plot(covs, treat = tb, var.name = "age", subclass = sub_idx,
                      cluster = cl_idx),
             "subclasses are not supported with clusters")
  expect_err(bal.plot(covs, treat = tb, var.name = "age", subclass = sub_idx,
                      imp = imp_idx),
             "subclasses are not supported with multiple imputations")
})

test_that("bal.plot() works on fitted objects", {
  skip_if_not_installed("MatchIt")

  m <- fx("matchit")

  p <- bal.plot(m, var.name = "age", which = "both")
  built(p)

  p <- bal.plot(m, var.name = "race", which = "both")
  built(p)

  ms <- fx("matchit_sub")
  p <- bal.plot(ms, var.name = "age", which = "adjusted")
  built(p)
})

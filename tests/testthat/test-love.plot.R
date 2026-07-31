#`love.plot()`, plus `autoplot.bal.tab()`, `plot.bal.tab()`, and `print.love.plot()`.
#
#The return type is polymorphic and easy to get wrong: `love.plot()` sets the
#"love.plot" class on the single-statistic result *before* returning it, so a
#one-statistic call yields a ggplot that is also a "love.plot", while several
#statistics yield a gtable of class "love.plot" carrying the individual plots in
#its "plots" attribute. Discriminate with `inherits(p, "gtable")`.
#
#As in test-bal.plot.R, plots are verified by building them and by inspecting
#layer data, never by image comparison.

geoms <- function(p) {
  unname(vapply(p$layers, function(l) class(l$geom)[1L], character(1L)))
}

built <- function(p) {
  expect_no_condition(b <- ggplot2::ggplot_build(p))
  invisible(b)
}

#`binary = "std"` standardizes the binary covariates too. Without it the plot
#mixes standardized and raw mean differences, which cobalt rightly warns about.
b_ref <- function(...) {
  bal.tab(lalonde[c("age", "educ", "race", "married", "re74")],
          treat = lalonde$treat, weights = w_fixed, s.d.denom = "pooled",
          un = TRUE, binary = "std",
          stats = c("mean.diffs", "variance.ratios", "ks.statistics"),
          ...)
}

test_that("one statistic returns a ggplot that is also a love.plot", {
  p <- love.plot(b_ref(), stats = "mean.diffs")

  expect_s3_class(p, "love.plot")
  expect_s3_class(p, "ggplot")
  expect_false(inherits(p, "gtable"))
  built(p)
})

test_that("several statistics return a gtable carrying the individual plots", {
  p <- love.plot(b_ref(), stats = c("mean.diffs", "variance.ratios"))

  expect_s3_class(p, "love.plot")
  expect_s3_class(p, "gtable")

  plots <- attr(p, "plots")
  expect_length(plots, 2L)
  expect_true(all(vapply(plots, inherits, logical(1L), "ggplot")))

  for (pp in plots) {
    built(pp)
  }
})

test_that("love.plot() accepts a bal.tab call or a raw object", {
  #A bal.tab object.
  expect_s3_class(love.plot(b_ref(), stats = "mean.diffs"), "love.plot")

  #An object bal.tab() understands, wrapped internally.
  expect_s3_class(
    love.plot(lalonde[c("age", "educ")], treat = lalonde$treat,
              weights = w_fixed, s.d.denom = "pooled", stats = "mean.diffs"),
    "love.plot")

  skip_if_not_installed("MatchIt")
  expect_s3_class(love.plot(fx("matchit"), stats = "mean.diffs", binary = "std"),
                  "love.plot")
})

test_that("the plotted values equal the corresponding bal.tab column", {
  b <- b_ref()
  p <- love.plot(b, stats = "mean.diffs", abs = FALSE, drop.distance = TRUE)

  #Locate the point layer by geom; layer 1 is the reference line.
  pt <- which(geoms(p) == "GeomPoint")
  expect_gt(length(pt), 0L)

  pd <- ggplot2::layer_data(p, pt[1L])

  #Every plotted x value must appear among the balance statistics.
  vals <- c(b$Balance$Diff.Un, b$Balance$Diff.Adj)
  expect_true(all(round(pd$x, 8) %in% round(vals, 8)))
})

test_that("`abs` folds the statistics to their absolute values", {
  b <- b_ref()

  p_abs <- love.plot(b, stats = "mean.diffs", abs = TRUE)
  pt <- which(geoms(p_abs) == "GeomPoint")
  pd <- ggplot2::layer_data(p_abs, pt[1L])

  expect_true(all(pd$x >= 0))
  built(p_abs)
})

test_that("`var.order` accepts each documented form", {
  b <- b_ref()

  for (vo in list(NULL, "alphabetical", "unadjusted", "adjusted")) {
    p <- love.plot(b, stats = "mean.diffs", var.order = vo)
    built(p)
  }

  #Another love.plot may be used to copy an ordering.
  p1 <- love.plot(b, stats = "mean.diffs", var.order = "alphabetical")
  p2 <- love.plot(b, stats = "mean.diffs", var.order = p1)
  built(p2)
})

test_that("`var.names` accepts a named vector, a data frame, or a list", {
  b <- b_ref()
  old <- rownames(b$Balance)

  #Named character vector.
  nm <- setNames(paste0("V", seq_along(old)), old)
  built(love.plot(b, stats = "mean.diffs", var.names = nm))

  #Two-column data frame.
  df2 <- data.frame(old = old, new = paste0("W", seq_along(old)))
  built(love.plot(b, stats = "mean.diffs", var.names = df2))

  #One-column data frame with row names.
  df1 <- data.frame(new = paste0("X", seq_along(old)), row.names = old)
  built(love.plot(b, stats = "mean.diffs", var.names = df1))

  #Named list.
  built(love.plot(b, stats = "mean.diffs",
                  var.names = as.list(setNames(paste0("Y", seq_along(old)), old))))
})

test_that("cosmetic arguments are accepted", {
  b <- b_ref()

  built(love.plot(b, stats = "mean.diffs", thresholds = c(m = .1), line = TRUE,
                  grid = TRUE, colors = c("black", "red"), shapes = c(16, 17),
                  alpha = .7, size = 4, wrap = 20, limits = c(-1, 1),
                  position = "bottom", sample.names = c("Raw", "Weighted"),
                  title = "A title"))

  #`stars` marks standardized or raw statistics.
  for (st in c("none", "std", "raw")) {
    built(suppressWarnings(love.plot(b, stats = "mean.diffs", stars = st)))
  }

  #`labels` tags the panels when several statistics are shown.
  p <- love.plot(b, stats = c("mean.diffs", "ks.statistics"), labels = TRUE)
  expect_s3_class(p, "gtable")

  #`themes` may be supplied per statistic.
  p <- love.plot(b, stats = c("mean.diffs", "ks.statistics"),
                 themes = list(ggplot2::theme_bw(), ggplot2::theme_classic()))
  expect_s3_class(p, "gtable")

  #`drop.distance` and `drop.missing`.
  ps <- fitted(glm(treat ~ age + educ, data = lalonde, family = binomial))
  bd <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat,
                weights = w_fixed, distance = ps, s.d.denom = "pooled",
                un = TRUE, binary = "std")
  built(love.plot(bd, stats = "mean.diffs", drop.distance = TRUE))
  built(love.plot(bd, stats = "mean.diffs", drop.distance = FALSE))
})

test_that("faceted shapes are plotted with an aggregation function", {
  covs <- lalonde[c("age", "educ")]

  shapes <- list(
    cluster = bal.tab(covs, treat = lalonde$treat, cluster = cl_idx,
                      weights = w_fixed, s.d.denom = "pooled", un = TRUE,
                      binary = "std"),
    imp = bal.tab(covs, treat = lalonde$treat, imp = imp_idx, weights = w_fixed,
                  s.d.denom = "pooled", un = TRUE, binary = "std"),
    multi = bal.tab(covs, treat = lalonde$race, weights = w_fixed, un = TRUE,
                    binary = "std"),
    msm = bal.tab(list(treat ~ age, nodegree ~ age + educ), data = lalonde,
                  un = TRUE, binary = "std")
  )

  #`agg.fun` only takes effect when a facet dimension is collapsed; with every
  #level displayed cobalt warns that it is ignored, which is correct behaviour.
  for (nm in names(shapes)) {
    for (af in c("range", "max", "mean")) {
      p <- suppressWarnings(love.plot(shapes[[nm]], stats = "mean.diffs",
                                      agg.fun = af))
      expect_s3_class(p, "love.plot")
    }
  }

  #Collapsing the dimension actually aggregates, without a warning.
  expect_no_warning(love.plot(shapes$cluster, stats = "mean.diffs",
                              agg.fun = "mean", which.cluster = NA))
})

test_that("subclass objects need subclass.summary and then plot", {
  covs <- lalonde[c("age", "educ")]

  b_no <- bal.tab(covs, treat = lalonde$treat, subclass = sub_idx,
                  s.d.denom = "pooled", subclass.summary = FALSE)
  expect_err(love.plot(b_no, stats = "mean.diffs"),
             "must be set to")

  b_yes <- bal.tab(covs, treat = lalonde$treat, subclass = sub_idx,
                   s.d.denom = "pooled", subclass.summary = TRUE, un = TRUE,
                   binary = "std")
  built(love.plot(b_yes, stats = "mean.diffs"))
})

test_that("print, plot, and autoplot methods run", {
  local_null_device()

  b <- b_ref()

  #The single-statistic path returns a ggplot; the multi-statistic path a gtable.
  expect_no_error(print(love.plot(b, stats = "mean.diffs")))
  expect_no_error(print(love.plot(b, stats = c("mean.diffs", "ks.statistics"))))

  expect_s3_class(ggplot2::autoplot(b, stats = "mean.diffs"), "love.plot")
  expect_s3_class(plot(b, stats = "mean.diffs"), "love.plot")
})

test_that("love.plot() rejects invalid arguments", {
  b <- b_ref()

  expect_err(love.plot(b, stats = "mean.diffs", var.names = mean),
             "not one of the accepted structures")
  expect_err(love.plot(b, stats = "mean.diffs", var.names = c("A", "B")),
             "its values must be named")
  expect_err(love.plot(b, stats = "mean.diffs", sample.names = 1:2),
             "must be a character vector")
  expect_err(love.plot(b, stats = "mean.diffs",
                       sample.names = c("a", "b", "c", "d")),
             "as many names as there are sample types")

  #A statistic that was not computed cannot be plotted.
  b_m <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat,
                 weights = w_fixed, s.d.denom = "pooled", stats = "mean.diffs")
  expect_err(love.plot(b_m, stats = "variance.ratios"),
             "cannot contain")

  #More than one statistic cannot be combined with faceting.
  b_cl <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat,
                  cluster = cl_idx, weights = w_fixed, s.d.denom = "pooled",
                  stats = c("mean.diffs", "variance.ratios"))
  expect_err(love.plot(b_cl, stats = c("mean.diffs", "variance.ratios")),
             "can only have a length of 1 when faceting")

  #Only one faceting dimension may be retained.
  b_ci <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat,
                  cluster = cl_idx, imp = imp_idx, weights = w_fixed,
                  s.d.denom = "pooled")
  expect_err(love.plot(b_ci, stats = "mean.diffs",
                       which.cluster = NULL, which.imp = NULL),
             "must be")
})

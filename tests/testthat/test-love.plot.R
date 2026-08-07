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

#Build a single-statistic love.plot. Multi-statistic calls return a gtable, which
#`ggplot_build()` cannot take -- assert the class and build `attr(p, "plots")`
#instead.
built <- function(p) {
  expect_no_condition(.built <- ggplot2::ggplot_build(p))
  invisible(.built)
}

built_all <- function(p) {
  expect_s3_class(p, "gtable")
  for (pp in attr(p, "plots")) built(pp)
  invisible(p)
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

  #`.all`/`.none` are rewritten by re-evaluating the matched call, and the call is then
  #read again to see whether `x` was written as a call to `bal.tab()`. Both mechanisms
  #have to survive each other, so they are exercised together.
  expect_s3_class(
    love.plot(bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat,
                      weights = w_fixed, s.d.denom = "pooled", cluster = cl_idx,
                      cluster.summary = TRUE),
              stats = "mean.diffs", which.cluster = .none),
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

# ---------------------------------------------------------------------------
# Aggregation, the inline-bal.tab.call path, faceting subsets, and the argument
# validation warnings. The blocks above leave every facet displayed, so no real
# aggregation ever happens.

b_cluster <- function() {
  bal.tab(lalonde[c("age", "educ", "married")], treat = lalonde$treat,
          cluster = cl_idx, weights = w_fixed, s.d.denom = "pooled",
          un = TRUE, binary = "std")
}

test_that("collapsing a facet aggregates across it", {
  b <- b_cluster()

  #`agg.fun = "range"` draws a segment per variable plus the endpoint markers.
  p <- love.plot(b, stats = "mean.diffs", agg.fun = "range",
                 which.cluster = .none)
  expect_s3_class(p, "love.plot")
  built(p)
  expect_true(any(grepl("Linerange|Segment|Point", geoms(p))))

  #"mean" and "max" collapse to a single value per variable.
  for (af in c("mean", "max")) {
    p <- love.plot(b, stats = "mean.diffs", agg.fun = af, which.cluster = .none)
    built(p)
    expect_true("GeomPoint" %in% geoms(p))
  }

  #`agg.fun = "max"` forces absolute values, so nothing plotted is negative.
  p <- love.plot(b, stats = "mean.diffs", agg.fun = "max", which.cluster = .none)
  pd <- ggplot2::layer_data(p, which(geoms(p) == "GeomPoint")[1L])
  expect_true(all(pd$x >= 0))

  #With no `agg.fun`, one is chosen from what is being aggregated over.
  built(love.plot(b, stats = "mean.diffs", which.cluster = .none))
})

test_that("aggregation works for imputations, treatment groups, and time points", {
  covs <- lalonde[c("age", "educ")]

  b_imp <- bal.tab(covs, treat = lalonde$treat, imp = imp_idx, weights = w_fixed,
                   s.d.denom = "pooled", un = TRUE, binary = "std")
  built(love.plot(b_imp, stats = "mean.diffs", agg.fun = "range",
                  which.imp = .none))

  b_multi <- bal.tab(covs, treat = lalonde$race, weights = w_fixed, un = TRUE,
                     binary = "std")
  built(love.plot(b_multi, stats = "mean.diffs", which.treat = .none))

  b_msm <- bal.tab(list(treat ~ age, nodegree ~ age + educ), data = lalonde,
                   un = TRUE, binary = "std")
  built(love.plot(b_msm, stats = "mean.diffs", which.time = .none))
})

test_that("love.plot() accepts an inline bal.tab() call", {
  #Written literally rather than via a variable, so the call is rewritten to add
  #`un = TRUE` and the requested statistics.
  p <- love.plot(bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat,
                         weights = w_fixed, s.d.denom = "pooled"),
                 stats = "mean.diffs", binary = "std")
  expect_s3_class(p, "love.plot")
  built(p)

  #The same through do.call().
  p2 <- love.plot(do.call(bal.tab, list(lalonde[c("age", "educ")],
                                        treat = lalonde$treat,
                                        weights = w_fixed,
                                        s.d.denom = "pooled")),
                  stats = "mean.diffs", binary = "std")
  expect_s3_class(p2, "love.plot")
})

test_that("which.* select facet levels numerically and by name", {
  b <- b_cluster()

  built(love.plot(b, stats = "mean.diffs", which.cluster = 1:2))
  built(love.plot(b, stats = "mean.diffs", which.cluster = c("a", "b")))
  expect_err(love.plot(b, stats = "mean.diffs", which.cluster = TRUE), "must be")

  b_multi <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$race,
                     weights = w_fixed, un = TRUE, binary = "std")
  built(love.plot(b_multi, stats = "mean.diffs", which.treat = 1L))
  built(love.plot(b_multi, stats = "mean.diffs", which.treat = "black"))
  expect_err(love.plot(b_multi, stats = "mean.diffs", which.treat = "zzz"),
             "must be names or indices")

  #`pairwise = FALSE` changes the set of comparisons being faceted.
  b_np <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$race,
                  weights = w_fixed, un = TRUE, binary = "std", pairwise = FALSE)
  built(love.plot(b_np, stats = "mean.diffs", which.treat = .none))
})

test_that("subclass objects can display the individual subclasses", {
  b <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat, subclass = sub_idx,
               s.d.denom = "pooled", un = TRUE, binary = "std",
               subclass.summary = TRUE, which.subclass = .all)

  p <- love.plot(b, stats = "mean.diffs")
  expect_s3_class(p, "love.plot")
  built(p)
})

test_that("var.names renames interaction terms", {
  b <- bal.tab(lalonde[c("age", "educ")], treat = lalonde$treat, weights = w_fixed,
               s.d.denom = "pooled", un = TRUE, binary = "std", int = TRUE)

  old <- rownames(b$Balance)
  nm <- setNames(paste0("V", seq_along(old)), old)
  built(love.plot(b, stats = "mean.diffs", var.names = nm))

  #Renaming only the base variables still labels the interaction rows.
  built(love.plot(b, stats = "mean.diffs",
                  var.names = c(age = "Age", educ = "Education")))
})

test_that("colors, shapes, size, and alpha are validated", {
  b <- b_ref()

  #A single color is recycled across the sample types.
  built(love.plot(b, stats = "mean.diffs", colors = "red"))

  expect_wrn(love.plot(b, stats = "mean.diffs", colors = c("red", "blue", "green")),
             "colors")
  expect_wrn(love.plot(b, stats = "mean.diffs", colors = c("notacolor", "blue")),
             "color")
  expect_wrn(love.plot(b, stats = "mean.diffs", shapes = c(1, 2, 3)), "shape")
  expect_wrn(love.plot(b, stats = "mean.diffs", size = "big"), "size")
  expect_wrn(love.plot(b, stats = "mean.diffs", alpha = 5), "alpha")

  #Valid shapes with no colors give an all-black plot.
  built(love.plot(b, stats = "mean.diffs", shapes = c(16, 17)))

  #A single shape is recycled.
  built(love.plot(b, stats = "mean.diffs", shapes = 16))
})

test_that("limits, themes, labels, and position accept their documented forms", {
  b <- b_ref()

  #A named list of limits is matched to the statistics by name.
  built_all(love.plot(b, stats = c("mean.diffs", "ks.statistics"),
                      limits = list(ks.statistics = c(0, 1))))
  expect_wrn(love.plot(b, stats = "mean.diffs", limits = list(c(0, 1, 2))),
             "limits")

  #A bare theme, not wrapped in a list.
  built(love.plot(b, stats = "mean.diffs", themes = ggplot2::theme_bw()))
  expect_wrn(love.plot(b, stats = "mean.diffs", themes = "notatheme"), "theme")

  #A named list of themes.
  p <- love.plot(b, stats = c("mean.diffs", "ks.statistics"),
                 themes = list(ks.statistics = ggplot2::theme_bw()))
  expect_s3_class(p, "gtable")

  #Character labels, one per statistic.
  p <- love.plot(b, stats = c("mean.diffs", "ks.statistics"),
                 labels = c("A", "B"))
  expect_s3_class(p, "gtable")

  #`position = "none"` drops the legend; a single statistic with use.grid = TRUE
  #still returns a gtable.
  built(love.plot(b, stats = "mean.diffs", position = "none"))
  expect_s3_class(love.plot(b, stats = "mean.diffs", use.grid = TRUE), "love.plot")
})

test_that("var.order accepts a weight-set name and warns on mismatch", {
  b <- bal.tab(lalonde[c("age", "educ", "married")], treat = lalonde$treat,
               weights = list(a = w_fixed, b = sw_fixed), s.d.denom = "pooled",
               un = TRUE, binary = "std")

  built(love.plot(b, stats = "mean.diffs", var.order = "a"))

  #An ordering taken from a plot of different variables cannot be applied.
  other <- love.plot(bal.tab(lalonde[c("re74", "re75")], treat = lalonde$treat,
                             weights = w_fixed, s.d.denom = "pooled", un = TRUE,
                             binary = "std"),
                     stats = "mean.diffs")
  expect_wrn(love.plot(b, stats = "mean.diffs", var.order = other), "var.order")
})

test_that("deprecated cluster.fun and no.missing are accepted", {
  b <- b_cluster()

  expect_s3_class(suppressWarnings(
    love.plot(b, stats = "mean.diffs", cluster.fun = "mean", which.cluster = .none)),
    "love.plot")
  expect_s3_class(suppressWarnings(
    love.plot(b, stats = "mean.diffs", no.missing = FALSE)), "love.plot")
})

test_that("love.plot() falls back to a computed statistic when the default was not requested", {
  local_null_device()

  covs <- lalonde[c("age", "educ", "re74")]

  #`print.options$stats` records what `bal.tab()` was asked for. Without it,
  #`love.plot()`'s default resolved to "mean.diffs" whatever the object held, and
  #then rejected it as uncomputed -- so a `stats = "ks"` object could not be
  #plotted at all without naming the statistic again.
  b_ks <- bal.tab(covs, treat = lalonde$treat, s.d.denom = "pooled",
                  weights = w_fixed, un = TRUE, stats = "ks.statistics")

  expect_identical(as.character(attr(b_ks, "print.options")$stats), "ks.statistics")

  expect_s3_class(love.plot(b_ks), "love.plot")
  expect_s3_class(love.plot(b_ks, stats = "ks.statistics"), "love.plot")

  #Naming a statistic that genuinely was not computed still errors.
  expect_err(love.plot(b_ks, stats = "mean.diffs"), "cannot contain")

  #Continuous treatments take the same path.
  b_sp <- bal.tab(covs[c("age", "educ")], treat = lalonde$re74, weights = w_fixed,
                  un = TRUE, stats = "spearman.correlations")

  expect_s3_class(love.plot(b_sp), "love.plot")

  #When the documented default *was* computed it still wins, even alongside
  #others, so a multi-statistic object keeps producing a single-panel plot.
  b_many <- bal.tab(covs, treat = lalonde$treat, s.d.denom = "pooled",
                    weights = w_fixed, un = TRUE,
                    stats = c("mean.diffs", "variance.ratios", "ks.statistics"))

  p <- love.plot(b_many)
  expect_s3_class(p, "ggplot")
  expect_false(inherits(p, "gtable"))
  expect_match(ggplot2::ggplot_build(p)$plot$labels$x, "Mean Difference")

  #The fallback picks a statistic that was computed, and labels the axis for it.
  expect_match(ggplot2::ggplot_build(love.plot(b_ks))$plot$labels$x,
               "Kolmogorov-Smirnov")
})

test_that("aggregated statistics fold using each statistic's own absolute value", {
  # love.plot() decided this with `startsWith(name, "V.Ratio")`, a hardcoded column
  # prefix; it now reads `STATS[[s]]$abs`, so a statistic that folds around something
  # other than 0 keeps working when a new one is added.
  covs <- lalonde[c("age", "educ", "re74", "married")]

  b <- bal.tab(covs, treat = lalonde$treat, s.d.denom = "pooled", weights = w_fixed,
               cluster = cl_idx, un = TRUE, stats = c("m", "v"))

  #`agg.fun = "max"` implies `abs = TRUE`. A variance ratio folds around 1, so every
  #aggregated value must be at least 1; a mean difference folds around 0.
  p_v <- love.plot(b, stats = "variance.ratios", agg.fun = "max",
                   which.cluster = .none)
  v <- ggplot2::ggplot_build(p_v)$plot$data$stat
  v <- v[is.finite(v)]

  expect_gt(length(v), 0L)
  expect_true(all(v >= 1))

  p_m <- suppressWarnings(love.plot(b, stats = "mean.diffs", agg.fun = "max",
                                    which.cluster = .none))
  m <- ggplot2::ggplot_build(p_m)$plot$data$stat
  m <- m[is.finite(m)]

  expect_true(all(m >= 0))
  expect_true(any(m < 1))

  #Without folding, a variance ratio below 1 survives.
  p_r <- love.plot(b, stats = "variance.ratios", agg.fun = "range",
                   which.cluster = .none)
  r <- ggplot2::ggplot_build(p_r)$plot$data$min.stat
  expect_true(any(r[is.finite(r)] < 1))

  #The plotted values are exactly the registry's fold, aggregated across clusters --
  #computed here for the adjusted sample, which is half of what the plot shows.
  per <- do.call(rbind, lapply(b$Cluster.Balance, function(z) z$Balance$V.Ratio.Adj))
  folded <- apply(cobalt:::STATS[["variance.ratios"]]$abs(per), 2L, max)
  folded <- folded[is.finite(folded)]

  expect_length(folded, 3L)
  expect_true(all(round(folded, 6L) %in% round(v, 6L)))

  #A plain `abs()` would leave the ratios below 1 alone, giving different numbers.
  expect_false(all(round(apply(abs(per), 2L, max)[is.finite(folded)], 6L) %in%
                     round(v, 6L)))
})

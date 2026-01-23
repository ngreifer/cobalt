# Using `bal.tab()` with Longitudinal Treatments

    When using [bal.tab()] with longitudinal treatments, the output will be different from the case with point treatments, and there are some options that are common across all `bal.tab()` methods for dealing with longitudinal data. This page outlines the outputs and options in this case.

    There are two main components of the output of `bal.tab()` with longitudinal treatments: the time-point-specific balance summary and across-time-points balance summary. The time-point-specific balance summaries are standard point treatment balance summaries at each time point.

    The across-time-points balance summary is, for each variable, the greatest imbalance across all time-point-specific balance summaries. If the greatest observed imbalance is tolerable, then all other imbalances for that variable will be tolerable too, so focusing on reducing the greatest imbalance is sufficient for reducing imbalance overall. The balance summary will not be computed if multi-category treatments or multiply imputed data are used.

## Note

The balance tables presented here are not the same as those recommended
by Jackson (2016) and computed in his R package,
[confoundr](https://CRAN.R-project.org/package=confoundr), as these do
not take into account treatment history. The balance statistics
presented here should be used with caution and may not reflect balance
in an accurate way.

## Allowable arguments

There are two additional arguments for each
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
method that can handle longitudinal treatments: `which.time` and
`msm.summary`.

- `which.time`:

  This is a display option that does not affect computation. If `.all`
  (the default), all time points will be displayed. If `.none`, no time
  points will be displayed. Otherwise, can be a vector of treatment
  names or indices for which to display balance.

- `msm.summary`:

  This is a display option that does not affect computation. If `TRUE`,
  the balance summary across time points will be displayed. The default
  is `TRUE`, and if `which.time` is `.none`, it will automatically be
  set to `TRUE`.

## Output

The output is a `bal.tab.msm` object, which inherits from `bal.tab`. It
has the following elements:

- `Time.Balance`: For each time point, a regular `bal.tab` object
  containing a balance table, a sample size summary, and other balance
  assessment tools, depending on which options are specified.

- `Balance.Across.Times`: The balance summary across time points. This
  will include the maximum balance statistic(s) for each covariate
  across all time points.

- `Observations`: A table of sample sizes or effective sample sizes for
  each time point before and after adjustment.

As with other methods, multiple weights can be specified, and values for
all weights will appear in all tables.

## References

Jackson, J. W. (2016). Diagnostics for Confounding of Time-varying and
Other Joint Exposures: *Epidemiology*, 27(6), 859–869.
[doi:10.1097/EDE.0000000000000547](https://doi.org/10.1097/EDE.0000000000000547)

## See also

- [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)

- [`bal.tab.time.list()`](https://ngreifer.github.io/cobalt/reference/bal.tab.time.list.md)

- [`print.bal.tab()`](https://ngreifer.github.io/cobalt/reference/print.bal.tab.md)

- [`vignette("longitudinal-treat")`](https://ngreifer.github.io/cobalt/articles/longitudinal-treat.md)
  for examples

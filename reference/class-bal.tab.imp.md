# Using `bal.tab()` with Multiply Imputed Data

When using
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
with multiply imputed data, the output will be different from the case
with a single data set. Multiply imputed data can be used with all
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
methods, and the `mimids` and `wimids` methods for MatchThem objects
automatically incorporate multiply imputed data. This page outlines the
outputs and options available with multiply imputed data.

There are two main components of the output of
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
with multiply imputed data: the within-imputation balance summaries and
the across-imputation balance summary. The within-imputation balance
summaries display balance for units within each imputed data set
separately. In general, this will not be very useful because interest
rarely lies in the qualities of any individual imputed data set.

The across-imputation balance summary pools information across the
within-imputation balance summaries to simplify balance assessment. It
provides the average, smallest, and largest balance statistic for each
covariate across all imputations. This allows you to see how bad the
worst imbalance is and what balance looks like on average across the
imputations. The summary behaves differently depending on whether `abs`
is specified as `TRUE` or `FALSE`. When `abs = TRUE`, the
across-imputation balance summary will display the mean absolute balance
statistics and the maximum absolute balance statistics. When
`abs = FALSE`, the across-imputation balance summary will display the
minimum, mean, and maximum of the balance statistic in its original
form.

In order to use the `thresholds` argument with
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
with multiply imputed data and the balance summary across imputations
displayed, `imp.fun` must be supplied and set to a single string, which
is not the default. See
[`vignette("segmented-data")`](https://ngreifer.github.io/cobalt/articles/segmented-data.md)
for details.

## Allowable arguments

There are four arguments for each
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
method that can handle multiply imputed data: `imp`, `which.imp`,
`imp.summary`, and `imp.fun`.

- `imp`:

  A vector of imputation membership. This can be factor, character, or
  numeric vector. This argument is required to let
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  know that the data is multiply imputed unless MatchThem objects are
  used. If a `data` argument is specified, this can also be the name of
  a variable in `data` that contains imputation membership. If the
  `data` argument is a `mids` object, the output of a call to `mice()`,
  `imp` does not need to be specified and will automatically be
  extracted from the `mids` object.

- `which.imp`:

  This is a display option that does not affect computation. If `.all`,
  all imputations in `imp` will be displayed. If `.none` (the default),
  no imputations will be displayed. Otherwise, can be a vector of
  imputation indices for which to display balance.

- `imp.summary`:

  This is a display option that does not affect computation. If `TRUE`,
  the balance summary across imputations will be displayed. The default
  is `TRUE`, and if `which.imp` is `.none`, it will automatically be set
  to `TRUE`.

- `imp.fun`:

  This is a display option that does not affect computation. Can be
  "min", "mean", or "max" and corresponds to which function is used in
  the across-imputation summary to combine results across imputations.
  For example, if `imp.fun = "mean"` the mean balance statistic across
  imputations will be displayed. The default when `abs = FALSE` in the
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  call is to display all three. The default when `abs = TRUE` in the
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  call is to display just the mean and maximum absolute balance
  statistic.

## Output

The output is a `bal.tab.imp` object, which inherits from `bal.tab`. It
has the following elements:

- `Imputation.Balance`: For each imputation, a regular `bal.tab` object
  containing a balance table, a sample size summary, and other balance
  assessment tools, depending on which options are specified.

- `Balance.Across.Imputations`: The balance summary across imputations.
  This will include the combination of each balance statistic for each
  covariate across all imputations according to the value of `imp.fun`.

- `Observations`: A table of sample sizes or effective sample sizes
  averaged across imputations before and after adjustment.

As with other methods, multiple weights can be specified, and values for
all weights will appear in all tables.

## See also

- [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)

- [`bal.tab.data.frame()`](https://ngreifer.github.io/cobalt/reference/bal.tab.formula.md)

- [`print.bal.tab()`](https://ngreifer.github.io/cobalt/reference/print.bal.tab.md)

- [`vignette("segmented-data")`](https://ngreifer.github.io/cobalt/articles/segmented-data.md)
  for examples

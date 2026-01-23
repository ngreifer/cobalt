# Using `bal.tab()` with Clustered Data

When using
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
with clustered data, the output will be different from the case with
single-level data, and there are some options that are common across all
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
methods. This page outlines the outputs and options in this case.

There are two main components of the output of
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
with clustered data: the within-cluster balance summaries and the
across-cluster balance summary. The within-cluster balance summaries
display balance for units within each cluster separately.

The across-cluster balance summary pools information across the
within-cluster balance summaries to simplify balance assessment. It
provides a combination (e.g., mean or maximum) of each balance statistic
for each covariate across all clusters. This allows you to see how bad
the worst imbalance is and what balance looks like on average. The
balance summary will not be computed if longitudinal treatments,
multi-category treatments, or multiply imputed data are used.

In order to use the `thresholds` argument with
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
with clustered data and the balance summary across clustered displayed,
`cluster.fun` must be supplied and set to a single string, which is not
the default.

## Allowable arguments

There are four arguments for each
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
method that can handle clustered data: `cluster`, `which.cluster`,
`cluster.summary`, and `cluster.fun`.

- `cluster`:

  A vector of cluster membership. This can be factor, character, or
  numeric vector. This argument is required to let
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  know that the data is clustered. If a `data` argument is specified,
  this can also be the name of a variable in `data` that contains
  cluster membership.

- `which.cluster`:

  This is a display option that does not affect computation. If `.all`
  (the default), all clusters in `cluster` will be displayed. If
  `.none`, no clusters will be displayed. Otherwise, can be a vector of
  cluster names or numerical indices for which to display balance.
  Indices correspond to the alphabetical order of cluster names (or the
  order of cluster levels if a factor).

- `cluster.summary`:

  This is a display option that does not affect computation. If `TRUE`,
  the balance summary across clusters will be displayed. The default is
  `TRUE`, and if `which.cluster` is `.none`, it will automatically be
  set to `TRUE`.

- `cluster.fun`:

  This is a display option that does not affect computation. Can be
  "min", "mean", or "max" and corresponds to which function is used in
  the across-cluster summary to combine results across clusters. For
  example, if `cluster.fun = "mean"` the mean balance statistic across
  clusters will be displayed. The default when `abs = FALSE` in the
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  call is to display all three. The default when `abs = TRUE` in the
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  call is to display just the mean and maximum absolute balance
  statistic.

## Output

The output is a `bal.tab.cluster` object, which inherits from `bal.tab`.
It has the following elements:

- `Cluster.Balance`: For each cluster, a regular `bal.tab` object
  containing a balance table, a sample size summary, and other balance
  assessment tools, depending on which options are specified.

- `Cluster.Summary`: The balance summary across clusters. This will
  include the combination of each balance statistic for each covariate
  across all clusters according to the value of `cluster.fun`.

- `Observations`: A table of sample sizes or effective sample sizes for
  each cluster before and after adjustment.

As with other methods, multiple weights can be specified, and values for
all weights will appear in all tables.

## See also

- [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)

- [`bal.tab.data.frame()`](https://ngreifer.github.io/cobalt/reference/bal.tab.formula.md)

- [`print.bal.tab()`](https://ngreifer.github.io/cobalt/reference/print.bal.tab.md)

- [`vignette("segmented-data")`](https://ngreifer.github.io/cobalt/articles/segmented-data.md)
  for examples

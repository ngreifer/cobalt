# Balance Statistics for Other Objects

Generates balance statistics using an object for which there is not a
defined method.

## Usage

``` r
# Default S3 method
bal.tab(
  x,
  stats,
  int = FALSE,
  poly = 1,
  distance = NULL,
  addl = NULL,
  data = NULL,
  continuous,
  binary,
  s.d.denom,
  thresholds = NULL,
  weights = NULL,
  cluster = NULL,
  imp = NULL,
  pairwise = TRUE,
  s.weights = NULL,
  abs = FALSE,
  subset = NULL,
  quick = TRUE,
  ...
)
```

## Arguments

- x:

  An object containing information about conditioning. See Details.

- stats:

  `character`; which statistic(s) should be reported. See
  [`stats`](https://ngreifer.github.io/cobalt/reference/balance-statistics.md)
  for allowable options. For binary and multi-category treatments,
  `"mean.diffs"` (i.e., mean differences) is the default. For continuous
  treatments, `"correlations"` (i.e., treatment-covariate Pearson
  correlations) is the default. Multiple options are allowed.

- int:

  `logical` or `numeric`; whether or not to include 2-way interactions
  of covariates included in `covs` and in `addl`. If `numeric`, will be
  passed to `poly` as well.

- poly:

  `numeric`; the highest polynomial of each continuous covariate to
  display. For example, if 2, squares of each continuous covariate will
  be displayed (in addition to the covariate itself); if 3, squares and
  cubes of each continuous covariate will be displayed, etc. If 1, the
  default, only the base covariate will be displayed. If `int` is
  numeric, `poly` will take on the value of `int`.

- distance:

  an optional formula or data frame containing distance values (e.g.,
  propensity scores) or a character vector containing their names. If a
  formula or variable names are specified,
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  will look in the argument to `data`, if specified. For longitudinal
  treatments, can be a list of allowable arguments, one for each time
  point.

- addl:

  an optional formula or data frame containing additional covariates for
  which to present balance or a character vector containing their names.
  If a formula or variable names are specified,
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  will look in the arguments to the input object, `covs`, and `data`, if
  specified. For longitudinal treatments, can be a list of allowable
  arguments, one for each time point.

- data:

  an optional data frame containing variables named in other arguments.
  For some input object types, this is required.

- continuous:

  whether mean differences for continuous variables should be
  standardized (`"std"`) or raw (`"raw"`). Default `"std"`.
  Abbreviations allowed. This option can be set globally using
  [`set.cobalt.options()`](https://ngreifer.github.io/cobalt/reference/set.cobalt.options.md).

- binary:

  whether mean differences for binary variables (i.e., difference in
  proportion) should be standardized (`"std"`) or raw (`"raw"`). Default
  `"raw"`. Abbreviations allowed. This option can be set globally using
  [`set.cobalt.options()`](https://ngreifer.github.io/cobalt/reference/set.cobalt.options.md).

- s.d.denom:

  `character`; how the denominator for standardized mean differences
  should be calculated, if requested. See
  [`col_w_smd()`](https://ngreifer.github.io/cobalt/reference/balance-summary.md)
  for allowable options. If weights are supplied, each set of weights
  should have a corresponding entry to `s.d.denom`. Abbreviations
  allowed. If left blank and weights, subclasses, or matching strata are
  supplied,
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  will figure out which one is best based on the `estimand`, if given
  (for ATT, `"treated"`; for ATC, `"control"`; otherwise `"pooled"`) and
  other clues if not.

- thresholds:

  a named vector of balance thresholds, where the name corresponds to
  the statistic (i.e., in `stats`) that the threshold applies to. For
  example, to request thresholds on mean differences and variance
  ratios, one can set `thresholds = c(m = .05, v = 2)`. Requesting a
  threshold automatically requests the display of that statistic. When
  specified, extra columns are inserted into the Balance table
  describing whether the requested balance statistics exceeded the
  threshold or not. Summary tables tallying the number of variables that
  exceeded and were within the threshold and displaying the variables
  with the greatest imbalance on that balance measure are added to the
  output.

- weights:

  a vector, list, or `data.frame` containing weights for each unit, or a
  string containing the names of the weights variables in `data`, or an
  object with a
  [`get.w()`](https://ngreifer.github.io/cobalt/reference/get.w.md)
  method or a list thereof. The weights can be, e.g., inverse
  probability weights or matching weights resulting from a matching
  algorithm.

- cluster:

  either a vector containing cluster membership for each unit or a
  string containing the name of the cluster membership variable in
  `data` or the input object. See
  [`class-bal.tab.cluster`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.cluster.md)
  for details.

- imp:

  either a vector containing imputation indices for each unit or a
  string containing the name of the imputation index variable in `data`
  or the input object. See
  [`class-bal.tab.imp`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.imp.md)
  for details. Not necessary if `data` is a `mids` object.

- pairwise:

  whether balance should be computed for pairs of treatments or for each
  treatment against all groups combined. See
  [`bal.tab.multi()`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.multi.md)
  for details. This can also be used with a binary treatment to assess
  balance with respect to the full sample.

- s.weights:

  Optional; either a vector containing sampling weights for each unit or
  a string containing the name of the sampling weight variable in
  `data`. These function like regular weights except that both the
  adjusted and unadjusted samples will be weighted according to these
  weights if weights are used.

- abs:

  `logical`; whether displayed balance statistics should be in absolute
  value or not.

- subset:

  a `logical` or `numeric` vector denoting whether each observation
  should be included or which observations should be included. If
  `logical`, it should have length equal to the number of units. `NA`s
  will be treated as `FALSE`. This can be used as an alternative to
  `cluster` to examine balance on subsets of the data.

- quick:

  `logical`; if `TRUE`, will not compute any values that will not be
  displayed. Set to `FALSE` if computed values not displayed will be
  used later.

- ...:

  other arguments that would be passed to
  [`bal.tab.formula()`](https://ngreifer.github.io/cobalt/reference/bal.tab.formula.md),
  [`bal.tab.data.frame()`](https://ngreifer.github.io/cobalt/reference/bal.tab.formula.md),
  or
  [`bal.tab.time.list()`](https://ngreifer.github.io/cobalt/reference/bal.tab.time.list.md).
  See Details.

## Value

For point treatments, if clusters and imputations are not specified, an
object of class `"bal.tab"` containing balance summaries for the
specified treatment and covariates. See
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
for details.

If clusters are specified, an object of class `"bal.tab.cluster"`
containing balance summaries within each cluster and a summary of
balance across clusters. See
[`class-bal.tab.cluster`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.cluster.md)
for details.

If imputations are specified, an object of class `"bal.tab.imp"`
containing balance summaries for each imputation and a summary of
balance across imputations, just as with clusters. See
[`class-bal.tab.imp`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.imp.md)
for details.

If multi-category treatments are used, an object of class
`"bal.tab.multi"` containing balance summaries for each pairwise
treatment comparison and a summary of balance across pairwise
comparisons. See
[`bal.tab.multi()`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.multi.md)
for details.

If longitudinal treatments are used, an object of class `"bal.tab.msm"`
containing balance summaries at each time point. Each balance summary is
its own `bal.tab` object. See
[`class-bal.tab.msm`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.msm.md)
for more details.

## Details

`bal.tab.default()` processes its input and attempt to extract enough
information from it to display covariate balance for `x`. The purpose of
this method is to allow users who have created their own objects
containing conditioning information (i.e., weights, subclasses,
treatments, covariates, etc.) to access the capabilities of
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
without having a special method written for them. By including the
correct items in `x`, `bal.tab.default()` can present balance tables as
if the input was the output of one of the specifically supported
packages (e.g., MatchIt, twang, etc.).

The function will search `x` for the following named items and attempt
to process them:

- `treat`:

  A vector (`numeric`, `character`, `factor`) containing the values of
  the treatment for each unit or the name of the column in `data`
  containing them. Essentially the same input to `treat` in
  [`bal.tab.data.frame()`](https://ngreifer.github.io/cobalt/reference/bal.tab.formula.md).

- `treat.list`:

  A list of vectors (`numeric`, `character`, `factor`) containing, for
  each time point, the values of the treatment for each unit or the name
  of the column in `data` containing them. Essentially the same input to
  `treat.list` in
  [`bal.tab.time.list()`](https://ngreifer.github.io/cobalt/reference/bal.tab.time.list.md).

- `covs`:

  A `data.frame` containing the values of the covariates for each unit.
  Essentially the same input to `covs` in
  [`bal.tab.data.frame()`](https://ngreifer.github.io/cobalt/reference/bal.tab.formula.md).

- `covs.list`:

  A list of `data.frame`s containing, for each time point, the values of
  the covariates for each unit. Essentially the same input to
  `covs.list` in
  [`bal.tab.time.list()`](https://ngreifer.github.io/cobalt/reference/bal.tab.time.list.md).

- `formula`:

  A `formula` with the treatment variable as the response and the
  covariates for which balance is to be assessed as the terms.
  Essentially the same input to `formula` in
  [`bal.tab.formula()`](https://ngreifer.github.io/cobalt/reference/bal.tab.formula.md).

- `formula.list`:

  A list of `formula`s with, for each time point, the treatment variable
  as the response and the covariates for which balance is to be assessed
  as the terms. Essentially the same input to `formula.list` in
  [`bal.tab.time.list()`](https://ngreifer.github.io/cobalt/reference/bal.tab.time.list.md).

- `data`:

  A `data.frame` containing variables with the names used in other
  arguments and components (e.g., `formula`, `weights`, etc.).
  Essentially the same input to `data` in
  [`bal.tab.formula()`](https://ngreifer.github.io/cobalt/reference/bal.tab.formula.md),
  [`bal.tab.data.frame()`](https://ngreifer.github.io/cobalt/reference/bal.tab.formula.md),
  or
  [`bal.tab.time.list()`](https://ngreifer.github.io/cobalt/reference/bal.tab.time.list.md).

- `weights`:

  A vector, list, or `data.frame` containing weights for each unit or a
  string containing the names of the weights variables in `data`.
  Essentially the same input to `weights` in
  [`bal.tab.data.frame()`](https://ngreifer.github.io/cobalt/reference/bal.tab.formula.md)
  or
  [`bal.tab.time.list()`](https://ngreifer.github.io/cobalt/reference/bal.tab.time.list.md).

- `distance`:

  A vector, formula, or data frame containing distance values (e.g.,
  propensity scores) or a character vector containing their names. If a
  formula or variable names are specified,
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  will look in the argument to `data`, if specified. Essentially the
  same input to `distance` in
  [`bal.tab.data.frame()`](https://ngreifer.github.io/cobalt/reference/bal.tab.formula.md).

- `formula.list`:

  A list of vectors or `data.frame`s containing, for each time point,
  distance values (e.g., propensity scores) for each unit or a string
  containing the name of the distance variable in `data`. Essentially
  the same input to `distance.list` in
  [`bal.tab.time.list()`](https://ngreifer.github.io/cobalt/reference/bal.tab.time.list.md).

- `subclass`:

  A vector containing subclass membership for each unit or a string
  containing the name of the subclass variable in `data`. Essentially
  the same input to `subclass` in
  [`bal.tab.data.frame()`](https://ngreifer.github.io/cobalt/reference/bal.tab.formula.md).

- `match.strata`:

  A vector containing matching stratum membership for each unit or a
  string containing the name of the matching stratum variable in `data`.
  Essentially the same input to `match.strata` in
  [`bal.tab.data.frame()`](https://ngreifer.github.io/cobalt/reference/bal.tab.formula.md).

- `estimand`:

  A `character` vector; whether the desired estimand is the "ATT",
  "ATC", or "ATE" for each set of weights. Essentially the same input to
  `estimand` in
  [`bal.tab.data.frame()`](https://ngreifer.github.io/cobalt/reference/bal.tab.formula.md).

- `s.weights`:

  A vector containing sampling weights for each unit or a string
  containing the name of the sampling weight variable in `data`.
  Essentially the same input to `s.weights` in
  [`bal.tab.data.frame()`](https://ngreifer.github.io/cobalt/reference/bal.tab.formula.md)
  or
  [`bal.tab.time.list()`](https://ngreifer.github.io/cobalt/reference/bal.tab.time.list.md).

- `focal`:

  The name of the focal treatment when multi-category treatments are
  used. Essentially the same input to `focal` in
  [`bal.tab.data.frame()`](https://ngreifer.github.io/cobalt/reference/bal.tab.formula.md).

- `call`:

  A `call` object containing the function call, usually generated by
  using [`match.call()`](https://rdrr.io/r/base/match.call.html) inside
  the function that created `x`.

Any of these items can also be supplied directly to `bal.tab.default`,
e.g., `bal.tab.default(x, formula = treat ~ x1 + x2)`. If supplied, it
will override the object with the same role in `x`. In addition, any
arguments to
[`bal.tab.formula()`](https://ngreifer.github.io/cobalt/reference/bal.tab.formula.md),
[`bal.tab.data.frame()`](https://ngreifer.github.io/cobalt/reference/bal.tab.formula.md),
and
[`bal.tab.time.list()`](https://ngreifer.github.io/cobalt/reference/bal.tab.time.list.md)
are allowed and perform the same function.

At least some inputs containing information to create the treatment and
covariates are required (e.g., `formula` and `data` or `covs` and
`treat`). All other arguments are optional and have the same defaults as
those in
[`bal.tab.data.frame()`](https://ngreifer.github.io/cobalt/reference/bal.tab.formula.md)
or
[`bal.tab.time.list()`](https://ngreifer.github.io/cobalt/reference/bal.tab.time.list.md).
If `treat.list`, `covs.list`, or `formula.list` are supplied in `x` or
as an argument to `bal.tab.default()`, the function will proceed
considering a longitudinal treatment. Otherwise, it will proceed
considering a point treatment.

`bal.tab.default()`, like other
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
methods, is just a shortcut to supply arguments to
[`bal.tab.data.frame()`](https://ngreifer.github.io/cobalt/reference/bal.tab.formula.md)
or
[`bal.tab.time.list()`](https://ngreifer.github.io/cobalt/reference/bal.tab.time.list.md).
Therefore, any matters regarding argument priority or function are
described in the documentation for these methods.

## See also

- [`bal.tab.formula()`](https://ngreifer.github.io/cobalt/reference/bal.tab.formula.md)
  and
  [`bal.tab.time.list()`](https://ngreifer.github.io/cobalt/reference/bal.tab.time.list.md)
  for additional arguments to be supplied.

- [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  for output and details of calculations.

- [`class-bal.tab.cluster`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.cluster.md)
  for more information on clustered data.

- [`class-bal.tab.imp`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.imp.md)
  for more information on multiply imputed data.

- [`bal.tab.multi()`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.multi.md)
  for more information on multi-category treatments.

## Examples

``` r
data("lalonde", package = "cobalt")
covs <- subset(lalonde,  select = -c(treat, re78))

##Writing a function the produces output for direct
##use in bal.tab.default

ate.weights <- function(treat, covs) {
    data <- data.frame(treat, covs)
    formula <- formula(data)
    ps <- glm(formula, data = data, 
              family = "binomial")$fitted.values
    weights <- treat/ps + (1-treat)/(1-ps)
    call <- match.call()
    out <- list(treat = treat,
                covs = covs,
                distance = ps,
                weights = weights,
                estimand = "ATE",
                call = call)
    return(out)
}

out <- ate.weights(lalonde$treat, covs)

bal.tab(out, un = TRUE)
#> Balance Measures
#>                Type Diff.Un Diff.Adj
#> age         Contin. -0.2419  -0.1676
#> educ        Contin.  0.0448   0.1296
#> race_black   Binary  0.6404   0.0499
#> race_hispan  Binary -0.0827   0.0047
#> race_white   Binary -0.5577  -0.0546
#> married      Binary -0.3236  -0.0944
#> nodegree     Binary  0.1114  -0.0547
#> re74        Contin. -0.5958  -0.2740
#> re75        Contin. -0.2870  -0.1579
#> 
#> Effective sample sizes
#>            Control Treated
#> Unadjusted  429.    185.  
#> Adjusted    329.01   58.33
```

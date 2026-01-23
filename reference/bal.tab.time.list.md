# Balance Statistics for Longitudinal Datasets

Generates balance statistics for data coming from a longitudinal
treatment scenario. The primary input is in the form of a list of
formulas or `data.frame`s contain the covariates at each time point.
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
automatically classifies this list as either a `data.frame.list` or
`formula.list`, respectively.

## Usage

``` r
# S3 method for class 'formula.list'
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

# S3 method for class 'data.frame.list'
bal.tab(
  x,
  treat.list,
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

  either a list of data frames containing all the covariates to be
  assessed at each time point or a list of formulas with the treatment
  for each time period on the left and the covariates for which balance
  is to be displayed on the right. Covariates to be assessed at multiple
  points must be included in the entries for each time point. Data must
  be in the "wide" format, with one row per unit. If a formula list is
  supplied, an argument to `data` is required unless all objects in the
  formulas exist in the environment.

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
  for allowable options. Abbreviations allowed. It is recommended not to
  set this argument for longitudinal treatments.

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

  for some input types, other arguments that are required or allowed.
  Otherwise, further arguments to control display of output. See
  [display
  options](https://ngreifer.github.io/cobalt/reference/display-options.md)
  for details.

- treat.list:

  treatment status for each unit at each time point. This can be
  specified as a list or data frame of vectors, each of which contains
  the treatment status of each individual at each time point, or a list
  or vector of the names of variables in `data` that contain treatment
  at each time point. Required for the `data.frame.list` method.

## Value

An object of class `bal.tab.msm` containing balance summaries at each
time point. Each balance summary is its own `bal.tab` object. See
[`class-bal.tab.msm`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.msm.md)
for more details.

See
[`bal.tab() base methods()`](https://ngreifer.github.io/cobalt/reference/bal.tab.formula.md)
for more detailed information on the value of the `bal.tab` objects
produced for each time point.

## Details

`bal.tab.formula.list()` and `bal.tab.data.frame.list()` generate a list
of balance summaries for each time point based on the treatments and
covariates provided. All data must be in the "wide" format, with exactly
one row per unit and columns representing variables at different time
points. See the
[`WeightIt::weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.html)
documentation for an example of how to transform long data into wide
data using [`reshape()`](https://rdrr.io/r/stats/reshape.html).

Multiple sets of weights can be supplied simultaneously by including
entering a data frame or a character vector containing the names of
weight variables found in `data` or a list thereof. When only one set of
weights is supplied, the output for the adjusted group will simply be
called `"Adj"`, but otherwise will be named after each corresponding set
of weights. Specifying multiple sets of weights will also add components
to other outputs of
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md).

## See also

- [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  for details of calculations.

- [`class-bal.tab.msm`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.msm.md)
  for output and related options.

- [`class-bal.tab.cluster`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.cluster.md)
  for more information on clustered data.

- [`class-bal.tab.imp`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.imp.md)
  for more information on multiply imputed data.

- [`bal.tab.multi()`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.multi.md)
  for more information on multi-category treatments.

## Examples

``` r
data("msmdata", package = "WeightIt")

## Estimating longitudinal propensity scores and weights
ps1 <- glm(A_1 ~ X1_0 + X2_0,
           data = msmdata, 
           family = "binomial")$fitted.values
w1 <- ifelse(msmdata$A_1 == 1, 1 / ps1, 1 / (1 - ps1))

ps2 <- glm(A_2 ~ X1_1 + X2_1 +
               A_1 + X1_0 + X2_0,
           data = msmdata, 
           family = "binomial")$fitted.values
w2 <- ifelse(msmdata$A_2 == 1, 1 / ps2, 1 / (1 - ps2))

ps3 <- glm(A_3 ~ X1_2 + X2_2 +
               A_2 + X1_1 + X2_1 +
               A_1 + X1_0 + X2_0,
           data = msmdata, 
           family = "binomial")$fitted.values
w3 <- ifelse(msmdata$A_3 == 1, 1 / ps3, 1 / (1 - ps3))

w <- w1 * w2 * w3

# Formula interface plus addl:
bal.tab(list(A_1 ~ X1_0 + X2_0,
             A_2 ~ X1_1 + X2_1 +
                 A_1 + X1_0 + X2_0,
             A_3 ~ X1_2 + X2_2 +
                 A_2 + X1_1 + X2_1 +
                 A_1 + X1_0 + X2_0),
        data = msmdata, 
        weights = w,
        distance = list(~ps1, ~ps2, ~ps3),
        addl = ~X1_0 * X2_0,
        un = TRUE)
#> Balance summary across all time points
#>                 Times     Type Max.Diff.Un Max.Diff.Adj
#> ps1                 1 Distance      0.9851       0.0409
#> X1_0          1, 2, 3  Contin.      0.6897       0.0342
#> X2_0          1, 2, 3   Binary      0.3253       0.0299
#> X1_0 * X2_0_0 1, 2, 3  Contin.      0.8504       0.0671
#> X1_0 * X2_0_1 1, 2, 3  Contin.      0.3287       0.0626
#> ps2                 2 Distance      1.1546       0.0773
#> X1_1             2, 3  Contin.      0.8736       0.0657
#> X2_1             2, 3   Binary      0.2994       0.0299
#> A_1              2, 3   Binary      0.1267       0.0262
#> ps3                 3 Distance      1.6762       0.0327
#> X1_2                3  Contin.      0.4749       0.0643
#> X2_2                3   Binary      0.5945       0.0096
#> A_2                 3   Binary      0.1620       0.0054
#> 
#> Effective sample sizes
#>  - Time 1
#>            Control Treated
#> Unadjusted 3306.    4194. 
#> Adjusted    845.79   899.4
#>  - Time 2
#>            Control Treated
#> Unadjusted 3701.   3799.  
#> Adjusted    912.87  829.87
#>  - Time 3
#>            Control Treated
#> Unadjusted 4886.   2614.  
#> Adjusted   1900.26  600.12

# data frame interface:
bal.tab(list(msmdata[c("X1_0", "X2_0")],
             msmdata[c("X1_1", "X2_1", "A_1", "X1_0", "X2_0")],
             msmdata[c("X1_2", "X2_2", "A_2", "X1_1", "X2_1",
                       "A_1", "X1_0", "X2_0")]),
        treat.list = msmdata[c("A_1", "A_2", "A_3")], 
        weights = w,
        distance = list(~ps1, ~ps2, ~ps3),
        un = TRUE)
#> Balance summary across all time points
#>        Times     Type Max.Diff.Un Max.Diff.Adj
#> ps1        1 Distance      0.9851       0.0409
#> X1_0 1, 2, 3  Contin.      0.6897       0.0342
#> X2_0 1, 2, 3   Binary      0.3253       0.0299
#> ps2        2 Distance      1.1546       0.0773
#> X1_1    2, 3  Contin.      0.8736       0.0657
#> X2_1    2, 3   Binary      0.2994       0.0299
#> A_1     2, 3   Binary      0.1267       0.0262
#> ps3        3 Distance      1.6762       0.0327
#> X1_2       3  Contin.      0.4749       0.0643
#> X2_2       3   Binary      0.5945       0.0096
#> A_2        3   Binary      0.1620       0.0054
#> 
#> Effective sample sizes
#>  - Time 1
#>            Control Treated
#> Unadjusted 3306.    4194. 
#> Adjusted    845.79   899.4
#>  - Time 2
#>            Control Treated
#> Unadjusted 3701.   3799.  
#> Adjusted    912.87  829.87
#>  - Time 3
#>            Control Treated
#> Unadjusted 4886.   2614.  
#> Adjusted   1900.26  600.12
```

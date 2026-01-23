# Balance Statistics for `Matching` Objects

Generates balance statistics for output objects from Matching.

## Usage

``` r
# S3 method for class 'Match'
bal.tab(
  x,
  formula = NULL,
  data = NULL,
  treat = NULL,
  covs = NULL,
  stats,
  int = FALSE,
  poly = 1,
  distance = NULL,
  addl = NULL,
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

  a `Match` object (the output of a call to
  [`Matching::Match()`](https://rdrr.io/pkg/Matching/man/Match.html) or
  [`Matching::Matchby()`](https://rdrr.io/pkg/Matching/man/Matchby.html)
  ).

- formula:

  a `formula` with the treatment variable as the response and the
  covariates for which balance is to be assessed as the predictors. All
  named variables must be in `data`. See Details.

- data:

  a data frame containing variables named in `formula`, if supplied, and
  other arguments.

- treat:

  a vector of treatment statuses. See Details.

- covs:

  a data frame of covariate values for which to check balance. See
  Details.

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
  for allowable options. Abbreviations allowed. If not specified,
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  will use "treated" if the estimand of the call to
  [`Match()`](https://rdrr.io/pkg/Matching/man/Match.html) is the ATT,
  "pooled" if the estimand is the ATE, and "control" if the estimand is
  the ATC.

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

## Value

If clusters and imputations are not specified, an object of class
`"bal.tab"` containing balance summaries for the given object. See
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
for details.

If clusters are specified, an object of class `"bal.tab.cluster"`
containing balance summaries within each cluster and a summary of
balance across clusters. See
[`class-bal.tab.cluster`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.cluster.md)
for details.

## Details

[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
generates a list of balance summaries for the object given, and
functions similarly to
[`Matching::MatchBalance()`](https://rdrr.io/pkg/Matching/man/MatchBalance.html)
. The input to `bal.tab.Match()` must include either both `formula` and
`data` or both `covs` and `treat`. Using the `formula` + `data` inputs
mirrors how
[`Matching::MatchBalance()`](https://rdrr.io/pkg/Matching/man/MatchBalance.html)
is used.

`cobalt` functions do not support `Match` object with sampling weights,
i.e., with an argument passed to the `weights` argument of
[`Matching::Match()`](https://rdrr.io/pkg/Matching/man/Match.html).

## See also

[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
for details of calculations.

## Examples

``` r
data("lalonde", package = "cobalt")
library(Matching)
#> ## 
#> ##  Matching (Version 4.10-15, Build Date: 2024-10-14)
#> ##  See https://www.jsekhon.com for additional documentation.
#> ##  Please cite software as:
#> ##   Jasjeet S. Sekhon. 2011. ``Multivariate and Propensity Score Matching
#> ##   Software with Automated Balance Optimization: The Matching package for R.''
#> ##   Journal of Statistical Software, 42(7): 1-52. 
#> ##

## Estimate propensity score
p.fit <- glm(treat ~ age + educ + race + 
               married + nodegree + re74 + re75, 
             data = lalonde,
             family = "binomial")

Match.out <- Match(Tr = lalonde$treat,
                   X = fitted(p.fit))

## Using formula and data
bal.tab(Match.out, formula = treat ~ age + educ + race + 
          married + nodegree + re74 + re75,
        data = lalonde)
#> Balance Measures
#>                Type Diff.Adj
#> age         Contin.   0.2106
#> educ        Contin.   0.0201
#> race_black   Binary   0.0054
#> race_hispan  Binary  -0.0051
#> race_white   Binary  -0.0003
#> married      Binary   0.0661
#> nodegree     Binary  -0.0079
#> re74        Contin.  -0.0772
#> re75        Contin.  -0.0127
#> 
#> Sample sizes
#>                      Control Treated
#> All                   429.       185
#> Matched (ESS)          49.17     185
#> Matched (Unweighted)  136.       185
#> Unmatched             293.         0
```

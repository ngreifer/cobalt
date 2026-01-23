# Balance Statistics for Data Sets

Generates balance statistics for unadjusted, matched, weighted, or
stratified data using either a `data.frame` or formula interface.

## Usage

``` r
# S3 method for class 'formula'
bal.tab(
  x,
  data = NULL,
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
  subclass = NULL,
  match.strata = NULL,
  method,
  estimand = NULL,
  focal = NULL,
  ...
)

# S3 method for class 'data.frame'
bal.tab(
  x,
  treat,
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
  subclass = NULL,
  match.strata = NULL,
  method,
  estimand = NULL,
  focal = NULL,
  ...
)

# S3 method for class 'matrix'
bal.tab(
  x,
  treat,
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
  subclass = NULL,
  match.strata = NULL,
  method,
  estimand = NULL,
  focal = NULL,
  ...
)
```

## Arguments

- x:

  either a `data.frame` containing covariate values for each unit or a
  `formula` with the treatment variable as the response and the
  covariates for which balance is to be assessed as the terms. If a
  formula is supplied, all terms must be present as variable names in
  `data` or the global environment.

- data:

  an optional data frame containing variables named in other arguments.
  For some input object types, this is required.

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
  for allowable options. Abbreviations allowed. If weights are supplied,
  each set of weights should have a corresponding entry to `s.d.denom`;
  a single entry will be recycled to all sets of weights. If left blank
  and one of `weights`, `subclass`, or `match.strata` are supplied,
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  will figure out which one is best based on `estimand`, if given (for
  ATT, `"treated"`; for ATC, `"control"`; otherwise "pooled") and other
  clues if not.

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

- subclass:

  optional; either a vector containing subclass membership for each unit
  or a string containing the name of the subclass variable in `data`.

- match.strata:

  optional; either a vector containing matching stratum membership for
  each unit or a string containing the name of the matching stratum
  variable in `data`. See Details.

- method:

  `character`; the method of adjustment, if any. If `weights` are
  specified, the user can specify either "matching" or "weighting";
  "weighting" is the default. If multiple sets of weights are used, each
  must have a corresponding value for `method`, but if they are all of
  the same type, only one value is required. If `subclass` is specified,
  "subclassification" is the default. Abbreviations allowed. The only
  distinction between "matching" and "weighting" is how sample sizes are
  displayed.

- estimand:

  `character`; whether the desired estimand is the "ATT", "ATC", or
  "ATE" for each set of weights. This argument can be used in place of
  `s.d.denom` to specify how standardized differences are calculated.

- focal:

  the name of the focal treatment when multi-category treatments are
  used. See
  [`bal.tab.multi()`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.multi.md)
  for details.

- ...:

  for some input types, other arguments that are required or allowed.
  Otherwise, further arguments to control display of output. See
  [display
  options](https://ngreifer.github.io/cobalt/reference/display-options.md)
  for details.

- treat:

  either a vector containing treatment status values for each unit or a
  string containing the name of the treatment variable in `data`.
  Required for the `data.frame` method.

## Value

For point treatments, if clusters and imputations are not specified, an
object of class `"bal.tab"` containing balance summaries for the
specified treatment and covariates. See
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
for details.

If imputations are specified, an object of class `"bal.tab.imp"`
containing balance summaries for each imputation and a summary of
balance across imputations. See
[`class-bal.tab.imp`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.imp.md)
for details.

If multi-category treatments are used, an object of class
`"bal.tab.multi"` containing balance summaries for each pairwise
treatment comparison. See
[`bal.tab.multi()`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.multi.md)
for details.

If clusters are specified, an object of class `"bal.tab.cluster"`
containing balance summaries within each cluster and a summary of
balance across clusters. See
[`class-bal.tab.cluster`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.cluster.md)
for details.

## Details

`bal.tab.data.frame()` generates a list of balance summaries for the
covariates and treatment status values given. `bal.tab.formula()` does
the same but uses a formula interface instead. When the formula
interface is used, the formula and data are reshaped into a treatment
vector and `data.frame` of covariates and then simply passed through the
`data.frame` method.

If `weights`, `subclass` and `match.strata` are all `NULL`, balance
information will be presented only for the unadjusted sample.

The argument to `match.strata` corresponds to a factor vector containing
the name or index of each pair/stratum for units conditioned through
matching, for example, using the optmatch package. If more than one of
`weights`, `subclass`, or `match.strata` are specified,
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
will attempt to figure out which one to apply. Currently only one of
these can be applied ta a time.
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
behaves differently depending on whether subclasses are used in
conditioning or not. If they are used,
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
creates balance statistics for each subclass and for the sample in
aggregate. See
[`class-bal.tab.subclass`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.subclass.md)
for more information.

Multiple sets of weights can be supplied simultaneously by entering a
`data.frame` or a character vector containing the names of weight
variables found in `data` or a list of weights vectors or names. The
arguments to `method`, `s.d.denom`, and `estimand`, if any, must be
either the same length as the number of sets of weights or of length
one, where the sole entry is applied to all sets. When standardized
differences are computed for the unadjusted group, they are done using
the first entry to `s.d.denom` or `estimand`. When only one set of
weights is supplied, the output for the adjusted group will simply be
called `"Adj"`, but otherwise will be named after each corresponding set
of weights. Specifying multiple sets of weights will also add components
to other outputs of
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md).

## See also

- [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  for details of calculations.

- [`class-bal.tab.cluster`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.cluster.md)
  for more information on clustered data.

- [`class-bal.tab.imp`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.imp.md)
  for more information on multiply imputed data.

- [`bal.tab.multi()`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.multi.md)
  for more information on multi-category treatments.

## Examples

``` r
data("lalonde", package = "cobalt")
lalonde$p.score <- glm(treat ~ age + educ + race, data = lalonde, 
                       family = "binomial")$fitted.values
covariates <- subset(lalonde, select = c(age, educ, race))

## Propensity score weighting using IPTW
lalonde$iptw.weights <- ifelse(lalonde$treat==1, 
                               1/lalonde$p.score, 
                               1/(1-lalonde$p.score))

# data frame interface:
bal.tab(covariates, treat = "treat", data = lalonde, 
        weights = "iptw.weights", s.d.denom = "pooled")
#> Balance Measures
#>                Type Diff.Adj
#> age         Contin.  -0.1242
#> educ        Contin.   0.0727
#> race_black   Binary   0.0053
#> race_hispan  Binary  -0.0025
#> race_white   Binary  -0.0029
#> 
#> Effective sample sizes
#>            Control Treated
#> Unadjusted  429.    185.  
#> Adjusted    344.33   65.47

# Formula interface:
bal.tab(treat ~ age + educ + race, data = lalonde, 
        weights = "iptw.weights", s.d.denom = "pooled")
#> Balance Measures
#>                Type Diff.Adj
#> age         Contin.  -0.1242
#> educ        Contin.   0.0727
#> race_black   Binary   0.0053
#> race_hispan  Binary  -0.0025
#> race_white   Binary  -0.0029
#> 
#> Effective sample sizes
#>            Control Treated
#> Unadjusted  429.    185.  
#> Adjusted    344.33   65.47

## Propensity score subclassification
lalonde$subclass <- findInterval(lalonde$p.score, 
                                 quantile(lalonde$p.score, 
                                          (0:6)/6), all.inside = TRUE)

# data frame interface:
bal.tab(covariates, treat = "treat", data = lalonde, 
        subclass = "subclass", disp.subclass = TRUE, 
        s.d.denom = "pooled")
#> Balance by subclass
#>  - - - Subclass 1 - - - 
#>                Type Diff.Adj
#> age         Contin.  -1.2029
#> educ        Contin.  -0.2551
#> race_black   Binary   0.0000
#> race_hispan  Binary   0.0000
#> race_white   Binary   0.0000
#> 
#>  - - - Subclass 2 - - - 
#>                Type Diff.Adj
#> age         Contin.   0.4108
#> educ        Contin.   0.3005
#> race_black   Binary   0.0000
#> race_hispan  Binary   0.0000
#> race_white   Binary   0.0000
#> 
#>  - - - Subclass 3 - - - 
#>                Type Diff.Adj
#> age         Contin.  -0.1400
#> educ        Contin.   0.0295
#> race_black   Binary   0.0000
#> race_hispan  Binary  -0.0833
#> race_white   Binary   0.0833
#> 
#>  - - - Subclass 4 - - - 
#>                Type Diff.Adj
#> age         Contin.   0.2294
#> educ        Contin.  -0.4409
#> race_black   Binary   0.3467
#> race_hispan  Binary  -0.3467
#> race_white   Binary   0.0000
#> 
#>  - - - Subclass 5 - - - 
#>                Type Diff.Adj
#> age         Contin.   0.4675
#> educ        Contin.   0.3427
#> race_black   Binary   0.0000
#> race_hispan  Binary   0.0000
#> race_white   Binary   0.0000
#> 
#>  - - - Subclass 6 - - - 
#>                Type Diff.Adj
#> age         Contin.   0.1293
#> educ        Contin.  -0.0838
#> race_black   Binary   0.0000
#> race_hispan  Binary   0.0000
#> race_white   Binary   0.0000
#> 

# Formula interface:
bal.tab(treat ~ age + educ + race, data = lalonde, 
        subclass = "subclass", disp.subclass = TRUE, 
        s.d.denom = "pooled")
#> Balance by subclass
#>  - - - Subclass 1 - - - 
#>                Type Diff.Adj
#> age         Contin.  -1.2029
#> educ        Contin.  -0.2551
#> race_black   Binary   0.0000
#> race_hispan  Binary   0.0000
#> race_white   Binary   0.0000
#> 
#>  - - - Subclass 2 - - - 
#>                Type Diff.Adj
#> age         Contin.   0.4108
#> educ        Contin.   0.3005
#> race_black   Binary   0.0000
#> race_hispan  Binary   0.0000
#> race_white   Binary   0.0000
#> 
#>  - - - Subclass 3 - - - 
#>                Type Diff.Adj
#> age         Contin.  -0.1400
#> educ        Contin.   0.0295
#> race_black   Binary   0.0000
#> race_hispan  Binary  -0.0833
#> race_white   Binary   0.0833
#> 
#>  - - - Subclass 4 - - - 
#>                Type Diff.Adj
#> age         Contin.   0.2294
#> educ        Contin.  -0.4409
#> race_black   Binary   0.3467
#> race_hispan  Binary  -0.3467
#> race_white   Binary   0.0000
#> 
#>  - - - Subclass 5 - - - 
#>                Type Diff.Adj
#> age         Contin.   0.4675
#> educ        Contin.   0.3427
#> race_black   Binary   0.0000
#> race_hispan  Binary   0.0000
#> race_white   Binary   0.0000
#> 
#>  - - - Subclass 6 - - - 
#>                Type Diff.Adj
#> age         Contin.   0.1293
#> educ        Contin.  -0.0838
#> race_black   Binary   0.0000
#> race_hispan  Binary   0.0000
#> race_white   Binary   0.0000
#> 
```

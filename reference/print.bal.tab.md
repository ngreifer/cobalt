# Print Results of a Call to `bal.tab()`

Prints
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
output in a clean way. Provides options for printing.

## Usage

``` r
# S3 method for class 'bal.tab'
print(
  x,
  imbalanced.only,
  un,
  disp.bal.tab,
  disp.call,
  stats,
  disp.thresholds,
  disp,
  which.subclass,
  subclass.summary,
  which.imp,
  imp.summary,
  imp.fun,
  which.treat,
  multi.summary,
  which.time,
  msm.summary,
  which.cluster,
  cluster.summary,
  cluster.fun,
  digits = max(3L, getOption("digits") - 3),
  ...
)
```

## Arguments

- x:

  a `bal.tab` object; the output of a call to
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md).

- imbalanced.only:

  `logical`; whether to display only the covariates that failed to meet
  at least one of balance thresholds. Depends only on whether threshold
  were initial set in the call to
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  and not on any arguments to
  [`print()`](https://rdrr.io/r/base/print.html) (except
  `disp.bal.tab`).

- un:

  `logical`; whether to display balance values for the unadjusted
  sample. Ignored (and set to `TRUE`) if no conditioning was performed.

- disp.bal.tab:

  `logical`; whether to display the table of balance statistics. If
  `FALSE`, only other values (e.g., the call, sample sizes, balance
  tallies, and maximum imbalances) will be presented.

- disp.call:

  `logical`; whether to display the function call for the input object,
  if any.

- stats:

  `character`; which statistic(s) should be reported. For binary or
  multi-category treatments, the options are "mean.diffs" for mean
  differences (standardized or not according the selected
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  options), "variance.ratios" for variance ratios, and "ks.statistics"
  for Kolmogorov-Smirnov statistics. "mean.diffs" is the default. For
  continuous treatments, the only option is "correlations" for
  treatment-covariate correlations. Multiple options are allowed.
  Abbreviations allowed. Statistics that weren't requested in the
  original call to
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  cannot be requested with
  [`print()`](https://rdrr.io/r/base/print.html) unless `quick = FALSE`
  in the original call.

- disp.thresholds:

  `logical`; whether to display thresholds for each statistic for which
  thresholds were originally requested in the call to
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md).
  Should be a named logical vector with names corresponding to the
  thresholds. For example, if thresholds for mean differences were
  requested in
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md),
  set `disp.thresholds = c(m = FALSE)` to prevent them from being
  printed. If a statistic was prevented from being displayed by another
  argument to [`print()`](https://rdrr.io/r/base/print.html), the
  thresholds will not be displayed.

- disp:

  `character`; which distribution summary statistics to display.
  Allowable options include "means" and "sds". Statistics that weren't
  requested in the original call to
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  cannot be requested with
  [`print()`](https://rdrr.io/r/base/print.html) unless `quick = FALSE`
  in the original call.

- which.subclass:

  when used with subclassification, which subclass(es) to display. If
  `NULL`, all subclasses will be displayed. If `NA`, no subclasses will
  be displayed. Otherwise, can be a vector of subclass indices for which
  to display balance. To display the subclasses requested in the
  original call to
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md),
  omit this argument. See
  [`class-bal.tab.subclass`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.subclass.md)
  for details.

- subclass.summary:

  `logical`; when used with subclassification, whether to display the
  subclass balance summary table. If `which.subclass` is `NA`,
  `subclass.summary` will be set to `TRUE`. See
  [`class-bal.tab.subclass`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.subclass.md)
  for details.

- which.imp:

  when used with multiply imputed data, which imputation(s) to display.
  If `NULL`, all imputations will be displayed. If `NA`, no imputations
  will be displayed. Otherwise, can be a vector of imputations numbers
  for which to display balance. To display the imputations requested in
  the original call to
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md),
  omit this argument. See
  [`class-bal.tab.imp`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.imp.md)
  for details.

- imp.summary:

  `logical`; when used with multiply imputed data, whether to display
  the imputation summary table. If `which.imp` is `NA`, `imp.summary`
  will be set to `TRUE`. See
  [`class-bal.tab.imp`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.imp.md)
  for details.

- imp.fun:

  `character`; when used with multiply imputed data, a character vector
  of functions of balance statistics to display when displaying balance
  across imputations. Can be "mean", "min", or "max". More than one are
  allowed. See
  [`class-bal.tab.imp`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.imp.md)
  for details.

- which.treat:

  when used with multi-category treatments, which treatments to display.
  See
  [`bal.tab.multi()`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.multi.md)
  for details.

- multi.summary:

  `logical`; when used with multi-category treatments, whether to
  display the balance summary table across pairwise comparisons. See
  [`bal.tab.multi()`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.multi.md)
  for details.

- which.time:

  when used with longitudinal treatments, which time periods to display
  if longitudinal treatments are used. See
  [`class-bal.tab.msm`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.msm.md)
  for details.

- msm.summary:

  `logical`; when used with longitudinal treatments, whether to display
  the balance summary table across time periods. See
  [`class-bal.tab.msm`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.msm.md)
  for details.

- which.cluster:

  when used with clustered data, which cluster(s) to display. If `NULL`,
  all clusters will be displayed. If `NA`, no clusters will be
  displayed. Otherwise, can be a vector of cluster names or numerical
  indices for which to display balance. Indices correspond to the
  alphabetical order of cluster names. To display the clusters requested
  in the original call to
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md),
  omit this argument. See
  [`class-bal.tab.cluster`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.cluster.md)
  for details.

- cluster.summary:

  `logical`; when used with clustered data, whether to display the
  cluster summary table. If `which.cluster` is `NA`, `cluster.summary`
  will be set to `TRUE`. See
  [`class-bal.tab.cluster`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.cluster.md)
  for details.

- cluster.fun:

  `character`; when used with clustered data, a character vector of
  functions of balance statistics to display when displaying balance
  across clusters. Can be "mean", "min", or "max". More than one are
  allowed. See
  [`class-bal.tab.cluster`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.cluster.md)
  for details.

- digits:

  the number of digits to display.

- ...:

  further arguments passed to or from other methods.

## Details

Simply calling
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
will print its results, but it can be useful to store the results into
an object and print them again later, possibly with different print
options specified. The [`print()`](https://rdrr.io/r/base/print.html)
function automatically dispatches the correct method for the `bal.tab`
object given.

Any parameter used in
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
for calculations, such as `int`, `addl`, or `distance`, cannot be used
with [`print()`](https://rdrr.io/r/base/print.html); only those
parameters listed above, those that solely determine printing options,
can be used. To change computation options, a new call to
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
must be performed.

Prior versions of [`print()`](https://rdrr.io/r/base/print.html) had
separate methods for each `bal.tab` class. Now they are dispatched
internally.

## Note

Unless `quick = FALSE` in the original call to
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
(which is not the default), some values may not be calculated, in which
case using [`print()`](https://rdrr.io/r/base/print.html) will not
display these values even when requested. For example, if `stats = "m"`
and `quick = TRUE` in the original call to
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
(the default for both), setting `stats = "ks"` in
[`print()`](https://rdrr.io/r/base/print.html) will not print the KS
statistics because they were not calculated.

## See also

[`print()`](https://rdrr.io/r/base/print.html),
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)

[`display-options`](https://ngreifer.github.io/cobalt/reference/display-options.md)
for further information on some of these options.

## Examples

``` r
data("lalonde", package = "cobalt")
library(WeightIt)

w.out <- weightit(treat ~ age + educ + married +
                    race + re74 + re75, 
                  data = lalonde)

b <- bal.tab(w.out, stats = c("m", "v", "ks"), 
             un = TRUE, v.threshold = 2)

print(b, un = FALSE, stats = c("m", "v"),
      disp.thresholds = c(v = FALSE))
#> Balance Measures
#>                 Type Diff.Adj V.Ratio.Adj
#> prop.score  Distance   0.1660      0.9426
#> age          Contin.  -0.1885      0.3730
#> educ         Contin.   0.0861      0.5164
#> married       Binary  -0.1043           .
#> race_black    Binary   0.0612           .
#> race_hispan   Binary   0.0104           .
#> race_white    Binary  -0.0715           .
#> re74         Contin.  -0.2825      0.8392
#> re75         Contin.  -0.1614      0.9440
#> 
#> Effective sample sizes
#>            Control Treated
#> Unadjusted  429.    185.  
#> Adjusted    330.88   65.26
```

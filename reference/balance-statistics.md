# Balance Statistics in `bal.tab` and `love.plot`

[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
and
[`love.plot()`](https://ngreifer.github.io/cobalt/reference/love.plot.md)
display balance statistics for the included covariates. The `stats`
argument in each of these functions controls which balance statistics
are to be displayed. The argument to `stats` should be a character
vector with the names of the desired balance statistics.

This page describes all of the available balance statistics and how to
request them. Abbreviations are allowed, so you can use the first few
letters of each balance statistics to request it instead of typing out
its whole name. That convention is used throughout the documentation.
For example, to request mean differences and variance ratios in
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md) or
[`love.plot()`](https://ngreifer.github.io/cobalt/reference/love.plot.md),
you could include `stats = c("m", "v")`. In addition, the `thresholds`
argument uses the same naming conventions and can be used to request
balance thresholds on each statistic. For example, to request a balance
threshold of .1 for mean differences, you could include
`thresholds = c(m = .1)`.

Below, each allowable entry to `stats` and `thresholds` are described,
along with other details or option that accompany them.

### Binary/Multi-Category Treatments

- `"mean.diffs"`:

  Mean differences as computed by
  [`col_w_smd()`](https://ngreifer.github.io/cobalt/reference/balance-summary.md).
  Can be abbreviated as `"m"`. Setting the arguments `continuous` and
  `binary` to either `"std"` or `"raw"` will determine whether
  standardized mean differences or raw mean differences are calculated
  for continuous and categorical variables, respectively. When
  standardized mean differences are requested, the `s.d.denom` argument
  controls how the standardization occurs. When `abs = TRUE`, negative
  values become positive. Mean differences are requested by default when
  no entry to `stats` is provided.

- `"variance.ratios"`:

  Variance ratios as computed by
  [`col_w_vr()`](https://ngreifer.github.io/cobalt/reference/balance-summary.md).
  Can be abbreviated as `"v"`. Will not be computed for binary
  variables. When `abs = TRUE`, values less than 1 will have their
  inverse taken. When used with `love.plot`, the x-axis scaled will be
  logged so that, e.g., .5 is as far away from 1 as 2 is.

- `"ks.statistics"`:

  Kolmogorov-Smirnov (KS) statistics as computed by
  [`col_w_ks()`](https://ngreifer.github.io/cobalt/reference/balance-summary.md).

- `"ovl.coefficients"`:

  Overlapping (OVL) statistics as computed by
  [`col_w_ovl()`](https://ngreifer.github.io/cobalt/reference/balance-summary.md).
  Can be abbreviated as `"ovl"`. Additional arguments passed to
  [`col_w_ovl()`](https://ngreifer.github.io/cobalt/reference/balance-summary.md),
  such as `integrate` or `bw`, can be supplied to
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  or
  [`love.plot()`](https://ngreifer.github.io/cobalt/reference/love.plot.md).

### Continuous Treatments

- `"correlations"`:

  Pearson correlations as computed by
  [`col_w_cov()`](https://ngreifer.github.io/cobalt/reference/balance-summary.md).
  Can be abbreviated as `"cor"`. Setting the arguments `continuous` and
  `binary` to either `"std"` or `"raw"` will determine whether
  correlations or covariances are calculated for continuous and
  categorical variables, respectively (they are both `"std"` by
  default). When correlations are requested, the `s.d.denom` argument
  controls how the standardization occurs. When `abs = TRUE`, negative
  values become positive. Pearson correlations are requested by default
  when no entry to `stats` is provided.

- `"spearman.correlations"`:

  Spearman correlations as computed by
  [`col_w_cov()`](https://ngreifer.github.io/cobalt/reference/balance-summary.md).
  Can be abbreviated as `"sp"`. All arguments are the same as those for
  `"correlations"`. When `abs = TRUE`, negative values become positive.

- `"distance.correlations"`:

  Distance correlations as computed by
  [`col_w_dcov()`](https://ngreifer.github.io/cobalt/reference/balance-summary.md).
  Can be abbreviated as `"dis"`. Setting the arguments `continuous` and
  `binary` to either `"std"` or `"raw"` will determine whether distance
  correlations or distance covariances are calculated for continuous and
  categorical variables, respectively (they are both `"std"` by
  default). When distance correlations are requested, the `s.d.denom`
  argument controls how the standardization occurs.

- `"mean.diffs.target"`:

  Mean differences computed between the weighted and unweighted sample
  to ensure the weighted sample is representative of the original
  population. Can be abbreviated as `"m"`. Setting the arguments
  `continuous` and `binary` to either `"std"` or `"raw"` will determine
  whether standardized mean differences or raw mean differences are
  calculated for continuous and categorical variables, respectively. The
  standardization factor will be computed in the unweighted sample. When
  `abs = TRUE`, negative values become positive. This statistic is only
  computed for the adjusted samples.

- `"ks.statistics.target"`:

  KS statistics computed between the weighted and unweighted sample to
  ensure the weighted sample is representative of the original
  population. Can be abbreviated as `"ks"`. This statistic is only
  computed for the adjusted samples.

- `"ovl.coefficients.target"`:

  Overlapping coefficients computed between the weighted and unweighted
  sample to ensure the weighted sample is representative of the original
  population. Can be abbreviated as `"ovl"`. This statistic is only
  computed for the adjusted samples.

If a statistic is requested in `thresholds`, it will automatically be
placed in `stats`. For example,
`bal.tab(..., stats = "m", thresholds = c(v = 2))` will display both
mean differences and variance ratios, and the variance ratios will have
a balance threshold set to 2.

## Examples

``` r
data(lalonde)

#Binary treatments
bal.tab(treat ~ age + educ + married + re74, data = lalonde,
        stats = c("m", "v", "ks"))
#> Error in bal.tab(treat ~ age + educ + married + re74, data = lalonde,     stats = c("m", "v", "ks")): The given response variable, `treat`, is not a variable in `data` or the global
#> environment.
love.plot(treat ~ age + educ + married + re74, data = lalonde,
          stats = c("m", "v", "ks"), binary = "std",
          thresholds = c(m = .1, v = 2))
#> Error in love.plot(treat ~ age + educ + married + re74, data = lalonde,     stats = c("m", "v", "ks"), binary = "std", thresholds = c(m = 0.1,         v = 2)): The given response variable, `treat`, is not a variable in `data` or the global
#> environment.

#Continuous treatments
bal.tab(re75 ~ age + educ + married + re74, data = lalonde,
        stats = c("cor", "sp"))
#> Error in bal.tab(re75 ~ age + educ + married + re74, data = lalonde, stats = c("cor",     "sp")): The variable "educ" cannot be found. Be sure it is entered correctly or supply
#> a dataset that contains this varialble to `data`.
love.plot(re75 ~ age + educ + married + re74, data = lalonde,
          thresholds = c(cor = .1, sp = .1))
#> Error in love.plot(re75 ~ age + educ + married + re74, data = lalonde,     thresholds = c(cor = 0.1, sp = 0.1)): The variable "educ" cannot be found. Be sure it is entered correctly or supply
#> a dataset that contains this varialble to `data`.
```

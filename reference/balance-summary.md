# Compute Balance and Summary Statistics for Covariates

These functions quickly compute balance statistics for the given
covariates. These functions are used in
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md),
but they are available for use in programming without having to call
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md) to
get them.

- `col_w_mean()` computes the (weighted) means for a set of covariates
  and weights and is essentially a weighted version of
  [`colMeans()`](https://rdrr.io/pkg/Matrix/man/colSums-methods.html).

- `col_w_sd()` computes the (weighted) standard deviations for a set of
  covariates and weights.

- `col_w_smd()` computes the (weighted) (absolute) (standardized)
  difference in means for a set of covariates, a binary treatment, and
  weights.

- `col_w_vr()` computes the (weighted) variance ratio for a set of
  covariates, a binary treatment, and weights.

- `col_w_ks()` computes the (weighted) Kolmogorov-Smirnov (KS) statistic
  for a set of covariates, a binary treatment, and weights.

- `col_w_ovl()` computes the complement of the (weighted) overlapping
  coefficient compliment for a set of covariates, a binary treatment,
  and weights (based on Franklin et al, 2014).

- `col_w_cov()` and `col_w_corr()` compute the (weighted) (absolute)
  treatment-covariate covariance or correlation for a set of covariates,
  a continuous treatment, and weights.

- `col_w_dcov()` and `col_w_dcorr()` compute the (weighted)
  treatment-covariate distance covariance or distance correlation for a
  set of covariates, a continuous treatment, and weights.

## Usage

``` r
col_w_mean(
  mat,
  weights = NULL,
  s.weights = NULL,
  subset = NULL,
  na.rm = TRUE,
  ...
)

col_w_sd(
  mat,
  weights = NULL,
  s.weights = NULL,
  bin.vars,
  subset = NULL,
  na.rm = TRUE,
  ...
)

col_w_smd(
  mat,
  treat,
  weights = NULL,
  std = TRUE,
  s.d.denom = "pooled",
  abs = FALSE,
  s.weights = NULL,
  bin.vars,
  subset = NULL,
  weighted.weights = weights,
  na.rm = TRUE,
  ...
)

col_w_vr(
  mat,
  treat,
  weights = NULL,
  abs = FALSE,
  s.weights = NULL,
  bin.vars,
  subset = NULL,
  na.rm = TRUE,
  ...
)

col_w_ks(
  mat,
  treat,
  weights = NULL,
  s.weights = NULL,
  bin.vars,
  subset = NULL,
  na.rm = TRUE,
  ...
)

col_w_ovl(
  mat,
  treat,
  weights = NULL,
  s.weights = NULL,
  bin.vars,
  subset = NULL,
  na.rm = TRUE,
  integrate = TRUE,
  steps = 1001L,
  ...
)

col_w_cov(
  mat,
  treat,
  weights = NULL,
  type = "pearson",
  std = FALSE,
  s.d.denom = "all",
  abs = FALSE,
  s.weights = NULL,
  bin.vars,
  subset = NULL,
  weighted.weights = weights,
  na.rm = TRUE,
  ...
)

col_w_corr(
  mat,
  treat,
  weights = NULL,
  type = "pearson",
  s.d.denom = "all",
  abs = FALSE,
  s.weights = NULL,
  bin.vars,
  subset = NULL,
  weighted.weights = weights,
  na.rm = TRUE,
  ...
)

col_w_dcov(
  mat,
  treat,
  weights = NULL,
  std = FALSE,
  s.d.denom = "all",
  s.weights = NULL,
  subset = NULL,
  weighted.weights = weights,
  na.rm = TRUE,
  ...
)

col_w_dcorr(
  mat,
  treat,
  weights = NULL,
  s.d.denom = "all",
  s.weights = NULL,
  subset = NULL,
  weighted.weights = weights,
  na.rm = TRUE,
  ...
)
```

## Arguments

- mat:

  a numeric matrix or a data frame containing the covariates for which
  the statistic is to be computed. If a data frame,
  [`splitfactor()`](https://ngreifer.github.io/cobalt/reference/splitfactor.md)
  with `drop.first = "if2"` will be called if any character or factor
  variables are present. This can slow down the function, so it's
  generally best to supply a numeric matrix. If a numeric vector is
  supplied, it will be converted to a 1-column matrix first.

- weights:

  `numeric`; an optional set of weights used to compute the weighted
  statistics. If sampling weights are supplied through `s.weights`, the
  `weights` should not incorporate these weights, as `weights` and
  `s.weights` will be multiplied together prior to computing the
  weighted statistics.

- s.weights:

  `numeric`; an optional set of sampling weights used to compute the
  weighted statistics. If weights are supplied through `weights`,
  `weights` and `s.weights` will be multiplied together prior to
  computing the weighted statistics. Some functions use `s.weights` in a
  particular way; for others, supplying `weights` and `s.weights` is
  equivalent to supplying their product to either `weights` or
  `s.weights`. See Details.

- subset:

  a `logical` vector with length equal to the number of rows of `mat`
  used to subset the data. See Details for notes on its use with
  `col_w_smd()`, `col_w_cov()`, and `col_w_corr()`.

- na.rm:

  `logical`; whether `NA`s should be ignored or not. If `FALSE`, any
  variable with any `NA`s will have its corresponding statistic returned
  as `NA`. If `TRUE`, any variable with any `NA`s will have its
  corresponding statistic computed as if the missing value were not
  there.

- ...:

  for all functions, additional arguments supplied to
  [`splitfactor()`](https://ngreifer.github.io/cobalt/reference/splitfactor.md)
  when `mat` is a data.frame. `data`, `var.name`, `drop.first`, and
  `drop.level` are ignored; `drop.first` is automatically set to
  `"if2"`. For `col_w_ovl()`, other arguments passed to
  [`density()`](https://rdrr.io/r/stats/density.html) besides `x` and
  `weights`. Note that the default value for `bw` when unspecified is
  `"nrd"` rather than the default in
  [`density()`](https://rdrr.io/r/stats/density.html), which is
  `"nrd0"`.

- bin.vars:

  a vector used to denote whether each variable is binary or not. Can be
  a `logical` vector with length equal to the number of columns of `mat`
  or a vector of numeric indices or character names of the binary
  variables. If missing (the default), the function will figure out
  which covariates are binary or not, which can increase computation
  time. If `NULL`, it will be assumed no variables are binary. All
  functions other than `col_w_mean()` treat binary variables different
  from continuous variables. If a factor or character variable is in
  `mat`, all the dummies created will automatically be marked as binary,
  but it should still receive an entry when `bin.vars` is supplied as
  `logical`.

- treat:

  a vector of treatment status for each individual. For `col_w_smd()`,
  `col_w_vr()`, `col_w_ks()`, and `col_w_ovl()`, `treat` should have
  exactly two unique values. For `col_w_cov()`, `col_w_corr()`,
  `col_w_dcov()`, and `col_w_dcorr()`, `treat` should be a many-valued
  numeric vector.

- std:

  `logical`; for `col_w_smd()`, whether the computed mean differences
  for each variable should be standardized; for `col_w_cov()`, whether
  treatment-covariate correlations should be computed (`TRUE`) rather
  than covariances (`FALSE`); for `col_w_dcov()`, whether
  treatment-covariate distance correlations should be computed (`TRUE`)
  rather than distance covariances (`FALSE`). Can have either length 1,
  whereby all variables will be standardized or not, or length equal to
  the number of columns of `mat`, whereby only variables with a value of
  `TRUE` will be standardized. See Details.

- s.d.denom:

  for `col_w_smd()`, `col_w_cov()`, and `col_w_dcov()` when `std` is
  `TRUE` for some variables, and for `col_w_corr()` and `col_w_dcorr()`,
  how the standardization factor should be computed. For `col_w_smd()`
  (i.e., when computing standardized mean differences), allowable
  options include

  - `"treated"` - uses the standard deviation of the variable in the
    treated group

  - `"control"` - uses the standard deviation of the variable in the
    control group

  - `"pooled"` - uses the square root of the average of the variances of
    the variable in the treated and control groups

  - `"all"` - uses the standard deviation of the variable in the full
    sample

  - `"weighted"` - uses the standard deviation of the variable in the
    full sample weighted by `weighted.weights`

  - `"hedges"` - uses the small-sample corrected version of Hedge's G
    described in the WWC Procedures Handbook (see References)

  - the name of one of the treatment values - uses the standard
    deviation of the variable in that treatment group.

  For `col_w_cov()`, `col_w_corr()`, `col_w_dcov()`, and
  `col_w_dcorr()`, only `"all"` and `"weighted"` are allowed.
  Abbreviations allowed. This can also be supplied as a numeric vector
  of standard deviations with length equal to the number of columns of
  `mat`; the values will be used as the standardization factors.

- abs:

  `logical`; for `col_w_smd()`, `col_w_cov()`, and `col_w_corr()`,
  whether the returned statistics should be in absolute value (`TRUE`)
  or not. For `col_w_vr()`, whether the ratio should always include the
  larger variance in the numerator, so that the ratio is always greater
  than or equal to 1. Default is `FALSE`.

- weighted.weights:

  for `col_w_smd()`, `col_w_cov()`, `col_w_corr()`, `col_w_dcov()`, and
  `col_w_dcorr()`, when `std = TRUE` and `s.d.denom = "weighted"`, a
  vector of weights to be applied to the computation of the denominator
  standard deviation. If not specified, will use the argument to
  `weights`. When `s.d.denom` is not `"weighted"`, this is ignored. The
  main purpose of this is to allow `weights` to be `NULL` while
  weighting the denominator standard deviations for assessing balance in
  the unweighted sample but using the standard deviations of the
  weighted sample.

- integrate:

  `logical`; for `col_w_ovl()`, whether to use
  [`integrate()`](https://rdrr.io/r/stats/integrate.html) to calculate
  the area of overlap for continuous variables. If `FALSE`, a midpoint
  Riemann sum will be used instead. The Riemann sum is a little slower
  and very slightly imprecise (unnoticibly in most contexts). When
  `TRUE`, [`integrate()`](https://rdrr.io/r/stats/integrate.html) will
  be tried, and if it fails, the Riemann sum will be used as a fallback.
  The default (`TRUE`) is to use
  [`integrate()`](https://rdrr.io/r/stats/integrate.html) when possible.

- steps:

  for `col_w_ovl()` when `integrate = FALSE`, the number of points to
  use to compute the Riemann sum to approximate the integral. Default is
  1001 for 1000 partitions.

- type:

  for `col_w_cov()` and `col_w_corr()`, the type of
  covariance/correlation to be computed. Allowable options include
  `"pearson"` and `"spearman"`. When `"spearman"` is requested, the
  covariates and treatment are first turned into ranks using
  [`rank()`](https://rdrr.io/r/base/rank.html) with `na.last = "keep"`.

## Value

A vector of balance statistics, one for each variable in `mat`. If `mat`
has column names, the output will be named as well.

## Details

`col_w_mean()` computes column weighted means for a matrix of variables.
It is similar to
[`colMeans()`](https://rdrr.io/pkg/Matrix/man/colSums-methods.html) but
(optionally) incorporates weights. `weights` and `s.weights` are
multiplied together prior to being used, and there is no distinction
between them. This could be used to compute the weighted means of each
covariate in the general population to examine the degree to which a
weighting method has left the weighted samples resembling the original
population.

`col_w_sd()` computes column weighted standard deviations for a matrix
of variables. `weights` and `s.weights` are multiplied together prior to
being used, and there is no distinction between them. The variance of
binary variables is computed as \\p(1-p)\\, where \\p\\ is the
(weighted) proportion of 1s, while the variance of continuous variables
is computed using the standard formula; the standard deviation is the
square root of this variance.

`col_w_smd()` computes the mean difference for each covariate between
treatment groups defined by `treat`. These mean differences can
optionally be weighted, standardized, and/or in absolute value. The
standardization factor is computed using the unweighted standard
deviation or variance when `s.weights` are absent, and is computed using
the `s.weights`-weighted standard deviation or variance when `s.weights`
are present, except when `s.d.denom = "weighted"`, in which case the
product of `weighted.weights` and `s.weights` (if present) are used to
weight the standardization factor. The standardization factor is
computed using the whole sample even when `subset` is used. Note that
unlike
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md),
`col_w_smd()` requires the user to specify whether each individual
variable should be standardized using `std` rather than relying on
`continuous` or `binary`. The weighted mean difference is computed using
the product of `weights` and `s.weights`, if specified. The variance of
binary variables is computed as \\p(1-p)\\, where \\p\\ is the
(weighted) proportion of 1s, while the variance of continuous variables
is computed using the standard formula.

`col_w_vr()` computes the variance ratio for each covariate between
treatment groups defined by `treat`. When `abs = TRUE`,
`pmax(out, 1/out)` is applied to the output so that the ratio is always
greater than or equal to 1. For binary variables, the variance is
computed as \\p(1-p)\\, where \\p\\ is the (weighted) proportion of 1s,
while the variance of continuous variables is computed using the
standard formula. Note that in
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md),
variance ratios are not computed for binary variables, while here, they
are (but likely should not be interpreted). `weights` and `s.weights`
are multiplied together prior to being used, and there is no distinction
between them. Because of how the weighted variance is computed, exactly
balanced groups may have variance ratios that differ slightly from 1.

`col_w_ks()` computes the KS statistic for each covariate using the
method implemented in twang. The KS statistics can optionally be
weighted. For binary variables, the KS statistic is just the difference
in proportions. `weights` and `s.weights` are multiplied together prior
to being used, and there is no distinction between them.

`col_w_ovl()` computes the complement of the overlapping coefficient as
described by Franklin et al. (2014). It does so by computing the density
of the covariate in the treated and control groups, then finding the
area where those density overlap, and subtracting that number from 1,
yielding a value between 0 and 1 where 1 indicates complete imbalance,
and 0 indicates perfect balance.
[`density()`](https://rdrr.io/r/stats/density.html) is used to model the
density in each group. The bandwidth of the covariate in the smaller
treatment group is used for both groups. The area of overlap can be
computed using `integrate`, which quickly and accurately computes the
integral, or using a midpoint Riemann sum with 1000 partitions, which
approximates the area more slowly. A reason to prefer the Riemann sum is
that `integrate` can fail for unknown reasons, though Riemann sums will
fail with some extreme distributions. When either method fails, the
resulting value will be `NA`. For binary variables, the complement of
the overlapping coefficient is just the difference in proportions.
`weights` and `s.weights` are multiplied together prior to being used,
and there is no distinction between them. The weights are used to
compute the weighted density by supplying them to the `weights` argument
of `density`.

`col_w_cov()` computes the covariances between a continuous treatment
and the covariates to assess balance for a continuous treatment as
recommended in Austin (2019). These covariances can optionally be
weighted or in absolute value or can be requested as correlations (i.e.,
standardized covariances). Each correlations is computed as the
covariance between the treatment and covariate divided by a
standardization factor, which is equal to the square root of the product
of the variance of treatment and the variance of the covariate. The
standardization factor is computed using the unweighted variances when
`s.weights` are absent, and is computed using the sampling weighted
variances when `s.weights` are present, except when
`s.d.denom = "weighted"`, in which case the product of
`weighted.weights` and `s.weights` (if present) are used to weight the
standardization factor. For this reason, the computed correlation can be
greater than 1 or less than -1. The standardization factor is always
computed using the whole sample even when `subset` is used. The
covariance is computed using the product of `weights` and `s.weights`,
if specified. The variance of binary variables is computed as
\\p(1-p)\\, where \\p\\ is the (weighted) proportion of 1s, while the
variance of continuous variables is computed using the standard formula.

`col_w_corr()` is a wrapper for `col_w_cov` with `std` set to `TRUE`.

`col_w_dcov()` computes the distance covariances between a continuous
treatment and the covariates to assess balance for a continuous
treatment. A multivariate version is described by Huling et al. (2023)
for computing a scalar value that represents the balance for all
covariates simultaneously; the statistic computed here is for one
covariate at a time. The distance covariances can optionally be weighted
or can be requested as distance correlations (i.e., standardized
distance covariances). The distance correlations are computed as the
distance covariance between the treatment and covariate divided by a
standardization factor, which is equal to the square root of the product
of the distance variance of treatment and the distance variance of the
covariate, where the distance variance is the distance covariance of a
variable with itself. The standardization factor is computed using the
unweighted distance variances when `s.weights` are absent, and is
computed using the sampling weighted distance variances when `s.weights`
are present, except when `s.d.denom = "weighted"`, in which case the
product of `weighted.weights` and `s.weights` (if present) are used to
weight the standardization factor. For this reason, the computed
distance correlation can be greater than 1. The standardization factor
is always computed using the whole sample even when `subset` is used.
The distance covariance is computed using the product of `weights` and
`s.weights`, if specified.

`col_w_dcorr()` is a wrapper for `col_w_dcov` with `std` set to `TRUE`.

## References

Austin, P. C. (2019). Assessing covariate balance when using the
generalized propensity score with quantitative or continuous exposures.
*Statistical Methods in Medical Research*, 28(5), 1365–1377.
[doi:10.1177/0962280218756159](https://doi.org/10.1177/0962280218756159)

Franklin, J. M., Rassen, J. A., Ackermann, D., Bartels, D. B., &
Schneeweiss, S. (2014). Metrics for covariate balance in cohort studies
of causal effects. *Statistics in Medicine*, 33(10), 1685–1699.
[doi:10.1002/sim.6058](https://doi.org/10.1002/sim.6058)

Huling, J. D., Greifer, N., & Chen, G. (2023). Independence Weights for
Causal Inference with Continuous Treatments. *Journal of the American
Statistical Association*, 0(0), 1–14.
[doi:10.1080/01621459.2023.2213485](https://doi.org/10.1080/01621459.2023.2213485)

What Works Clearinghouse. (2020). WWC Procedures Handbook (Version 4.1).
Retrieved from <https://ies.ed.gov/ncee/wwc/Handbooks>

## See also

- [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)

- [`bal.compute()`](https://ngreifer.github.io/cobalt/reference/bal.compute.md)

- [balance-statistics](https://ngreifer.github.io/cobalt/reference/balance-statistics.md)

## Examples

``` r
data("lalonde", package = "cobalt")

treat <- lalonde$treat
covs <- subset(lalonde, select = -c(treat, re78))
covs0 <- splitfactor(covs, drop.first = "if2")
bin.vars <- c(FALSE, FALSE, TRUE, TRUE, TRUE,
              TRUE, TRUE, FALSE, FALSE)
W <- WeightIt::weightit(treat ~ covs, method = "glm", 
                        estimand = "ATE")
weights <- W$weights

round(data.frame(
    m0 = col_w_mean(covs0, weights = weights, subset = treat == 0),
    sd0 = col_w_sd(covs0, weights = weights,
                   bin.vars = bin.vars, subset = treat == 0),
    m1 = col_w_mean(covs0, weights = weights, subset = treat == 1),
    sd1 = col_w_sd(covs0, weights = weights,
                   bin.vars = bin.vars, subset = treat == 1),
    smd = col_w_smd(covs0, treat = treat, weights = weights,
                    std = TRUE, bin.vars = bin.vars),
    vr = col_w_vr(covs0, treat = treat, weights = weights,
                  bin.vars = bin.vars),
    ks = col_w_ks(covs0, treat = treat, weights = weights,
                  bin.vars = bin.vars),
    ovl = col_w_ovl(covs0, treat = treat, weights = weights,
                  bin.vars = bin.vars),
    row.names = colnames(covs0)
), 4)
#>                    m0       sd0        m1       sd1     smd     vr     ks
#> age           27.1000   10.8071   25.5663    6.5640 -0.1676 0.3689 0.1912
#> educ          10.2863    2.7430   10.6064    2.0631  0.1296 0.5657 0.0768
#> race_black     0.3979    0.4895    0.4478    0.4973  0.1302 1.0322 0.0499
#> race_hispan    0.1170    0.3215    0.1217    0.3269  0.0156 1.0344 0.0047
#> race_white     0.4851    0.4998    0.4305    0.4951 -0.1378 0.9815 0.0546
#> married        0.4089    0.4916    0.3146    0.4643 -0.2102 0.8920 0.0944
#> nodegree       0.6250    0.4841    0.5702    0.4950 -0.1157 1.0456 0.0547
#> re74        4552.7364 6339.3397 2932.1845 5743.4197 -0.2740 0.8208 0.3121
#> re75        2172.0386 3161.2645 1658.0651 3091.1829 -0.1579 0.9562 0.1526
#>                ovl
#> age         0.2509
#> educ        0.1075
#> race_black  0.0499
#> race_hispan 0.0047
#> race_white  0.0546
#> married     0.0944
#> nodegree    0.0547
#> re74        0.3240
#> re75        0.1209

# Compare to bal.tab():
bal.tab(covs, treat = treat, weights = weights,
        disp = c("m", "sd"),
        stats = c("m", "v", "ks", "ovl"),
        estimand = "ATE", method = "weighting",
        binary = "std")
#> Balance Measures
#>                Type   M.0.Adj  SD.0.Adj   M.1.Adj  SD.1.Adj Diff.Adj
#> age         Contin.   27.1000   10.8071   25.5663    6.5640  -0.1676
#> educ        Contin.   10.2863    2.7430   10.6064    2.0631   0.1296
#> race_black   Binary    0.3979    0.4895    0.4478    0.4973   0.1302
#> race_hispan  Binary    0.1170    0.3215    0.1217    0.3269   0.0156
#> race_white   Binary    0.4851    0.4998    0.4305    0.4951  -0.1378
#> married      Binary    0.4089    0.4916    0.3146    0.4643  -0.2102
#> nodegree     Binary    0.6250    0.4841    0.5702    0.4950  -0.1157
#> re74        Contin. 4552.7364 6339.3397 2932.1845 5743.4197  -0.2740
#> re75        Contin. 2172.0386 3161.2645 1658.0651 3091.1829  -0.1579
#>             V.Ratio.Adj KS.Adj OVL.Adj
#> age              0.3689 0.1912  0.2509
#> educ             0.5657 0.0768  0.1075
#> race_black            . 0.0499  0.0499
#> race_hispan           . 0.0047  0.0047
#> race_white            . 0.0546  0.0546
#> married               . 0.0944  0.0944
#> nodegree              . 0.0547  0.0547
#> re74             0.8208 0.3121  0.3240
#> re75             0.9562 0.1526  0.1209
#> 
#> Effective sample sizes
#>            Control Treated
#> Unadjusted  429.    185.  
#> Adjusted    329.01   58.33
```

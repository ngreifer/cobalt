# cobalt: Covariate Balance Tables and Plots

A set of tools for assessing covariate balance in observational studies
numerically and graphically. The functions provide integration with the
major R packages used for balancing covariates, including MatchIt,
WeightIt, twang, CBPS, and many others, and support objects not made
using these packages. They support binary, multi-category and continuous
treatments, point and longitudinal treatments, and clustered and
multiply imputed data.

The main functions of cobalt are the following:

- [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md) -
  generate tables of balance statistics before and after matching,
  weighting, or subclassification

- [`bal.plot()`](https://ngreifer.github.io/cobalt/reference/bal.plot.md) -
  generate plots to assess balance visually on one covariate at a time

- [`love.plot()`](https://ngreifer.github.io/cobalt/reference/love.plot.md) -
  generate plots to summarize and report balance statistics

Other functions include
[`get.w()`](https://ngreifer.github.io/cobalt/reference/get.w.md) for
extracting weights from objects produced by other packages,
[`col_w_smd()`](https://ngreifer.github.io/cobalt/reference/balance-summary.md)
(and friends documented on the same page) for computing (weighted)
balance statistics outside of
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md),
[`bal.compute()`](https://ngreifer.github.io/cobalt/reference/bal.compute.md)
for computing scalar balance statistics efficiently, and
[`splitfactor()`](https://ngreifer.github.io/cobalt/reference/splitfactor.md)
for splitting factor variables in a dataset into dummy variables.

cobalt has several vignettes, which can be accessed using
`vignette(package = "cobalt")` or visiting the website at
<https://ngreifer.github.io/cobalt/>.

## See also

Useful links:

- <https://ngreifer.github.io/cobalt/>

- <https://github.com/ngreifer/cobalt>

- Report bugs at <https://github.com/ngreifer/cobalt/issues>

## Author

**Maintainer**: Noah Greifer <noah.greifer@gmail.com>
([ORCID](https://orcid.org/0000-0003-3067-7154))

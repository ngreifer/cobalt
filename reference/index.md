# Package index

## Primary Balance Assessment Functions

- [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  : Display Balance Statistics in a Table
- [`bal.plot()`](https://ngreifer.github.io/cobalt/reference/bal.plot.md)
  : Visualize Distributional Balance
- [`love.plot()`](https://ngreifer.github.io/cobalt/reference/love.plot.md)
  : Display Balance Statistics in a Love Plot

## bal.tab() Methods

- [`bal.tab(`*`<CBPS>`*`)`](https://ngreifer.github.io/cobalt/reference/bal.tab.CBPS.md)
  :

  Balance Statistics for `CBPS` Objects

- [`bal.tab(`*`<Match>`*`)`](https://ngreifer.github.io/cobalt/reference/bal.tab.Match.md)
  :

  Balance Statistics for `Matching` Objects

- [`bal.tab(`*`<cem.match>`*`)`](https://ngreifer.github.io/cobalt/reference/bal.tab.cem.match.md)
  :

  Balance Statistics for `cem` Objects

- [`bal.tab(`*`<default>`*`)`](https://ngreifer.github.io/cobalt/reference/bal.tab.default.md)
  : Balance Statistics for Other Objects

- [`bal.tab(`*`<designmatch>`*`)`](https://ngreifer.github.io/cobalt/reference/bal.tab.designmatch.md)
  :

  Balance Statistics for `designmatch` Objects

- [`bal.tab(`*`<ebalance>`*`)`](https://ngreifer.github.io/cobalt/reference/bal.tab.ebalance.md)
  :

  Balance Statistics for `ebalance` Objects

- [`bal.tab(`*`<formula>`*`)`](https://ngreifer.github.io/cobalt/reference/bal.tab.formula.md)
  [`bal.tab(`*`<data.frame>`*`)`](https://ngreifer.github.io/cobalt/reference/bal.tab.formula.md)
  [`bal.tab(`*`<matrix>`*`)`](https://ngreifer.github.io/cobalt/reference/bal.tab.formula.md)
  : Balance Statistics for Data Sets

- [`bal.tab(`*`<matchit>`*`)`](https://ngreifer.github.io/cobalt/reference/bal.tab.matchit.md)
  :

  Balance Statistics for `MatchIt` Objects

- [`bal.tab(`*`<mimids>`*`)`](https://ngreifer.github.io/cobalt/reference/bal.tab.mimids.md)
  :

  Balance Statistics for `MatchThem` Objects

- [`bal.tab(`*`<optmatch>`*`)`](https://ngreifer.github.io/cobalt/reference/bal.tab.optmatch.md)
  :

  Balance Statistics for `optmatch` Objects

- [`bal.tab(`*`<ps>`*`)`](https://ngreifer.github.io/cobalt/reference/bal.tab.ps.md)
  :

  Balance Statistics for `twang` Objects

- [`bal.tab(`*`<sbwcau>`*`)`](https://ngreifer.github.io/cobalt/reference/bal.tab.sbwcau.md)
  :

  Balance Statistics for `sbw` Objects

- [`bal.tab(`*`<formula.list>`*`)`](https://ngreifer.github.io/cobalt/reference/bal.tab.time.list.md)
  [`bal.tab(`*`<data.frame.list>`*`)`](https://ngreifer.github.io/cobalt/reference/bal.tab.time.list.md)
  : Balance Statistics for Longitudinal Datasets

- [`bal.tab(`*`<weightit>`*`)`](https://ngreifer.github.io/cobalt/reference/bal.tab.weightit.md)
  :

  Balance Statistics for `WeightIt` Objects

## bal.tab() Classes

- [`class-bal.tab.cluster`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.cluster.md)
  :

  Using
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  with Clustered Data

- [`class-bal.tab.imp`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.imp.md)
  :

  Using
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  with Multiply Imputed Data

- [`class-bal.tab.msm`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.msm.md)
  :

  Using
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  with Longitudinal Treatments

- [`class-bal.tab.multi`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.multi.md)
  :

  Using
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  with Multi-Category Treatments

- [`class-bal.tab.subclass`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.subclass.md)
  :

  Using
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  with Subclassified Data

## Balance and Display Options

- [`balance-statistics`](https://ngreifer.github.io/cobalt/reference/balance-statistics.md)
  :

  Balance Statistics in `bal.tab` and `love.plot`

- [`display-options`](https://ngreifer.github.io/cobalt/reference/display-options.md)
  :

  Options for Displaying
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  Output

## Helper Functions

- [`col_w_mean()`](https://ngreifer.github.io/cobalt/reference/balance-summary.md)
  [`col_w_sd()`](https://ngreifer.github.io/cobalt/reference/balance-summary.md)
  [`col_w_smd()`](https://ngreifer.github.io/cobalt/reference/balance-summary.md)
  [`col_w_vr()`](https://ngreifer.github.io/cobalt/reference/balance-summary.md)
  [`col_w_ks()`](https://ngreifer.github.io/cobalt/reference/balance-summary.md)
  [`col_w_ovl()`](https://ngreifer.github.io/cobalt/reference/balance-summary.md)
  [`col_w_cov()`](https://ngreifer.github.io/cobalt/reference/balance-summary.md)
  [`col_w_corr()`](https://ngreifer.github.io/cobalt/reference/balance-summary.md)
  [`col_w_dcov()`](https://ngreifer.github.io/cobalt/reference/balance-summary.md)
  [`col_w_dcorr()`](https://ngreifer.github.io/cobalt/reference/balance-summary.md)
  : Compute Balance and Summary Statistics for Covariates

- [`bal.compute()`](https://ngreifer.github.io/cobalt/reference/bal.compute.md)
  [`bal.init()`](https://ngreifer.github.io/cobalt/reference/bal.compute.md)
  [`available.stats()`](https://ngreifer.github.io/cobalt/reference/bal.compute.md)
  : Efficiently compute scalar balance statistics

- [`f.build()`](https://ngreifer.github.io/cobalt/reference/f.build.md)
  : Convenient Formula Generation

- [`get.w()`](https://ngreifer.github.io/cobalt/reference/get.w.md) :
  Extract Weights from Preprocessing Objects

- [`print(`*`<bal.tab>`*`)`](https://ngreifer.github.io/cobalt/reference/print.bal.tab.md)
  :

  Print Results of a Call to
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)

- [`set.cobalt.options()`](https://ngreifer.github.io/cobalt/reference/set.cobalt.options.md)
  [`get.cobalt.options()`](https://ngreifer.github.io/cobalt/reference/set.cobalt.options.md)
  :

  Set and Get Options in `cobalt`

- [`splitfactor()`](https://ngreifer.github.io/cobalt/reference/splitfactor.md)
  [`unsplitfactor()`](https://ngreifer.github.io/cobalt/reference/splitfactor.md)
  : Split and Unsplit Factors into Dummy Variables

- [`var.names()`](https://ngreifer.github.io/cobalt/reference/var.names.md)
  :

  Extract Variable Names from `bal.tab` Objects

## Data

- [`lalonde`](https://ngreifer.github.io/cobalt/reference/lalonde.md)
  [`lalonde_mis`](https://ngreifer.github.io/cobalt/reference/lalonde.md)
  : Data from National Supported Work Demonstration and PSID, as
  analyzed by Dehejia and Wahba (1999).

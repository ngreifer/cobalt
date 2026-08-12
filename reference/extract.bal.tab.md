# Extract a Balance Table for Further Use

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) extracts
the balance statistics computed by
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md) as
a tidy data frame, one row per covariate, sample, and statistic.
[`format()`](https://rdrr.io/r/base/format.html) returns the balance
table exactly as
[`print.bal.tab()`](https://ngreifer.github.io/cobalt/reference/print.bal.tab.md)
displays it, as a data frame of formatted strings, ready to be passed to
a table-rendering function.

Together they are meant to remove the need to pick apart a `bal.tab`
object by hand when reporting balance in a document.

## Usage

``` r
# S3 method for class 'bal.tab'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  var.names = NULL,
  wide = FALSE
)

# S3 method for class 'bal.tab'
format(
  x,
  ...,
  var.names = NULL,
  digits = max(3L, getOption("digits") - 3L),
  component = "balance"
)
```

## Arguments

- x:

  a `bal.tab` object; the output of a call to
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md).

- row.names, optional:

  ignored; present for consistency with the
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  generic.

- ...:

  arguments passed to
  [`print.bal.tab()`](https://ngreifer.github.io/cobalt/reference/print.bal.tab.md)
  to control which statistics, samples, and covariates are included,
  e.g. `stats`, `disp`, `un`, `imbalanced.only`, or `disp.thresholds`.
  The `.all` and `.none` shorthands are accepted as they are by
  [`print()`](https://rdrr.io/r/base/print.html). Arguments that were
  not computed in the original call to
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  cannot be requested here, for the same reason they cannot be requested
  in [`print()`](https://rdrr.io/r/base/print.html).

- var.names:

  an optional object providing alternate names for the variables, which
  will otherwise be returned as they are stored. Entries given here add
  to those given in the original call to
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  and replace any entry they name. See
  [`display-options`](https://ngreifer.github.io/cobalt/reference/display-options.md)
  for how to specify it.

- wide:

  `logical`; for
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html),
  whether to return the table in the layout
  [`print()`](https://rdrr.io/r/base/print.html) uses, with one column
  per sample and statistic, rather than the default tidy layout. Default
  is `FALSE`.

- digits:

  for [`format()`](https://rdrr.io/r/base/format.html), the number of
  significant digits to display. Default is the same as for
  [`print()`](https://rdrr.io/r/base/print.html).

- component:

  for [`format()`](https://rdrr.io/r/base/format.html), which table to
  return: `"balance"` (the default) for the balance table, or
  `"observations"` for the sample size table.

## Value

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns a
data frame with one row per covariate, sample, and statistic, and these
columns:

- `variable`:

  the covariate name, as it appears in the row names of the balance
  table.

- `type`:

  the covariate type: `"Binary"`, `"Contin."`, or `"Distance"`.

- `sample`:

  `"Unadjusted"`, or the name of the set of weights.

- `stat`:

  the name of the statistic, using the same names as
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)'s
  `stats` argument (e.g., `"mean.diffs"`), with `"mean"` and `"sd"` for
  the distribution summary statistics.

- `group`:

  for `"mean"` and `"sd"`, the name of the treatment group the value
  describes, as
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  names it: the two group names for a binary treatment, the level for a
  multi-category one, `"Uncensored"` or `"Full"` for a censoring
  indicator, and `"All"` for a continuous treatment, which has no
  groups, or for the full sample when `pairwise = FALSE`. `NA` for a
  statistic that contrasts two groups, which belongs to neither of them.

- `estimate`:

  the value of the statistic.

- `threshold`:

  the balance verdict, e.g. `"Balanced, <0.1"`, when a threshold was
  requested for that statistic; `NA` for a row it does not cover, such
  as a mean.

- `threshold.value`:

  the numeric threshold, on the same rows.

The two threshold columns are present only when a threshold is on
display for at least one statistic; when none is, they would be empty
throughout and are omitted. The rest of the columns are always present.

When the data are segmented – by cluster, imputation, treatment pair,
time point, or subclass – one further column per level of segmentation
identifies it. Segmentation is always a column, never a nested list, so
the result is a single rectangle whatever the shape of the input.

A multi-category treatment is reported one pair of groups at a time, but
a mean or a standard deviation belongs to a group rather than to a
comparison, and is the same in every pair that group appears in. Such a
row therefore appears once, with `pair` set to `NA`; only the statistics
that contrast two groups carry a `pair`. The same applies to the full
sample's own means when `pairwise = FALSE`, which would otherwise be
repeated against every group.

With `wide = TRUE`, the columns are those
[`print()`](https://rdrr.io/r/base/print.html) displays, with the
covariate names moved from the row names into a `variable` column and
any segmentation columns retained.

[`format()`](https://rdrr.io/r/base/format.html) returns a data frame of
character vectors with the covariate names as row names, formatted
exactly as [`print()`](https://rdrr.io/r/base/print.html) displays them:
rounded to `digits`, padded to a common number of decimal places, with
`NA` shown as `"."`.

## Details

Both functions accept every argument
[`print()`](https://rdrr.io/r/base/print.html) accepts, and resolve them
the same way, so the same warnings are raised when a requested value was
not computed because `quick = TRUE` in the original call to
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md).

The `variable` column and the row names carry the covariates' display
names, so `var.names` – given here or in the original call to
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md) –
changes them. What the covariates are stored under is unaffected. See
[`display-options`](https://ngreifer.github.io/cobalt/reference/display-options.md).

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the balance statistics themselves, from each innermost balance table. It
does not return the summaries across clusters, imputations, treatment
pairs, or time points, which are aggregates of those statistics;
[`format()`](https://rdrr.io/r/base/format.html) returns the summary
when that is what [`print()`](https://rdrr.io/r/base/print.html)
displays.

## See also

- [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)

- [`print.bal.tab()`](https://ngreifer.github.io/cobalt/reference/print.bal.tab.md)

- [`love.plot()`](https://ngreifer.github.io/cobalt/reference/love.plot.md)
  for a graphical alternative

- [`display-options`](https://ngreifer.github.io/cobalt/reference/display-options.md)
  for `var.names`

## Examples

``` r
data("lalonde", package = "cobalt")

b <- bal.tab(treat ~ age + educ + race + re74, data = lalonde,
             s.d.denom = "pooled", stats = c("m", "ks"),
             thresholds = c(m = .1), un = TRUE)

#Tidy: one row per covariate, sample, and statistic
head(as.data.frame(b))
#>     variable    type     sample          stat group    estimate
#> 1        age Contin. Unadjusted    mean.diffs  <NA> -0.24190362
#> 2        age Contin. Unadjusted ks.statistics  <NA>  0.15772696
#> 3       educ Contin. Unadjusted    mean.diffs  <NA>  0.04475509
#> 4       educ Contin. Unadjusted ks.statistics  <NA>  0.11137151
#> 5 race_black  Binary Unadjusted    mean.diffs  <NA>  0.64044604
#> 6 race_black  Binary Unadjusted ks.statistics  <NA>  0.64044604
#>            threshold threshold.value
#> 1 Not Balanced, >0.1             0.1
#> 2               <NA>              NA
#> 3     Balanced, <0.1             0.1
#> 4               <NA>              NA
#> 5 Not Balanced, >0.1             0.1
#> 6               <NA>              NA

#The layout print() shows
as.data.frame(b, wide = TRUE)
#>      variable    Type     Diff.Un     M.Threshold.Un      KS.Un
#> 1         age Contin. -0.24190362 Not Balanced, >0.1 0.15772696
#> 2        educ Contin.  0.04475509     Balanced, <0.1 0.11137151
#> 3  race_black  Binary  0.64044604 Not Balanced, >0.1 0.64044604
#> 4 race_hispan  Binary -0.08273168     Balanced, <0.1 0.08273168
#> 5  race_white  Binary -0.55771436 Not Balanced, >0.1 0.55771436
#> 6        re74 Contin. -0.59575159 Not Balanced, >0.1 0.44703585

#Ready for knitr::kable() or any other table renderer
format(b)
#>                Type Diff.Un     M.Threshold.Un  KS.Un
#> age         Contin. -0.2419 Not Balanced, >0.1 0.1577
#> educ        Contin.  0.0448     Balanced, <0.1 0.1114
#> race_black   Binary  0.6404 Not Balanced, >0.1 0.6404
#> race_hispan  Binary -0.0827     Balanced, <0.1 0.0827
#> race_white   Binary -0.5577 Not Balanced, >0.1 0.5577
#> re74        Contin. -0.5958 Not Balanced, >0.1 0.4470

format(b, component = "observations")
#>     Control Treated
#> All     429     185

#print()'s arguments work here too
as.data.frame(b, stats = "ks", un = FALSE)
#>      variable    type     sample          stat group   estimate
#> 1         age Contin. Unadjusted ks.statistics  <NA> 0.15772696
#> 2        educ Contin. Unadjusted ks.statistics  <NA> 0.11137151
#> 3  race_black  Binary Unadjusted ks.statistics  <NA> 0.64044604
#> 4 race_hispan  Binary Unadjusted ks.statistics  <NA> 0.08273168
#> 5  race_white  Binary Unadjusted ks.statistics  <NA> 0.55771436
#> 6        re74 Contin. Unadjusted ks.statistics  <NA> 0.44703585
```

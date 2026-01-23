# Using `bal.tab()` with Multi-Category Treatments

When using
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
with multi-category treatments, the output will be different from the
case with binary or continuous treatments, and there are some options
that are common across all
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
methods. This page outlines the outputs and options in this case.

There are two main components of the output of
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
with multi-category treatments: the two-group treatment comparisons and
the balance summary. The two-group treatment comparisons are standard
binary treatment comparison either for pairs of groups (e.g., for
treatments A, B, and C, "A vs. B", "A vs. C", and "B vs. C") or each
group against all the groups (i.e., the entire sample).

The balance summary is, for each variable, the greatest imbalance across
all two-group comparisons. So, for variable X1, if "A vs. B" had a
standardized mean difference of 0.52, "A vs. C" had a standardized mean
difference of .17, and "B vs. C" had a standardized mean difference of
.35, the balance summary would have 0.52 for the value of the
standardized mean difference for X1. The same goes for other variables
and other measures of balance. If the greatest observed imbalance is
tolerable, then all other imbalances for that variable will be tolerable
too, so focusing on reducing the greatest imbalance is sufficient for
reducing imbalance overall. (Note that when `s.d.denom = "pooled"`,
i.e., when the estimand is the ATE, the pooled standard deviation in the
denominator will be the average of the standard deviations across all
treatment groups, not just those used in the pairwise comparison.) The
balance summary will not be computed if multiply imputed data are used.

## Note

In versions 4.3.1 and earlier, setting `pairwise = FALSE` would compare
each group to the full adjusted sample. Now, each group is compared to
the full *un*adjusted sample (unadjusted except for `s.weights`, if
supplied).

In versions 4.3.1 and earlier, `pairwise` was ignored with non-`NULL`
`focal` and was automatically set to `FALSE`. `pairwise` can be
specified and its default is now `TRUE`, so balance between all
treatment groups will be computed by default rather than only between
each non-group and the focal group. To recover previous functionality,
set `pairwise = FALSE` with non-`NULL` `focal`.

## Allowable arguments

There are four arguments for each
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
method that can handle multi-category treatments: `pairwise`, `focal`,
`which.treat`, and `multi.summary`.

- `pairwise`:

  Whether to compute the two-group comparisons pairwise or not. If
  `TRUE`,
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  will compute comparisons for each pair of treatments. This can be
  valuable if treatments are to be compared with one another (which is
  often the case). If `FALSE`,
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  will compute balance for each treatment group against the full
  unadjusted sample when `focal` is `NULL` and for each non-focal group
  against the focal group otherwise.

- `focal`:

  When one group is to be compared to multiple control groups in an ATT
  analysis, the group considered "treated" is the focal group. By
  specifying the name or index of the treatment condition considered
  focal,
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  will only compute and display pairwise balance for treatment
  comparisons that include the focal group when `pairwise = FALSE`.

- `which.treat`:

  This is a display option that does not affect computation. When
  displaying the `bal.tab` output, which treatments should be displayed?
  If a vector of length 1 is entered, all comparisons involving that
  treatment group will be displayed. If a vector of length 2 or more is
  entered, all comparisons involving treatments that both appear in the
  input will be displayed. For example, inputting `"A"` will display "A
  vs. B" and "A vs. C", while entering `c("A", "B")` will only display
  "A vs. B". `.none` indicates no treatment comparisons will be
  displayed, and `.all` indicates all treatment comparisons will be
  displayed. `.none` is the default.

- `multi.summary`:

  If `TRUE`, the balance summary across all comparisons will be computed
  and displayed. This includes one row for each covariate with maximum
  balance statistic across all pairwise comparisons. Note that, if
  variance ratios or KS statistics are requested in addition to mean
  differences, the displayed values may not come from the same pairwise
  comparisons; that is, the greatest standardized mean difference and
  the greatest variance ratio may not come from the same comparison. The
  default is `TRUE`, and if `which.treat` is `.none`, it will
  automatically be set to `TRUE`.

## Output

The output is a `bal.tab.multi` object, which inherits from `bal.tab`.
It has the following elements:

- `Pair.Balance`:For each pair of treatment groups, a regular `bal.tab`
  object containing a balance table, a sample size summary, and other
  balance assessment tools, depending on which options are specified. If
  `pairwise` is `FALSE`, the comparisons will be between each group and
  the groups combined (labeled "All") when `focal` is `NULL` and between
  each non-focal group and the focal group otherwise.

- `Balance.Across.Pairs`: The balance summary across two-group
  comparisons. This will include the greatest (i.e., maximum) absolute
  balance statistics(s) for each covariate across all comparisons
  computed. Thresholds can be requested for each balance measure as with
  binary treatments.

- `Observations`: A table of sample sizes or effective sample sizes for
  each treatment group before and after adjustment.

As with other methods, multiple weights can be specified, and values for
all weights will appear in all tables.

## See also

- [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)

- [`bal.tab.data.frame()`](https://ngreifer.github.io/cobalt/reference/bal.tab.formula.md)

- [`print.bal.tab()`](https://ngreifer.github.io/cobalt/reference/print.bal.tab.md)

- [`vignette("segmented-data")`](https://ngreifer.github.io/cobalt/articles/segmented-data.md)
  for examples

# Using `bal.tab()` with Subclassified Data

When using
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
with subclassified data, i.e., data split into subclasses where balance
may hold, the output will be different from the standard,
non-subclassified case, and there is an additional option for
controlling display. This page outlines the outputs and options in this
case.

There are two main components of the output of
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
with subclassified data: the balance within subclasses and the balance
summary across subclasses. The within-subclass balance displays
essentially are standard balance displays for each subclass, except that
only "adjusted" values are available, because the subclassification
itself is the adjustment.

The balance summary is, for each variable, like a weighted average of
the balance statistics across subclasses. This is computed internally by
assigning each individual a weight based on their subclass and treatment
group membership and then computing weighted balance statistics as usual
with these weights. This summary is the same one would get if subclasses
were supplied to the `match.strata` argument rather than to `subclass`.
Because the means and mean differences are additive, their computed
values will be weighted averages of the subclass-specific values, but
for other statistics, the computed values will not be.

## Allowable arguments

There are three arguments for
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
that relate to subclasses: `subclass`, `which.subclass`, and
`subclass.summary`.

- `subclass`:

  For the `data.frame` and formula methods of
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md),
  a vector of subclass membership or the name of the variable in `data`
  containing subclass membership. When using subclassification with a
  function compatible with cobalt, such as
  [`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.html)
  in MatchIt, this argument can be omitted because the subclasses are in
  the output object.

- `which.subclass`:

  This is a display option that does not affect computation. If `.all`,
  all subclasses in `subclass` will be displayed. If `.none` (the
  default), no subclasses will be displayed. Otherwise, can be a vector
  of subclass indices for which to display balance.

- `subclass.summary`:

  This is a display option that does not affect computation. If `TRUE`,
  the balance summary across subclasses will be displayed. The default
  is `TRUE`, and if `which.subclass` is `.none`, it will automatically
  be set to `TRUE`.

## Output

The output is a `bal.tab.subclass` object, which inherits from
`bal.tab`. It has the following elements:

- `Subclass.Balance`: A list of data frames containing balance
  information for each covariate in each subclass.

- `Balance.Across.Subclass`: A data frame containing balance statistics
  for each covariate aggregated across subclasses and for the original
  sample (i.e., unadjusted). See
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  for details on what this includes.

- `Observations`: A table of sample sizes in each subclass and overall.

## See also

- [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)

- [`bal.tab.data.frame()`](https://ngreifer.github.io/cobalt/reference/bal.tab.formula.md)

- [`print.bal.tab()`](https://ngreifer.github.io/cobalt/reference/print.bal.tab.md)

# Using `bal.tab()` with Longitudinal Treatments

When using
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
with longitudinal treatments, the output will be different from the case
with point treatments, and there are some options that are common across
all
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
methods for dealing with longitudinal data. This page outlines the
outputs and options in this case.

There are two main components of the output of
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
with longitudinal treatments: the time-point-specific balance summary
and across-time-points balance summary. The time-point-specific balance
summaries are standard point treatment balance summaries at each time
point.

The across-time-points balance summary is, for each variable, the
greatest imbalance across all time-point-specific balance summaries. If
the greatest observed imbalance is tolerable, then all other imbalances
for that variable will be tolerable too, so focusing on reducing the
greatest imbalance is sufficient for reducing imbalance overall. The
balance summary will not be computed if multi-category treatments or
multiply imputed data are used, or if the time points are not all of the
same kind (see below).

## Note

The balance tables presented here are not the same as those recommended
by Jackson (2016) and computed in his R package,
[confoundr](https://CRAN.R-project.org/package=confoundr), as these do
not take into account treatment history. The balance statistics
presented here should be used with caution and may not reflect balance
in an accurate way.

## Mixing treatments and censoring

One entry of the list may be a censoring indicator marked with
[`.cens()`](https://ngreifer.github.io/cobalt/reference/cens.md) rather
than a treatment, as in `list(A1 ~ x, .cens(C1) ~ x, A2 ~ x)`, which is
how a joint treatment-and-censoring model is written for
[`WeightIt::weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.html).
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
produces one table per entry, of whichever kind that entry is: an
ordinary balance table for a treatment, and the censoring balance table
described at
[`class-bal.tab.cens`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.cens.md)
for an indicator.

Each entry is assessed among the units still under observation entering
it. A unit is under observation until a censoring indicator earlier in
the list marks it censored; the risk set is accumulated from the
indicators themselves rather than read off the treatments, so it makes
no difference whether the data records a treatment for a unit that has
already dropped out or leaves it missing. A missing treatment for a unit
that is still under observation is an error naming the time point it
appeared in. A censoring entry's own comparison uses the risk set as it
stands entering it, so it includes the units it is about to remove;
those are the full sample it compares against.

A censoring balance table and a treatment balance table say different
things about different samples, so a list that mixes them gets no
balance summary across time points, exactly as a list mixing continuous
and binary treatments does not. Each time point's own table still
reports its sample sizes. A list in which every entry is a censoring
indicator is not a mixture and is summarized as usual.

Each time point's default `s.d.denom` is the one its own kind of model
implies – `"pooled"` for a binary or multi-category treatment, `"all"`
for a continuous one, and `"full"` for a censoring indicator – so that a
model gives the same numbers in a list as it would on its own. A value
supplied to `s.d.denom` is shared by every time point, which is one
reason it is recommended not to set it for longitudinal treatments.

## Naming the time points

Each time point is named for its position in the list, whether the model
there is a treatment or a censoring model, and the variable it is about,
as in `1. Treatment: A_1` or `2. Censoring: C_1`. The number is the
position in the list rather than the treatment period, so that it is the
number `which.time` takes and a censoring model has one too;
`which.time` also accepts the variable name on its own (e.g., `"C_1"`).

## Allowable arguments

There are two additional arguments for each
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
method that can handle longitudinal treatments: `which.time` and
`msm.summary`.

- `which.time`:

  This is a display option that does not affect computation. If `.all`
  (the default), all time points will be displayed. If `.none`, no time
  points will be displayed. Otherwise, can be a vector of treatment
  names or indices for which to display balance.

- `msm.summary`:

  This is a display option that does not affect computation. If `TRUE`,
  the balance summary across time points will be displayed. The default
  is `TRUE`, and if `which.time` is `.none`, it will automatically be
  set to `TRUE`.

## Output

The output is a `bal.tab.msm` object, which inherits from `bal.tab`. It
has the following elements:

- `Time.Balance`: For each time point, a regular `bal.tab` object
  containing a balance table, a sample size summary, and other balance
  assessment tools, depending on which options are specified.

- `Balance.Across.Times`: The balance summary across time points. This
  will include the maximum balance statistic(s) for each covariate
  across all time points. Absent when the time points are not all of the
  same kind.

- `Observations`: A list with a table of sample sizes or effective
  sample sizes for each time point, before and after adjustment. Always
  present, but displayed only alongside the balance summary across time
  points: it gathers in one place what each time point's own table has
  already reported, so without the summary there is nothing it adds.

As with other methods, multiple weights can be specified, and values for
all weights will appear in all tables.

## References

Jackson, J. W. (2016). Diagnostics for Confounding of Time-varying and
Other Joint Exposures: *Epidemiology*, 27(6), 859–869.
[doi:10.1097/EDE.0000000000000547](https://doi.org/10.1097/EDE.0000000000000547)

## See also

- [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)

- [`bal.tab.time.list()`](https://ngreifer.github.io/cobalt/reference/bal.tab.time.list.md)

- [`print.bal.tab()`](https://ngreifer.github.io/cobalt/reference/print.bal.tab.md)

- [`.cens()`](https://ngreifer.github.io/cobalt/reference/cens.md) and
  [`class-bal.tab.cens`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.cens.md)
  for a censoring indicator among the time points

- [`vignette("longitudinal-treat")`](https://ngreifer.github.io/cobalt/articles/longitudinal-treat.md)
  for examples

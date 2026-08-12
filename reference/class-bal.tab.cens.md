# Using `bal.tab()` with a Censoring Indicator

When
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md) is
given a censoring indicator marked with
[`.cens()`](https://ngreifer.github.io/cobalt/reference/cens.md) rather
than a treatment, its target is the full at-risk sample, and what
matters is whether the units still under observation resemble that
sample once weighted. This page outlines the output in this case.

Balance is therefore assessed between two samples built from the same
units:

- **Uncensored**: the units with `C == 0`, carrying the censoring
  weights.

- **Full**: every at-risk unit, i.e., every unit with a non-missing `C`,
  carrying a weight of 1.

Setting `un = TRUE` adds the same comparison with the uncensored units
unweighted, which is the imbalance the weights were estimated to remove.

The two samples are compared using the same statistics available for
binary treatments (`mean.diffs`, `variance.ratios`, `ks.statistics`, and
`ovl.coefficients`), because internally they are stacked into a binary
comparison. The two samples name their own columns rather than taking a
binary treatment's positional `0` and `1`, so the balance table has
`M.Uncensored` and `SD.Uncensored` for the uncensored sample and
`M.Full` and `SD.Full` for the full one, and `Diff` is the difference
between the full sample and the uncensored one.

[`bal.plot()`](https://ngreifer.github.io/cobalt/reference/bal.plot.md)
works too, showing the weighted uncensored sample against the unweighted
full sample.

## Allowable arguments

Every argument to
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
applies as it does with a binary treatment, with the following
exceptions.

- `s.d.denom`:

  allowable values are `"full"` (the default), which uses the standard
  deviation of the covariate in the full at-risk sample; `"uncensored"`,
  which uses that of the uncensored sample; and `"pooled"`, `"all"`,
  `"weighted"`, and `"hedges"`, which behave as they do for a binary
  treatment.

- `estimand` and `focal`:

  do not apply and are ignored with a warning. The target is the full
  at-risk sample.

- `subclass`:

  applies: subclassifying is itself a way of solving a censoring
  problem, in that within each subclass the uncensored units should
  resemble every at-risk unit in it. See "With subclasses" below.

- `match.strata`:

  does not apply.

`cluster` and `imp` apply as usual, and produce a `bal.tab.cluster` or
`bal.tab.imp` object whose per-cluster or per-imputation components are
the censoring balance tables described here.

## Among longitudinal treatments

A censoring indicator can appear among a list of longitudinal
treatments, as in `list(A1 ~ x, .cens(C1) ~ x, A2 ~ x)`, which is how a
joint treatment-and-censoring model is written for
[`WeightIt::weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.html);
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
also accepts such a `weightitMSM` object directly. Each entry of the
list gets a table of its own kind, and each is assessed among the units
still under observation entering it, so the full sample a censoring
indicator is compared against is the risk set at that time point rather
than the original cohort, and a treatment after it is assessed only
among the units it did not censor. A list that mixes censoring with
treatment gets no balance summary across time points. See
[`class-bal.tab.msm`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.msm.md).

## With subclasses

Subclassification is an alternative to weighting for solving a censoring
problem: within each subclass, the units still under observation should
resemble every at-risk unit in that subclass. Supplying `subclass`
therefore produces a `bal.tab.cens` object that also inherits from
`bal.tab.subclass`, with a balance table for each subclass and a summary
across them, as described at
[`class-bal.tab.subclass`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.subclass.md).

The summary across subclasses is subclassification expressed as
censoring weights: a unit still under observation in subclass \\k\\
receives \\n_k / n\_{k1}\\, where \\n_k\\ is the number of at-risk units
in the subclass and \\n\_{k1}\\ the number of them still under
observation, and the full sample is left unweighted. This is the same
summary one would get from supplying those weights to a censoring model
directly.

The sample sizes have one column per subclass, with `Full`,
`Uncensored`, and `Censored` rows.

## Output

The output is a `bal.tab.cens` object, which inherits from `bal.tab.bin`
and `bal.tab` (or from `bal.tab.subclass` when `subclass` is supplied),
and has the same elements as an ordinary binary-treatment `bal.tab`
object. Its `Observations` component differs: it has a single `Total`
column, because the target sample is the same in every row, and the rows

- `Full`: the number (or effective number) of at-risk units,

- `Uncensored`: the number of units still under observation,

- `Adjusted` (or one row per set of weights): their effective number
  once weighted,

- `Censored`: the number of units censored, omitted when there are none.

`Uncensored` and `Censored` sum to `Full`.

## See also

- [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)

- [`.cens()`](https://ngreifer.github.io/cobalt/reference/cens.md) for
  marking an indicator as censoring

- [`class-bal.tab.cluster`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.cluster.md),
  [`class-bal.tab.imp`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.imp.md),
  [`class-bal.tab.subclass`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.subclass.md),
  and
  [`class-bal.tab.msm`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.msm.md)
  for the segmented and longitudinal cases

- [`treat`](https://ngreifer.github.io/cobalt/reference/treat-class.md)
  for the attributes that decide how a treatment is compared and what
  its groups are called

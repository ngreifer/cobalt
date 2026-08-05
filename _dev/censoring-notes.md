# Censoring balance: the design, and why

Branch `censoring-balance`, off `master` at `ea39e53`. Phase 1: a censoring model on its
own. Combining censoring with treatment balance across time points is not done.

## What is being assessed

A censoring model asks a different question from a treatment model. There is no second
treatment group to compare against — the target is the **full at-risk sample**, and the
question is whether the units still under observation, once weighted, look like it.

So balance is assessed between two samples built from the same data:

| | rows | weight |
|---|---|---|
| **Uncensored** | the units with `C == 0` | the censoring weights |
| **Full** | every at-risk unit (`C` not missing) | 1 |

`un = TRUE` gives the same comparison with the uncensored block unweighted, which is
exactly "how far off was the observed sample before weighting".

## The one structural decision

**The two samples are stacked into a genuine binary pseudo-treatment, and everything
downstream is the ordinary binary machinery.**

`X$treat` holds the censoring indicator all the way through `x2base()` and
`.assign_X_class()`, so `subset_X()`, `process_cluster()`, and the imputation expansion
all still see one row per unit. The stacking happens once, at the leaf, in
`base.bal.tab.cens()`, which then hands a `type = "bin"` `X` to `base.bal.tab.base()`.

Consequences, all of them wanted:

- **Every binary statistic works unchanged** — `mean.diffs`, `variance.ratios`,
  `ks.statistics`, `ovl.coefficients` — because the stacked treatment really is binary.
  No new `STATS` entries, no new `type`, and no third branch in `balance_table()`.
- `print()`, `love.plot()`, `as.data.frame()`, and `format()` need no new methods:
  `bal.tab_print` finds no `bal.tab.cens` method and falls through to `bal.tab.bin`'s.
- `cluster` and `imp` compose for free, because `cens` is a *leaf* class ranked below
  them in `.assign_X_class()`.

The alternative — a `type = "cens"` with its own column spec so the moment columns could
read `M.Uncens`/`M.Full` instead of `M.0`/`M.1` — was rejected. It buys two nicer column
names and costs a third branch in `balance_table()` and `balance_table_subclass()`, new
`STATS` types, and a fork in every place that switches on `type`. The correspondence
(`M.0` is the uncensored sample, `M.1` the full one) is documented on
`?class-bal.tab.cens` instead.

**The `.target` statistics are not used and were not touched.** They do their own 2n
doubling inside `$fun` for *continuous* treatments and hardcode `s.d.denom <- "1"`;
nothing about them fits here.

## Which group is which, and why it matters

`0` (control) is the **uncensored** sample; `1` (treated) is the **full** sample. So
`Diff` is `full − uncensored`, and `s.d.denom` defaults to the full sample.

This is not arbitrary — it is exactly the stacking WeightIt's `?.cens` documents as the
manual workaround:

```r
u <- which(W$treat == 0)
bal.tab(rbind(covs[u, ], covs),
        treat = c(rep(0L, length(u)), rep(1L, nrow(covs))),
        weights = c(W$weights[u], rep(1, nrow(covs))),
        estimand = "ATT", s.d.denom = "treated")
```

Matching it means the numbers agree with what users have been told to compute, and it
gives the test suite an **independent oracle**: `test-censoring.R` defines
`cens_manual()`, which is that recipe run through cobalt's ordinary binary path with no
censoring code involved, and asserts `expect_equal(b$Balance, manual$Balance)` over all
four statistics and both moments, for the plain case, with `s.weights`, with `NA`s in
`C`, per cluster, per imputation, and on a `weightit` object. Every number in the feature
is pinned to a statement of what the comparison *is*, not to what this code produced.

It also lines up with cobalt's ATT convention — the focal/target group is the treated
one, the denominator comes from it, and the reweighted group is the control.

## Things to remember

1. **`subset_X()` is what does the stacking.** It takes indices and nothing stops an
   index from repeating, so `subset_X(X, c(uncensored, at.risk))` stacks every per-unit
   slot the same way, including slots added later. Do not replace this with hand-written
   `rbind()`s per slot; that is how `distance`, `addl`, `cluster`, and `discarded` get
   silently left out of step.
2. **`.cens_s.d.denom()` must be idempotent.** A wrapper (`cluster`, `imp`) resolves
   `s.d.denom` once via `.resolve_s.d.denom()` and hands the result down to every child,
   so whatever it returns has to survive being resolved again. That is why it stays in
   the censoring vocabulary (`"full"`, `"uncensored"`, …) and the leaf translates to
   `"treated"`/`"control"` only at the point of use. Returning `"treated"` from the
   wrapper would make the child's own `match.arg()` reject it.
3. **cobalt does not define `[.treat`.** `.cens()` returns class `treat` to match
   WeightIt's contract, but registering an S3 method on that class in both packages
   would collide. `process_treat()` converts the tag to a `processed.treat` immediately
   and `subset_processed.treat()` carries it — with an early return, because the ordinary
   path narrows `treat_names` and re-runs `assign.treat.type()`, which would turn a
   censoring indicator back into a binary treatment.
4. **A censoring `processed.treat` has no `treat_names`/`treat_vals`**, like a continuous
   one. The two named samples belong to the *stacked* pseudo-treatment, built by
   `.cens_pseudo_treat()`.
5. **`process_stats()` maps censoring onto `type = "bin"`.** That one line is what makes
   the binary statistic set and the binary column layout apply.
6. **`.is_cens_formula()` matches the `::` call form too**, so a formula written as
   `WeightIt::.cens(C) ~ x` works. The marker is stripped, not evaluated, in both
   `get_treat_from_formula()` and `get_covs_from_formula()`, so the treatment keeps the
   indicator's own name (`C`, not `.cens(C)`).
7. **`base.bal.tab.base()` gained a `.obs` argument.** The censoring leaf reshapes `X`
   before handing it over, so the sample sizes cannot be counted from what arrives; they
   describe the units, not the stack. Computing and discarding a wrong table would have
   worked but is exactly the kind of thing that later reads as intentional.
8. **`.cluster_check()` needed its own censoring branch.** The binary check requires
   every treatment level in every cluster, which for censoring would reject a cluster in
   which nobody happened to be censored. The real requirement is at least one
   *uncensored* unit per cluster.
9. **`C` may be `NA`**, and this is the only treatment type for which cobalt allows a
   missing treatment value. In a longitudinal model a unit censored at an earlier time
   has no later indicator. Such units are in neither sample: the at-risk sample is the
   units with a non-missing `C`.

## Sample sizes

One column, because the target is the same in every row:

```
Sample sizes
            Total
Full          614
Uncensored    322
Adjusted    308.3
Censored      292
```

`Uncensored + Censored == Full`. One row per set of weights when there is more than one.
The `Censored` row is dropped when nothing is censored, following the `Discarded`
precedent.

## Gates

The suite is **2,517 passed / 0 failed / 12 skipped** (2,429 before this work; 88 new
assertions in `test-censoring.R`). The golden set is **129 identical / 0 differing** and
`check_col_spec()` **323 tables / 0 mismatched**, with seven new cells: `cens_default`,
`cens_all_stats`, `cens_s_weights`, `cens_cluster`, `cens_imp`, and the two
`obj_weightit_cens*` cells the per-fixture loop generates from the new
`weightit_cens` fixture.

`R CMD check --as-cran` is clean. Note from the previous branch: **`cem` must not be
installed** or it hangs the build; see `_dev/refactor-notes.md` § "`R CMD check`".

## Not in this phase

- **Censoring combined with treatment across time points**, via the MSM architecture.
  The shape is different from what is here: a joint treatment-and-censoring
  `weightitMSM` wants ordinary treatment balance at each time point restricted to that
  time's `at.risk` column — which is what removes the missing treatments that currently
  error — *plus* a censoring table at each censoring time point. `.assign_X_class()`
  ranks `msm` above the leaves, so the per-time-point `X` is where `cens` would appear;
  `base.bal.tab.msm()` already rebuilds one `X` per time point and would need to choose
  the leaf per time point rather than assuming one type throughout.
- **`bal.plot()`**, which errors with a message saying so. The two samples share units —
  the uncensored ones appear in both — and nothing in its faceting is built for that.
- **Subclassification**, rejected: it is not a way of estimating censoring weights.

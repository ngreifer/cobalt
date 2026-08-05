# Censoring balance: the design, and why

Branch `censoring-balance`, off `master` at `ea39e53`. Phase 1: a censoring model on its
own, with or without subclassification, in `bal.tab()` and `bal.plot()`. Combining
censoring with treatment balance across time points is not done.

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
`.stack_cens_X()`, and the result is handed to `base.bal.tab.base()` or
`base.bal.tab.subclass()` with `type = "bin"`.

Consequences, all of them wanted:

- **Every binary statistic works unchanged** — `mean.diffs`, `variance.ratios`,
  `ks.statistics`, `ovl.coefficients` — because the stacked treatment really is binary.
  No new `STATS` entries, no new `type`, and no third branch in `balance_table()`.
- `print()`, `love.plot()`, `as.data.frame()`, and `format()` need no new methods:
  `bal.tab_print` finds no `bal.tab.cens` method and falls through to `bal.tab.bin`'s.
- `cluster`, `imp`, and `subclass` compose for free, because `cens` is a *leaf* class
  ranked below them in `.assign_X_class()`.
- `bal.plot()` needed one branch: stack, then plot the two samples as any binary
  treatment.

A `type = "cens"` with its own column spec was considered and rejected as the way to get
readable moment columns. It would have cost a third branch in `balance_table()` and
`balance_table_subclass()`, new `STATS` types, and a fork in every place that switches on
`type`. The names are instead a property of the treatment — see `group_labels` below — so
the columns read `M.Uncensored`/`M.Full` with no new `type`.

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
four statistics and both moments — for the plain case, with `s.weights`, with `NA`s in
`C`, per cluster, per imputation, per subclass, and on a `weightit` object. Every number
in the feature is pinned to a statement of what the comparison *is*, not to what this code
produced.

It also lines up with cobalt's ATT convention — the focal/target group is the treated one,
the denominator comes from it, and the reweighted group is the control.

## The `treat` class, and where the column names come from

A processed treatment is an object of class `treat` (`R/treat.R`, `?treat-class`), which
is also the class `WeightIt::.cens()` returns. It carries `treat.type`, `treat.name`,
`treat_names`, `treat_vals`, and `group_labels`, and `[` preserves all of them. This
replaced the old `processed.treat` class and gathered `assign.treat.type()`,
`get.treat.type()`, the accessors, `process_treat()`, and `subset_treat()` into one file
from three.

`group_labels` is what makes the censoring columns readable. It is how a treatment's
groups are labelled in a balance table's column names, and it defaults to `c("0", "1")`
because that is cobalt's positional convention for a binary treatment — `M.0` has always
meant the control group's mean, and changing that would break every existing table. The
stacked pseudo-treatment sets it to `c("Uncensored", "Full")` instead, and the label flows
from `.bal_tab_col_spec(group.labels =)` through `print.options$group.labels` so that
`print()` can rebuild the layout of the table in front of it. `treat_labels()` does the
same job for `bal.plot()`'s legend.

Adding `group.labels` to `print.options` changed every object cobalt produces, so the
whole golden set was regenerated. Before doing that, the harness was run with the new
field scrubbed out of both sides: **126 of 129 cells were unchanged**, and the 3 that
differed were exactly the censoring cells that display moment columns.

## Subclassification

Subclassification solves a censoring problem, not only a treatment one: within each
subclass, the units still under observation should resemble every at-risk unit in it.
`base.bal.tab.subclass.cens()` therefore stacks as the plain leaf does and hands the
result to `base.bal.tab.subclass(type = "bin")`.

The one thing it adds is `X$estimand <- "ATT"`. `base.bal.tab.subclass()` builds its
across-subclass summary with `strata2weights()`, which for an ATT gives the focal group a
weight of 1 and every other unit `n_{k,focal} / n_{k,t}`. With the full sample as the
focal group that is `n_k / n_{k,uncensored}` for the uncensored units and 1 for the full
sample — subclassification expressed as censoring weights, which is exactly the summary
wanted. Nothing else was needed, and the test asserts it both against the manual stacking
and against those weights written out by hand.

## `s.d.denom`

Accepts `"full"` (the default), `"uncensored"`, and the four generic values
(`"pooled"`, `"all"`, `"weighted"`, `"hedges"`). `"full"`/`"uncensored"` map onto the
stacked treatment's treated/control groups. Nothing is inferred, so no `note:` is
printed — unlike a treatment model, a censoring model's target is defined by the design.

## Sample sizes

One column, because the target is the same in every row:

```
Sample sizes
            Total
Full          614
Uncensored    492
Adjusted    491.9
Censored      122
```

`Uncensored + Censored == Full`. One row per set of weights when there is more than one.
The `Censored` row is dropped when nothing is censored, following the `Discarded`
precedent. With subclasses the table is transposed — one column per subclass, plus `All` —
with `Full`, `Uncensored`, and `Censored` rows.

## Things to remember

1. **`[.cobalt.treat`, not `[.treat`.** WeightIt 2.0.0 registers its own `[.treat`. Two
   packages registering the same method means whichever loads second overwrites the other
   *and R says so on the console* — on every `library(WeightIt)`, since it imports cobalt.
   Dispatching one class earlier avoids the collision without giving up the shared `treat`
   class; when WeightIt adopts this class the extra layer can go. Relatedly, the class has
   to come **first** in `class(x)`: a binary treatment is a factor underneath, and
   `[.factor` would otherwise win and drop every attribute.
2. **`subset_X()` is what does the stacking.** It takes indices and nothing stops an index
   from repeating, so `subset_X(X, c(uncensored, at.risk))` stacks every per-unit slot the
   same way, including slots added later. Do not replace this with hand-written `rbind()`s
   per slot; that is how `distance`, `addl`, `cluster`, `subclass`, and `discarded` get
   silently left out of step.
3. **`.cens_s.d.denom()` must be idempotent.** A wrapper (`cluster`, `imp`) resolves
   `s.d.denom` once via `.resolve_s.d.denom()` and hands the result down to every child,
   so whatever it returns has to survive being resolved again. That is why it stays in the
   censoring vocabulary (`"full"`, `"uncensored"`, …) and the leaf translates to
   `"treated"`/`"control"` only at the point of use. Returning `"treated"` from the
   wrapper would make the child's own `match.arg()` reject it.
4. **`subset_treat()` needs an early return for censoring.** Its ordinary path narrows
   `treat_names` to the levels still present and re-runs `assign.treat.type()`, which
   would turn a subset with nobody censored in it back into a binary treatment — or error,
   since a single remaining value is not a treatment.
5. **A censoring treatment is a plain 0/1 vector**, not a factor, so that `C == 0` keeps
   working. Its groups are named `Uncensored`/`Censored`; the *stacked* pseudo-treatment
   built by `.cens_pseudo_treat()` is the one whose groups are `Uncensored`/`Full`.
6. **`process_stats()` maps censoring onto `type = "bin"`.** That one line is what makes
   the binary statistic set and the binary column layout apply.
7. **`.is_cens_formula()` matches the `::` call form too**, so a formula written as
   `WeightIt::.cens(C) ~ x` works. The marker is stripped, not evaluated, in both
   `get_treat_from_formula()` and `get_covs_from_formula()`, so the treatment keeps the
   indicator's own name (`C`, not `.cens(C)`).
8. **`base.bal.tab.base()` and `base.bal.tab.subclass()` gained a `.obs` argument.** The
   censoring leaves reshape `X` before handing it over, so the sample sizes cannot be
   counted from what arrives; they describe the units, not the stack. Computing and
   discarding a wrong table would have worked but is exactly the kind of thing that later
   reads as intentional.
9. **`.cluster_check()` needed its own censoring branch.** The binary check requires every
   treatment level in every cluster, which for censoring would reject a cluster in which
   nobody happened to be censored. The real requirement is at least one *uncensored* unit
   per cluster.
10. **`C` may be `NA`**, and this is the only treatment type for which cobalt allows a
    missing treatment value. In a longitudinal model a unit censored at an earlier time
    has no later indicator. Such units are in neither sample: the at-risk sample is the
    units with a non-missing `C`.
11. **`bal.plot()` keeps its data in the layers, one per group**, not in `p$data`, and the
    column naming the panel is `which` for samples and `subclass` for subclasses. The test
    helper `plot_data()` unions them; anything asserting on what a plot shows has to go
    through something like it.

## Gates

The suite is **2,560 passed / 0 failed / 12 skipped** (2,429 before this work; 131 new
assertions in `test-censoring.R`). The golden set is **131 identical / 0 differing** and
`check_col_spec()` **343 tables / 0 mismatched**, with nine new cells: `cens_default`,
`cens_all_stats`, `cens_s_weights`, `cens_cluster`, `cens_imp`, `cens_subclass`,
`cens_subclass_cluster`, and the two `obj_weightit_cens*` cells the per-fixture loop
generates from the new `weightit_cens` fixture.

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
  the leaf per time point rather than assuming one type throughout. `bal.plot()` rejects
  a censoring indicator among longitudinal treatments for the same reason.
- **`match.strata`**, rejected: it turns strata into weights before the two samples exist,
  and the same strata say the same thing given to `subclass`.

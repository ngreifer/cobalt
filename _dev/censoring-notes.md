# Censoring balance: the design, and why

Branch `censoring-balance`, off `master` at `ea39e53`. A censoring model on its own, with
or without subclassification, and interleaved with treatments across time points, in
`bal.tab()` and `bal.plot()`.

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
groups are labeled in a balance table's column names, and it defaults to `c("0", "1")`
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

The table carries no `ss.type`, so no row is starred. That attribute marks which rows of a
*treatment's* table are effective sample sizes, and it pays for itself there because such
a table can hold several sets of weights fit by different methods, only some of them
weighting; the star is what the heading cannot say in that case. Here the only row that
could be an effective sample size is the adjusted one, the heading says so, and the rest
are counts of units — so the star was pure noise, and a treatment's table does not carry
one either.

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

## Across time points

A joint treatment-and-censoring model interleaves the two — `A1 ~ x`, `.cens(C1) ~ x`,
`A2 ~ x` — and the rule is **one table per entry, of whichever kind that entry is**. Most
of that already worked: `.assign_X_class()` is called per time point, so a censoring entry
already reached the censoring leaf, and `print.bal.tab.msm` calls a bare `print()` on each
child, so heterogeneous children already rendered. What was missing was the plumbing.

**The risk set is accumulated from the censoring indicators, not read off the treatments.**
Everyone starts at risk; each censoring entry removes the units censored there from every
entry after it. That is `weightitMSM()`'s own rule, it needs nothing but `treat.list`, and
— the reason it is worth stating — it does not depend on a data convention. `weightitMSM()`
leaves a censored unit's later treatments `NA`, so an `NA`-based rule would look right on
its output and silently include censored units for a user who records those treatments
instead. `NA`s are therefore a *check*: after the restriction, a remaining `NA` in an
ordinary treatment is a unit at risk whose treatment is unknown, and that is an error. A
`test-censoring.R` case asserts that blanking post-censoring treatments and recording them
give the same tables, which is exactly what the discarded rule would have got wrong.

A censoring entry's own risk set still holds the units it is about to remove — they are the
full sample it compares against — so `at.risk[[i]]` is the set *entering* `i`.

**No summary across a mixture.** A censoring table and a treatment table have neither the
same columns nor the same meaning, so a list mixing them gets no `Balance.Across.Times`,
exactly as a list mixing continuous and binary treatments does not. The existing
`all_the_same(treat.types)` gate already produced that and was left alone; an all-censoring
list is not a mixture and is still summarized.

`Observations` is now always *computed*, where it used to come only with the summary, but it
is still only *printed* with it. The collected table says what each time point's own table
already says, gathered in one place, so it belongs with the one other thing that gathers the
time points together; a mixed model has the numbers in its per-time-point tables and does
not need them twice. (An earlier pass here printed it unconditionally, which is what showing
it to Noah corrected.)

Things worth remembering:

12. **`subset_X()` is what restricts a time point**, with the same argument the censoring
    leaf stacks with. When nothing is censored the risk set is all `TRUE` and `subset_X()`
    returns `X` untouched, so a longitudinal model without censoring is bit-identical —
    which is what keeps the existing golden cells from moving.
13. **`process_treat()` gained `.missing.okay`**, and `process_treat.list()` passes `TRUE`
    when any entry is a censoring indicator. Without an indicator in the list a missing
    treatment stays an error. `subset_X()` passes `anyNA(x)`: subsetting never introduces
    a missing value, so a subset may have one exactly when its source did, which is what
    lets `cluster` and `imp` wrap a censored longitudinal model.
14. **The per-time-point `s.d.denom` default reads the treatment type, not `X.class`.** A
    time point wrapped in clusters or imputations has the class of the wrapper, and the
    vocabulary `s.d.denom` is written in belongs to the treatment. The old `X.class`
    switch gave `"pooled"` to a clustered *continuous* longitudinal treatment, which
    `.get_s.d.denom.cont()` rejects outright — that combination simply errored.
15. **A `weightitMSM` comes in two shapes, and both have to work.** WeightIt 2.1.0 treats
    censoring as a treatment type, so `treat.list`/`covs.list` hold every model in the
    order it was fit and there is nothing to do. The development versions between 2.0.0
    and 2.1.0 instead segregated them: `treat.list`/`covs.list` held only the treatment
    models, `cens.list`/`cens.covs.list` only the censoring ones, and `cens.time` recorded
    where each censoring model belonged. `.weightitMSM_models()` puts the segregated shape
    back in order and passes the interleaved one through, so everything downstream sees
    one shape. Two traps in the old one: the covariate lists are unnamed, and the
    censoring models' names live on `cens.time` — as the deparsed left side of the
    formula, marker and all (`WeightIt::.cens(C)`), which `.uncens_name()` strips so the
    time point goes by the indicator's own name. Both that function and this branch can go
    once no such object is left in the wild; a `test-censoring.R` case converts whichever
    shape the fixture is built with into the other and asserts the two agree, so it keeps
    testing both whichever WeightIt is installed.
16. **`bal.plot()`'s longitudinal branch is index-driven now.** It used to replicate every
    per-unit vector `sum(appears.in.time)` times, which assumes each time point
    contributes exactly `n` rows; a censoring time point contributes
    `n_uncensored + n_atrisk` and a restricted treatment time point fewer than `n`. It
    now builds one row-index block per displayed time point and indexes everything by the
    concatenation. `.cens_stack_index()` is that block for a censoring time point, and is
    also what `.stack_cens_X()` stacks with, so there is one definition.
17. **`treat.names` in `bal.plot()` is now one name per time point, not per *displayed*
    time point.** Three lines downstream (`which.time %in% treat.names[appears.in.time]`,
    `treat.names[X[["time"]]]`) index it as if it were full-length, so a variable missing
    from some time point made `which.time` by name select the wrong one.
18. **A time point is named `1. Treatment: A_1` / `2. Censoring: C_1`**, built by
    `.msm_time_labels()` and carried in `print.options$time.labels`. This is the form
    `WeightIt::summary.weightitMSM()` uses, so a model goes by the same name in both
    packages, and it is what the section heading and the collected sample sizes use.
    `bal.plot()` builds the same labels for its facet strips but still matches
    `which.time` against the bare names, so relabeling changes only what is displayed.
    `love.plot()` is left alone: it already faceted by the variable name, and its
    `which.time` matches character input against exactly that facet value.

## Gates

The suite is **2,619 passed / 0 failed / 12 skipped** (2,429 before this work; 190 new
assertions in `test-censoring.R`). The golden set is **137 identical / 0 differing** and
`check_col_spec()` **367 tables / 0 mismatched**, with fifteen new cells: `cens_default`,
`cens_all_stats`, `cens_s_weights`, `cens_cluster`, `cens_imp`, `cens_subclass`,
`cens_subclass_cluster`, `msm_cens_mixed`, `msm_cens_mixed_disp`, `msm_cens_only`,
`msm_cens_cluster`, and the two `obj_weightit_cens*` and two `obj_weightitmsm_cens*` cells
the per-fixture loop generates from the new `weightit_cens` and `weightitmsm_cens`
fixtures. **No existing cell moved**, at any stage of this work.

`R CMD check --as-cran` reports what the branch point reports. Two caveats about this
machine: **`cem` must not be installed** or it hangs the build (see
`_dev/refactor-notes.md` § "`R CMD check`"), and under the agent sandbox the Intel OpenMP
runtime writes `OMP: Warning #179 ... /tmp` to stderr, which `R CMD check` reports as a
`checking Rd files ... WARNING` with that line as its whole content. A check of the
unmodified branch point in the same sandbox produces it identically, so it is the
environment, not the package. Build with vignettes: `--no-build-vignettes` adds two
unrelated `inst/doc` warnings of its own.

## Not done

- **`match.strata`**, rejected: it turns strata into weights before the two samples exist,
  and the same strata say the same thing given to `subclass`.
- **Per-time-point weight sets.** `weightitMSM` returns one weight vector, the product
  across all its models, so balance at an early time point is assessed with weights that
  include later censoring. That is existing cobalt behavior for MSMs and matches
  WeightIt's own documented recipe; it is now said out loud in `?class-bal.tab.msm`.
- **Per-time-point headings.** They stay `Time: <index>`. In a mixed model the index is
  not the study time point, but each child's own sample-size table identifies its kind and
  the summaries' `Times` columns give the indices — not worth churning every existing msm
  golden cell over.
- ~~**Multiply imputed data in *stacked* form** (an `imp` longer than the data) produces
  an empty balance table.~~ Fixed in a follow-up commit, and it was never about censoring:
  `length_imp_process()` replicates with `[`, which drops a matrix's `co.names`, so a
  covariate matrix supplied for one imputation arrived at `.get_C2()` unable to say what
  any of its columns were and came back with none. `.copy_attrs()` puts back what `[`
  dropped, in both `length_imp_process()` and `subset_X()`; only the attributes the result
  does not already have are copied, which is what keeps the row-shaped ones (`dim`,
  `dimnames`, `names`, `row.names`) from being carried over stale. The same commit stopped
  `subset_X()` reading a zero-*column* matrix as absent — `is_null()` is a length test, and
  `.get_C2()` returns an n × 0 matrix by design when nothing is left to assess.

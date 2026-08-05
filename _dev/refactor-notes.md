# The `refactor-bal.tab` branch: what changed, and what to know before changing it again

53 commits on top of `dafd24d`, the tip of `master` when this started. `R/` went
**20,271 → 18,463 lines (−1,808)**, and that is *after* adding 350 lines of new feature and
implementing continuous subclassification, so the refactor itself removed about 2,200.

Line counts elsewhere in this file are measured against `48a3094` — two commits earlier,
the tip of `tests-and-bug-fixes`, where the measuring started — which is 15 lines smaller.
`x2base.R` 4,262 → 2,586 (−39%); `print.bal.tab.R` 1,487 → 1,146 (−23%).

Every commit was verified against three gates. **Read the next section before touching
anything** — the gates are the reason this was safe, and they are cheap to run.

---

## 1. The gates

### `_dev/refactor-golden.R` — the golden-output harness

122 `bal.tab()` cells covering every fixture in `tests/testthat/helper-fixtures.R`
crossed with the option grid, plus nested shapes. Each cell stores the whole object, the
captured stdout of `print()` under **10** argument sets, the conditions those raise
(whitespace-collapsed, so cli's line wrapping is not a diff), and `love.plot()`'s
`$plot$data` for 10 cells. That is ~1,220 captured outputs.

```r
devtools::load_all(); source("_dev/refactor-golden.R")
compare_golden()    # objects via all.equal(), printed output via identical()
check_col_spec()    # 311 tables: columns re-derived from print.options must match
golden_build(filter = "regex")   # regenerate only cells a change intentionally alters
```

A cell whose fixture cannot build is recorded `unavailable = TRUE`, and `compare_golden()`
reports an "availability changed" line. **That field earned its place repeatedly** — it is
how a deleted function and a mislabelled `.msm` flag were both caught, because the symptom
was an error, not a wrong number.

`check_col_spec()` is deliberately *not* tautological: the table builders call
`.bal_tab_col_spec()`, so comparing them to it proves nothing. It instead re-derives each
table's columns from that object's own `print.options` and requires an exact match — the
property `print()`'s column selection and `as.data.frame()` both depend on.

### The test suite

```bash
NOT_CRAN=true COBALT_SLOW_TESTS=true Rscript -e 'devtools::load_all(quiet=TRUE); r <- testthat::test_local(reporter="silent"); d <- as.data.frame(r); cat(sum(d$passed), sum(d$failed), sum(d$skipped), "\n")'
```

Now **2,429 passed / 0 failed / 12 skipped** (was 1,653 at the start of this work). A
changed *skip* count means a fixture broke, not that a refactor worked.

### `R CMD check`

Three invocations wedged before this one was found. What works:

```bash
mkdir -p /tmp/chk && cd /tmp/chk
TMPDIR=/tmp R CMD build /path/to/cobalt
TMPDIR=/tmp NOT_CRAN=true _R_CHECK_FORCE_SUGGESTS_=false _R_CHECK_CRAN_INCOMING_=false \
  R CMD check --as-cran --no-manual cobalt_*.tar.gz
```

Every part matters. `devtools::check()` and a `pueue`-launched check both hang
indefinitely (at "installing the package" and "checking use of S3 registration"
respectively) unless `TMPDIR` is set explicitly — the OMP `Can't set size of /tmp file`
warning is the clue. `_R_CHECK_CRAN_INCOMING_=false` skips a CRAN network call that took
73 s once and hung the next time. Do **not** pass `--no-build-vignettes`: it produces two
spurious warnings, and building them is the only check that they still knit.

**`cem` must not be installed when checking this package.** It Depends on `tcltk`, which
on macOS needs X11, and `bal.tab.cem.match.Rd` makes `prepare_Rd` load it. Sometimes that
only warns — three checks then report `couldn't connect to display …org.xquartz:0` as
2 warnings and a note whose body is otherwise empty. But it can also **hang `R CMD build`
indefinitely at "installing the package"**, which is what two of the earlier wedges
actually were, not the `TMPDIR` problem. The diagnosis is
`lsof -p <pid> | grep -i 'tcltk\|libX11'`: if `tcltk.so`, `libtk8.6`, and `libX11` are
mapped into a build that is producing no output, that is it. `R CMD REMOVE cem` fixes it,
costs nothing — the `cem` example is guarded by `@examplesIf rlang::is_installed("cem") &&
FALSE` and so never runs, and the `cem` test fixture is a stored
`tests/testthat/fixtures/cem_match.rds` — and makes the check come out clean.

**Current result, with `cem` removed: `Status: OK`. 0 errors, 0 warnings, 0 notes.**

---

## 2. What changed

### New seams, in the order a new feature would reach for them

| | |
|---|---|
| `R/bal.tab-columns.R` | `.bal_tab_col_spec()` — **the column layout of a balance table as a value**: one row per column with `name`, `quantity`, `stat`, `sample`, `agg.fun`, `group`. `.p.ops_col_spec()` derives it from a printed object's options; `.keep_bal_cols()` decides which to show; `.stat_cols()`/`.stat_of_col()` serve `love.plot()`. |
| `R/bal.tab-dispatch.R` | `.bal.tab_dispatch()` — all 14 `bal.tab()` methods are now one line. |
| `R/extract.bal.tab.R` | `as.data.frame()` and `format()` for `bal.tab`, plus `.bal.tab_leaves()`, the traversal that yields every innermost balance table with the segmentation identifying it. **Feature 1 will want that traversal.** |
| `x2base.R` | `.x2base_data()`, `.reject_args()`, `.finish_X()` — the head, guards, and tail every adapter shares. |
| `functions_for_processing.R` | `process_stats_and_thresholds()`, `.get_std_defaults()`, `.resolve_s.d.denom()`, `.threshold_label()`/`.threshold_verdicts()`, `.bal.tab_prepare()`, `.bal.tab_summarize()`, `.C_keep()`/`.C_rename()`. |
| `print.bal.tab.R` | `.rewrite_all_none()`, `.display_args()`, `.resolve_p.ops()` (the argument-resolution seam `print()`, `format()`, and `as.data.frame()` share), `.print_observations()`, and seven `.cat_*()` block helpers. |

### The `STATS` registry is now the single extension point

Adding a balance statistic used to mean editing **8 files**, because stat-name literals
lived outside the registry — `c("correlations", "spearman.correlations",
"distance.correlations")` appeared 4×, and `"mean.diffs" %in% X$stats` 5×. Two new fields
killed both families: `needs_s.d.denom` and `agg_fun`. `love.plot()`'s
`startsWith(name, "V.Ratio")` — the last hardcoded prefix — now reads `STATS[[s]]$abs`.

`agg_fun` is stored **by name, not as a function**: `STATS` is built at load time, before
`.geom_mean` exists.

### Twelve bugs fixed

Nine were found while planning; three more surfaced only when merging near-duplicates
revealed they had diverged. All twelve have a regression test and a NEWS entry. The three
worth remembering as a *class*:

- **`print()` on a multi-category `bal.tab` with `quick = FALSE`** errored
  (`'names' attribute [7] must be the same length as the vector [3]`).
- **Narrowing `cluster.fun`/`imp.fun` in `print()`** errored (`undefined columns
  selected`) when `bal.tab()` had computed several.
- **`cluster.summary = TRUE` with subclassified data** errored (`missing value where
  TRUE/FALSE needed`).

The first two are now *unrepresentable*, not merely fixed: they were length mismatches
between a positionally-built display index and the actual columns, and there is now
exactly one spec row per column.

### The four open behaviour questions, resolved

`_dev/stale-code-candidates.md` is gone: it held four decisions rather than refactors, and
all four were made. Recorded here because each changed behaviour.

- **A multi-category ATC is an ATT.** With more than two groups there is no single control
  group for an ATC to name, so `focal` identifies the reference group either way and is
  required either way. The two used to be separate branches in
  `process_focal_and_estimand()` that could disagree; reachable only from `bal.init()`,
  since `bal.tab()` resolves `focal`/`estimand` in each `x2base()` method.
- **`s.d.denom` is honoured per time point.** `base.bal.tab.msm()` used to *overwrite* it
  with `switch(X.class, cont = "all", "pooled")`. That value is still the default — a
  longitudinal treatment targets the ATE and those are the ATE's denominators — but a
  supplied value now survives, and an unusable one is now rejected instead of ignored.
- **`bal.plot(var.name = )` still refuses a split dummy**, and now says which variable to
  supply instead. A factor's dummies exist so `bal.tab()` can summarize the factor one
  level at a time; `bal.plot()` plots the factor.
- **Subclassification with a continuous treatment is implemented** (see below).

### Subclassification with continuous treatments

`base.bal.tab.subclass.cont()` was a stub, and `base.bal.tab.subclass()` already took a
`type`. What was missing was the summary across subclasses, and the reason it was missing
matters: **a binary treatment can express subclassification as weights** — that is exactly
what `strata2weights()` does, and the summary is then an ordinary `balance_table()` on the
reweighted sample. **A continuous treatment cannot.** No set of unit weights makes a
continuous treatment independent of the covariates within a subclass.

So `balance_table_across_subclass()` combines the subclass-specific statistics instead,
weighting each subclass by its share of the subclassified units — the same share
`strata2weights()` gives it, which is worth checking if the aggregation is ever revisited:
for the ATE that function assigns `n_k / n_{k,t}`, so subclass *k*'s total weight is
proportional to `n_k`, **not** `n_k²`. `s.weights` make the share a population share, which
is what keeps the aggregated means equal to the unadjusted ones. Standard deviations are
combined in quadrature, so the summary value is the pooled within-subclass SD.

Two things to know if this is touched again:

- **The layout must be the one a `balance_table()` with something to adjust would produce.**
  My first version passed the default `threshold.samples`, so the table carried an
  `R.Threshold.Un` column that `print()`'s spec does not predict — and since
  `.keep_bal_cols()` returns a logical vector indexed against the table, R *recycled* it and
  silently displayed the wrong columns. `check_col_spec()` is the gate for this; the two new
  `subclass_cont*` golden cells put continuous subclassification under it.
- `balance_table_subclass()`'s continuous SD branch read `subset = treat == in.subclass`
  where it meant `subset = in.subclass`. Nothing had ever run it.

The abandoned `balance_table_across_subclass_cont()` (marked `# !!! NEEDS TO BE UPDATED !!!`)
is deleted. It happened to reach the same aggregation rule, but it was called with
`subclass.obs = out[["Observations"]]` at a point where `Observations` is still NULL, so it
could not have run.

### Feature 2

`as.data.frame.bal.tab()` (tidy: one row per covariate × sample × statistic, with each
level of segmentation as a **column**, never a nested list) and `format.bal.tab()` (the
printed table as formatted strings, so `knitr::kable(format(b))` is publication-ready).
Both accept every `print()` argument and resolve them identically.

The plan had routed this through a new `bal.tab_blocks()` generic. **That turned out to be
unnecessary** — once the column layout and its selection were values, the only remaining
shared piece was three lines of argument resolution. The six `bal.tab_print.*` methods were
left untouched, which also preserved their per-block `cat()` conventions.

---

## 3. Things to remember — most cost a broken build

### R semantics

1. **A named element of `...` matches a formal of the same name even when the caller
   supplied that formal positionally.** A helper formal called `treat` therefore received
   the user's raw `treat` argument. **Dot-prefix any helper formal whose name is also a
   `bal.tab()` argument** — `treat`, `data`, `imp`, `weights`, `stats`, `focal`,
   `estimand`, `subclass`, `cluster`, `s.weights`. Hence `.treat`, `.x`, `.args`, `.imp`,
   `.length`, `.call`, `.msm`, `.env`.
2. **`paste()` treats a zero-length argument as `""`, it does not drop it.**
   `paste(NULL, "Diff", "Un", sep = ".")` is `".Diff.Un"`. Use `.paste_col()`.
3. **`rep.int()` strips names; plain `[` indexing does not.** Passing `seps["int"]` where
   the original passed `rep.int(seps["int"], 1L)` leaked an `"int"` name into every
   interaction term's `component` vector, and thence into the returned object.
4. **`anyNA(NULL)` is `FALSE`, not `TRUE`.** Guard registry lookups with `is_not_null()`.
5. **`[[` on a list with an absent *character* name returns NULL, like `$`** — it errors
   only for out-of-range *numeric* indices. So `$` → `[[` is safe, and `[[` also avoids
   partial matching. (`[.data.frame` *does* preserve non-standard attributes; `[[` on a
   list does not preserve the list's own.)
6. **cli's `{?s}` marks the plural.** `contain{?s}` read "contain" for one and "contains"
   for several — backwards. Use `{?s/}` on a verb. Worth grepping for others.
7. **cli does not re-render markup inside an interpolated value**, and `{.name}` starting
   with a dot is a reserved style. Never paste a source string into `expect_error()`; match
   the *rendered* text (that is what `expect_err()`/`expect_wrn()`/`expect_msg()` in
   `tests/testthat/helpers.R` are for).
8. **A replacement function is called as `f(x) <- v`.** Searching for the literal `f<-`
   finds only its own definition — deleting `treat_names<-` and `treat_vals<-` on that
   basis broke 110 golden cells. An unreferenced-function search must look for both forms,
   **and must exclude comments**: `get_ints_from_co.names()` looked live because its only
   mention was inside a commented-out block.

### Process

9. **When splitting a function, re-thread its arguments.** The `.get_s.d.denom()` rewrite
   dropped `quietly` going into the estimand branch, so a wrapper asking for silence would
   have warned. No golden cell combines an invalid `estimand` with `quietly`; only the new
   unit tests caught it.
10. **Before deleting a local, check for reads that precede any assignment.** Removing a
    guard also removed `focal <- ...get("focal")`, and eight methods read that local
    afterwards.
11. **Do not delimit a function by "the next definition matching a pattern".** A script
    using "up to the next `x2base.*`" swallowed two helper functions that are not methods —
    deleting one outright and giving the other its neighbour's `.msm = TRUE`. Brace-match,
    and audit every generated call site against the previous revision.
12. **Tests before a rewrite, and they must pin the artifact, not the surface.** This
    worked three times. For `get_covs_from_formula()` the artifact is the `co.names`
    *decomposition* — a rewrite could get every column name right and every label wrong.

### Architecture invariants that will bite

13. **`X`'s slot names are the contract.** `.finish_X()` fills each slot from the
    like-named local in the calling `x2base()` method (`X_SLOTS` / `X_MSM_SLOTS`). A slot
    never assigned stays NULL, which is what every consumer expects of an absent one —
    which is why the reflection loop was *kept* rather than replaced with explicit
    arguments. `names(X)` must stay the full slot set including NULLs: that is the only
    thing stopping a user's `weights=`/`distance=`/`addl=` from reaching
    `balance_table()` as a duplicate formal.
14. **`.assign_X_class()` ranks `cluster` above `msm`**, so **cluster is the only shape
    that can receive a longitudinal `X`**. The `if (is_null(X[["covs.list"]]))` guard in
    the cluster wrapper is deliberate; do not generalise it to imp/multi. Relatedly,
    `.cluster_check()` needs `X[["treat"]] %or% X[["treat.list"]]` — it used to reach the
    list only by `$` partial matching.
15. **`s.d.denom` nullification differs three ways and must stay visible.** A leaf clears
    it when nothing is standardized; cluster and imp keep the user's value
    (`%or% X[["s.d.denom"]]`); the multi-category wrapper precomputes one denominator per
    weight set into `s.d.denom.list` for its children to share. Hiding this behind a shared
    default makes each per-stratum child re-derive its own — and a covariate constant within
    one cluster can pick a *different* denominator, plus emit a spurious note.
16. **`print_process()` returns a hand-written subset of `print.options`, not all of it.**
    `weight.names`, `quick`, and `stats` were computed and discarded; anything reading
    `p.ops` at print time silently got NULL. Adding a field to `print.options` is not
    enough — pass it through too.
17. **Do not attach column metadata to `$Balance`.** `[` on a data frame preserves
    non-standard attributes, and `test-bal.tab.data.frame.R` compares a `quick = FALSE`
    object's columns against a `quick = TRUE` one. The column spec must stay a pure
    function of `p.ops`.
18. **`balance_table_subclass()`'s `compute` is deliberately not intersected with `stats`**,
    unlike `balance_table()`'s. `test-print.bal.tab.R` depends on it. `compute` is also the
    one attribute `balance_table()` returns *untouched* — `disp` and `thresholds` both lose
    entries whose columns came out non-finite — which is why
    `balance_table_across_subclass()` can rebuild the layout from it.
19. **A balance table with a column its `print.options` do not predict displays the wrong
    columns, silently.** `.keep_bal_cols()` returns a logical vector that is indexed against
    the table, so an extra column makes R recycle it rather than error. `check_col_spec()`
    exists for exactly this; add a golden cell for any new table shape.
20. **The subclass `print_process()` must *remove* `nweights`, not set it to NULL.** For a
    nested object `p.ops` is `c(parent, child)` and `$` takes the first match, so a
    present-but-NULL entry masks the parent's count.
21. **Separators are substituted at the very end of `.get_C2()`.** Everything above
    compares names built with the separators the pieces arrived with; hoisting the
    substitution changes which columns count as duplicates.
22. **`.get_C2()`'s pieces carry their own `co.names`.** `.C_keep()` is the only thing that
    removes a column, so the matrix and its parsed names cannot desynchronise. Keep it that
    way.
23. **`x2base.ps` probes a positional `stop.method` from raw `...`** via
    `...names()`/`...elt()`, reachable only from `bal.plot()`. And **`bal.plot()` passes its
    own lazy `...` through**, so converting `...get("foo")` to `A[["foo"]]` forces those
    promises and changes when an erroring unused argument fails. `.finish_X()` forwards
    `...` unforced for this reason.

### Fixture hazards

24. `lalonde` is sorted treated-first (185 of 614), so `rep(..., length.out = n)` indices
    are collinear with treatment at some periods. The `mids` fixture needs `seed = 5678L`
    or `mice()` is nondeterministic across sessions. `covr` needs `NOT_CRAN` in the
    **shell** environment, not `Sys.setenv()` in-R.

---

## 4. What could be done next

**Ranked by value, with the reason each was not done.**

1. **Feature 1: censoring / target balance.** Designed in the plan file; the next piece of
   work. No new statistics are needed — the three `STATS[["*.target"]]` entries already do
   the doubling internally. It needs a new `X.class = "target"` as a *leaf* peer to
   `binary`/`cont` at the bottom of `.assign_X_class()`, so `cluster` and `imp` compose for
   free. `.bal.tab_leaves()` and `.bal_tab_col_spec()` are the seams it should build on.

   **Correction to the plan file, verified 2026-08-05.** The plan says to discard
   `_dev/Under_construction/target.bal.tab.R` because it "uses an n0+n stacking that
   disagrees for `ovl.*`/`energy.dist`". **That is wrong on both counts.** Read the two
   side by side:

   - `STATS[["mean.diffs.target"]]$fun` does `C <- rbind(C, C)`,
     `treat <- rep(c(0, 1), each = n)`, `weights <- c(weights, rep.int(1, n))`.
   - The prototype (`target.bal.tab.R:114-121`) does `covs <- rbind(covs, covs)`,
     `weights <- rbind(weights, <all 1s>)`, and stacks `s.weights`, `distance`, `addl`,
     `discarded` to match.

   Both are the **same 2n doubling**: the weighted sample against the same units
   unweighted. There is no n0+n stacking anywhere and no demonstrated disagreement.

   What they genuinely differ on is the *grouping*, and that is the real design question:

   - The registry compares the **whole weighted sample** against the whole unweighted
     sample — two groups.
   - The prototype builds `k + 1` levels (each treatment group plus a `target` level) and
     compares **each treatment group** against the unweighted sample pairwise, with
     `s.d.denom` hardcoded to `"treated"`.

   The prototype's framing is a superset and answers the more useful question for a
   censoring model. So it is worth *reading* rather than discarding — the parts to leave
   behind are its hardcoded `s.d.denom`, and its `stop()`s for `cluster` and `subclass`,
   which the new `X.class = "target"` leaf makes unnecessary.
2. **`get_covs_from_formula()` is still 345 lines** and has no duplication left that I could
   find; the length is genuine surface area (`.` expansion, nested data frames and matrices,
   backtick-quoting, missingness indicators, single-level factors, the `co.names`
   construction). It now has `test-get-covs-from-formula.R` as a gate, so a future attempt
   is safe. **`get_covs_and_treat_from_formula2()` in WeightIt cannot replace it** — checked
   against WeightIt 2.0.0: circular dependency (WeightIt Imports cobalt), a different return
   contract, and decisively no `co.names`. Do not revisit without new information.
3. **`.use_tc_fd()` (93 lines)** resolves the treat/covs versus formula/data conventions. A
   decision table over genuinely different inputs; long but not repetitive.
4. **The four `bal.tab()` wrapper *heads* are still four functions.** Deliberate: they
   differ on eight axes, and the msm one does not subset at all (it reshapes `X`). Only the
   shared tail was extracted. A 9-field spec with two callback slots would read worse than
   what it replaced.
5. **`base.bal.tab()`'s hand-rolled `switch()` on `attr(X, "X.class")` stays.** 12 clear
   lines, and a string key is a better dispatch key than an S3 class here.
6. **`_dev/Under_construction/_as.data.frame.bal.tab.R` is superseded** by
   `R/extract.bal.tab.R`. It was a 21-line `reshape(direction = "long")` sketch; the
   shipped version is long-by-default with segmentation as columns, and `format()` covers
   the printed-table case. `bal.tab2tableone.R` was assessed and rejected in the plan: the
   prototype *fabricates* `median`/`p25`/`p75`/`skew`/`kurt` that `bal.tab` never computed.
   Neither is referenced by anything.

7. **Tests still to write:** `bal.plot()` has the thinnest coverage of the exported
   surface. And `love.plot()`'s aggregation path is pinned by exactly one test — the golden
   set captures `love.plot()` data for 10 cells, **none of which aggregate**.

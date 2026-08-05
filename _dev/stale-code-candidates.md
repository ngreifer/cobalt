# Stale / unreachable code candidates

Compiled while building out the test suite. Each entry is code that could not be
reached from any public entry point, or that is provably dead. Nothing here has
been deleted — this is input to the planned refactor, not a change log.

Line numbers are as of the commit that added this file; re-grep before acting.

## Status

The **Provably dead** entries below were removed on the `refactor-bal.tab` branch
(commit "Remove provably dead code"), along with the never-reassigned
`bad.s.d.denom`/`bad.estimand` flags and their three warning branches, the thirteen
`set_class()` calls on the `X` list that nothing read, and the `target.summary`
option that nothing read. The `treatATT` alias was repaired rather than removed.

The **Unreachable defensive guards** are deliberately kept: they cost nothing and
document the invariant they check. The **Behaviour worth a decision** items are
still open.

## Provably dead

- **`R/functions_for_processing.R:1351–1353`** — a repeat of the check at `:1342`.
  Both test `any(abs(subset) > n)`; in between, `subset` is only *filtered*, never
  grown, so the second test cannot be true if the first was false. Verified:
  `bal.tab(covs, treat = t, subset = c(1L, 10000L))` errors from `:1343`.

- **`R/balance-summary.R:422`** — `all_pos_w <- all(weights >= 0)` is always `TRUE`,
  because `arg::arg_gte(weights, 0)` at `:410` has already rejected negatives. Its
  only use is `if (all(all_pos_w) && ...)` at `:464`, so the false branch is dead.

- **`R/x2base.R:3584–3587`** — guards `P[["distance.list"]]`, but `P` is built with
  `make_list(names(Q))` and `Q` has no `distance.list` entry, so that slot is
  always `NULL`. (A component *named* `distance.list` in the input is picked up
  under `P[["distance"]]`, since `distance.list` is one of the accepted aliases for
  it. `.x2base_default_msm()` was reading the always-`NULL` slot; that has been
  fixed to read `obj[["distance"]]`.)

- **Commented-out `arg::err()` calls**: `R/x2base.R:916`, `R/utils.R:971`.

- **`R/bal.plot.R:962`** — the facet-ordering fallback does
  `facets[order(facet$length, facet$facet, ...)]`, but `facet` is a character
  vector, so `facet$length` would error. The branch appears unreachable because
  `facet` always contains `"which"` or is exactly `"subclass"`.

- **`R/class-bal.tab.subclass.R:242`** — `base.bal.tab.subclass.cont()` is a stub
  whose entire body is `arg::err("subclasses are not yet compatible with
  continuous treatments")`; the real implementation is commented out below it. The
  error is reachable (`bal.tab(covs, treat = <continuous>, subclass = s)`), but the
  method is a placeholder rather than an implementation.

## Unreachable defensive guards

These cannot fire even when the enclosing internal function is called directly
with adversarial input, because an upstream `match_arg()`/type check has already
normalized the value. They are left in place as guards and are simply never
covered.

- `R/functions_for_processing.R:936` — `s.d.denom` is `match_arg`'d in
  `.get_s.d.denom()` and in `col_w_smd()` before reaching `.compute_s.d.denom()`.
- `R/functions_for_processing.R:1008` — `.get_length_X()`'s final `else`; by the
  time `.assign_X_class()` has run, `X` always has one of the four components.
- `R/functions_for_processing.R:3031`, `:3125` — `samplesize()` requires
  `method == "subclassification"`, which is *derived from* a non-`NULL` `subclass`,
  so `is_null(subclass)` is impossible.
- `R/utils.R:1005` — `probably.a.bug()`, and its sole call site at
  `R/functions_for_processing.R:997`: `get.treat.type()` returns one of three
  values, all handled above.
- `R/base.bal.tab.R:41` and `R/class-bal.tab.subclass.R:57` — both check
  `type == "bin" && get.treat.type(X$treat) != "binary"`, but these functions are
  only dispatched when `.assign_X_class()` set the class *because* the treatment is
  binary.
- `R/balance-summary.R:419` — a non-negative, not-all-zero weight vector cannot sum
  to zero, and both conditions are guaranteed by `:410` and `:412–415`.

## Dead alias

- **`treatATT` as an alias for `focal` in `x2base.default()`** (`R/x2base.R:3505`).
  `x2base.default()` does `names(x) <- tolower(names(x))` at `:3511`, then looks up
  `x[["treatATT"]]`, which can never match. Verified: neither `treatATT` nor
  `treatatt` sets `focal` — both fall through to the no-focal path. Either lowercase
  the alias in `Q` or drop it. (`x2base.mnps()` handles twang's `treatATT`
  separately, so real `mnps` objects are unaffected.)

## Behaviour worth a decision

- **`estimand = "ATC"` with a multi-category treatment and no `focal`** does not
  error from the `data.frame` or `formula` entry points, though
  `R/functions_for_processing.R:1666` intends it to. Reachable only via
  fitted-object paths, if at all.

- **On the cluster + longitudinal path**, `s.d.denom` comes from
  `base.bal.tab.msm()`'s hardcoded `"pooled"`/`"all"` per time point rather than a
  shared full-sample denominator, so a user-supplied `s.d.denom` is ignored. This
  predates the fix that made cluster + longitudinal work at all.

- **`var.name` does not accept split-dummy names in `bal.plot()`.**
  `bal.plot(covs, treat = t, var.name = "race_black")` errors with "is not the name
  of an available covariate", and so does the same call on a `bal.tab` object, even
  though `race_black` is a row name in the balance table. Only the parent factor
  (`"race"`) works. Worth confirming whether the documented `unsplitfactor()`
  round-trip is meant to make the dummy names resolvable.

## Assessed for replacement, deliberately left (2026-08-05)

A pass over every function in `R/functions_for_processing.R`, asking what it does, who
calls it, and whether something simpler would do. Six were replaced and one collapsed
(see `git log`). These were examined and kept, with reasons, so the next pass does not
have to re-derive them.

**Two of the three big ones are done** (`.get_s.d.denom()` and `.get_C2()`; see
`git log`). One remains:

- **`get_covs_from_formula()` (385 lines)** — the formula parser. Handles `.`, `poly()`,
  transformations, factor splitting, `offset()`, and per-time-point naming. It is long
  because the surface is large, not because it repeats itself; there is no obvious
  duplication to fold. Any rewrite should start by pinning its behaviour with a
  dedicated test file, since it currently has no direct tests — only what reaches it
  through `bal.tab()`.

**Kept as they are, for reasons that will not change:**

- `model.frame2()` — is a `stats::model.frame()` wrapper whose whole point is turning
  "object 'x' not found" into a message naming the variable and suggesting `data`.
- `find_perfect_col()`, `co.cbind()`, `.mids_complete()`, `strata2weights()`,
  `intapprox()` — each does one thing at about the length that thing takes.
- `.use_tc_fd()` (93 lines) — resolves the treat/covs versus formula/data calling
  conventions. Long, but it is a decision table over genuinely different inputs.
- `.baltal()` recovers the numeric threshold by regex from the label
  `"Balanced, <0.1"` that `.threshold_label()` generated. Parsing a string the package
  itself just produced is backwards, and every caller has the threshold to hand. Worth
  changing; not done here because it means touching `.baltal()`'s signature and its
  four call sites in the same commit as nothing else.

## Still open

All three functions flagged for rewrite are done: `.get_s.d.denom()`, `.get_C2()`, and
`get_covs_from_formula()`. Each has direct unit tests written *before* the rewrite; that
sequence is what made all three safe and is worth repeating.

**`get_covs_and_treat_from_formula2()` in WeightIt cannot replace
`get_covs_from_formula()`**, checked 2026-08-05 against WeightIt 2.0.0. Three
independent reasons, so do not revisit this without new information:

1. WeightIt Imports cobalt, so cobalt depending on WeightIt would be circular.
2. Different contract. WeightIt's returns a four-element list -- `reported.covs`,
   `model.covs`, `simple.covs`, `treat` -- for a model-fitting pipeline; cobalt's returns
   one matrix and takes `factor_sep`/`int_sep` separately rather than a single `sep`.
3. The decisive one: WeightIt's attaches **no `co.names`**. Its `model.covs` produces the
   same *column names* as cobalt's, but not the decomposition of each name into
   components and types (`base`/`fsep`/`level`/`isep`/`na`) that `.get_C2()`,
   `love.plot()`, `var.names()`, and `.baltal()` all read. Building that structure is
   most of what makes cobalt's version longer, so the extra length is earning its keep.

**What is left in `get_covs_from_formula()` (345 lines).** No duplication remains that I
can see. The length is the surface: `.` expansion, data frames and matrices nested in
`data`, backtick-quoting of names that need it, missingness indicators inserted after
the variable they belong to, single-level factors, and the `co.names` construction
itself. It is now covered by `test-get-covs-from-formula.R`, so a future attempt has a
gate.

**Three lessons from the rewrites**, none of which the golden set would have caught alone:

- Splitting a function into helpers means re-threading its arguments. The
  `.get_s.d.denom()` rewrite dropped `quietly` going into the estimand branch, so a
  wrapper asking for silence would have warned. Only the new unit tests caught it.
- `rep.int()` strips names; plain `[` indexing of a named vector does not. Passing
  `seps["int"]` where the original passed `rep.int(seps["int"], 1L)` leaked an `"int"`
  name into every interaction term's `component` vector. Two golden cells caught it.
- cli's `{?s}` appends "s" in the *plural*. `contain{?s}` therefore read "contain" for
  one variable and "contains" for several -- backwards. Worth grepping for other
  `{?s}` markers attached to a verb rather than a noun.

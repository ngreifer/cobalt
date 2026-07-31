# Stale / unreachable code candidates

Compiled while building out the test suite. Each entry is code that could not be
reached from any public entry point, or that is provably dead. Nothing here has
been deleted — this is input to the planned refactor, not a change log.

Line numbers are as of the commit that added this file; re-grep before acting.

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

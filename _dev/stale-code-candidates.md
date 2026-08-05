# Open questions about cobalt's behaviour

What this file used to be: a list of unreachable and provably dead code found while
building out the test suite. All of that has been dealt with — see
`_dev/refactor-notes.md` for what happened to each category, including the guards that
were examined and deliberately kept, so they do not get re-flagged.

What is left here is the residue: four things where the code does something defensible
but not obviously intended, and where the right answer is a decision rather than a
refactor. Each was re-verified on 2026-08-05 against the current branch.

## 1. `estimand = "ATC"` with a multi-category treatment and no `focal`

`bal.tab(covs, treat = <3 levels>, estimand = "ATC", weights = w)` does not error. It
notes `s.d.denom` not specified; assuming "pooled"` and proceeds.

`.s.d.denom_from_estimand()` returns NULL for an ATT or ATC on a non-binary treatment,
deliberately, so that the focal group can decide instead — and with no focal group it
falls through to inferring from the weights. That is reasonable. But
`process_focal_and_estimand()` contains a branch that intends to *error* on this
combination, which never fires from the `data.frame` or `formula` entry points.

**Decision needed:** either drop the unreachable branch, or make this an error. Silently
computing a pooled denominator for a requested ATC is defensible only because the note
says what happened.

## 2. `s.d.denom` is ignored on the cluster + longitudinal path

Verified: `bal.tab(list(f1, f2), data = d, weights = "w", cluster = cl,
s.d.denom = "treated")` gives numerically identical results to the same call with
`s.d.denom = "pooled"`.

`base.bal.tab.msm()` overwrites `X[["s.d.denom"]]` per time point with a hardcoded
`switch(X.class, cont = "all", "pooled")`, so whatever the user asked for is discarded.
This predates the fix that made clustered longitudinal data work at all.

**Decision needed:** honour the user's value when supplied and fall back to the hardcoded
default otherwise. The hardcoded value exists because each time point's `X` is rebuilt
from `covs.list`/`treat.list` and never went through `.get_s.d.denom()`; the fix is to
let it, which is small but changes numbers for anyone who supplied `s.d.denom` with a
longitudinal treatment.

## 3. `bal.plot(var.name = )` does not accept split-dummy names

`bal.plot(lalonde[c("age", "race")], treat = t, var.name = "race_black")` errors with
`"race_black" is not the name of an available covariate`, though `race_black` is a row
name in the corresponding balance table. Only the parent factor (`"race"`) works, and it
plots the whole factor.

The previous version of this note also claimed the same call fails on a `bal.tab` object.
That is wrong, and the reason is worth recording: `bal.plot()` does not accept `bal.tab`
objects *at all* — its `x` is "any object for which there is support in `bal.tab()`", and
a `bal.tab` is not one of those. Such a call reaches `x2base.default()` and fails with an
unrelated `treat.list` message. Confirmed identical at `48a3094`, so this is long-standing
and not something the refactor changed.

**Decision needed:** whether `var.name` should resolve a dummy's name by splitting the
parent factor, which is what a user reading the balance table would expect.

## 4. `base.bal.tab.subclass.cont()` is a stub

Subclassification with a continuous treatment errors. The machinery half-exists:
`base.bal.tab.subclass()` takes `type` and would accept `"cont"`, `samplesize()` has a
continuous subclassification branch, and `balance_table_across_subclass_cont()` exists —
but none of the three has ever been exercised, and the last is marked
`# !!! NEEDS TO BE UPDATED !!!`.

**Decision needed:** implement it (which means checking those three paths produce correct
numbers) or delete the unreachable scaffolding and keep the error. Leaving it half-built
is the worst of the three.

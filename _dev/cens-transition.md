# The `.cens()` transition: where it stands

One function and one `[` method, shared between cobalt and WeightIt, instead of a copy in
each. Two of the three steps have shipped. This records what is done, the one correction
to the original plan, and the ordering of what is left.

## Done

**WeightIt 2.1.0 (on CRAN).** `.cens()` is cobalt's, re-exported (`importFrom(cobalt,.cens)`
plus `export(.cens)`), so `WeightIt::.cens` and `cobalt::.cens` are the same object and
`base::library()`'s conflict check — which filters identical objects — no longer reports a
masked `.cens` when both packages are attached. WeightIt registers no `[` method; it tags
its processed treatments `c("cobalt.treat", "treat")`, exactly as `process_treat()` does
here, so cobalt's method is what subsets them.

**cobalt 5.0.0.9000 (this version).** The method is registered on both `treat` and
`cobalt.treat`, the function renamed to `[.treat`:

```r
#' @exportS3Method `[` treat
#' @rawNamespace S3method("[",cobalt.treat,"[.treat")
```

roxygen will not take two `@exportS3Method` tags in one block, and the second registration
needs the three-argument form anyway because the function's name no longer matches the
class, so it goes through `@rawNamespace`. The documentation that described cobalt as
having defined its own `.cens()` was updated to say WeightIt re-exports this one.

## The correction

The original plan here was to **rename** `[.cobalt.treat` to `[.treat` and drop the alias
and the class in one step, on the reasoning that WeightIt 2.0.0 leaving CRAN removed the
collision. That is wrong, and the way it is wrong is worth keeping.

The collision does not disappear when 2.0.0 leaves CRAN; it disappears when 2.0.0 leaves
people's libraries, which is a different and much later date. For a user who has cobalt-new
and WeightIt 2.0.0 installed, R reports the overwrite and WeightIt's method wins — and that
method copies `treat.type`, `treat.name`, `treated`, and `control`, but not `treat_names`,
`treat_vals`, or `group_labels`.

The console message is unavoidable the moment cobalt claims `treat` at all. What the alias
buys is everything else: cobalt's own objects carry `cobalt.treat` first, find the alias
first, and never reach the overwritten method. Measured with stub packages standing in for
each shape:

| cobalt shape | overwrite message | cobalt's own objects keep their attributes |
|---|---|---|
| `[.treat` only, class `treat` | yes | **no** |
| both registrations, class `cobalt.treat` first | yes | yes |

So registering both dominates renaming: same message, no breakage. The rename is the step
that looks like cleanup and is actually the one that breaks silently.

## What is left, in order

1. **WeightIt drops the alias from its own objects.** `.treat_classes <- "treat"` in
   `R/treat.R`, and `Imports: cobalt (>= <the version carrying `[.treat`>)`. Blocked until
   that cobalt is on CRAN: against cobalt 5.0.0 a treatment classed only `treat` finds no
   `[` method at all and loses every attribute on the first subset.
2. **cobalt drops the alias.** Delete the `@rawNamespace` line and `"cobalt.treat"` from the
   two `set_class()` calls (`process_treat()` in `R/treat.R`, `.cens()` in `R/cens.R`).
   Safe once the WeightIt from step 1 is the floor, because from then on nothing sets the
   alias and nothing registers a competing `[.treat`.

Each step is safe against the released version of the other package, which is the property
to preserve: at no point does either package require an unreleased version of the other.

`test-censoring.R` pins the part that would otherwise break quietly — it subsets a processed
treatment under all three class shapes (`c("cobalt.treat", "treat")`, `treat` alone,
`cobalt.treat` alone) and asserts every attribute present survives. The `treat`-alone case
is step 1's shape, so it is already tested here before WeightIt moves to it.

## Two things that stay duplicated, on purpose

`.make_cens_treat()` and the `.cens()`-formula helpers exist in both packages. cobalt cannot
depend on WeightIt, WeightIt should not reach into cobalt with `:::`, and each validates its
own inputs at its own call sites. Leave them.

The two `.make_cens_treat()` implementations differ only in the wording of the `i` bullet in
the error and in where the all-`NA` check sits (`allNA(C)` here, `nunique(C) == 0L` inside
WeightIt's `assign_treat_type()`). Same error, same text, different file. Change one and
change the other — an all-`NA` indicator should be rejected in the same words whichever
package the user reached it through. Both packages' tests match on `"only the values"`, so
neither is pinned to the other's exact phrasing.

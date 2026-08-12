# Mark a Censoring Indicator

`.cens()` marks a variable as a censoring indicator rather than a
treatment, so that
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
assesses the balance of the units still under observation against the
full at-risk sample. It is most often used on the left side of a
formula, as in `.cens(C) ~ x1 + x2`, but it can also be called directly
and supplied to
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)'s
`treat` argument.

## Usage

``` r
.cens(x)
```

## Arguments

- x:

  a censoring indicator: a numeric variable taking only the values 0
  (still under observation) and 1 (censored), a logical variable, or a
  factor with levels `0`/`1` or `FALSE`/`TRUE`. Missing values are
  allowed and preserved.

## Value

`x` coerced to a 0/1 numeric vector of class `treat` (see
[`treat`](https://ngreifer.github.io/cobalt/reference/treat-class.md))
with a `"treat.type"` attribute of `"censoring"`. Any value other than
0, 1, or `NA` throws an error.

Inside a formula the marker is stripped before the formula is processed,
so `.cens()` is not actually evaluated there and the treatment name
remains that of the indicator itself (e.g., `C` rather than `.cens(C)`).

## Details

Censoring is considered its own treatment type in cobalt, distinct from
binary, multi-category, and continuous treatments. What matters is
whether the units still under observation resemble the full sample once
weighted. See
[`class-bal.tab.cens`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.cens.md)
for the output this produces and the arguments that control it.

This function is deliberately identical in name and contract to
[`WeightIt::.cens()`](https://ngreifer.github.io/WeightIt/reference/dot-cens.html),
so that the same code works whichever package is attached; cobalt
defines its own only to avoid depending on WeightIt.
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
recognizes an indicator tagged by either.

## Note on the coding convention

The survival convention is used: `1` means the unit is *censored* (drops
out of observation) and `0` means it remains under observation. This is
the opposite of an "observed" or "event" indicator.

Missing values are permitted because a unit censored at an earlier time
point has no later indicator. Units with a missing indicator are not at
risk, so they are in neither sample.

## See also

- [`class-bal.tab.cens`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.cens.md)
  for the output of
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  with a censoring indicator

- [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)

## Examples

``` r
data("lalonde", package = "cobalt")

# A censoring indicator: 1 = lost to follow-up
set.seed(1234)
lalonde$C <- rbinom(nrow(lalonde), 1,
                    prob = plogis(-1.5 + 0.05 * lalonde$age))

# Inverse probability of censoring weights
W <- WeightIt::weightit(.cens(C) ~ age + educ + race,
                        data = lalonde, method = "glm")

# Balance of the weighted uncensored units against the
# full at-risk sample
bal.tab(W, un = TRUE)
#> Balance Measures
#>                 Type Diff.Un Diff.Adj
#> prop.score  Distance  0.1816  -0.0024
#> age          Contin.  0.1713  -0.0028
#> educ         Contin. -0.0277  -0.0130
#> race_black    Binary -0.0266   0.0041
#> race_hispan   Binary -0.0132  -0.0012
#> race_white    Binary  0.0398  -0.0029
#> 
#> Effective sample sizes
#>             Total
#> Full       614.  
#> Uncensored 322.  
#> Adjusted   308.33
#> Censored   292.  

# The same thing without WeightIt, using the weights directly
bal.tab(.cens(C) ~ age + educ + race, data = lalonde,
        weights = W$weights, un = TRUE)
#> Balance Measures
#>                Type Diff.Un Diff.Adj
#> age         Contin.  0.1713  -0.0028
#> educ        Contin. -0.0277  -0.0130
#> race_black   Binary -0.0266   0.0041
#> race_hispan  Binary -0.0132  -0.0012
#> race_white   Binary  0.0398  -0.0029
#> 
#> Effective sample sizes
#>             Total
#> Full       614.  
#> Uncensored 322.  
#> Adjusted   308.33
#> Censored   292.  
```

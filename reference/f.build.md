# Convenient Formula Generation

`f.build()` returns a [`formula`](https://rdrr.io/r/stats/formula.html)
of the form `y ~ x1 + x2 + ...` from a data frame input. It can be much
quicker to use `f.build()` than to hand-write the precise formula, which
may contain errors. It can be used in place of a formula in, for
example, [`glm()`](https://rdrr.io/r/stats/glm.html),
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.html),
or
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md).
It provides similar functionality to
[`reformulate()`](https://rdrr.io/r/stats/delete.response.html).

## Usage

``` r
f.build(y = NULL, rhs = NULL)
```

## Arguments

- y:

  the quoted name of the response (left hand side) variable in the
  formula. Only one variable is supported. If missing, `NULL`, or the
  empty string (`""`), the formula will have no response variable. If
  `rhs` is not supplied, `y` will replace `rhs` and `y` will be set to
  `""`.

- rhs:

  a data frame whose variable names will be the terms on the right hand
  side of the formula, or a character vector whose values will be the
  terms on the right hand side of the formula. If missing, the argument
  to `y` will replace `rhs` and `y` will be set to `""`; in essence,
  `f.build("x")` is the same as `f.build("", "x")`, both producing
  `~ x`.

## Value

a `formula` object.

## See also

[`reformulate()`](https://rdrr.io/r/stats/delete.response.html)

## Examples

``` r
data(lalonde)
covs <- subset(lalonde, select = -c(treat, re78))
#> Error in eval(substitute(select), nl, parent.frame()): object 'treat' not found
lm(f.build("treat", covs), data = lalonde)
#> Error in eval(mf, parent.frame()): object 'covs' not found
```

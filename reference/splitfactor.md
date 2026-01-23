# Split and Unsplit Factors into Dummy Variables

`splitfactor()` splits factor variables into dummy (0/1) variables. This
can be useful when functions do not process factor variables well or
require numeric matrices to operate. `unsplitfactor()` combines dummy
variables into factor variables, undoing the operation of
`splitfactor()`.

## Usage

``` r
splitfactor(
  data,
  var.name,
  drop.level = NULL,
  drop.first = TRUE,
  drop.singleton = FALSE,
  drop.na = TRUE,
  sep = "_",
  replace = TRUE,
  split.with = NULL,
  check = TRUE
)

unsplitfactor(
  data,
  var.name,
  dropped.level = NULL,
  dropped.na = TRUE,
  sep = "_",
  replace = TRUE
)
```

## Arguments

- data:

  A `data.frame` containing the variables to be split or unsplit. In
  `splitfactor()`, can be a factor variable to be split.

- var.name:

  For `splitfactor()`, the names of the factor variables to split. If
  not specified, will split all factor variables in `data`. If `data` is
  a factor, the stem for each of the new variables to be created. For
  `unsplitfactor()`, the name of the previously split factor. If not
  specified and `data` is the output of a call to `splitfactor()`, all
  previously split variables will be unsplit.

- drop.level:

  The name of a level of `var.name` for which to drop the dummy
  variable. Only works if there is only one variable to be split.

- drop.first:

  Whether to drop the first dummy created for each factor. If `"if2"`,
  will only drop the first category if the factor has exactly two
  levels. The default is to always drop the first dummy (`TRUE`).

- drop.singleton:

  Whether to drop a factor variable if it only has one level.

- drop.na:

  If `NA`s are present in the variable, how to handle them. If `TRUE`,
  no new dummy will be created for `NA` values, but all created dummies
  will have `NA` where the original variable was `NA`. If `FALSE`, `NA`
  will be treated like any other factor level, given its own column, and
  the other dummies will have a value of 0 where the original variable
  is `NA`.

- sep:

  A character separating the the stem from the value of the variable for
  each dummy. For example, for `"race_black"`, `sep = "_"`.

- replace:

  Whether to replace the original variable(s) with the new variable(s)
  (`TRUE`) or the append the newly created variable(s) to the end of the
  data set (`FALSE`).

- split.with:

  A list of vectors or factors with lengths equal to the number of
  columns of `data` that are to be split in the same way `data` is. See
  Details.

- check:

  Whether to make sure the variables specified in `var.name` are
  actually factor (or character) variables. If splitting non-factor (or
  non-character) variables into dummies, set `check = FALSE`. If
  `check = FALSE` and `data` is a `data.frame`, an argument to
  `var.name` must be specified.

- dropped.level:

  The value of each original factor variable whose dummy was dropped
  when the variable was split. If left empty and a dummy was dropped,
  the resulting factor will have the value `NA` instead of the dropped
  value. There should be one entry per variable to unsplit. If no dummy
  was dropped for a variable, an entry is still required, but it will be
  ignored.

- dropped.na:

  If `TRUE`, will assume that `NA`s in the variables to be unsplit
  correspond to `NA` in the unsplit factor (i.e., that `drop.na = TRUE`
  was specified in `split.factor()`). If `FALSE`, will assume there is a
  dummy called "var.name_stem_NA" (e.g., "x_NA") that contains 1s where
  the unsplit factor should be `NA` (i.e., that `drop.na = FALSE` was
  specified in `split.factor()`. If `NA`s are stored in a different
  column with the same stem, e.g., "x_miss", that name (e.g., "miss")
  can be entered instead.

## Value

For `splitfactor()`, a `data.frame` containing the original data set
with the newly created dummies. For `unsplitfactor()`. a `data.frame`
containing the original data set with the newly created factor
variables.

## Details

If there are `NA`s in the variable to be split, the new variables
created by `splitfactor()` will have `NA` where the original variable is
`NA`.

When using `unsplitfactor()` on a `data.frame` that was generated with
`splitfactor()`, the arguments `dropped.na`, and `sep` are unnecessary.

If `split.with` is supplied, the elements will be split in the same way
`data` is. For example, if `data` contained a 4-level factor that was to
be split, the entries of `split.with` at the same index as the factor
and would be duplicated so that resulting entries will have the same
length as the number of columns of `data` after being split. The
resulting values are stored in the `"split.with"` attribute of the
output object. See Examples.

## See also

[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)

## Examples

``` r
data("lalonde", package = "cobalt")

lalonde.split <- splitfactor(lalonde, "race", 
                             replace = TRUE, 
                             drop.first = TRUE)
# A data set with "race_hispan" and "race_white" instead
# of "race".

lalonde.unsplit <- unsplitfactor(lalonde.split, "race", 
                                 replace = TRUE,
                                 dropped.level = "black")

all.equal(lalonde, lalonde.unsplit) #TRUE
#> [1] TRUE

# Demonstrating the use of split.with:
to.split <- list(letters[1:ncol(lalonde)],
                 1:ncol(lalonde))

lalonde.split <- splitfactor(lalonde, split.with = to.split,
                             drop.first = FALSE)
attr(lalonde.split, "split.with")
#> [[1]]
#>       treat         age        educ  race_black race_hispan  race_white 
#>         "a"         "b"         "c"         "d"         "d"         "d" 
#>     married    nodegree        re74        re75        re78 
#>         "e"         "f"         "g"         "h"         "i" 
#> 
#> [[2]]
#>       treat         age        educ  race_black race_hispan  race_white 
#>           1           2           3           4           4           4 
#>     married    nodegree        re74        re75        re78 
#>           5           6           7           8           9 
#> 

```

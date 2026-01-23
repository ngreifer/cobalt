# Set and Get Options in `cobalt`

Makes it easier to set cobalt options. `set.cobalt.options()` is
essentially a wrapper for
[`options()`](https://rdrr.io/r/base/options.html) but performs several
checks, and `get.cobalt.options()` is essentially a wrapper for
[`getOption()`](https://rdrr.io/r/base/options.html).

## Usage

``` r
set.cobalt.options(..., default = FALSE)

get.cobalt.options(...)
```

## Arguments

- ...:

  For `set.cobalt.options()`,
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  parameters and the values they should take. These should be the name
  of the parameter in
  [`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
  without `"cobalt_"` preceding them. See examples. If any values are
  `NULL`, the corresponding options will be set back to their defaults.

  For `get.cobalt.options()`, one or more strings containing the name of
  a parameter option to be retrieved. See examples. If empty, all
  available options and their values will be returned.

- default:

  if `TRUE`, sets all cobalt options not named in `...` to their default
  values.

## Details

When an option is set to `NULL`, it is set to its default value. The
defaults are not displayed but are listed on the help pages where they
appear. Most options correspond to display options, which can be
accessed
[here](https://ngreifer.github.io/cobalt/reference/display-options.md).
Some others (e.g., `continuous` and `binary`) are described on the
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
help page.

## See also

- [`options()`](https://rdrr.io/r/base/options.html)

- [`display-options`](https://ngreifer.github.io/cobalt/reference/display-options.md)
  for some arguments that can be set via options.

## Examples

``` r
# Set un to be TRUE to always display unadjusted 
# balance measures and set binary to "std" to 
# produce standardized mean differences for 
# binary variables.

set.cobalt.options(un = TRUE, binary = "std")

# Note: the above is equivalent to:
# options(cobalt_un = TRUE, cobalt_binary = "std")
# but performs some additional checks

get.cobalt.options("un", "binary")
#> $un
#> [1] TRUE
#> 
#> $binary
#> [1] "std"
#> 

# Note: the above is equivalent to:
# getOption("cobalt_un")
# getOption("cobalt_binary")

# Return all cobalt options to their defaults

set.cobalt.options(default = TRUE)

# View all available options
get.cobalt.options()
#> $stats
#> NULL
#> 
#> $un
#> NULL
#> 
#> $continuous
#> NULL
#> 
#> $binary
#> NULL
#> 
#> $imbalanced.only
#> NULL
#> 
#> $disp
#> NULL
#> 
#> $disp.means
#> NULL
#> 
#> $disp.sds
#> NULL
#> 
#> $disp.v.ratio
#> NULL
#> 
#> $disp.ks
#> NULL
#> 
#> $disp.subclass
#> NULL
#> 
#> $disp.bal.tab
#> NULL
#> 
#> $cluster.summary
#> NULL
#> 
#> $cluster.fun
#> NULL
#> 
#> $imp.summary
#> NULL
#> 
#> $imp.fun
#> NULL
#> 
#> $multi.summary
#> NULL
#> 
#> $msm.summary
#> NULL
#> 
#> $target.summary
#> NULL
#> 
#> $subclass.summary
#> NULL
#> 
#> $int_sep
#> NULL
#> 
#> $factor_sep
#> NULL
#> 
#> $center
#> NULL
#> 
#> $orth
#> NULL
#> 
#> $remove_perfect_col
#> NULL
#> 
#> $disp.call
#> NULL
#> 
```

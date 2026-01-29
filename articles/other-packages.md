# Using cobalt with Other Preprocessing Packages

This is an appendix to the main vignette, “Covariate Balance Tables and
Plots: A Guide to the cobalt Package”, accessible at
[`vignette("cobalt")`](https://ngreifer.github.io/cobalt/articles/cobalt.md).
It contains descriptions and demonstrations of several utility functions
in *cobalt* and the use of
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
with *twang*, *Matching*, *optmatch*, *CBPS*, *ebal*, *designmatch*,
*sbw*, *MatchThem*, and *cem*. Note that *MatchIt* can perform most of
the functions that *Matching*, *optmatch*, and *cem* can, and *WeightIt*
can perform most of the functions that *twang*, *CBPS*, *ebal*, and
*sbw* can. Because *cobalt* has been optimized to work with *MatchIt*
and *WeightIt*, it is recommended to use those packages to simplify
preprocessing and balance assessment, but we recognize users may prefer
to use the packages described in this vignette.

## Utilities

In addition to its main balance assessment functions, *cobalt* contains
several utility functions. These are meant to reduce the typing and
programming burden that often accompany the use of R with a diverse set
of packages.

### `splitfactor()` and `unsplitfactor()`

Some functions (outside of *cobalt*) are not friendly to factor or
character variables, and require numeric variables to operate correctly.
For example, some regression-style functions, such as `ebalance()` in
*ebal*, can only take in non-singular numeric matrices. Other functions
will process factor variables, but will return output in terms of dummy
coded version of the factors. For example,
[`lm()`](https://rdrr.io/r/stats/lm.html) will create dummy variables
out of a factor and drop the reference category to create regression
coefficients.

To prepare data sets for use in functions that do not allow factors or
to mimic the output of functions that split factor variables, users can
use
[`splitfactor()`](https://ngreifer.github.io/cobalt/reference/splitfactor.md),
which takes in a data set and the names of variables to split, and
outputs a new data set with newly created dummy variables. Below is an
example splitting the `race` variable in the Lalonde data set into
dummies, eliminating the reference category (`"black"`):

``` r
head(lalonde)
```

    ##   treat age educ   race married nodegree re74 re75       re78
    ## 1     1  37   11  black       1        1    0    0  9930.0460
    ## 2     1  22    9 hispan       0        1    0    0  3595.8940
    ## 3     1  30   12  black       0        0    0    0 24909.4500
    ## 4     1  27   11  black       0        1    0    0  7506.1460
    ## 5     1  33    8  black       0        1    0    0   289.7899
    ## 6     1  22    9  black       0        1    0    0  4056.4940

``` r
lalonde.split <- splitfactor(lalonde, "race")
head(lalonde.split)
```

    ##   treat age educ race_hispan race_white married nodegree re74 re75       re78
    ## 1     1  37   11           0          0       1        1    0    0  9930.0460
    ## 2     1  22    9           1          0       0        1    0    0  3595.8940
    ## 3     1  30   12           0          0       0        0    0    0 24909.4500
    ## 4     1  27   11           0          0       0        1    0    0  7506.1460
    ## 5     1  33    8           0          0       0        1    0    0   289.7899
    ## 6     1  22    9           0          0       0        1    0    0  4056.4940

It is possible to undo the action of
[`splitfactor()`](https://ngreifer.github.io/cobalt/reference/splitfactor.md)
with
[`unsplitfactor()`](https://ngreifer.github.io/cobalt/reference/splitfactor.md),
which takes in a data set with dummy variables formed from
[`splitfactor()`](https://ngreifer.github.io/cobalt/reference/splitfactor.md)
or otherwise and recreates the original factor variable. If the
reference category was dropped, its value needs to be supplied.

``` r
lalonde.unsplit <- unsplitfactor(lalonde.split, "race", 
                                 dropped.level = "black")
head(lalonde.unsplit)
```

    ##   treat age educ   race married nodegree re74 re75       re78
    ## 1     1  37   11  black       1        1    0    0  9930.0460
    ## 2     1  22    9 hispan       0        1    0    0  3595.8940
    ## 3     1  30   12  black       0        0    0    0 24909.4500
    ## 4     1  27   11  black       0        1    0    0  7506.1460
    ## 5     1  33    8  black       0        1    0    0   289.7899
    ## 6     1  22    9  black       0        1    0    0  4056.4940

Notice the original data set and the unsplit data set look identical. If
the input to
[`unsplitfactor()`](https://ngreifer.github.io/cobalt/reference/splitfactor.md)
is the output of a call to
[`splitfactor()`](https://ngreifer.github.io/cobalt/reference/splitfactor.md)
(as it was here), you don’t need to tell
[`unsplitfactor()`](https://ngreifer.github.io/cobalt/reference/splitfactor.md)
the name of the split variable or the value of the dropped level. It was
done here for illustration purposes.

### `get.w()`

[`get.w()`](https://ngreifer.github.io/cobalt/reference/get.w.md) allows
users to extract weights from the output of a call to a preprocessing
function in one of the supported packages. Because each package stores
weights in different ways, it can be helpful to have a single function
that applies equally to all outputs. *twang* has a function called
[`get.weights()`](https://rdrr.io/pkg/twang/man/get.weights.html) that
performs the same functions with slightly finer control for the output
of a call to [`ps()`](https://rdrr.io/pkg/twang/man/ps.html).

## `bal.tab()`

The next sections describe the use of
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
with packages other than those described in the main vignette. Even if
you are using
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
with one of these packages, it may be useful to read the main vignette
at
[`vignette("cobalt")`](https://ngreifer.github.io/cobalt/articles/cobalt.md)
to understand
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)’s
main options, which are not detailed here.

### Using `bal.tab()` with *twang*

Generalized boosted modeling (GBM), as implemented in *twang*, can be an
effective way to generate propensity scores and weights for use in
propensity score weighting.
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
functions similarly to the functions
[`bal.table()`](https://rdrr.io/pkg/twang/man/bal.table.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) when used with GBM in
*twang*. Below is a simple example of its use:

``` r
#GBM PS weighting for the ATT
data("lalonde", package = "cobalt") ##If not yet loaded
covs0 <- subset(lalonde, select = -c(treat, re78))
f <- reformulate(names(covs0), "treat")

ps.out <- twang::ps(f, data = lalonde, 
                    stop.method = c("es.mean", "es.max"), 
                    estimand = "ATT", n.trees = 1000,
                    verbose = FALSE)

bal.tab(ps.out, stop.method = "es.mean")
```

    ## Balance Measures
    ##                 Type Diff.Adj
    ## prop.score  Distance   0.5189
    ## age          Contin.   0.0400
    ## educ         Contin.  -0.0819
    ## race_black    Binary   0.0250
    ## race_hispan   Binary  -0.0008
    ## race_white    Binary  -0.0242
    ## married       Binary  -0.0116
    ## nodegree      Binary   0.0864
    ## re74         Contin.   0.0691
    ## re75         Contin.   0.0953
    ## 
    ## Effective sample sizes
    ##            Control Treated
    ## Unadjusted  429.       185
    ## Adjusted     33.03     185

The output looks a bit different from *twang*’s
[`bal.table()`](https://rdrr.io/pkg/twang/man/bal.table.html) output.
First is the original call to
[`ps()`](https://rdrr.io/pkg/twang/man/ps.html). Next is the balance
table containing mean differences for the covariates included in the
input to [`ps()`](https://rdrr.io/pkg/twang/man/ps.html). Last is a
table displaying sample size information, similar to what would be
generated using *twang*’s
[`summary()`](https://rdrr.io/r/base/summary.html) function. The
“effective” sample size is displayed when weighting is used; it is
calculated as is done in *twang*. See the *twang* documentation,
[`?bal.tab`](https://ngreifer.github.io/cobalt/reference/bal.tab.md), or
“Details on Calculations” at
[`vignette("cobalt")`](https://ngreifer.github.io/cobalt/articles/cobalt.md)
for details on this calculation.

When using
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
with *twang*, the user must specify the `ps` object, the output of a
call to [`ps()`](https://rdrr.io/pkg/twang/man/ps.html), as the first
argument. The second argument, `stop.method`, is the name of the stop
method(s) for which balance is to be assessed, since a `ps` object may
contain more than one if so specified.
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
can display the balance for more than one stop method at a time by
specifying a vector of stop method names. If this argument is left empty
or if the argument to `stop.method` does not correspond to any of the
stop methods in the `ps` object,
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
will default to displaying balance for all stop methods available.
Abbreviations are allowed for the stop method, which is not case
sensitive.

The other arguments to
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
when using it with *twang* have the same form and function as those
given when using it without a conditioning package, except for
`s.d.denom`. If the estimand of the stop method used is the ATT,
`s.d.denom` will default to `"treated"` if not specified, and if the
estimand is the ATE, `s.d.denom` will default to `"pooled"`, mimicking
the behavior of *twang*. The user can specify their own argument to
`s.d.denom`, but using the defaults is advised.

If sampling weights are used in the call to
[`ps()`](https://rdrr.io/pkg/twang/man/ps.html), they will be
automatically incorporated into the
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
calculations for both the adjusted and unadjusted samples, just as
*twang* does.

`mnps` objects resulting from fitting models in *twang* with
multi-category treatments are also compatible with *cobalt*. See the
section “Using *cobalt* with multi-category treatments” at
[`vignette("cobalt")`](https://ngreifer.github.io/cobalt/articles/cobalt.md).
`iptw` objects resulting from fitting models in *twang* with
longitudinal treatments are also compatible with *cobalt*. See
[`vignette("longitudinal-treat")`](https://ngreifer.github.io/cobalt/articles/longitudinal-treat.md).
`ps.cont` objects resulting from using `ps.cont()` in `twangContinuous`,
which implements GBM for continuous treatments, are also compatible. See
the section “Using *cobalt* with continuous treatments” at
[`vignette("cobalt")`](https://ngreifer.github.io/cobalt/articles/cobalt.md).

### Using `bal.tab()` with *Matching*

The *Matching* package is used for propensity score matching, and was
also the first package to implement genetic matching. *MatchIt* calls
*Matching* to use genetic matching and can accomplish many of the
matching methods *Matching* can, but *Matching* is still a widely used
package with its own strengths.
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
functions similarly to *Matching*’s
[`MatchBalance()`](https://rdrr.io/pkg/Matching/man/MatchBalance.html)
command, which yields a thorough presentation of balance. Below is a
simple example of the use of
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
with *Matching*:

``` r
#1:1 NN PS matching w/ replacement
data("lalonde", package = "cobalt") #If not yet loaded
covs0 <- subset(lalonde, select = -c(treat, re78))
f <- reformulate(names(covs0), "treat")

fit <- glm(f, data = lalonde, family = binomial)
p.score <- fit$fitted.values
match.out <- Matching::Match(Tr = lalonde$treat, X = p.score,
                             estimand = "ATT")

bal.tab(match.out, formula = f, data = lalonde,
        distance = ~ p.score)
```

    ## Balance Measures
    ##                 Type Diff.Adj
    ## p.score     Distance   0.0043
    ## age          Contin.   0.2106
    ## educ         Contin.   0.0201
    ## race_black    Binary   0.0054
    ## race_hispan   Binary  -0.0051
    ## race_white    Binary  -0.0003
    ## married       Binary   0.0661
    ## nodegree      Binary  -0.0079
    ## re74         Contin.  -0.0772
    ## re75         Contin.  -0.0127
    ## 
    ## Sample sizes
    ##                      Control Treated
    ## All                   429.       185
    ## Matched (ESS)          49.17     185
    ## Matched (Unweighted)  136.       185
    ## Unmatched             293.         0

The output looks quite different from *Matching*’s
[`MatchBalance()`](https://rdrr.io/pkg/Matching/man/MatchBalance.html)
output. Rather than being stacked vertically, balance statistics are
arranged horizontally in a table format, allowing for quick balance
checking. Below the balance table is a summary of the sample size before
and after matching, similar to what *Matching*’s
[`summary()`](https://rdrr.io/r/base/summary.html) command would
display. The sample size can include an “ESS” and “unweighted” value;
the “ESS” value is the effective sample size resulting from the matching
weights, while the “unweighted” is the count of units with nonzero
matching weights.

The input to
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md) is
similar to that given to
[`MatchBalance()`](https://rdrr.io/pkg/Matching/man/MatchBalance.html):
the `Match` object resulting from the call to
[`Match()`](https://rdrr.io/pkg/Matching/man/Match.html), a formula
relating treatment to the covariates for which balance is to be
assessed, and the original data set. This is not the only way to call
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md):
instead of a formula and a data set, one can also input a data frame of
covariates and a vector of treatment status indicators, just as when
using
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
without a conditioning package. For example, the code below will yield
the same results as the call to
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
above:

``` r
bal.tab(match.out, treat = lalonde$treat, covs = covs0,
        distance = ~ p.score)
```

The other arguments to
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
when using it with *Matching* have the same form and function as those
given when using it without a conditioning package, except for
`s.d.denom`. If the estimand of the original call to
[`Match()`](https://rdrr.io/pkg/Matching/man/Match.html) is the ATT,
`s.d.denom` will default to `"treated"` if not specified; if the
estimand is the ATE, `s.d.denom` will default to `"pooled"`; if the
estimand is the ATC, `s.d.denom` will default to `"control"`. The user
can specify their own argument to `s.d.denom`, but using the defaults is
advisable. In addition, the use of the `addl` argument is unnecessary
because the covariates are entered manually as arguments, so all
covariates for which balance is to be assessed can be entered through
the `formula` or `covs` argument. If the covariates are stored in two
separate data frames, it may be useful to include one in `formula` or
`covs` and the other in `addl`.

### Using `bal.tab()` with *optmatch*

The *optmatch* package is useful for performing optimal pairwise or full
matching. Most functions in *optmatch* are subsumed in *MatchIt*, but
*optmatch* sees use from those who want finer control of the matching
process than *MatchIt* allows. The output of calls to functions in
*optmatch* is an *optmatch* object, which contains matching stratum
membership for each unit in the given data set. Units that are matched
with each other are assigned the same matching stratum. The user guide
for *optmatch* recommends using the *RItools* package for balance
assessment, but below is an example of how to use
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
for the same purpose. Note that some results will differ between
*cobalt* and *RItools* because of differences in how balance is
calculated in each.

``` r
#Optimal full matching on the propensity score
data("lalonde", package = "cobalt") #If not yet loaded
covs0 <- subset(lalonde, select = -c(treat, re78))
f <- reformulate(names(covs0), "treat")

fit <- glm(f, data = lalonde, family = binomial)
p.score <- fit$fitted.values #get the propensity score
fm <- optmatch::fullmatch(treat ~ p.score, data = lalonde)

bal.tab(fm, covs = covs0, distance = ~ p.score)
```

    ## Balance Measures
    ##                 Type Diff.Adj
    ## p.score     Distance   0.0053
    ## age          Contin.   0.1479
    ## educ         Contin.  -0.0158
    ## race_black    Binary   0.0086
    ## race_hispan   Binary  -0.0062
    ## race_white    Binary  -0.0024
    ## married       Binary   0.0573
    ## nodegree      Binary   0.0037
    ## re74         Contin.  -0.0696
    ## re75         Contin.  -0.0175
    ## 
    ## Sample sizes
    ##                      Control Treated
    ## All                   429.       185
    ## Matched (ESS)          51.53     185
    ## Matched (Unweighted)  429.       185

Most details for the use of
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
with *optmatch* are similar to those when using
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
with *Matching*. Users can enter either a formula and a data set or a
vector of treatment status and a set of covariates. Unlike with
*Matching*, entering the treatment variable is optional as it is already
stored in the *optmatch* object.
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md) is
compatible with both `pairmatch()` and `fullmatch()` output.

### Using `bal.tab()` with *CBPS*

The *CBPS* (Covariate Balancing Propensity Score) package is a great
tool for generating covariate balancing propensity scores, a class of
propensity scores that are quite effective at balancing covariates among
groups. *CBPS* includes functions for estimating propensity scores for
binary, multi-category, and continuous treatments.
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
functions similarly to *CBPS*’s
[`balance()`](https://rdrr.io/pkg/CBPS/man/balance.html) command. Below
is a simple example of its use with a binary treatment:

``` r
#CBPS weighting
data("lalonde", package = "cobalt") #If not yet loaded
covs0 <- subset(lalonde, select = -c(treat, re78))
f <- reformulate(names(covs0), "treat")

#Generating covariate balancing propensity score weights for ATT
cbps.out <- CBPS::CBPS(f, data = lalonde)
```

    ## [1] "Finding ATT with T=1 as the treatment.  Set ATT=2 to find ATT with T=0 as the treatment"

``` r
bal.tab(cbps.out)
```

    ## Balance Measures
    ##                 Type Diff.Adj
    ## prop.score  Distance   0.0244
    ## age          Contin.   0.1055
    ## educ         Contin.  -0.0246
    ## race_black    Binary   0.0126
    ## race_hispan   Binary  -0.0036
    ## race_white    Binary  -0.0090
    ## married       Binary   0.0114
    ## nodegree      Binary   0.0204
    ## re74         Contin.  -0.0145
    ## re75         Contin.   0.0043
    ## 
    ## Effective sample sizes
    ##            Control Treated
    ## Unadjusted  429.       185
    ## Adjusted    104.16     185

First is the original call to
[`CBPS()`](https://rdrr.io/pkg/CBPS/man/CBPS.html). Next is the balance
table containing mean differences for the covariates included in the
input to [`CBPS()`](https://rdrr.io/pkg/CBPS/man/CBPS.html). Last is a
table displaying sample size information. The “effective” sample size is
displayed when weighting (rather than matching or subclassification) is
used; it is calculated as is done in *twang*. See the *twang*
documentation,
[`?bal.tab`](https://ngreifer.github.io/cobalt/reference/bal.tab.md), or
“Details on Calculations” at
[`vignette("cobalt")`](https://ngreifer.github.io/cobalt/articles/cobalt.md)
for details on this calculation.

The other arguments to
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
when using it with *CBPS* have the same form and function as those given
when using it without a conditioning package, except for `s.d.denom`. If
the estimand of the original call to
[`CBPS()`](https://rdrr.io/pkg/CBPS/man/CBPS.html) is the ATT,
`s.d.denom` will default to `"treated"` if not specified, and if the
estimand is the ATE, `s.d.denom` will default to `"pooled"`. The user
can specify their own argument to `s.d.denom`, but using the defaults is
advisable.

`CBPSContinuous` objects resulting from fitting models in *CBPS* with
continuous treatments are also compatible with *cobalt*. See the section
“Using *cobalt* with continuous treatments” at
[`vignette("cobalt")`](https://ngreifer.github.io/cobalt/articles/cobalt.md).
*CBPS* objects resulting from fitting models in *CBPS* with
multi-category treatments are also compatible with *cobalt*. See the
section “Using *cobalt* with multi-category treatments” at
[`vignette("cobalt")`](https://ngreifer.github.io/cobalt/articles/cobalt.md).
`CBMSM` objects resulting from fitting models in *CBPS* with
longitudinal treatments are also compatible with *cobalt*. See
[`vignette("longitudinal-treat")`](https://ngreifer.github.io/cobalt/articles/longitudinal-treat.md).

### Using `bal.tab()` with *ebal*

The *ebal* package implements entropy balancing, a method of weighting
for the ATT that yields perfect balance on all desired moments of the
covariate distributions between groups. Rather than estimate a
propensity score, entropy balancing generates weights directly that
satisfy a user-defined moment condition, specifying which moments are to
be balanced. Note that all the functionality of *ebal* is contained
within *WeightIt*. *ebal* does not have its own balance assessment
function; thus, *cobalt* is the only way to assess balance without
programming, which the *ebal* documentation instructs. Below is a simple
example of using
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
with *ebal*:

``` r
#Entropy balancing
data("lalonde", package = "cobalt") #If not yet loaded
covs0 <- subset(lalonde, select = -c(treat, re78, race))

#Generating entropy balancing weights
e.out <- ebal::ebalance(lalonde$treat, covs0)
```

    ## Converged within tolerance

``` r
bal.tab(e.out, treat = lalonde$treat, covs = covs0)
```

    ## Balance Measures
    ##             Type Diff.Adj
    ## age      Contin.       -0
    ## educ     Contin.       -0
    ## married   Binary       -0
    ## nodegree  Binary        0
    ## re74     Contin.       -0
    ## re75     Contin.       -0
    ## 
    ## Effective sample sizes
    ##            Control Treated
    ## Unadjusted  429.       185
    ## Adjusted    247.64     185

First is the balance table containing mean differences for covariates
included in the original call to `ebalance()`. In general, these will
all be very close to 0. Next is a table displaying effective sample size
information. See
[`?bal.tab`](https://ngreifer.github.io/cobalt/reference/bal.tab.md) or
“Details on Calculations” at
[`vignette("cobalt")`](https://ngreifer.github.io/cobalt/articles/cobalt.md)
for details on this calculation. A common issue when using entropy
balancing is small effective sample size, which can yield low precision
in effect estimation when using weighted regression, so it is important
that users pay attention to this measure.

The input is similar to that for using
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
with *Matching* or *optmatch*. In addition to the `ebalance` object, one
must specify either both a formula and a data set or both a treatment
vector and a data frame of covariates.

### Using `bal.tab()` with *designmatch*

The *designmatch* package implements various matching methods that use
optimization to find matches that satisfy certain balance constraints.
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
functions similarly to *designmatch*’s `meantab()` command but provides
additional flexibility and convenience. Below is a simple example of
using
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
with *designmatch*:

``` r
#Mixed integer programming matching
library("designmatch")
data("lalonde", package = "cobalt") #If not yet loaded
covs0 <- subset(lalonde, select = -c(treat, re78, race))

#Matching for balance on covariates
dmout <- bmatch(lalonde$treat,
                dist_mat = NULL,
                subset_weight = NULL,
                mom = list(covs = covs0,
                           tols = absstddif(covs0, lalonde$treat, .05)),
                n_controls = 1,
                total_groups = 185)

bal.tab(dmout, treat = lalonde$treat, covs = covs0)
```

The input is similar to that for using
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
with *Matching* or *optmatch*. In addition to the `bmatch()` output
object, one must specify either both a formula and a data set or both a
treatment vector and a data frame of covariates. The output is similar
to that of *optmatch*.

### Using `bal.tab()` with *sbw*

The *sbw* package implements optimization-based weighting to estimate
weights that satisfy certain balance constraints and have minimal
variance.
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
functions similarly to *sbw*’s
[`summarize()`](https://rdrr.io/pkg/sbw/man/summarize.html) function but
provides additional flexibility and convenience. Below is a simple
example of using
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
with *sbw*:

``` r
#Optimization-based weighting
data("lalonde", package = "cobalt") #If not yet loaded
lalonde_split <- splitfactor(lalonde, drop.first = "if2")
cov.names <- setdiff(names(lalonde_split), c("treat", "re78"))

#Estimating balancing weights for the ATT
sbw.out <- sbw::sbw(lalonde_split,
                    ind = "treat",
                    bal = list(bal_cov = cov.names,
                               bal_alg = FALSE, 
                               bal_tol = .001),
                    par = list(par_est = "att"))
```

    ##   Building the weighting problem... 
    ##   quadprog optimizer is open... 
    ##   Finding the optimal weights... 
    ##   Optimal weights found.

``` r
bal.tab(sbw.out, un = TRUE, disp.means = TRUE)
```

    ## Balance Measures
    ##                Type    M.0.Un    M.1.Un Diff.Un   M.0.Adj   M.1.Adj Diff.Adj
    ## age         Contin.   28.0303   25.8162 -0.3094   25.8054   25.8162   0.0015
    ## educ        Contin.   10.2354   10.3459  0.0550   10.3431   10.3459   0.0014
    ## race_black   Binary    0.2028    0.8432  0.6404    0.8428    0.8432   0.0004
    ## race_hispan  Binary    0.1422    0.0595 -0.0827    0.0594    0.0595   0.0001
    ## race_white   Binary    0.6550    0.0973 -0.5577    0.0978    0.0973  -0.0005
    ## married      Binary    0.5128    0.1892 -0.3236    0.1897    0.1892  -0.0005
    ## nodegree     Binary    0.5967    0.7081  0.1114    0.7076    0.7081   0.0005
    ## re74        Contin. 5619.2365 2095.5737 -0.7211 2102.3624 2095.5737  -0.0014
    ## re75        Contin. 2466.4844 1532.0553 -0.2903 1528.7633 1532.0553   0.0010
    ## 
    ## Effective sample sizes
    ##            Control Treated
    ## Unadjusted  429.       185
    ## Adjusted    108.99     185

The output is similar to the output of a call to
[`summarize()`](https://rdrr.io/pkg/sbw/man/summarize.html). Rather than
stack several balance tables vertically, each with their own balance
summary, here they are displayed horizontally. Note that due to
differences in how *sbw* and *cobalt* compute the standardization factor
in the standardized mean difference, values may not be identical between
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
and [`summarize()`](https://rdrr.io/pkg/sbw/man/summarize.html). Also
note that
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)’s
default is to display raw rather than standardized mean differences for
binary variables.

### Using `bal.tab()` with *MatchThem*

The *MatchThem* package is essentially a wrapper for
[`matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.html)
from *MatchIt* and
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.html)
from *WeightIt* but for use with multiply imputed data. Using
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md) on
`mimids` or `wimids` objects from *MatchThem* activates the features
that accompany multiply imputed data; balance is assessed within each
imputed dataset and aggregated across imputations. See `?bal.tab.imp` or
[`vignette("segmented-data")`](https://ngreifer.github.io/cobalt/articles/segmented-data.md)
for more information about using *cobalt* with multiply imputed data.
Below is a simple example of using
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
with *MatchThem*:

``` r
#PS weighting on multiply imputed data
data("lalonde_mis", package = "cobalt")

#Generate imputed data sets
m <- 10 #number of imputed data sets
imp.out <- mice::mice(lalonde_mis, m = m, print = FALSE) 

#Matching for balance on covariates
mt.out <- MatchThem::matchthem(treat ~ age + educ + married +
                                   race + re74 + re75, 
                               datasets = imp.out,
                               approach = "within", 
                               method = "nearest",
                               estimand = "ATT")

bal.tab(mt.out)
```

    ## Balance summary across all imputations
    ##                 Type Min.Diff.Adj Mean.Diff.Adj Max.Diff.Adj
    ## distance    Distance       0.9503        0.9579       0.9644
    ## age          Contin.      -0.0423        0.0199       0.0884
    ## educ         Contin.      -0.1855       -0.1449      -0.1048
    ## married       Binary      -0.0270       -0.0151       0.0000
    ## race_black    Binary       0.3730        0.3730       0.3730
    ## race_hispan   Binary      -0.1892       -0.1768      -0.1568
    ## race_white    Binary      -0.2162       -0.1962      -0.1838
    ## re74         Contin.      -0.0799       -0.0600      -0.0355
    ## re75         Contin.      -0.0960       -0.0576      -0.0370
    ## 
    ## Average sample sizes across imputations
    ##             0   1
    ## All       429 185
    ## Matched   185 185
    ## Unmatched 244   0

``` r
#Weighting for balance on covariates
wt.out <- MatchThem::weightthem(treat ~ age + educ + married +
                                    race + re74 + re75, 
                                datasets = imp.out,
                                approach = "within", 
                                method = "glm",
                                estimand = "ATE")

bal.tab(wt.out)
```

    ## Balance summary across all imputations
    ##                 Type Min.Diff.Adj Mean.Diff.Adj Max.Diff.Adj
    ## prop.score  Distance       0.1430        0.1547       0.1618
    ## age          Contin.      -0.1939       -0.1886      -0.1830
    ## educ         Contin.       0.0769        0.0834       0.0913
    ## married       Binary      -0.1052       -0.0992      -0.0892
    ## race_black    Binary       0.0514        0.0563       0.0611
    ## race_hispan   Binary       0.0048        0.0102       0.0125
    ## race_white    Binary      -0.0733       -0.0665      -0.0596
    ## re74         Contin.      -0.3094       -0.2839      -0.2585
    ## re75         Contin.      -0.1680       -0.1600      -0.1486
    ## 
    ## Average effective sample sizes across imputations
    ##                 0      1
    ## Unadjusted 429.   185.  
    ## Adjusted   330.76  65.58

The input is similar to that for using
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
with *MatchIt* or *WeightIt*.

### Using `bal.tab()` with *cem*

The *cem* package implements coarsened exact matching for binary and
multi-category treatments.
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
functions similarly to `cems`’s
[`imbalance()`](https://rdrr.io/pkg/cem/man/imbalance.html). Below is a
simple example of using
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
with *cem*:

``` r
#Coarsened exact matching
data("lalonde", package = "cobalt") #If not yet loaded

#Matching for balance on covariates
cem.out <- cem::cem("treat", data = lalonde, drop = "re78")

bal.tab(cem.out, data = lalonde, stats = c("m", "ks"))
```

The input is similar to that for using
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
with *Matching* or *optmatch*. In addition to the
[`cem()`](https://rdrr.io/pkg/cem/man/cem.html) output object, one must
specify either both a formula and a data set or both a treatment vector
and a data frame of covariates. Unlike with *Matching*, entering the
treatment variable is optional as it is already stored in the output
object. The output is similar to that of *optmatch*.

When using [`cem()`](https://rdrr.io/pkg/cem/man/cem.html) with multiply
imputed data (i.e., by supplying a list of data.frames to the `datalist`
argument in [`cem()`](https://rdrr.io/pkg/cem/man/cem.html)), an
argument to `imp` should be specified to
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md) or
a `mids` object from the *mice* package should be given as the argument
to `data`. See `?bal.tab.imp` or
[`vignette("segmented-data")`](https://ngreifer.github.io/cobalt/articles/segmented-data.md)
for more information about using *cobalt* with multiply imputed data.
Below is an example of using *cem* with multiply imputed data from
*mice*:

``` r
#Coarsened exact matching on multiply imputed data
data("lalonde_mis", package = "cobalt")

#Generate imputed data sets
m <- 10 #number of imputed data sets
imp.out <- mice::mice(lalonde_mis, m = m, print = FALSE) 
imp.data.list <- mice::complete(imp.out, "all")

#Match within each imputed dataset
cem.out.imp <- cem::cem("treat", datalist = imp.data.list,
                        drop = "re78")

bal.tab(cem.out.imp, data = imp.out)
```

### Using `bal.tab()` with other packages

It is possible to use
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
with objects that don’t come from these packages using the `default`
method. If an object that doesn’t correspond to the output from one of
the specifically supported packages is passed as the first argument to
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md),
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
will do its best to process that object as if it did come from a
supported package. It will search through the components of the object
for items with names like `"treat"`, `"covs"`, `"data"`, `"weights"`,
etc., that have the correct object types. Any additional arguments can
be specified by the user.

The goal of the `default` method is to allow package authors to rely on
*cobalt* as a substitute for any balancing function they might otherwise
write. By ensuring compatibility with the `default` method, package
authors can have their users simply supply the output of a compatible
function into *cobalt* functions without having to write a specific
method in *cobalt*. A package author would need to make sure the output
of their package contained enough information with correctly named
components; if so, *cobalt* functions can be used as conveniently with
the output as it is with specifically supported packages.

Below, we demonstrate this capability with the output of *optweight*,
which performs a version of propensity score weighting using
optimization, similar to *sbw*. No
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
method has been written with *optweight* output in mind; rather,
*optweight* was written to have output compatible with the `default`
method of
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md).

``` r
#Optimization-based weighting
data("lalonde", package = "cobalt")

#Estimate the weights using optimization
ow.out <- optweight::optweight(treat ~ age + educ + married + race +
                                 re74 + re75,
                               data = lalonde,
                               estimand = "ATE",
                               tols = .01)

#Note the contents of the output object:
names(ow.out)
```

    ##  [1] "weights"     "treat"       "covs"        "s.weights"   "b.weights"  
    ##  [6] "estimand"    "focal"       "norm"        "call"        "tols"       
    ## [11] "target.tols" "duals"       "info"        "solver"

``` r
#Use bal.tab() directly on the output
bal.tab(ow.out)
```

    ## Balance Measures
    ##                Type Diff.Adj
    ## age         Contin.    -0.01
    ## educ        Contin.     0.01
    ## married      Binary    -0.01
    ## race_black   Binary     0.01
    ## race_hispan  Binary    -0.00
    ## race_white   Binary    -0.01
    ## re74        Contin.    -0.01
    ## re75        Contin.     0.01
    ## 
    ## Effective sample sizes
    ##            Control Treated
    ## Unadjusted  429.    185.  
    ## Adjusted    349.26   52.04

The output is treated as output from a specifically supported package.
See
[`?bal.tab.default`](https://ngreifer.github.io/cobalt/reference/bal.tab.default.md)
for more details and another example.

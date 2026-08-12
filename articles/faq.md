# Frequently Asked Questions

### How are standardized mean differences computed in *cobalt*?

Most questions below are related to this one, so here I will try to
explain in complete detail how standardized mean differences (SMDs) are
computed.

First, it is important to know that by default, **mean differences for
binary covariates are not standardized**. That means in the `Diff.Adj`
column, etc., what you are seeing the is *raw difference in proportion*
for binary variables. (By raw, I mean unstandardized, but weights may be
applied if relevant.) To request SMDs for binary covariates, set
`binary = "std"` in the call to
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md) or
[`love.plot()`](https://ngreifer.github.io/cobalt/reference/love.plot.md).
See also the question below.

For continuous covariates, the standardized mean difference is computed
as \\ \text{SMD} = \frac{\bar{x}\_1 - \bar{x}\_0}{s^\*} \\ where
\\\bar{x}\_1\\ is the mean of the covariate in the treated group,
\\\bar{x}\_0\\ is the mean of the covariate in the control group, and
\\s^\*\\ is a standardization factor (not necessarily a standard
deviation!). After matching or weighting, the weighted standardized mean
difference is computed as \\ \text{SMD}^w = \frac{\bar{x}\_{1w} -
\bar{x}\_{0w}}{s^\*} \\ where \\\bar{x}\_{1w}\\ is the weighted mean of
the covariate in the treated group, i.e., \\\bar{x}\_{1w} =
\frac{1}{\sum\_{i:A_i = 1}{w_i}}\sum\_{i:A_i=1}{w_ix_i}\\, and similarly
for the control group. Critically, the standardization factor \\s^\*\\
is the same before and after weighting. I will repeat, **the
standardization factor** \\s^\*\\ **is the same before and after
weighting**. I don’t mean it has the same formula, it mean it is
literally the same value. I explain in more detail in a question below
why this is the case.

How is the standardization factor computed? This depends on the argument
to `s.d.denom` supplied to
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md) or
[`love.plot()`](https://ngreifer.github.io/cobalt/reference/love.plot.md).
When `s.d.denom` is not supplied, this is determined by the argument
supplied to `estimand`, and when that is not supplied, the estimand is
guessed based on the form of the weights, if any. By default, with no
weights supplied and no argument to `s.d.denom` or `estimand`,
`s.d.denom` is set to `"pooled"`, and a note will appear saying so. That
note doesn’t appear if weights are supplied or balance is assessed on
the output of a another package, as the estimand, and therefore
`s.d.denom`, can be determined automatically.

Below are the formulas for the standardization factor corresponding to
each value of `s.d.denom`:

- `"pooled"`: \\s^\* = \sqrt{\frac{s_1^2 + s_0^2}{2}}\\
- `"treated"`: \\s^\* = s_1\\
- `"control"`: \\s^\* = s_0\\
- `"all"`: \\s^\* = s\\
- `"weighted"`: \\s^\* = s_w\\
- `"hedges"`: \\s^\* = \frac{1}{1 - \frac{3}{4(n - 2) -
  1}}\sqrt{\frac{(n_1 - 1)s_1^2 + (n_0 - 1)s_0^2}{n - 2}}\\

where \\s_1\\ is the standard deviation of the treated group, \\s_0\\ is
the standard deviation of the control group, \\s\\ is the standard
deviation of the whole sample ignoring treatment group membership,
\\s_w\\ is the weighted standard deviation of the whole sample, \\n_1\\
and \\n_0\\ are sizes of the treated and control groups, respectively,
and \\n = n_1 + n_0\\.[^1] For continuous covariates, the unweighted
standard deviation is computed as usual, i.e., as \\ s =
\sqrt{\frac{1}{n-1}\sum_i{(x_i - \bar{x})^2}} \\ and the weighted
standard deviation is computed as \\ s = \sqrt{\frac{\sum\_{i}
w\_{i}}{(\sum\_{i} w\_{i})^2 - \sum\_{i=1}^{n} w^2\_{i}}\sum_i{w_i(x_i -
\bar{x}\_w)^2}} \\ For binary covariates, the unweighted standard
deviation is computed as \\s = \sqrt{\bar{x}(1-\bar{x})}\\ and the
weighted standard deviation is computed as \\s =
\sqrt{\bar{x}\_w(1-\bar{x}\_w)}\\.

When sampling weights are supplied, all standard deviations in the
standardization factor are computed incorporating the sampling weights.
When `s.d.denom = "weighted"`, the standardization factor is computed
using the weights used to balance the sample (i.e., the matching or
weighting weights), even for the unadjusted sample. Remember, the
standardization factor is ALWAYS the same before and after adjustment.

I know some of these formulas seem overly complicated for such simple
statistics, but they are required to keep things consistent and not
dependent on the scale of the weights.

### Why are mean differences not standardized for binary covariates?

Ultimately, bias in the treatment effect estimate is a function of
imbalance. That bias is indifferent to whether you measure that
imbalance using a standardized or unstandardized mean difference. The
reason we use SMDs is that covariates naturally are on a variety of
different scales, and when trying to quickly assess whether a sample is
balanced, it is productive to unify the scales of the covariates. That
way, balance on a covariate measured with large numbers (e.g., days in
hospital or prior earnings in dollars) can be assessed alongside balance
on a covariate measured with small numbers (e.g., number of
comorbidities or years of education).

With binary covariates, though, they are already on a comprehensible
scale, so there is no need to standardize. In addition, the scale is
intuitive for people; a difference in proportion of .1 when both groups
have 100 people means that there is an imbalance of 10 people on the
covariate. Are 10 people being different enough to cause bias in the
estimate? That can be assessed substantively without needing to take the
additional step of translating the variable’s scale into something
meaningful.

Another important reason why mean difference are not standardized is
that it is possible for two covariates with the same imbalance to have
vastly different mean differences. For example, consider the following
dataset. `X1` and `X2` both have a mean difference of .1; if they both
affected the outcome equally, then each would contribute to the bias in
the estimate to the same extent.

``` r

treat <- rep(1:0, each = 20)
X1 <- c(rep(0:1, c(1, 19)), rep(0:1, c(3, 17)))
X2 <- c(rep(0:1, c(9, 11)), rep(0:1, c(11, 9)))

bal.tab(treat ~ X1 + X2,
        binary = "raw",
        disp = "means",
        s.d.denom = "treated")
#> Balance Measures
#>      Type M.0.Un M.1.Un Diff.Un
#> X1 Binary   0.85   0.95     0.1
#> X2 Binary   0.45   0.55     0.1
#> 
#> Sample sizes
#>     Control Treated
#> All      20      20
```

But if we standardized the mean differences, not only do we move away
from an actually interpretable statistic (i.e., what does it mean to
divide by the standard deviation of a binary variable?), we see that the
standardized mean differences vary by a huge amount, with `X1` having
twice the imbalance of `X2`.

``` r

bal.tab(treat ~ X1 + X2,
        binary = "std",
        s.d.denom = "treated")
#> Balance Measures
#>      Type Diff.Un
#> X1 Binary  0.4588
#> X2 Binary  0.2010
#> 
#> Sample sizes
#>     Control Treated
#> All      20      20
```

Why does this happen? The standard deviation of a binary variable is a
function of its mean (in particular, it is \\s = \sqrt{p(1-p)}\\) where
\\p\\ is the mean of the variable). That means information about the
mean of the variable, which is unrelated to imbalance, contaminates the
standardized mean difference, which is supposed to measure imbalance. In
this case, standardizing the mean difference only adds confusion and
reduces interpretability. That is why mean differences for binary
variables are unstandardized by default.

You can always change this by setting `binary = "std"` in the call to
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md) or
setting `set.cobalt.options(binary = "std")` to change the option for
the whole session. One advantage of using standardized mean differences
for binary variables is that they are always larger than the raw mean
difference (because the standardization factor is always less than 1),
which means if you use the standardized mean difference as your balance
criterion, you will always seek better balance than using the raw mean
differences. The balance statistics computed by
[`bal.compute()`](https://ngreifer.github.io/cobalt/reference/bal.compute.md)
that involve the standardized mean difference standardize all variables,
including binary variables.

### Why do you use the same standardization factor before and after adjustment?

It is important to remember that bias is a function of the difference in
means of a covariate, and standardization is a just tool to aid in
balance assessment. As a tool, it should reflect imbalance accurately
(i.e., without incorporating extraneous information), but there is no
statistical “truth” about the nature of the standardization factor. I
use the same standardization factor before and after adjustment as
recommended by Stuart
([2008](#ref-stuartDevelopingPracticalRecommendations2008)). The
rationale is that by isolating the SMD to reflect changes in the
difference in means, one can more accurately assess improvement in
balance rather than combining information about the difference in means
with information about the variability of the covariate, which may
change in a variety of ways after adjustment. I describe a specific
example of how allowing the standardization factor to change can cause
problems [here](https://stats.stackexchange.com/a/565705/116195).

### How do I extract the balance tables from the `bal.tab()` object?

Use [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) to
get the balance statistics as data, or
[`format()`](https://rdrr.io/r/base/format.html) to get the table
exactly as [`print()`](https://rdrr.io/r/base/print.html) displays it.

``` r

data("lalonde")

b <- bal.tab(treat ~ age + educ + race + married + re74,
             data = lalonde, s.d.denom = "treated",
             disp = "means", stats = c("m", "v"),
             thresholds = c(m = .1))

head(as.data.frame(b))
#>   variable    type     sample            stat   group   estimate
#> 1      age Contin. Unadjusted            mean Control 28.0303030
#> 2      age Contin. Unadjusted            mean Treated 25.8162162
#> 3      age Contin. Unadjusted      mean.diffs    <NA> -0.3094453
#> 4      age Contin. Unadjusted variance.ratios    <NA>  0.4399955
#> 5     educ Contin. Unadjusted            mean Control 10.2354312
#> 6     educ Contin. Unadjusted            mean Treated 10.3459459
#>            threshold threshold.value
#> 1               <NA>              NA
#> 2               <NA>              NA
#> 3 Not Balanced, >0.1             0.1
#> 4               <NA>              NA
#> 5               <NA>              NA
#> 6               <NA>              NA
```

The result is tidy: one row per covariate, sample, and statistic, with
the balance verdict and the threshold alongside. Those last two columns
appear only when a threshold is on display, since they would otherwise
be empty throughout; the rest are always there. That shape is the same
whatever the input, which makes it convenient for plotting or for your
own aggregation. When the data are segmented—by cluster, imputation,
treatment pair, time point, or subclass—each level of segmentation
becomes a column rather than a nested list, so the result is always a
single rectangle.

Set `wide = TRUE` to get the layout
[`print()`](https://rdrr.io/r/base/print.html) uses instead, with one
column per sample and statistic:

``` r

as.data.frame(b, wide = TRUE)
#>      variable    Type       M.0.Un       M.1.Un     Diff.Un     M.Threshold.Un
#> 1         age Contin.   28.0303030 2.581622e+01 -0.30944526 Not Balanced, >0.1
#> 2        educ Contin.   10.2354312 1.034595e+01  0.05496466     Balanced, <0.1
#> 3  race_black  Binary    0.2027972 8.432432e-01  0.64044604 Not Balanced, >0.1
#> 4 race_hispan  Binary    0.1421911 5.945946e-02 -0.08273168     Balanced, <0.1
#> 5  race_white  Binary    0.6550117 9.729730e-02 -0.55771436 Not Balanced, >0.1
#> 6     married  Binary    0.5128205 1.891892e-01 -0.32363132 Not Balanced, >0.1
#> 7        re74 Contin. 5619.2365064 2.095574e+03 -0.72108381 Not Balanced, >0.1
#>   V.Ratio.Un
#> 1  0.4399955
#> 2  0.4958934
#> 3         NA
#> 4         NA
#> 5         NA
#> 6         NA
#> 7  0.5181285
```

[`format()`](https://rdrr.io/r/base/format.html) returns the printed
table as a data frame of formatted strings—rounded, decimal-aligned,
with missing values shown as `.`—which is what you want when the
destination is a document rather than more code:

``` r

format(b)
#>                Type    M.0.Un    M.1.Un Diff.Un     M.Threshold.Un V.Ratio.Un
#> age         Contin.   28.0303   25.8162 -0.3094 Not Balanced, >0.1     0.4400
#> educ        Contin.   10.2354   10.3459  0.0550     Balanced, <0.1     0.4959
#> race_black   Binary    0.2028    0.8432  0.6404 Not Balanced, >0.1          .
#> race_hispan  Binary    0.1422    0.0595 -0.0827     Balanced, <0.1          .
#> race_white   Binary    0.6550    0.0973 -0.5577 Not Balanced, >0.1          .
#> married      Binary    0.5128    0.1892 -0.3236 Not Balanced, >0.1          .
#> re74        Contin. 5619.2365 2095.5737 -0.7211 Not Balanced, >0.1     0.5181
```

Because it is an ordinary data frame, it goes straight into any table
renderer, so a balance table in a Quarto or R Markdown document is one
line:

``` r

knitr::kable(format(b))
```

Use `component = "observations"` for the sample size table.

Both functions accept every argument
[`print()`](https://rdrr.io/r/base/print.html) accepts—`stats`, `disp`,
`un`, `imbalanced.only`, `disp.thresholds`, the `which.*` arguments—and
resolve them the same way, so you can select what to report without
recomputing anything:

``` r

as.data.frame(b, stats = "v", un = FALSE)
#>       variable    type     sample            stat   group     estimate
#> 1          age Contin. Unadjusted            mean Control 2.803030e+01
#> 2          age Contin. Unadjusted            mean Treated 2.581622e+01
#> 3          age Contin. Unadjusted variance.ratios    <NA> 4.399955e-01
#> 4         educ Contin. Unadjusted            mean Control 1.023543e+01
#> 5         educ Contin. Unadjusted            mean Treated 1.034595e+01
#> 6         educ Contin. Unadjusted variance.ratios    <NA> 4.958934e-01
#> 7   race_black  Binary Unadjusted            mean Control 2.027972e-01
#> 8   race_black  Binary Unadjusted            mean Treated 8.432432e-01
#> 9   race_black  Binary Unadjusted variance.ratios    <NA>           NA
#> 10 race_hispan  Binary Unadjusted            mean Control 1.421911e-01
#> 11 race_hispan  Binary Unadjusted            mean Treated 5.945946e-02
#> 12 race_hispan  Binary Unadjusted variance.ratios    <NA>           NA
#> 13  race_white  Binary Unadjusted            mean Control 6.550117e-01
#> 14  race_white  Binary Unadjusted            mean Treated 9.729730e-02
#> 15  race_white  Binary Unadjusted variance.ratios    <NA>           NA
#> 16     married  Binary Unadjusted            mean Control 5.128205e-01
#> 17     married  Binary Unadjusted            mean Treated 1.891892e-01
#> 18     married  Binary Unadjusted variance.ratios    <NA>           NA
#> 19        re74 Contin. Unadjusted            mean Control 5.619237e+03
#> 20        re74 Contin. Unadjusted            mean Treated 2.095574e+03
#> 21        re74 Contin. Unadjusted variance.ratios    <NA> 5.181285e-01
```

They also take `var.names`, which gives the covariates the names you
want to report them under. Give it to
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
and it applies wherever the object is displayed afterwards, including
[`love.plot()`](https://ngreifer.github.io/cobalt/reference/love.plot.md);
give it here and it adds to that:

``` r

format(b, var.names = c(age = "Age", race = "Race/Ethnicity",
                        re74 = "Earnings, 1974"))
#>                          Type    M.0.Un    M.1.Un Diff.Un     M.Threshold.Un
#> Age                   Contin.   28.0303   25.8162 -0.3094 Not Balanced, >0.1
#> educ                  Contin.   10.2354   10.3459  0.0550     Balanced, <0.1
#> Race/Ethnicity_black   Binary    0.2028    0.8432  0.6404 Not Balanced, >0.1
#> Race/Ethnicity_hispan  Binary    0.1422    0.0595 -0.0827     Balanced, <0.1
#> Race/Ethnicity_white   Binary    0.6550    0.0973 -0.5577 Not Balanced, >0.1
#> married                Binary    0.5128    0.1892 -0.3236 Not Balanced, >0.1
#> Earnings, 1974        Contin. 5619.2365 2095.5737 -0.7211 Not Balanced, >0.1
#>                       V.Ratio.Un
#> Age                       0.4400
#> educ                      0.4959
#> Race/Ethnicity_black           .
#> Race/Ethnicity_hispan          .
#> Race/Ethnicity_white           .
#> married                        .
#> Earnings, 1974            0.5181
```

A name given for a base variable reaches every name it appears in, which
is why one entry for `race` renames all three of its levels. Only the
displayed names change—`b$Balance` still has the stored ones—so nothing
downstream has to know about them. See `?display-options` for the other
ways to specify it.

The underlying components are still there if you want them: the balance
table lives in `b$Balance` for an unsegmented object, and in a
differently named component otherwise (`b$Cluster.Balance` and so on).
`str(b, give.attr = FALSE)` shows the structure.

### How are balance statistics computed when using subclassification?

Subclassification involves creating strata (usually based on the
propensity score), within which covariates are ideally balanced.
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
lets you assess balance both within and across subclasses.

One must always remember that the standardized mean difference uses the
standardization factor computed in the original sample, i.e., prior to
subclassification. Let’s take a look below using *MatchIt*:

``` r

# PS Subclassification
msub <- MatchIt::matchit(treat ~ age + educ + race + married + re74,
                         data = lalonde, method = "subclass",
                         estimand = "ATE", min.n = 4)

# Balance in the first subclass
bal.tab(msub, which.sub = 1, binary = "std")
#> Balance by subclass
#> 
#> ─── Subclass 1 ──────────────
#> 
#>                 Type Diff.Adj
#> distance    Distance   0.1574
#> age          Contin.  -1.0433
#> educ         Contin.  -0.2759
#> race_black    Binary   0.0000
#> race_hispan   Binary   0.0000
#> race_white    Binary   0.0000
#> married       Binary  -1.1135
#> re74         Contin.  -1.8353
```

Let’s see where the number `-1.0433` came from (the standardized mean
difference for `age`). We compute the mean of `age` in each treatment
group in subclass 1, and then divide it by the pooled standard deviation
of age (because we requested the ATE) in *the original sample*.

``` r

m0 <- mean(lalonde$age[lalonde$treat == 0 & msub$subclass == 1])
m1 <- mean(lalonde$age[lalonde$treat == 1 & msub$subclass == 1])

s0 <- sd(lalonde$age[lalonde$treat == 0])
s1 <- sd(lalonde$age[lalonde$treat == 1])

(m1 - m0) / sqrt((s1^2 + s0^2) / 2)
#> [1] -1.043294
```

A common mistake is to compute the standard deviation within each
subclass. There are a few reasons why this is bad: 1) it suffers from
the same problem that changing the standardization factor does with
matching or weighting, i.e., that balance can appear to be worse because
the standardization factor shrank even as the means got closer together;
2) when there is no or little variation of a covariate within a
subclass, which is desirable, the standardization factor will be tiny,
making the SMD potentially appear huge; and 3) the same variable will
use different standardization factors across subclasses, which means the
same difference in means, which contribute to bias equally, will have
different balance statistics.

A related question is how the balance statistics are computed across
subclasses to compute an overall balance statistic for the sample. For
(standardized) mean differences, it is as easy as computing the average
of the statistic across subclasses, where the statistics are weighted
corresponding to the number of units in the subclass in the target group
(e.g., the treated units for the ATT, all units for the ATE, etc.).
Below I’ll demonstrate how to do that manually for the `age` covariate:

``` r

# SMDs across subclasses for age
smds <- sapply(1:6, function(s) {
    m0 <- mean(lalonde$age[lalonde$treat == 0 & msub$subclass == s])
    m1 <- mean(lalonde$age[lalonde$treat == 1 & msub$subclass == s])
    
    s0 <- sd(lalonde$age[lalonde$treat == 0])
    s1 <- sd(lalonde$age[lalonde$treat == 1])
    
    (m1 - m0) / sqrt((s1^2 + s0^2) / 2)
})

# Sample size in each subclass
ns <- table(msub$subclass)

# Summary SMD for age
weighted.mean(smds, ns)
#> [1] -0.2354095

bal.tab(msub)
#> 
#> Balance measures across subclasses
#>                 Type Diff.Adj
#> distance    Distance   0.1081
#> age          Contin.  -0.2354
#> educ         Contin.   0.0075
#> race_black    Binary   0.0535
#> race_hispan   Binary  -0.0420
#> race_white    Binary  -0.0115
#> married       Binary  -0.1160
#> re74         Contin.  -0.3200
#> Sample sizes by subclass
#>           1   2  3   4   5   6 All
#> Control 102 100 88  72  39  28 429
#> Treated   4   4  9  30  62  76 185
#> Total   106 104 97 102 101 104 614
```

This works for mean differences but not other statistics. So the way
*cobalt* actually does this is compute stratification weights, and then
compute the balance statistics using the stratification weights in the
full sample. Stratification weights are first computed by computing the
proportion of treated units in each sample, and then using the formulas
to compute propensity score weights from propensity scores. Here’s how I
do that manually for `age`:

``` r

# Compute proportion of treated units in each subclass
prop1 <- sapply(1:6, function(s) mean(lalonde$treat[msub$subclass == s]))

# Assign to each unit
ps <- prop1[msub$subclass]

# Compute ATE weights
w <- ifelse(lalonde$treat == 1, 1 / ps, 1 / (1 - ps))

# Compute weighted KS statistic
col_w_ks(lalonde$age, treat = lalonde$treat,
         weights = w)
#> [1] 0.1658923

bal.tab(msub, stats = "ks")
#> 
#> Balance measures across subclasses
#>                 Type KS.Adj
#> distance    Distance 0.2187
#> age          Contin. 0.1659
#> educ         Contin. 0.0627
#> race_black    Binary 0.0535
#> race_hispan   Binary 0.0420
#> race_white    Binary 0.0115
#> married       Binary 0.1160
#> re74         Contin. 0.3038
#> Sample sizes by subclass
#>           1   2  3   4   5   6 All
#> Control 102 100 88  72  39  28 429
#> Treated   4   4  9  30  62  76 185
#> Total   106 104 97 102 101 104 614
```

### Why don’t I get the same balance statistics when using *cobalt* as I do when using *tableone*?

*tableone* is another package that provides tools for balance
assessment. One strength that the package has is its beautiful,
publication-ready tables that include summary statistics for the
covariates, clean variable names, and clean headings. But it does not
incorporate best practices in balance assessment in favor of
transparency. This differs from the ethos of *cobalt*, which is to
provide highly customizable balance statistics that reflect best
practices and use well-reasoned decisions. This is not an insult to
*tableone* but is meant to reflect the different purposes *cobalt* and
*tableone* have. They should not be used interchangeably or expect to
yield identical results because they use different formulas for
computing certain statistics, most notably the SMD.

Below are some of the reasons why SMDs might differ between *tableone*
and *cobalt*:

- *tableone* always uses the pooled standard deviation (i.e., the
  standardizaton factor setting `s.d.denom = "pooled"`) as the
  standardization factor, while *cobalt* determines the standardization
  factor based on the estimand (though by default or when the ATE is the
  estimand, the two should be aligned).
- *tableone* uses the weighted standardization factor in the SMD,
  whereas *cobalt* always uses the standardization factor computed in
  the unadjusted sample. For matching, this means *tableone* computes
  the standardization factor in the matched sample, while *cobalt* uses
  the original sample.
- *tableone* uses
  [`survey::svyvar()`](https://rdrr.io/pkg/survey/man/surveysummary.html)
  to compute weighted variances, whereas *cobalt* uses the formula
  described previously (and implemented in
  [`col_w_sd()`](https://ngreifer.github.io/cobalt/reference/balance-summary.md)).
  These values will differ by small amounts when the weights are not
  constant.
- For multi-category covariates, *tableone* uses a single statistic
  described by Yang and Dalton
  ([2012](#ref-yangUnifiedApproachMeasuring2012)) to summarize balance,
  whereas *cobalt* provides a balance statistic for each level of the
  covariate. There is no reason to prefer the statistic used by
  *tableone*; it does not have any relationship to the bias of the
  estimate and can mask large differences in some categories when there
  are many categories. See
  [here](https://stats.stackexchange.com/a/496608/116195) for a more
  detailed answer.

In practice, these differences will be small. Obviously, I recommend
using *cobalt* instead for balance assessment, and I recommend reporting
the balance statistics *cobalt* produces. That said, if you understand
what *tableone* is doing and are okay with the choices it makes, there
is no denying that it can produce beautiful tables.

### Why doesn’t `thresholds` work with `bal.tab()` with multiply imputed or clustered data?

This question was asked
[here](https://github.com/ngreifer/cobalt/issues/90) and
[here](https://stackoverflow.com/q/79562364/6348551). With multiply
imputed data, the default output of
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md) is
the balance summary across imputations, which contains, for each balance
statistic and for each covariate, the minimum, mean, and maximum value
of that balance statistic for that covariate across imputations. When
you request a balance threshold using `thresholds`, it isn’t clear to
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
which of those summaries the threshold is to be applied to. To get
thresholds to appear, supply an argument to `imp.fun` to request just
one summary, e.g., `imp.fun = "mean"`, and the thresholds will be
applied to that summary.

For clustered data, the same is true, but the across-cluster balance
summary is not displayed by default. To request a single summary, use
`cluster.fun`.

Note this does not apply to
[`love.plot()`](https://ngreifer.github.io/cobalt/reference/love.plot.md),
which will produce thresholds even when the default `agg.fun` (`"range"`
for multiply imputed data) is requested.

## References

Stuart, Elizabeth A. 2008. “Developing Practical Recommendations for the
Use of Propensity Scores: Discussion of ‘A Critical Appraisal of
Propensity Score Matching in the Medical Literature Between 1996 and
2003’ by Peter Austin, Statistics in Medicine.” *Statistics in Medicine*
27 (12): 2062–65. <https://doi.org/10.1002/sim.3207>.

Yang, Dongsheng, and Jarrod E Dalton. 2012. *A Unified Approach to
Measuring the Effect Size Between Two Groups Using SAS®*. 6.

[^1]: For multi-category treatments, all standardization factors are
    computed using the full data, not just the groups being compared.
    For example, the pooled standard deviation involves computing the
    mean of all the group-specific variances, not just the two being
    compared. Similarly, in the `"hedges"` formula, \\n-2\\ is replaced
    with \\n-k\\, where \\k\\ is the number of treatment groups.

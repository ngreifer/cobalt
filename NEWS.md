`cobalt` News and Updates
======

# cobalt 4.5.5

* Minor updates to `bal.plot()` to prevent warnings due to `ggplot2` 3.5.0.

* Improved processing when no covariates are specified.

* Fixed a bug when multiple weights are specified, `s.d.denom` is specified, and either all variables are continuous and `continuous = "raw"`, all variables are binary and `binary = "raw"`, or both `continuous = "raw"` and `binary = "raw"`.

* Documentation updates

# cobalt 4.5.4

* Minor update to accommodate `ggplot2` 3.5.0. Thanks to @teunbrand. (#80)

* `bal.tab()` no longer throws a note about `s.d.denom` when `binary = "raw"` and `continuous = "raw"` (i.e., only raw mean differences are requested).

* `col_w_smd()` now correctly includes all data in computing the standardization factor for standardized mean differences when `subset` is supplied, consistent with the documentation. Previously, only the subsetted units were included in the standardization factor. This does not affect any results from `bal.tab()` or `bal.compute()`, which already used the correct units.

* Added a new vignette for frequently asked questions, which describes in further detail why some choices were made. See `vignette("faq")`.

# cobalt 4.5.3

* Fixed a bug when missing values were present in continuous covariates. Thanks to @vnusinfo. (#76)

* Fixed a bug when using `bal.tab()` with the `cluster` argument supplied with the `caret` package loaded. Thanks to @BorgeJorge. (#77)

* When `cluster` is specified, categorical variables that perfectly coincide with the cluster variable are now correctly removed.

* Perfectly colinear variables are no longer removed (unless they are binary variables split from the same factor). This should speed up evaluation and reduce the probability of false positives being removed.

* Variables with a single value are now more reliably categorized as "binary" in tables and calculations.

# cobalt 4.5.2

* Fixed a bug when using `bal.compute()` with a treatment variable with levels named "treated" and "control".

* Fixed a bug when using `bal.tab()` with `mnps` objects from `twang`. Thanks to @sherwinkuah. (#74)

* Fixed a bug when using `addl` without a dataset supplied. (#71)

* Fixed a bug when using `subset` to remove clusters lacking full representation in all treatment groups when `cluster` is specified. (#70)

# cobalt 4.5.1

* Added a new function `available.stats()` which lists the available balance statistics for use with `bal.init()` and `bal.compute()`.

* The interfaces to `bal.compute()` and `bal.init()` have changed slightly. The arguments have a slightly different order to match other `cobalt` functions. `bal.compute()` now is a generic function with a method for `bal.init` objects and a default method. The default method accepts the same arguments as `bal.init()` (and optionally an additional `weights` argument) and computes the (weighted) balance statistic directly. For most uses, the `bal.init() |> bal.compute()` workflow should be preferred.

* Added a new multivariate balance statistic, the kernel distance as described by Zhu, Savage, and Ghosh (2018). This can be requested in `bal.compute()` and `bal.init()` by setting `stat = "kernel.dist"`. In most cases, this will perform similarly to the energy distance.

* `s.weights` are more compatible with matching methods and can now be supplied to `bal.tab()` with `mimids` and `wimids` objects from `MatchThem`. Thanks to Helen Wright for pointing out this issue.

* Fixed a bug when using `bal.init()` with non-`NULL` `s.weights`.

* Fixed bugs when using `bal.init()` and `bal.compute()` with multi-category treatments.

* Fixed a bug when using `col_w_ovl()` with missing data.

* Fixed a bug when using `bal.tab()` with `stat = "spearman.correlations"`.

* Fixed a bug when using `bal.plot()` with longitudinal treatments.

* Fixed a bug in while the display options `factor_sep` and `int_sep` were not functioning correctly.

* Fixed a bug when no covariates are supplied.

* Improved some errors all around, and particularly in `col_w_smd()` and friends, `bal.init()`, and `var.names()`.

# cobalt 4.5.0

* Added new functions `bal.compute()` and `bal.init()`, which are used for compute scalar balance statistics efficiently for use in optimizing balance. A new vignette, `vignette("optimizing-balance")` is available as well.

* When `focal` is specified with multi-category treatments (by the user or implicitly by the supplied object), `pairwise` can be set to `TRUE` to request balance between each pair of treatment groups and to `FALSE` to request balance only between each non-focal group and the focal group. Previously only the behavior of setting `pairwise` to `FALSE` was supported. Now the default is for `pairwise` to be `TRUE`. To recover balance results for version prior to this one, set `pairwise = FALSE` with non-`NULL` `focal`.

* With `optmatch` objects, the `estimand` argument can now be supplied to `bal.tab()`, etc., to control how the matching weights are computed from the subclass/pair membership. This is consistent with how `get.w()` uses the same argument.

* Fixed a bug in which using `.` in formulas incorrectly included the treatment among the covariates.

* Fixed a bug in which formulas supplied as character strings were not correctly interpreted as formulas. Thanks to @istallworthy.

* Fixed a bug in which a spurious warning about dropping weights would occur when using `bal.plot()` with a density.

* Documentation updates, including some new pages and the use of `roxygen2`.

# cobalt 4.4.1

* Fixed a bug when covariates with nonstandard names are extracted from model objects (#63). Thanks to @markdanese.

* Fixed a bug when "0" and "1" are the names of two of the treatment levels in a multinomial treatment.

* Fixed a bug with the default method of `bal.tab()` which was ignoring components of the supplied object.

* Fixed a bug where `bal.plot()` would ignore `s.weights`. They are now included correctly.

* The call to the original balancing function is now hidden by default. To request it be displayed, set `disp.call = TRUE` in the call to `bal.tab()` or `print.bal.tab()` or use `set.cobalt.options(disp.call = TRUE)` to display it for the session.

# cobalt 4.4.0

* Added support in `bal.plot()` for negative weights with `type = "density"`.

* Added support for `ps.cont()` objects from the `twangContinuous` package. `ps.cont` objects from `WeightIt` are no longer supported.

* Major documentation overhaul. More arguments are explained at `help("bal.tab")` and a new package help page can be found at `help("cobalt-package")`.

* The function call is no longer included in the `bal.tab()` results for objects from `twang`.

* Fixed a bug when some predictors were binary in some clusters and continuous in others. Variables now have a stable type across partitions.

* Fixed a bug where binary variables were not being correctly processed when using the `formula` interface.

* When using `poly`, orthogonal polynomials can be requested by setting `orth = TRUE`.

* Improved appearance of conditional examples in `pkgdown` site.

* Removed `mlogit` from Suggests.

* Returned `sbw` to Suggests.

* Updated the logo, thanks to [Ben Stillerman](https://stillben.com).

# cobalt 4.3.2

* When `pairwise = FALSE` with binary or multi-category treatments, the balance statistics now refer to the difference between each group and the original full sample, unadjusted except possibly by `s.weights`. Previously, they referred to the difference between each group and the combined adjusted sample.

* When subclassification is used and some units are discarded, `bal.tab()` now reports the number of discarded units along with with the number of units in each subclass in the sample sizes table.(#59)

* Fixed several bugs when using `love.plot()` with subclassification that were caused by the last update. Thanks to Mario Lawes for pointing them out.

* Fixed a bug in how `get.w()` computed weights for `Match` objects resulting from `Matching::Match()` with `estimand = "ATE"`. Results now agree with `Matching::MatchBalance()`.

* Fixed a bug that would occur when using `cobalt` functions without attaching the package (e.g., `cobalt::bal.tab()`). (#53)

* Fixed a bug that would occur with ordinal treatments.

* Added better support for negative weights.

* Fixed typos (#54, many identified and fixed by @jessecambon).

# cobalt 4.3.1

* Added support for objects from the new version of `MatchThem`.

* Fixed a bug and improved speed when using `match.strata`.

# cobalt 4.3.0

* Returned `cem` to Suggests.

* Added ability to display threshold summaries with multiply imputed datasets, clustered datasets, multi-category treatments, and longitudinal treatments.

* Added `pairwise` argument for binary treatments. When set to `FALSE`, `bal.tab()` will display balance between each treatment group and the full sample (i.e., the target population). This functionality already existed for multi-category treatments; indeed, for binary treatments, it works by treating the treatment as multi-category.

* Added two new `stats` options in `bal.tab()` and `love.plot()` for continuous treatments: `"mean.diffs.target"` (abbreviated as `"m"`) and `"ks.statistics.target"` (abbreviated as `"ks"`). These compute (standardized) mean differences and KS statistics between the weighted and unweighted samples to ensure the weighted sample is representative of the original population. These statistics are only computed for the adjusted sample (i.e., they will not appear in the absence of adjustment).

* With subclassification methods, the arguments `which.subclass` and `subclass.summary` have been added to display balance on individual subclasses and control output of the balance across subclasses summary. These arguments replace the `disp.subclass` argument, which can still be used.

* When using `bal.plot()` with clustered or multiply imputed data, the `which.cluster` and `which.imp` arguments can be set to `.none` to display balance ignoring cluster membership and combining across imputations.

* Changed processing of the `print()` method. Now there is only one `print()` method (`print.bal.tab()`) for all `bal.tab` objects. Processing is a little smoother and some printing bugs have been fixed. `"as.is"` can no longer be supplied to keep the print setting as-is; simply omit the corresponding argument to use the options as specified in the call to `bal.tab()`.

* An additional argument, `disp.call`, can be supplied to `bal.tab()` and `print.bal.tab()` to control printing of the `call` component of the input object, which contains the original function call. Set to `FALSE` to hide the call. This option is documented in `?display_options` and can also be set using `set.cobalt.options()`.

* The balance table component of `bal.tab` objects is smaller because some extraneous columns are no longer produced. In particular, if no threshold is requested, no threshold columns will be produced. This does not affect display, but makes it easier to extract balance statistics from `bal.tab` objects (e.g., for exporting as a table). This does mean that previously saved `bal.tab` objects produced by earlier versions of `cobalt` will not be able to be printed correctly.

* Fixed a bug where `bal.plot()` would incorrectly process 2-level factor variables (#48).

* Fixed a bug where `love.plot()` would not display variables in the correct order when using aggregation and setting `var.order = NULL`. Thanks to Florian Kaiser.

* Fixed a bug in `love.plot()` where the color of points could be incorrect.

* Fixed a bug in `love.plot()` where samples were not always displayed in the right order. Now they are displayed in the same order they are in `bal.tab()`.

* Fixed a bug in `love.plot()` when the weight names had spaces in them.

* Added an error message when not all clusters contain all treatment levels. Thanks to Rachel Visontay.

* Fixed a bug when supplying the `weights` argument as a list of supported objects (e.g., `weightit` objects) if they were unnamed. Samples are more conveniently named.

* Fixed a bug in `col_w_mean()`, `col_w_smd()`, and friends that occurred when few nonzero weights were present. Now an informative error is thrown.

* Updates to documentation.

# cobalt 4.2.4

* Sampling weights now function correctly with subclassification.

* Fixed a bug in `print.bal.tab()` when no units were unmatched but some were discarded.

* Fixed an issue with the version number for `gridExtra` in `DESCRIPTION`. (#47)

* `cem` removed from Suggests because it has been removed from CRAN.

* Updated to support `MatchIt` 4.0.0, which includes sampling weights and improved processing of the covariates.

# cobalt 4.2.3

* Fixed bugs in processing functions in formulas, including `rms` functions and `poly()`. (#40)

* Fixed a bug in how KS statistics were computed with `col_w_ks()`. Results now agree with those from `MatchIt` and `twang`.

* Fixed bugs in processing small and partially empty subclasses. 

* In functions that compute weights from matching strata (e.g., `get.w()` for some types of objects), an `estimand` argument can be supplied to choose which formula is used to compute the weights. Subclass propensity scores are computed as the number of treated units in each subclass, and then stabilized weights are computed from those propensity scores using the standard formulas.

* Effective sample sizes now print only up to two digits (believe me, you don't need three) and print more cleanly with whole numbers.

# cobalt 4.2.2

* Fixed a bug due to new version of `sbw`.

* Minor improvements to error messages and documentation.

# cobalt 4.2.1

* Fixed a bug where `int` and `poly` were ignored with binary and continuous treatments.

* Fixed a bug where subclass balance statistics were incorrectly computed. Thanks to Mario Lawes.

* Improved processing of inappropriately given S4 objects.

* Removed `bal.tab` methods for atomic vectors (which were undocumented). The errors they would provide when inappropriately supplied were unhelpful.

* Fixed a bug with `backports` 1.1.7 not running correctly.

* Fixed a bug with `str2expression` for R versions below 3.6.0. Thanks to @kthohr and @jimmyg909.

* When data is segmented (i.e., with a multi-category or longitudinal treatment or when clusters or multiple imputations are specified), the balance summary across segments will not be computed or displayed when individual segment balance is requested. See `?display_options` to see the defaults for the different segment types, some of which have changed. 

* Updated some warnings.

# cobalt 4.2.0

* Added support for `Matchby` objects resulting from a call to `Matchby()` in the `Matching` package. These function identically to `Match` objects.

* When using `formula` inputs, interaction terms (e.g., `X1 * X2`) will now correctly be resolved and displayed as an interaction term. This makes it easier to check balance on specific interactions rather than having to set `int = TRUE` or create a separate interaction variable in the data. Interaction terms specified in this way will be ignored by `int` and `poly` when they are used. When using `var.names` with `love.plot()`, changing the names of the base components of the interaction will also change their name in the interaction term, consistent with `int` behavior. If a `formula` was supplied in the input object to `bal.tab()` or other functions, terms in that formula will now also be included in balance reports.

* Arguments to `addl` can now be specified as a one-sided formula (e.g., `~ X1 + X2 * X3`). This makes it easy to take advantage of the above changes to the formula interface to add additional interaction terms. The formula will look at all available datasets in the conditioning object or supplied to `bal.tab()` and at the global environment. If supplying a single variable that exists in the global environment, it makes sense to supply it as a formula (e.g., `addl = ~ X1`) rather than as just the variable (e.g., `addl = X1`). Doing the former will retain the name of the variable. The same can be done with `distance`. If variables in `addl` are perfectly correlated with or have the same name as supplied covariates, those variables will be removed from `addl`.

* If only one argument is provided to `f.build()` (e.g., `f.build("x")`), it will be treated as the right-hand-side of the formula with no left-hand-side (e.g., the above will evaluate to `~ x`).

* Fixed bug that caused `match.strata` input to be ignored.

* Improved processing and error reporting when using the default `bal.tab()` method.

* Speed improvements due to changes in how formulas are processed (now using `model.matrix()` directly rather than `splitfactor()` to process factors) and other small fixes. This is what enables the above changes to the formula capabilities.

* Improved documentation for `weightitMSM` objects from `WeightIt` and `CBMSM` objects from `CBPS`.

* Fixed bug when printing `bal.tab` objects with continuous treatments.

* Fixed bug when using multi-category treatments with numbers as the level names.

* Fixed bug when using `mnps` objects from `twang` with multiple stop methods.

* Fixed bug when requesting means or standard deviations with segmented data.

* Fixed bug where the x-axis in `love.plot` was always "Standardized Mean Differences" even when it wasn't supposed to be for `stats = "mean.diffs"`.

* `rlang` is now in IMPORTS. 

* General speed and stability improvements.

# cobalt 4.1.0

* Added support for `sbwcau` objects from `sbw`. See Appendix 1 or `?bal.tab.sbw` for an example.

* Added support for `cem.match` objects from `cem`. See Appendix 1 or `?bal.tab.cem.match` for an example.

* Added `stats` argument to `bal.tab()` and `print()` to replace `disp.v.ratio` and `disp.ks`. This argument functions similarly to how it does in `love.plot()`; for example, to request mean differences and variance ratios, one can enter `stats = c("m", "v")`. One consequence of this is that it is possible to request statistics that don't include mean differences. See `?display_options` for more details. The old arguments still work (and probably always will) but you should use `stats` instead. The goal here was to unify syntax across `bal.tab()`, `print()`, and `love.plot()`. A new help page specifically for the `stats` argument can be viewed at `?balance.stats`.

* Added `thresholds` argument to `bal.tab()` to replace `m.threshold`, `v.threshold`, etc. This argument functions similarly to how it does in `love.plot()`; for example, to request thresholds for mean differences and variance ratios, one can enter `thresholds = c(m = .1, v = 2)`. The old arguments still work (and probably always will) but you should use `thresholds` instead. The goal here was to unify syntax across `bal.tab()`, `print()`, and `love.plot()`.

* Added `disp.means` option to `bal.plot` to display the mean of the covariate as a line on density plots and histograms.

* Added `"hedges"` as an option to `s.d.denom`. This will compute the standardized mean difference using the formula for the small sample-corrected Hedge's G as described in the What Works Clearinghouse Procedures Handbook.

* With multi-category treatments when `pairwise = FALSE`, rather than computing balance between each treatment group and the other treatment groups, balance is now computed between each treatment group and the entire sample.

* In `print()`, the arguments `disp.m.threshold`, `disp.v.threshold`, `disp.ks.threshold`, and `disp.r.threshold`, which could be set to `FALSE` to prevent the corresponding balance thresholds and summaries from being printed, have been replaced with `disp.thresholds`. Named entries can be set to `FALSE`. The goal here was to unify syntax across `bal.tab()` and `print()`.

* A new balance statistic, the overlapping coefficient (OVL), is allowed with binary and multi-category treatments. This is described in Belitser et al. (2011) and Franklin et al. (2014) for assessing balance. Generally, for each covariate, the overlapping coefficient is the area of the probability density functions for each sample that overlap. Here I follow Franklin et al. (2014) and report 1 - (OVL) so that values close to zero indicate good balance (i.e., completely overlapping distributions) and values close to 1 indicate poor balance (i.e., completely non-overlapping distributions). To estimate and display the OVL, set include `"ovl"` in the `stats` argument in a call to `bal.tab()` or `love.plot()` (or you can use the old syntax by setting `disp.ovl = TRUE`). The balance threshold can be requested by including `"ovl"` in the `thresholds` argument (or you can use the old syntax by using the `ovl.threshold` argument).

* Spearman correlations can be requested for continuous treatments by adding `"sp"` to the `stats` argument.

* The argument `weights` can now be supplied to any `bal.tab()` call to request balance on additional weights beyond the weights from the object on which `bal.tab()` is called. This argument takes a named list, where each element is a vector of weights, the name of a variable containing weights in an available dataset, or an object with a `get.w()` method (e.g., the output of another preprocessing function). This should make it easier to compare balancing methods without having to specify the covariates and treatment using the `formula` or `data.frame` methods.

* `ggplot2` version 3.3.0 is required, which removes some warnings and makes it so `ggstance` doesn't need to be imported.

* When there are more than 900 variables to compute balance statistics on in `bal.tab` (which can happen quickly when `int = TRUE` and categorical variables have many categories), to avoid major slowdowns, checks for redundancy of variables are forgone. This will dramatically increase the speed of `bal.tab` in these scenarios. This option can be changed with the `cobalt` option `"remove_perfect_col"` which can be set to `TRUE` or or `FALSE`. Set to `FALSE` to improve speed at the expense of possibly having redundant variables appear.

* Fixed a bug when using the default `bal.tab` method with objects containing longitudinal treatments.

* Fixed a bug when using `bal.tab` with continuous treatments and clusters.

* Fixed a bug in `love.plot()` when using subclassification.

* Fixed a bug when using `bal.tab` with longitudinal treatments and multiple sets of weights.

* Fixed a bug when using `col_w_ovl()`. OVL values are now more accurate.

* Speedups and other small fixes.

# cobalt 4.0.0

**Major Updates**

* Added support for `mimids` and `wimids` objects from `MatchThem`.

* Major restructuring so that clusters, longitudinal treatments, multi-category treatments, and multiply imputed data can all be used with each other. These are layers in the following order: clusters, time points, treatment categories, and imputations. Summaries across these layers are handled slightly differently from how they used to be; importantly, summaries are not nested, and only the lowest layer present can have a summary. For example, if multiply imputed data is used with multi-category treatments, there will be a summary across imputations (the lowest layer) but not across treatment pairs. `love.plot` allows multiple forms of faceting and aggregating and is extremely flexible in this regard.

* Major changes to appearance of `bal.plot()` to be more in line with `love.plot()`, including new `grid` and `position` options to control the presence of the grid and the position of the legend.

* Formula interfaces now accept `poly(x, .)` and other matrix-generating functions of variables, including the `rms`-class-generating functions from the `rms` package (e.g., `pol()`, `rcs()`, etc.) (the `rms` package must be loaded to use these latter ones) and the `basis`-class-generating functions from the `splines` package (i.e., `bs()` and `ns()`). A bug in an early version of this was found by @ahinton-mmc.

**Minor Updates and Bug Fixes**

***`s.d.denom` and `estimand`***

* `s.d.denom` can now use the name of a treatment rather than just `"treated"` or `"control"`. In addition, `s.d.denom` can be `"weighted"` to use the weighted sample's standardization factors, an option available for continuous treatments, too.

* Improved guessing of the estimand when not provided.

* Estimands besides ATT can now be used with subclasses. The estimand can be inferred from the provided subclasses. Works with `match.strata` as well, which function like subclasses. In addition, it is not always assumed that `MatchIt` objects are targeting the ATT, for example, with subclassification or calipers.

***`bal.plot`***

* Added `sample.names` argument in `bal.plot` in response to [this post](https://stackoverflow.com/questions/57970679/change-name-of-groups-in-bal-plot) on Cross Validated.

* Added functionality to the `which` argument in `bal.plot`, allowing more specificity when multiple sets of weights are used.

* Added `type = "ecdf"` option to `bal.plot` for categorical treatments with continuous covariates to display empirical cumulative density plots as an alternative to density plots.

* When using `bal.plot` with continuous treatments and continuous covariates, the points are shaded based on their weights; this behavior is controlled by the new `alpha.weight` argument, which replaces the functionality of `size.weight` (which was kind of ugly and not very informative) and is `TRUE` by default. Now it's more apparent which points are influential in the weighted sample. In addition, a line illustrating the unweighted covariate mean is present.

* The default of the `grid` argument is now `FALSE` in `bal.plot()` and `love.plot()`. Previously it was `TRUE`. This make the plots cleaner at the outset.

***Other improvements***

* Added new function `col_w_cov()` to compute treatment-covariate covariances (i.e., unstandardized correlations) for continuous treatments. `continuous` and `binary` can be set to `"raw"` in `bal.tab()` and `std` can be set to `FALSE` in `col_w_cov()` to request treatment-covariate covariances instead of correlations. `col_w_corr()` is now a wrapper for `col_w_cov()` with `std = TRUE`. To get more functionality out of the `std` argument (e.g., to standardize the covariances for some covariates but not others), use `col_w_cov()`. 

* Balance summary functions (e.g., `col_w_sd()`, `col_w_smd()`, etc.) process binary variables slightly differently. If `bin.vars` is missing, the function will figure out which variables are binary. If `NULL`, it will be assumed no variables are binary. Entering values for `bin.vars` can be done more flexibly. When a factor variable is supplied as part of `mat` and is split internally by `splitfactor()`, extra values will be automatically added to `bin.vars` with the newly created dummies considered binary variables.

* Bug fixes when binary factor treatments are used, thanks to Moaath Mustafa Ali. 

* `bal.tab()` no longer tells you whether it assumes matching or weighting when certain non-package-related methods are used.

* Improvements to assessment of subclass balance. For binary treatments, balance statistics other than mean differences can now be requested. The across-subclass balance summary uses subclassification weights (processed in the same way `match.strata` is) instead of simply taking a weighted average across subclasses (which is not valid for non-additive statistics like variance ratios or KS statistics). For continuous treatments, a balance summary across subclasses can now be produced. This uses a weighted average of the subclass-specific balance statistics.

* The default in `love.plot()` for `abs` is now to be whatever it is in the (implicit) call to `bal.tab()`, which is usually `FALSE`. Previously `abs` was not aligned between `love.plot()` and `bal.tab()`.

* `s.weights` can now be manually supplied to methods that usually come with their own sampling weights, such as `twang` and `WeightIt`.

* Speedup of `splitfactor()`.

* `splitfactor()` now has a `split.with` option to split one or more vectors in concert with the data set being split.

* `splitfactor()` and `unsplitfactor()` are a little smarter and more in sync.

* All functions work better inside other functions like `lapply()` or `purrr::map()`, thanks to @the-Zian.

* Updates to the vignettes; Appendix 2 is particularly different.

* Other bug fixes and performance improvements here and there.

# cobalt 3.9.0

* Added vignette for use of `love.plot`.

* Changed `grid` version requirement.

* Updated README.

* Fixed bugs that would occur when using `love.plot()` with various combinations of `var.order`, multiple `stats`, and `agg.fun = "range"`.

* Fixed bugs that would occur when using `bal.tab()` with objects from the `Matching` package. Calculated statistics are now the same as those generated using `Matching::MatchBalance`. Changes based on updates to `get.w.Match()`.

* Added balance summary functions `col_w_mean`, `col_w_sd`, `col_w_smd`, `col_w_vr`, `col_w_ks`, `col_w_ovl`, and `col_w_corr`. These make it easier to get quick, simple summaries of balance without calling `bal.tab`, for example, for use in programming other functions. Some of these are now used inside `bal.tab` to increase speed and simplify internal syntax.

* Other small bug fixes.

# cobalt 3.8.0

* Added the ability to display balance on multiple measures (e.g., mean differences, variance ratios, KS statistics) at the same time with `love.plot()`.

* Bug fixes that make `bal.tab()` and `love.plot()` more usable within other functions and especially when called with `do.call()`.

* Made it easier to get proper `bal.tab` output when using `matchit()` with an argument to `distance` (in the call to `matchit()`). Include the original dataset in the `data` argument of `bal.tab()` to get the variables to display correctly.

* Changed the default shape in `love.plot()` to `"circle"`, which is a solid circle. I found this a prettier alternative to the open circle, especially on Windows. To get back open circles you set `shapes = "circle filled"` (yes, that is a bit confusing).

* Added ability to hide the gridlines easily in `love.plot()`.

* Changed the calculation of standard deviations (and standardized differences in proportion) for binary variables to be more in line with recommendations, as noted by @mbloechl05. Note this will make these values different from those in `MatchIt::summary` by a small amount.

* The KS statistic is now computed for binary variables. It is simply the difference in proportion.

* Allowed some methods to accept `mids` objects (the output of a call to `mice::mice()`) in the `data` argument to supply multiply imputed data. This essentially replaces `data = complete(imp.out, "long"), imp = ".imp"` with `data = imp.put`, assuming `imp.out` is a `mids` object.

* Other bug fixes and improvements.

# cobalt 3.7.0

* Changes to some `bal.tab` defaults: `quick` is now set to `TRUE` by default. Adjusted and unadjusted means, standard deviations, and mean differences will always be computed, regardless of `quick`. Variance ratios and KS statistics will only be computed if `quick = FALSE` or `disp.v.ratio` or `disp.ks`, respectively, are `TRUE`.

* Variance ratios now respond to `abs`. When `abs = FALSE`, the default in `bal.tab`, the variance is ratio is the variance of the treated (1) divided by the variance of the control (0). When `abs = TRUE`, the numerator of the variance ratio is the larger variance and the denominator is the smaller variance, which was the old behavior. `v.threshold` still responds as if `abs` was set to `TRUE`, just like with mean differences. Any time variance ratios are aggregated (e.g., across imputations or clusters), the "mean" variance ratio is the geometric mean to account for the asymmetry in the ratios.

* `love.plot` has several changes that make it much more user-friendly. First, rather than supplying a `bal.tab` object to `love.plot`, you can simply supply the arguments that would have gone into the `bal.tab()` call straight into `love.plot()`. Second, if `quick = TRUE` (the new default) and the first argument to `love.plot()` is a call to `bal.tab()` (or arguments provided to `bal.tab()`) and `stat` is set to `"variance.ratios"` or `"ks.statistics"`, `bal.tab()` will be re-called with the corresponding `disp` argument set to `TRUE` so that `love.plot()` will display those statistics regardless of `quick`. This will not work if the argument supplied to `love.plot()` is a `bal.tab` object. Third, because unadjusted mean differences are computed regardless of `quick`, there will never be a circumstance in which only adjusted values will be displayed. If `quick = TRUE`, `un = FALSE`, and `stat` is `"variance.ratios"` or `"ks.statistics"`, `un` will automatically be set to `TRUE` in the `bal.tab()` re-call.

* When using `which.` arguments (e.g., `which.cluster`, `which.imp`, etc.), instead of supplying `NULL` and `NA`, you can supply `.all` and `.none` (not in quotes). This should make them easier to use. Note that these new inputs are not variables; they are keywords and are evaluated using nonstandard evaluation. If you actually have objects with those names, they will be ignored.

* Bugs in scoping related to the formula interface have been solved, in particularly making `bal.tab()` more usable within other functions.

* Fixed bug occurring when using `matchit` objects having set `discard` to something other than `NULL` and `reestimate = TRUE` in the call to `matchit()`. Thank you to Weiyi Xie for finding this bug.

* Fixed bug occurring when using balance thresholds with subclassification.

* Fixed bug occurring when printing `bal.tab` output for continuous treatments with clusters.

* Fixed bug occurring when using `bal.tab()` on `mnps` objects with multiple stop methods.

# cobalt 3.6.1

* Fixed bug when installed version of R was earlier than 3.5.0.

# cobalt 3.6.0

* Added `poly` argument to `bal.tab()` to display polynomials of continuous covariates (e.g., squares, cubes, etc.). This used to only be available with the `int` argument, which also displayed all interactions. Now, the polynomials can be requested separately. When `int = TRUE`, squares of the covariates will no longer be displayed; to replicate the old behavior, set `int = 2`, which is equivalent to `int = TRUE, poly = 2`.

* Fixed a bug where using `subset` would produce an error.

* Fixed a bug when using multiply imputed data with binary treatments that were factors or characters.

* Updated the `bal.tab` documentation to make it easier to navigate to the right page.

* Small documentation and syntax updates.

* Added the hidden and undocumented argument `center` to `bal.tab`, which, when set to `TRUE`, centers the covariates at the mean of the entire unadjusted sample prior to computing interactions and polynomials.

* Added `set.cobalt.options` function to more easily set the global options that can be used as defaults to some arguments. For example, `set.global.options(binary = "std")` makes it so that standardized mean difference are always displayed for binary covariates (in the present R session). The options can be retrieved with `get.cobalt.options`.

# cobalt 3.5.0

* Several changes to `bal.tab()` display options (i.e., `imbalanced.only`, `un`, `disp.means`, `disp.v.ratio`, `disp.ks`, `disp.bal.tab`, `disp.subclass`, and parameters related to the display of balance tables with multinomial treatments, clusters, multiple imputations, and longitudinal treatments). First, the named arguments have been removed from the method-specific functions in order to clean them up and make it easier to add new functions, but they are still available to be specified. Second, a help page devoted just to these functions has been created, which can be accessed with `?options-display`. Third, global options for these arguments can be set with `options()` so they don't need to be typed each time. For example, if you wanted `un = TRUE` all the time, you could set `options(cobalt_un = TRUE)` once and not have to include it in the call to `bal.tab()`.

* Added `disp.sds` option to display standard deviations for each group in `bal.tab()`. This works in all the same places `disp.means` does.

* Added `cluster.fun` and `imp.fun` options to request that only certain functions (e.g., mean or maximum) of the balance statistics are displayed in the summary across clusters/imputations. Previously this option was only available by call `print()`. These parameters are part of the display options described above, so they are documented in `?options-display` and not in the `bal.tab` help files.

* Added `factor_sep` and `int_sep` options to change the separators between variable names when factor variables and interactions are displayed. This functionality had been available since version 3.4.0 but was not documented. It is now documented in the new `display_options` help page.

* In `bal.tab()`, `continuous` and `binary` can be specified with the global options `"cobalt_continuous"` and `"cobalt_binary"`, respectively, so that a global setting (e.g., to set `binary = "std"` to view standardized mean difference rather than raw differences in proportion for binary variables) can be used instead of specifying the argument each time in the call to `bal.tab()`.

* Minor updates to `f.build()` to process inputs more flexibly. The left hand side can now be empty, and the variables on the right hand side can now contain spaces.

* Fixed a bug when logical treatments were used. Thanks to @victorn1.

* Fixed a bug that would occur when a variable had only one value. Thanks to @victorn1.

* Made it so the names of 0/1 and logical variables are not printed with `"_1"` appended to them. Thanks to @victorn1 for the suggestion.

* Major updates to the organization of the code and help files. Certain functions have simplified syntax, relying more on `...`, and help pages have been shorted and consolidated for some methods. In particular, the code and help documents for the `Matching`, `optmatch`, `ebal`, and `designmatch` methods of `bal.tab()` have been consolidated since they all rely on exactly the same syntax.

# cobalt 3.4.1

* Fixed a bug that would occur when `imabalanced.only = TRUE` in `bal.tab()` but all variables were balanced.

* Fixed a bug where the mean of a binary variable would be displayed as 1 minus its mean.

* Fixed a bug that would occur when missingness patterns were the same for multiple variables.

* Fixed a bug that would occur when a distance measure was to be assessed with `bal.tab()` and there were missing values in the covariates (thanks to Laura Helmkamp).

* Fixed a bug that would occur when `estimand` was supplied by the user when using the `default` method of `bal.tab()`.

* Fixed a bug where non-standard variable names (like `"I(age^2)"`) would cause an error.

* Fixed a bug where treatment levels that had different numbers of characters would yield an error.

* Added `disp.means` option to `bal.tab` with continuous treatments.

# cobalt 3.4.0

* Added `default` method for `bal.tab` so it can be used with specially formatted output from other packages (e.g., from `optweight`). `bal.plot` should work with these outputs too. This, of course, will never be completely bug-free because infinite inputs are possible and cannot all be processed perfectly. Don't try to break this function :)

* Fixed some bugs occurring when standardized mean differences are not finite, thanks to Noémie Kiefer.

* Speed improvements in `bal.plot`, especially with multiple facets, and in `bal.tab`.

* Added new options to `bal.plot`, including the ability to display histograms rather than densities and mirrored rather than overlapping plots. This makes it possible to make the popular mirrored histogram plot for propensity scores. In addition, it's now easier to change the colors of the components of the plots.

* Made behavior around binary variables with interactions more like documentation, where interactions with both levels of the variable are present (thanks to @victorn1). Also, replaced `_` with ` * ` as the delimiter between variable names in interactions. For the old behavior, use `int_sep = "_"` in `bal.tab`.

* Expanded the flexibility of `var.names` in `love.plot` so that replacing the name of a variable will replace it everywhere it appears, including interactions. Thanks to @victorn1 for the suggestion.

* Added `var.names` function to extract and save variable names from `bal.tab` objects. This makes it a lot easier to create replacement names for use in `love.plot`. Thanks to @victorn1 for the suggestion.

* When weighted correlations are computed for continuous treatments, the denominator of the correlation now uses the unweighted standard deviations. See `?bal.tab` for the rationale.

# cobalt 3.3.0

* Added methods for objects from the `designmatch` package.

* Added methods for `ps.cont` objects from the `WeightIt` package.

* Fixed bugs resulting form changes to how formula inputs are handled.

* Cleaned up some internal functions, also fixing some related bugs

* Added `subset` option in all `bal.tab()` methods (and consequently in `bal.plot()`) that allows users to specify a subset of the data to assess balance on (i.e., instead of the whole data set). This provides a workaround for methods were the `cluster` option isn't allowed (e.g., longitudinal treatments) but balance is desired on subsets of the data. However, in most cases, `cluster` with `which.cluster` specified makes more sense.

* Updated help files, in particular, more clearly documenting methods for `iptw` objects from `twang` and `CBMSM` objects from `CBPS`.

* Added pretty printing with `crayon`, inspired by Jacob Long's `jtools` package

* Added `abs` option to `bal.tab` to display absolute values of statistics, which can be especially helpful for aggregated output. This also affects how `love.plot()` handles aggregated balance statistics.

# cobalt 3.2.3

* Added support for data with missing covariates. `bal.tab()` will produce balance statistics for the non-missing values and will automatically create a new variable indicating whether the variable is missing or not and produce balance statistics on this variable as well. 

* Fixed a bug when displaying maximum imbalances with subclassification.

* Fixed a bug where the unadjusted statistics were not displayed when using `love.plot()` with subclasses. (Thanks to Megha Joshi.)

* Add the ability to display individual subclass balance using `love.plot()` with subclasses.

* Under-the-hood changes to how `weightit` objects are handled.

* Objects in the environment are now handled better by `bal.tab()` with the formula interface. The `data` argument is now optional if all variables in the formula exist in the environment.

# cobalt 3.2.2

* Fixed a bug when using `get.w()` (and `bal.tab()`) with `mnps` objects from `twang` with only one stop method.

* Fixed a bug when using `bal.tab()` with `twang` objects that contained missing covariate values.

* Fixed a bug when using `int = TRUE` in `bal.tab()` with few covariates.

* Fixed a bug when variable names had special characters.

* Added ability to check higher order polynomials by setting `int` to a number.

* Changed behavior of `bal.tab()` with multinomial treatments and `s.d.denom = "pooled"` to use the pooled standard deviation from the entire sample, not just the paired treatments.

* Restored some vignettes that required `WeightIt`.

# cobalt 3.2.1

* Edits to vignettes and help files to respond to missing packages. Some vignette items may not display if packages are (temporarily) unavailable.

* Fixed issue with sampling weights in `CBPS` objects. (Thanks to @kkranker on Github.)

* Added more support for sampling weights in `get.w()` and help files.

# cobalt 3.2.0

* Added support for longitudinal treatments in `bal.tab()`, `bal.plot()`, and `love.plot()`, including output from `iptw()` in `twang`, `CBMSM()` from `CBPS`, and `weightitMSM()` from `WeightIt`.

* Added a vignette to explain use with longitudinal treatments.

* Edits to help files.

* Added ability to change density options in `bal.plot()`.

* Added support for `imp` in `bal.tab()` for `weightit` objects.

* Fixed bugs when limited variables were present. (One found and fixed by @sumtxt on Github.)

* Fixed bug with multiple methods when weights were entered as a list.

# cobalt 3.1.0

* Added full support for tibbles.

* Examples for `weightit` methods in documentation and vignette now work.

* Improved speed and performance.

* Added `pairwise` option for `bal.tab()` with multinomial treatments.

* Increased flexibility for displaying balance using `love.plot()` with clustered or multiply imputed data.

* Added `imbalanced.only` and `disp.bal.tab` options to `bal.tab()`.

* Fixes to the vignettes. Also, creation of a new vignette to simplify the main one.

# cobalt 3.0.0

* Added support for multinomial treatments in `bal.tab()`, including output from `CBPS` and `twang`.

* Added support for `weightit` objects from `WeightIt`, including for multinomial treatments.

* Added support for `ebalance.trim` objects from `ebal`.

* Fixes to the vignette.

* Fixes to `splitfactor()` to handle tibbles better. 

* Fixed bug when using `bal.tab()` with multiply imputed data without adjustment. Fixed bug when using `s.weights` with the `formula` method of `bal.tab()`.

# cobalt 2.2.0

* Added `disp.ks` and `ks.threshold` options to `bal.tab()` to display Kolmogorov-Smirnov statistics before and after preprocessing.

* Added support for sampling weights, which are applied to both control and treated units, using option `s.weights` in `bal.tab()`. Sampling weights are also now compatible with the sampling weights in `ps` objects from `twang`; the default is to apply the sampling weights before and after adjustment, mimicking the behavior of `bal.table()` in `twang`.

* Changed behavior of `bal.tab()` for `ps` objects to allow displaying balance for more than one stop method at a time, and to default to displaying balance for all available stop methods. The `full.stop.method` argument in `bal.tab()` has been renamed `stop.method`, but `full.stop.method` still works. `get.w()` for `ps` objects has also gone through some changes to be more like `twang`'s `get.weights()`.

* Added support in `bal.tab()` and `bal.plot()` for subclassification with continuous treatments.

* Added support in `splitfactor()` and `unsplitfactor()` for `NA` values

* Fixed a bug in `love.plot()` caused when `var.order` was specified to be a sample that was not present.

# cobalt 2.1.0

* Added support in `bal.tab()`, `bal.plot()`, and `love.plot()` for examining balance on multiple weight specifications at a time

* Added new utilities `splitfactor()`, `unsplitfactor()`, and `get.w()`

* Added option in `bal.plot()` to display points sized by weights when treatment and covariate are continuous

* Added `which = "both"` option in `bal.plot()` to simultaneously display plots for both adjusted and unadjusted samples; changed argument syntax to accommodate

* Allowed `bal.plot()` to display balance for multiple clusters and imputations simultaneously

* Allowed `bal.plot()` to display balance for multiple subclasses simultaneously with `which.sub`

* Fixes to `love.plot()` to ensure adjusted points are in front of unadjusted points; changed colors and shape defaults and allowable values

* Fixed bug where `s.d.denom` and `estimand` were not functioning correctly in `bal.tab()`

* `distance`, `addl`, and `weights` can now be specified as lists of the usual arguments

# cobalt 2.0.0

* Added support for matching using the `optmatch` package or by specifying matching strata.

* Added full support (`bal.tab()`, `love.plot()`, and `bal.plot()`) for multiply imputed data, including for clustered data sets.

* Added support for multiple distance measures, including special treatment in `love.plot()`

* Adjusted specifications in `love.plot()` for color and shape of points, and added option to generate a line connecting the points.

* Adjusted `love.plot()` display to perform better on Windows.

* Added capabilities for `love.plot()` and `bal.plot()` to display plots for multiple groups at a time

* Added flexibility to `f.build()`.

* Updated `bal.plot()`, giving the capability to view multiple plots for subclassified or clustered data. Multinomial treatments are also supported.

* Created a new vignette for clustered and multiply imputed data

* Speed improvements

* Fixed a bug causing mislabelling of categorical variables

* Changed calculation of weighted variance to be in line with recommendations; `CBPS` can now be used with standardized weights

# cobalt 1.3.1

* Added support for entropy balancing through the `ebal` package.

* Changed default color scheme of `love.plot()` to be black and white and added options for color, shape, and size of points.

* Added sample size calculations for continuous treatments.

* Edits to the vignette.
    
# cobalt 1.3.0

* Increased capabilities for cluster balance in `bal.tab()` and `love.plot()`

* Increased information and decreased redundancy when assessing balance on interactions

* Added `quick` option for `bal.tab()` to increase speed

* Added options for `print()`

* Bug fixes

* Speed improvements

* Edits to the vignette

# cobalt 1.2.0

* Added support for continuous treatment variables in `bal.tab()`, `bal.plot()`, and `love.plot()`

* Added balance assessment within and across clusters

* Other small performance changes to minimize errors and be more intuitive

* Major revisions and adjustments to the vignette

# cobalt 1.1.0

* Added a vignette.

* Fixed error in `bal.tab.Match()` that caused wrong values and and warning messages when used.

* Added new capabilities to `bal.plot()`, including the ability to view unadjusted sample distributions, categorical variables as such, and the distance measure. Also updated documentation to reflect these changes and make `which.sub` more focal.

* Allowed subclasses to be different from simply 1:S by treating them like factors once input is numerical

* Changed column names in Balance table output to fit more compactly, and updated documentation to reflect these changes.

* Other small performance changes to minimize errors and be more intuitive.
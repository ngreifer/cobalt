cobalt News and Updates
======
Version 3.1.0

* Added full support for tibbles.

* Examples for `weightit` methods in documentation and vignette now work.

* Improved speed and performance.

* Added `pairwise` option for `bal.tab()` with multinomial treatments.

* Increased flexibility for displaying balance using `love.plot()` with clustered or multiply imputed data.

* Added `imbalanced.only` and `disp.bal.tab` options to `bal.tab()`.

* Fixes to the vignettes. Also, creation of a new vignette to simplify the main one.

Version 3.0.0

* Added support for multinomial treatments in `bal.tab()`, including output from `CBPS` and `twang`.

* Added support for `weightit` objects from `WeightIt`, including for multinomial treatments.

* Added support for `ebalance.trim` objects from `ebal`.

* Fixes to the vignette.

* Fixes to `splitfactor()` to handle tibbles better. 

* Fixed bug when using `bal.tab()` with multiply imputed data without adjustment. Fixed bug when using `s.weights` with the `formula` method of `bal.tab()`.

Version 2.2.0

* Added `disp.ks` and `ks.threshold` options to `bal.tab()` to display Kolmogorov-Smirnov statistics before and after preprocessing.

* Added support for sampling weights, which are applied to both control and treated units, using option `s.weights` in `bal.tab()`. Sampling weights are also now compatible with the sampling weights in `ps` objects from `twang`; the default is to apply the sampling weights before and after adjustment, mimicking the behavior of `bal.table()` in `twang`.

* Changed behavior of `bal.tab()` for `ps` objects to allow displaying balance for more than one stop method at a time, and to default to displaying balance for all available stop methods. The `full.stop.method` argument in `bal.tab()` has been renamed `stop.method`, but `full.stop.method` still works. `get.w()` for `ps` objects has also gone through some changes to be more like `twang`'s `get.weights()`.

* Added support in `bal.tab()` and `bal.plot()` for subclassification with continuous treatments.

* Added support in `splitfactor()` and `unsplitfactor()` for `NA` values

* Fixed a bug in `love.plot()` caused when `var.order` was specified to be a sample that was not present.

Version 2.1.0

* Added support in `bal.tab()`, `bal.plot()`, and `love.plot()` for examining balance on multiple weight specifications at a time

* Added new utilities `splitfactor()`, `unsplitfactor()`, and `get.w()`

* Added option in `bal.plot()` to display points sized by weights when treatment and covariate are continuous

* Added `which = "both"` option in `bal.plot()` to simultaneously display plots for both adjusted and unadjusted samples; changed argument syntax to accommodate

* Allowed `bal.plot()` to display balance for mutliple clusters and imputations simultaneously

* Allowed `bal.plot()` to display balance for mutliple subclasses simultaneously with `which.sub`

* Fixes to `love.plot()` to ensure adjusted points are in front of unadjusted points; changed colors and shape defaults and allowable values

* Fixed bug where `s.d.denom` and `estimand` were not functioning correctly in `bal.tab()`

* `distance`, `addl`, and `weights` can now be specified as lists of the usual arguments

Version 2.0.0

* Added support for matching using the `optmatch` package or by specifying matching strata.

* Added full support (`bal.tab()`, `love.plot()`, and `bal.plot()`) for multiply imputed data, including for clustered data sets.

* Added support for multiple distance measures, including special treatment in `love.plot()`

* Adjusted specifications in `love.plot()` for color and shape of points, and added option to generate a line connecting the points.

* Adjusted `love.plot()` display to perform better on Windows.

* Added capabilties for `love.plot()` and `bal.plot()` to display plots for multiple groups at a time

* Added flexibility to `f.build()`.

* Updated `bal.plot()`, giving the capability to view multiple plots for subclassified or clustered data. Multinomial treatments are also supported.

* Created of a new vignette for clustered and multiply imputed data

* Speed improvements

* Fixed a bug causing mislabelling of categorical variables

* Changed calculation of weighted variance to be in line with recommendations; CBPS can now be used with standardized weights

Version 1.3.1

* Added support for entropy balancing through the ebal package.

* Changed default color scheme of `love.plot()` to be black and white and added options for color, shape, and size of points.

* Added sample size calculations for continuous treatments.

* Edits to the vignette.
    
Version 1.3.0

* Increased capabilities for cluster balance in `bal.tab()` and `love.plot()`

* Increased information and decreased redundancy when assessing balance on interactions

* Added "quick" option for `bal.tab()` to increase speed

* Added options for `print()`

* Bug fixes

* Speed improvements

* Edits to the vignette

Version 1.2.0

* Added support for continous treatment variables in `bal.tab()`, `bal.plot()`, and `love.plot()`

* Added balance assessment within and across clusters

* Other small performance changes to minimize errors and be more intuitive

* Major revisions and adjustments to the vignette

Version 1.1.0

* Added a vignette.

* Fixed error in bal.tab.Match that caused wrong values and and warning messages when used.

* Added new capabilities to bal.plot, including the ability to view unadjusted sample distributions, categorical variables as such, and the distance measure. Also updated documentation to reflect these changes and make which.sub more focal.

* Allowed subclasses to be different from simply 1:S by treating them like factors once input is numerical

* Changed column names in Balance table output to fit more compactly, and updated documentation to reflect these changes.

* Other small performance changes to minimize errors and be more intuitive.
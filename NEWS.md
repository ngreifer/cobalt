cobalt News and Updates
======
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
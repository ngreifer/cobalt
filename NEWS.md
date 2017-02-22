cobalt News and Updates
======
Version 1.3.2

*Adjusted specifications in `love.plot()` for color and shape of poins, and added option to generate a line

*Adjusted `love.plot()` display to perform better on Windows

Version 1.3.1

*Added support for entropy balancing through the ebal package.

*Changed default color scheme of `love.plot()` to be black and white and added options for color, shape, and size of points

*Added sample size calculations for continuous treatments

*Edits to the vignette
    
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
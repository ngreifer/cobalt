Things to Do
======

* Add default method that takes in an arbitrary object and other inputs and checks for treatment, covariates, and a conditioning method, etc. Add in documentation and vignette that programmers of other packages can rely on it by ensuring the inputs match.

* Add method for multinomial treatments with multiply imputed data

* Add method to add own balance function. Could be supplied as a name (e.g., skew.diff) that takes in several arguments: treat, covs (or maybe one cov), weights, and maybe other arguments.

* Add ability to submit multiple objects at the same time rather than by using df/formula and providing weights.

* Add ability to examine pair distances for matched sets. Maybe mean absolute and max absolute distance between pairs for each covariate. For binary, proportion of units with the same values.

* Allow users to enter matched data set and unmatched data set for comparison. Could have an additional argument in df/formula method for `matched.data`, for exampe.
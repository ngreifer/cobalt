Things to Do
======

* Add method for multinomial treatments with multiply imputed data

* Add method to add own balance function. Could be supplied as a name (e.g., skew.diff) that takes in several arguments: treat, covs (or maybe one cov), weights, and maybe other arguments.

* Add ability to submit multiple objects at the same time rather than by using df/formula and providing weights.

* Add ability to examine pair distances for matched sets. Maybe mean absolute and max absolute distance between pairs for each covariate. For binary, proportion of units with the same values.

* Allow users to enter matched data set and unmatched data set for comparison. Could have an additional argument in df/formula method for `matched.data`, for exampe.

* Add vignette specifically about `love.plot` and its capabilities
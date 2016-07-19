
<!-- README.md is generated from README.Rmd. Please edit that file -->
cobalt
======

Welcome to `cobalt`, which stands for **Co**variate **Bal**ance **T**ables (and Plots). `cobalt` allows users to assess balance on covariate distributions in preprocessed groups generated through weighting, matching, or subclassification, such as by using the propensity score. `cobalt`'s primary function is `bal.tab()`, which stands for "balance table", and essentially repalces (or supplements) the balance assessment tool found in the R packages twang, MatchIt, CBPS, and Matching. To examine how `bal.tab()` integrates with these packages, see the help file for `bal.tab()` with `?bal.tab`, which links to the methods used for each package. Each page has examples of how `bal.tab()` is used with the package.

Also included are two plotting functions, `bal.plot()` and `love.plot()`, and a utility function, `f.build()`. See the help files for these functions to learn what they do and how to use them.

Please remember to cite this package when using it to analyze data. For example, in a manuscript, write: "Matching was performed using Matching (Sekhon, 2011), and covariate balance was assessed using cobalt (Greifer, 2016) in R (R Core team, 2016)." Use `citation("cobalt")` to generate a bibliographic reference for the `cobalt` package. Thank you very much!

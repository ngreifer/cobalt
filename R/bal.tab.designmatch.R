#' Balance Statistics for `designmatch` Objects
#' 
#' @description Generates balance statistics for output objects from \pkg{designmatch}.
#' 
#' @inheritParams bal.tab.Match
#' @param x the output of a call to \pkgfun{designmatch}{bmatch} or related wrapper functions from the \pkg{designmatch} package.
#' @param s.d.denom `character`; how the denominator for standardized mean differences should be calculated, if requested. See [col_w_smd()] for allowable options. Abbreviations allowed. If not specified, will be set to `"treated"`.
#' 
#' @inherit bal.tab.Match return
#' 
#' @details `bal.tab()` generates a list of balance summaries for the object given, and functions similarly to \pkgfun{designmatch}{meantab}. Note that output objects from \pkg{designmatch} do not have their own class; `bal.tab()` first checks whether the object meets the criteria to be treated as a `designmatch` object before dispatching the correct method. Renaming or removing items from the output object can create unintended consequences.
#' 
#' The input to `bal.tab.designmatch()` must include either both `formula` and `data` or both `covs` and `treat`. Using the `covs` + `treat` input mirrors how \pkgfun{designmatch}{meantab} is used (note that to see identical results to `meantab()`, `s.d.denom` must be set to `"pooled"`).
#' 
#' @inherit bal.tab.Match seealso
#' 
#' @examplesIf (requireNamespace("designmatch", quietly = TRUE))
#' \donttest{
#' data("lalonde", package = "cobalt")
#' 
#' library(designmatch)
#' covariates <- as.matrix(lalonde[c("age", "educ", "re74", "re75")])
#' treat <- lalonde$treat
#' dmout <- bmatch(treat,
#'                 total_groups = sum(treat == 1),
#'                 mom = list(covs = covariates,
#'                            tols = absstddif(covariates, 
#'                                             treat, .05))
#' )
#' 
#' ## Using treat and covs
#' bal.tab(dmout, treat = treat, covs = covariates)
#' }

#' @exportS3Method bal.tab designmatch
bal.tab.designmatch <- bal.tab.Match

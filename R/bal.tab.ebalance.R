#' Balance Statistics for `ebalance` Objects
#' 
#' @description Generates balance statistics for output objects from \pkg{ebal}.
#' 
#' @inheritParams bal.tab.Match
#' @param x an `ebalance` object (the output of a call to \pkgfun{ebal}{ebalance} or \pkgfun{ebal}{ebalance.trim}).
#' @param s.d.denom `character`; how the denominator for standardized mean differences should be calculated, if requested. See [col_w_smd()] for allowable options. Abbreviations allowed. If not specified, will be set to `"treated"`.
#' 
#' @inherit bal.tab.Match return
#' 
#' @details
#' `bal.tab()` generates a list of balance summaries for the object given. The input to `bal.tab.ebalance()` must include either both `formula` and `data` or both `covs` and `treat`.
#' 
#' @inherit bal.tab.Match seealso
#' 
#' @examplesIf rlang::is_installed("ebal")
#' data("lalonde", package = "cobalt")
#' 
#' covs <- subset(lalonde, select = -c(re78, treat))
#' covs0 <- splitfactor(covs)
#' 
#' e.out <- ebal::ebalance(lalonde$treat, covs0)
#' 
#' ## Using formula and data
#' bal.tab(e.out, formula = treat ~ age + educ + race +
#'             married + nodegree + re74 + re75,
#'         data = lalonde)
#' 
#' ## Using treat and covs
#' bal.tab(e.out, treat = lalonde$treat, covs = covs)

#' @exportS3Method bal.tab ebalance
bal.tab.ebalance <- bal.tab.Match

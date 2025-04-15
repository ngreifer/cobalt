#' @title Balance Statistics for `sbw` Objects
#' 
#' @description
#' Generates balance statistics for `sbwcau` objects from \pkg{sbw}.
#' 
#' @inheritParams bal.tab
#' @param x an `sbwcau` object; the output of a call to \pkgfun{sbw}{sbw}.
#' @param s.d.denom `character`; how the denominator for standardized mean differences should be calculated, if requested. See [col_w_smd()] for allowable options. Abbreviations allowed. If not specified, `bal.tab()` will figure out which one is best based on the `par` component of the `sbwcau` object: if "att", `"treated"`; if "atc", `"control"`; otherwise `"pooled"`.
#' 
#' @returns
#' If clusters are not specified, an object of class `"bal.tab"` containing balance summaries for the `sbwcau` object. See [bal.tab()] for details.
#' 
#' If clusters are specified, an object of class `"bal.tab.cluster"` containing balance summaries within each cluster and a summary of balance across clusters. See [`class-bal.tab.cluster`] for details.
#' 
#' @details
#' `bal.tab.sbwcau()` generates a list of balance summaries for the `sbwcau` object given, and functions similarly to \pkgfun{sbw}{summarize}.
#' 
#' @seealso
#' * [bal.tab()] for details of calculations.
#' 
#' @examplesIf requireNamespace("sbw", quietly = TRUE)
#' library(sbw)
#' data("lalonde", package = "cobalt")
#' 
#' ## Stable balancing weights for the ATT
#' sbw.out <- sbw(splitfactor(lalonde, drop.first = "if2"),
#'                ind = "treat",
#'                bal = list(bal_cov = c("age", "educ", "race_black", 
#'                                       "race_hispan", "race_white", 
#'                                       "married", "nodegree", 
#'                                       "re74", "re75"),
#'                           bal_alg = FALSE, 
#'                           bal_tol = .001),
#'                par = list(par_est = "att"))
#' 
#' bal.tab(sbw.out, un = TRUE, poly = 2)

#' @exportS3Method bal.tab sbwcau
bal.tab.sbwcau <- function(x,
                           stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, s.weights = NULL, abs = FALSE, subset = NULL, quick = TRUE,
                           ...) {
  
  args <- try_chk(c(as.list(environment()), list(...))[-1L])
  
  #Adjustments to arguments
  
  args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
  args[lengths(args) == 0L & names(args) %nin% names(match.call())[-1L]] <- NULL
  
  #Initializing variables
  X <- do.call("x2base.sbwcau", c(list(x), args), quote = TRUE)
  
  args[names(X)] <- NULL
  
  X <- .assign_X_class(X)
  
  do.call("base.bal.tab", c(list(X), args),
          quote = TRUE)
}

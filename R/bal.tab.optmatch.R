#' @title Balance Statistics for `optmatch` Objects
#' @description Generates balance statistics for output objects from \pkg{optmatch}.
#' 
#' @inheritParams bal.tab.Match
#' @param x an `optmatch` object (the output of a call to \pkgfun{optmatch}{pairmatch} or \pkgfun{optmatch}{fullmatch}).
#' @param estimand `character`; whether the desired estimand is the "ATT", "ATC", or "ATE". Default is "ATT".
#' @param s.d.denom `character`; how the denominator for standardized mean differences should be calculated, if requested. See [col_w_smd()] for allowable options. Abbreviations allowed. If unspecified, `bal.tab()` will figure out which one is best based on the `estimand`, if given (for ATT, `"treated"`; for ATC, `"control"`; otherwise `"pooled"`) and other clues if not.
#' 
#' @inherit bal.tab.Match return
#' 
#' @details `bal.tab()` generates a list of balance summaries for the object given. The input to `bal.tab.optmatch()` must include either both `formula` and `data` or just `covs` (`treat` is not necessary).
#' 
#' @inherit bal.tab.Match seealso
#' 
#' @examplesIf requireNamespace("optmatch", quietly = TRUE)
#' data("lalonde", package = "cobalt")
#' 
#' lalonde$prop.score <- glm(treat ~ age + educ + race + 
#'                               married + nodegree + re74 + re75, 
#'                           data = lalonde, family = binomial)$fitted.values
#' pm <- optmatch::pairmatch(treat ~ prop.score, data = lalonde)
#' 
#' ## Using formula and data; LHS of formula not required
#' bal.tab(pm, formula = ~ age + educ + race +
#'             married + nodegree + re74 + re75,
#'         data = lalonde)
#' 
#' ## Using covs
#' covs <- subset(lalonde, select = -c(re78, treat))
#' bal.tab(pm, covs = covs)

#' @exportS3Method bal.tab optmatch
bal.tab.optmatch <- function(x, formula = NULL, data = NULL, treat = NULL, covs = NULL, estimand = NULL,
                             stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, s.weights = NULL, abs = FALSE, subset = NULL, quick = TRUE,
                             ...) {
  
  args <- try_chk(c(as.list(environment()), list(...))[-1L])
  
  #Adjustments to arguments
  
  args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
  args[lengths(args) == 0L & names(args) %nin% names(match.call())[-1L]] <- NULL
  
  #Initializing variables
  X <- do.call("x2base", c(list(x), args), quote = TRUE) 
  
  args[names(X)] <- NULL
  
  X <- .assign_X_class(X)
  
  do.call("base.bal.tab", c(list(X), args),
          quote = TRUE)
}

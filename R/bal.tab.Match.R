#' @title Balance Statistics for `Matching` Objects
#' @description Generates balance statistics for output objects from \pkg{Matching}.
#' 
#' @inheritParams bal.tab
#' @param x a `Match` object (the output of a call to \pkgfun{Matching}{Match} or \pkgfun{Matching}{Matchby}).
#' @param formula a `formula` with the treatment variable as the response and the covariates for which balance is to be assessed as the predictors. All named variables must be in `data`. See Details.
#' @param data a data frame containing variables named in `formula`, if supplied, and other arguments.
#' @param treat a vector of treatment statuses. See Details.
#' @param covs a data frame of covariate values for which to check balance. See Details.
#' @param s.d.denom `character`; how the denominator for standardized mean differences should be calculated, if requested. See [col_w_smd()] for allowable options. Abbreviations allowed. If not specified, `bal.tab()` will use "treated" if the estimand of the call to `Match()` is the ATT, "pooled" if the estimand is the ATE, and "control" if the estimand is the ATC.
#' 
#' @returns If clusters and imputations are not specified, an object of class `"bal.tab"` containing balance summaries for the given object. See [bal.tab()] for details.
#' 
#' If clusters are specified, an object of class `"bal.tab.cluster"` containing balance summaries within each cluster and a summary of balance across clusters. See [`class-bal.tab.cluster`] for details.
#' 
#' @details `bal.tab()` generates a list of balance summaries for the object given, and functions similarly to \pkgfun{Matching}{MatchBalance}. The input to `bal.tab.Match()` must include either both `formula` and `data` or both `covs` and `treat`. Using the `formula` + `data` inputs mirrors how \pkgfun{Matching}{MatchBalance} is used.
#' 
#' `cobalt` functions do not support `Match` object with sampling weights, i.e., with an argument passed to the `weights` argument of `Matching::Match()`.
#' 
#' @seealso [bal.tab()] for details of calculations.
#' 
#' @examplesIf requireNamespace("Matching", quietly = TRUE)
#' library(Matching); data("lalonde", package = "cobalt")
#' 
#' p.score <- glm(treat ~ age + educ + race + 
#'                    married + nodegree + re74 + re75, 
#'                data = lalonde, family = "binomial")$fitted.values
#' Match.out <- Match(Tr = lalonde$treat, X = p.score)
#' 
#' ## Using formula and data
#' bal.tab(Match.out, formula = treat ~ age + educ + race + 
#'             married + nodegree + re74 + re75, data = lalonde)

#' @exportS3Method bal.tab Match
bal.tab.Match <-      function(x, formula = NULL, data = NULL, treat = NULL, covs = NULL,
                               stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, s.weights = NULL, abs = FALSE, subset = NULL, quick = TRUE,
                               ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) .err(conditionMessage(e)))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    args[vapply(args, is_null, logical(1L)) & names(args) %nin% names(match.call())[-1]] <- NULL
    
    #Initializing variables
    X <- do.call("x2base", c(list(x), args), quote = TRUE) 
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- .assign_X_class(X)
    
    do.call("base.bal.tab", c(list(X), args),
            quote = TRUE)
}

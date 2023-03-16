#' @title Balance Statistics for \code{cem} Objects
#' 
#' @description
#' Generates balance statistics for \code{cem.match} objects from \pkg{cem}.
#' 
#' @inheritParams bal.tab
#' @param x a \code{cem.match} or \code{cem.match.list} object; the output of a call to \pkgfun{cem}{cem}.
#' @param data a data frame containing variables named in other arguments. An argument to \code{data} is \strong{required}. It must be the same data used in the call to \code{cem()} or a \code{mids} object from which the data supplied to \code{datalist} in the `cem()` call originated.
#' @param s.d.denom \code{character}; how the denominator for standardized mean differences should be calculated, if requested. See [col_w_smd()] for allowable options. Abbreviations allowed. If not specified, will be set to `"treated"`, where the treated group corresponds to the \code{baseline.group} in the call to \code{cem()}.
#' 
#' @return
#' If clusters and imputations are not specified, an object of class \code{"bal.tab"} containing balance summaries for the \code{cem.match} object. See \fun{bal.tab} for details.
#' 
#' If imputations are specified, an object of class \code{"bal.tab.imp"} containing balance summaries for each imputation and a summary of balance across imputations. See \code{\link[=class-bal.tab.imp]{bal.tab.imp}} for details.
#' 
#' If \code{cem()} is used with multi-category treatments, an object of class \code{"bal.tab.multi"} containing balance summaries for each pairwise treatment comparison. See \code{\link[=class-bal.tab.multi]{bal.tab.multi}} for details.
#' 
#' If clusters are specified, an object of class \code{"bal.tab.cluster"} containing balance summaries within each cluster and a summary of balance across clusters. See \code{\link[=class-bal.tab.cluster]{bal.tab.cluster}} for details.
#' 
#' @details
#' \code{bal.tab.cem.match()} generates a list of balance summaries for the \code{cem.match} object given, and functions similarly to \pkgfun{cem}{imbalance}.
#' 
#' @seealso [bal.tab()] for details of calculations.
#' 
#' @examplesIf requireNamespace("cem", quietly = TRUE)
#' data("lalonde", package = "cobalt")
#' 
#' ## Coarsened exact matching
#' cem.out <- cem::cem("treat", data = lalonde, drop = "re78")
#' 
#' bal.tab(cem.out, data = lalonde, un = TRUE, 
#'         stats = c("m", "k"))

#' @exportS3Method bal.tab cem.match
bal.tab.cem.match <-  function(x, data,
                               stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, s.weights = NULL, abs = FALSE, subset = NULL, quick = TRUE,
                               ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) .err(conditionMessage(e)))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    args[vapply(args, is_null, logical(1L)) & names(args) %nin% names(match.call())[-1]] <- NULL
    
    #Initializing variables
    X <- do.call("x2base.cem.match", c(list(x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    return(out)
}


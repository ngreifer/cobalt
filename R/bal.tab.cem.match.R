#' @title Balance Statistics for `cem` Objects
#' 
#' @description
#' Generates balance statistics for `cem.match` objects from \pkg{cem}.
#' 
#' @inheritParams bal.tab
#' @param x a `cem.match` or `cem.match.list` object; the output of a call to \pkgfun{cem}{cem}.
#' @param data a data frame containing variables named in other arguments. An argument to `data` is **required**. It must be the same data used in the call to `cem()` or a `mids` object from which the data supplied to `datalist` in the `cem()` call originated.
#' @param s.d.denom `character`; how the denominator for standardized mean differences should be calculated, if requested. See [col_w_smd()] for allowable options. Abbreviations allowed. If not specified, will be set to `"treated"`, where the treated group corresponds to the `baseline.group` in the call to `cem()`.
#' 
#' @returns
#' If clusters and imputations are not specified, an object of class `"bal.tab"` containing balance summaries for the `cem.match` object. See [bal.tab()] for details.
#' 
#' If imputations are specified, an object of class `"bal.tab.imp"` containing balance summaries for each imputation and a summary of balance across imputations. See [`class-bal.tab.imp`] for details.
#' 
#' If `cem()` is used with multi-category treatments, an object of class `"bal.tab.multi"` containing balance summaries for each pairwise treatment comparison. See [`bal.tab.multi()`][class-bal.tab.multi] for details.
#' 
#' If clusters are specified, an object of class `"bal.tab.cluster"` containing balance summaries within each cluster and a summary of balance across clusters. See [`class-bal.tab.cluster`] for details.
#' 
#' @details
#' `bal.tab.cem.match()` generates a list of balance summaries for the `cem.match` object given, and functions similarly to \pkgfun{cem}{imbalance}.
#' 
#' @seealso [bal.tab()] for details of calculations.
#' 
#' @examplesIf requireNamespace("cem", quietly = TRUE) && FALSE
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
  
  args <- try_chk(c(as.list(environment()), list(...))[-1L])
  
  #Adjustments to arguments
  
  args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
  args[lengths(args) == 0L & names(args) %nin% names(match.call())[-1L]] <- NULL
  
  #Initializing variables
  X <- do.call("x2base.cem.match", c(list(x), args), quote = TRUE)
  
  args[names(X)] <- NULL
  
  X <- .assign_X_class(X)
  
  do.call("base.bal.tab", c(list(X), args),
          quote = TRUE)
}


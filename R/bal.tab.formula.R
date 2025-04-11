#' @title Balance Statistics for Data Sets
#' 
#' @description Generates balance statistics for unadjusted, matched, weighted, or stratified data using either a `data.frame` or formula interface.
#' 
#' @inheritParams bal.tab
#' @param x either a `data.frame` containing covariate values for each unit or a `formula` with the treatment variable as the response and the covariates for which balance is to be assessed as the terms. If a formula is supplied, all terms must be present as variable names in `data` or the global environment.
#' @param treat either a vector containing treatment status values for each unit or a string containing the name of the treatment variable in `data`. Required for the `data.frame` method.
#' @param s.d.denom `character`; how the denominator for standardized mean differences should be calculated, if requested. See [col_w_smd()] for allowable options. Abbreviations allowed. If weights are supplied, each set of weights should have a corresponding entry to `s.d.denom`; a single entry will be recycled to all sets of weights. If left blank and one of `weights`, `subclass`, or `match.strata` are supplied, `bal.tab()` will figure out which one is best based on `estimand`, if given (for ATT, `"treated"`; for ATC, `"control"`; otherwise "pooled") and other clues if not.
#' @param subclass optional; either a vector containing subclass membership for each unit or a string containing the name of the subclass variable in `data`.
#' @param match.strata optional; either a vector containing matching stratum membership for each unit or a string containing the name of the matching stratum variable in `data`. See Details.
#' @param method `character`; the method of adjustment, if any. If `weights` are specified, the user can specify either "matching" or "weighting"; "weighting" is the default. If multiple sets of weights are used, each must have a corresponding value for `method`, but if they are all of the same type, only one value is required. If `subclass` is specified, "subclassification" is the default. Abbreviations allowed. The only distinction between "matching" and "weighting" is how sample sizes are displayed.
#' @param estimand `character`; whether the desired estimand is the "ATT", "ATC", or "ATE" for each set of weights. This argument can be used in place of `s.d.denom` to specify how standardized differences are calculated.
#' @param focal the name of the focal treatment when multi-category treatments are used. See [`bal.tab.multi()`][class-bal.tab.multi] for details.
#' 
#' @returns 
#' For point treatments, if clusters and imputations are not specified, an object of class `"bal.tab"` containing balance summaries for the specified treatment and covariates. See [bal.tab()] for details.
#' 
#' If imputations are specified, an object of class `"bal.tab.imp"` containing balance summaries for each imputation and a summary of balance across imputations. See [`class-bal.tab.imp`] for details.
#' 
#' If multi-category treatments are used, an object of class `"bal.tab.multi"` containing balance summaries for each pairwise treatment comparison. See [`bal.tab.multi()`][class-bal.tab.multi] for details.
#' 
#' If clusters are specified, an object of class `"bal.tab.cluster"` containing balance summaries within each cluster and a summary of balance across clusters. See [`class-bal.tab.cluster`] for details.
#' 
#' @details 
#' `bal.tab.data.frame()` generates a list of balance summaries for the covariates and treatment status values given. `bal.tab.formula()` does the same but uses a formula interface instead.  When the formula interface is used, the formula and data are reshaped into a treatment vector and `data.frame` of covariates and then simply passed through the `data.frame` method.  
#' 
#' If `weights`, `subclass` and `match.strata` are all `NULL`, balance information will be presented only for the unadjusted sample.
#' 
#' The argument to `match.strata` corresponds to a factor vector containing the name or index of each pair/stratum for units conditioned through matching, for example, using the \pkg{optmatch} package. If more than one of `weights`, `subclass`, or `match.strata` are specified, `bal.tab()` will attempt to figure out which one to apply. Currently only one of these can be applied ta a time. `bal.tab()` behaves differently depending on whether subclasses are used in conditioning or not. If they are used, `bal.tab()` creates balance statistics for each subclass and for the sample in aggregate. See [`class-bal.tab.subclass`] for more information.
#' 
#' Multiple sets of weights can be supplied simultaneously by entering a `data.frame` or a character vector containing the names of weight variables found in `data` or a list of weights vectors or names. The arguments to `method`, `s.d.denom`, and `estimand`, if any, must be either the same length as the number of sets of weights or of length one, where the sole entry is applied to all sets. When standardized differences are computed for the unadjusted group, they are done using the first entry to `s.d.denom` or `estimand`. When only one set of weights is supplied, the output for the adjusted group will simply be called `"Adj"`, but otherwise will be named after each corresponding set of weights. Specifying multiple sets of weights will also add components to other outputs of `bal.tab()`.
#' 
#' @seealso 
#' * [bal.tab()] for details of calculations.
#' * [`class-bal.tab.cluster`] for more information on clustered data.
#' * [`class-bal.tab.imp`] for more information on multiply imputed data.
#' * [`bal.tab.multi()`][class-bal.tab.multi] for more information on multi-category treatments.
#' 
#' @examples
#' data("lalonde", package = "cobalt")
#' lalonde$p.score <- glm(treat ~ age + educ + race, data = lalonde, 
#'                        family = "binomial")$fitted.values
#' covariates <- subset(lalonde, select = c(age, educ, race))
#' 
#' ## Propensity score weighting using IPTW
#' lalonde$iptw.weights <- ifelse(lalonde$treat==1, 
#'                                1/lalonde$p.score, 
#'                                1/(1-lalonde$p.score))
#' 
#' # data frame interface:
#' bal.tab(covariates, treat = "treat", data = lalonde, 
#'         weights = "iptw.weights", s.d.denom = "pooled")
#' 
#' # Formula interface:
#' bal.tab(treat ~ age + educ + race, data = lalonde, 
#'         weights = "iptw.weights", s.d.denom = "pooled")
#' 
#' ## Propensity score subclassification
#' lalonde$subclass <- findInterval(lalonde$p.score, 
#'                                  quantile(lalonde$p.score, 
#'                                           (0:6)/6), all.inside = TRUE)
#' 
#' # data frame interface:
#' bal.tab(covariates, treat = "treat", data = lalonde, 
#'         subclass = "subclass", disp.subclass = TRUE, 
#'         s.d.denom = "pooled")
#' 
#' # Formula interface:
#' bal.tab(treat ~ age + educ + race, data = lalonde, 
#'         subclass = "subclass", disp.subclass = TRUE, 
#'         s.d.denom = "pooled")

#' @exportS3Method bal.tab formula
#' @rdname bal.tab.formula
bal.tab.formula <-    function(x, data = NULL,
                               stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, s.weights = NULL, abs = FALSE, subset = NULL, quick = TRUE,
                               subclass = NULL, match.strata = NULL, method, estimand = NULL, focal = NULL, ...) {
  
  args <- try_chk(c(as.list(environment()), list(...))[-1L])
  
  #Adjustments to arguments
  
  args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
  args[lengths(args) == 0L & names(args) %nin% names(match.call())[-1L]] <- NULL
  
  #Initializing variables
  X <- do.call("x2base.formula", c(list(x = x), args), quote = TRUE)
  
  args[names(X)] <- NULL
  
  X <- .assign_X_class(X)
  
  do.call("base.bal.tab", c(list(X), args),
          quote = TRUE)
}
#' @exportS3Method bal.tab data.frame
#' @rdname bal.tab.formula
bal.tab.data.frame <- function(x, treat,
                               stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, s.weights = NULL, abs = FALSE, subset = NULL, quick = TRUE,
                               subclass = NULL, match.strata = NULL, method, estimand = NULL, focal = NULL, ...) {
  
  args <- try_chk(c(as.list(environment()), list(...))[-1L])
  
  #Adjustments to arguments
  
  args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
  args[lengths(args) == 0L & names(args) %nin% names(match.call())[-1L]] <- NULL
  
  #Initializing variables
  X <- do.call("x2base.data.frame", c(x = list(x), args), quote = TRUE)
  
  args[names(X)] <- NULL
  
  X <- .assign_X_class(X)
  
  do.call("base.bal.tab", c(list(X), args),
          quote = TRUE)
}

#' @exportS3Method bal.tab matrix
#' @rdname bal.tab.formula
bal.tab.matrix <- bal.tab.data.frame

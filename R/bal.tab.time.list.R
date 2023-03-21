#' @title Balance Statistics for Longitudinal Datasets
#' 
#' @description
#' Generates balance statistics for data coming from a longitudinal treatment scenario. The primary input is in the form of a list of formulas or `data.frame`s contain the covariates at each time point. `bal.tab()` automatically classifies this list as either a `data.frame.list` or `formula.list`, respectively.
#' 
#' @inheritParams bal.tab
#' @param x either a list of data frames containing all the covariates to be assessed at each time point or a list of formulas with the treatment for each time period on the left and the covariates for which balance is to be displayed on the right. Covariates to be assessed at multiple points must be included in the entries for each time point. Data must be in the "wide" format, with one row per unit. If a formula list is supplied, an argument to `data` is required unless all objects in the formulas exist in the environment.
#' @param treat.list treatment status for each unit at each time point. This can be specified as a list or data frame of vectors, each of which contains the treatment status of each individual at each time point, or a list or vector of the names of variables in `data` that contain treatment at each time point. Required for the `data.frame.list` method.
#' @param s.d.denom `character`; how the denominator for standardized mean differences should be calculated, if requested. See [col_w_smd()] for allowable options. Abbreviations allowed. It is recommended not to set this argument for longitudinal treatments.
#' 
#' @returns 
#' An object of class `bal.tab.msm` containing balance summaries at each time point. Each balance summary is its own `bal.tab` object. See [`bal.tab.msm()`][class-bal.tab.msm] for more details.
#' 
#' See [`bal.tab() base methods()`][bal.tab.formula] for more detailed information on the value of the `bal.tab` objects produced for each time point.
#' 
#' @details 
#' `bal.tab.formula.list()` and `bal.tab.data.frame.list()` generate a list of balance summaries for each time point based on the treatments and covariates provided. All data must be in the "wide" format, with exactly one row per unit and columns representing variables at different time points. See the \pkgfun{WeightIt}{weightitMSM} documentation for an example of how to transform long data into wide data using [reshape()]. 
#' 
#' Multiple sets of weights can be supplied simultaneously by including entering a data frame or a character vector containing the names of weight variables found in `data` or a list thereof. When only one set of weights is supplied, the output for the adjusted group will simply be called `"Adj"`, but otherwise will be named after each corresponding set of weights. Specifying multiple sets of weights will also add components to other outputs of `bal.tab()`.
#' 
#' @seealso 
#' * [bal.tab()] for details of calculations.
#' * [`bal.tab.msm()`][class-bal.tab.msm] for output and related options.
#' * [`bal.tab.cluster()`][class-bal.tab.cluster] for more information on clustered data.
#' * [`bal.tab.imp()`][class-bal.tab.imp] for more information on multiply imputed data.
#' * [`bal.tab.multi()`][class-bal.tab.multi] for more information on multi-category treatments.
#' 
#' @examplesIf requireNamespace("twang", quietly = TRUE)
#' data("iptwExWide", package = "twang")
#' library("cobalt")
#' 
#' ## Estimating longitudinal propensity scores and weights
#' ps1 <- glm(tx1 ~ age + gender + use0,
#'            data = iptwExWide, 
#'            family = "binomial")$fitted.values
#' w1 <- ifelse(iptwExWide$tx1 == 1, 1/ps1, 1/(1-ps1))
#' ps2 <- glm(tx2 ~ age + gender + use0 + tx1 + use1,
#'            data = iptwExWide, 
#'            family = "binomial")$fitted.values
#' w2 <- ifelse(iptwExWide$tx2 == 1, 1/ps2, 1/(1-ps2))
#' ps3 <- glm(tx3 ~ age + gender + use0 + tx1 + use1 + tx2 + use2,
#'            data = iptwExWide, 
#'            family = "binomial")$fitted.values
#' w3 <- ifelse(iptwExWide$tx3 == 1, 1/ps3, 1/(1-ps3))
#' 
#' w <- w1*w2*w3
#' 
#' # Formula interface plus addl:
#' bal.tab(list(tx1 ~ use0 + gender,
#'              tx2 ~ use0 + gender + use1 + tx1,
#'              tx3 ~ use0 + gender + use1 + tx1 + use2 + tx2),
#'         data = iptwExWide, 
#'         weights = w,
#'         distance = list(~ps1, ~ps2, ~ps3),
#'         addl = ~age*gender,
#'         un = TRUE)
#' 
#' # data frame interface:
#' bal.tab(list(iptwExWide[c("use0", "gender")],
#'              iptwExWide[c("use0", "gender", "use1", "tx1")],
#'              iptwExWide[c("use0", "gender", "use1", "tx1", "use2", "tx2")]),
#'         treat.list = iptwExWide[c("tx1", "tx2", "tx3")], 
#'         weights = w,
#'         distance = list(~ps1, ~ps2, ~ps3),
#'         un = TRUE)

#' @exportS3Method bal.tab formula.list
#' @name bal.tab.time.list
bal.tab.formula.list <- function(x,
                                 stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, s.weights = NULL, abs = FALSE, subset = NULL, quick = TRUE,
                                 ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) .err(conditionMessage(e)))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    args[vapply(args, is_null, logical(1L)) & names(args) %nin% names(match.call())[-1]] <- NULL
    
    X <- do.call("x2base.formula.list", c(list(formula.list = x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    return(out)
}

#' @exportS3Method bal.tab data.frame.list
#' @rdname bal.tab.time.list
bal.tab.data.frame.list <- function(x, treat.list, 
                                    stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, s.weights = NULL, abs = FALSE, subset = NULL, quick = TRUE,
                                    ...) {
    
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) .err(conditionMessage(e)))
    
    #Adjustments to arguments
    
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    args[vapply(args, is_null, logical(1L)) & names(args) %nin% names(match.call())[-1]] <- NULL
    
    X <- do.call("x2base.data.frame.list", c(list(covs.list = x), args), quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    return(out)
}

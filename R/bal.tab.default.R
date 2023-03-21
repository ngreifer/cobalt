#' @title Balance Statistics for Other Objects
#' 
#' @description Generates balance statistics using an object for which there is not a defined method. 
#' 
#' @inheritParams bal.tab
#' @param x An object containing information about conditioning. See Details.
#' @param ... other arguments that would be passed to [bal.tab.formula()], [bal.tab.data.frame()], or [bal.tab.time.list()]. See Details.
#' 
#' @returns 
#' For point treatments, if clusters and imputations are not specified, an object of class `"bal.tab"` containing balance summaries for the specified treatment and covariates. See [bal.tab()] for details.
#' 
#' If clusters are specified, an object of class `"bal.tab.cluster"` containing balance summaries within each cluster and a summary of balance across clusters. See [`bal.tab.cluster()`][class-bal.tab.cluster] for details.
#' 
#' If imputations are specified, an object of class `"bal.tab.imp"` containing balance summaries for each imputation and a summary of balance across imputations, just as with clusters. See [`bal.tab.imp()`][class-bal.tab.imp] for details.
#' 
#' If multi-category treatments are used, an object of class `"bal.tab.multi"` containing balance summaries for each pairwise treatment comparison and a summary of balance across pairwise comparisons. See [`bal.tab.multi()`][class-bal.tab.multi] for details.
#' 
#' If longitudinal treatments are used, an object of class `"bal.tab.msm"` containing balance summaries at each time point. Each balance summary is its own `bal.tab` object. See [`bal.tab.msm()`][class-bal.tab.msm] for more details.
#' 
#' @details 
#' `bal.tab.default()` processes its input and attempt to extract enough information from it to display covariate balance for `x`. The purpose of this method is to allow users who have created their own objects containing conditioning information (i.e., weights, subclasses, treatments, covariates, etc.) to access the capabilities of `bal.tab()` without having a special method written for them. By including the correct items in `x`, `bal.tab.default()` can present balance tables as if the input was the output of one of the specifically supported packages (e.g., \pkg{MatchIt}, \pkg{twang}, etc.).
#' 
#' The function will search `x` for the following named items and attempt to process them:
#'     \describe{
#'         \item{`treat`}{A vector (`numeric`, `character`, `factor`) containing the values of the treatment for each unit or the name of the column in `data` containing them. Essentially the same input to `treat` in [bal.tab.data.frame()].
#'         }
#'         \item{`treat.list`}{A list of vectors (`numeric`, `character`, `factor`) containing, for each time point, the values of the treatment for each unit or the name of the column in `data` containing them. Essentially the same input to `treat.list` in [bal.tab.time.list()].
#'         }
#'         \item{`covs`}{A `data.frame` containing the values of the covariates for each unit. Essentially the same input to `covs` in [bal.tab.data.frame()].
#'         }
#'         \item{`covs.list`}{A list of `data.frame`s containing, for each time point, the values of the covariates for each unit. Essentially the same input to `covs.list` in [bal.tab.time.list()].
#'         }
#'         \item{`formula`}{A `formula` with the treatment variable as the response and the covariates for which balance is to be assessed as the terms. Essentially the same input to `formula` in [bal.tab.formula()].
#'         }
#'         \item{`formula.list`}{A list of `formula`s with, for each time point, the treatment variable as the response and the covariates for which balance is to be assessed as the terms. Essentially the same input to `formula.list` in [bal.tab.time.list()].
#'         }
#'         \item{`data`}{A `data.frame` containing variables with the names used in other arguments and components (e.g., `formula`, `weights`, etc.). Essentially the same input to `data` in [bal.tab.formula()], [bal.tab.data.frame()], or [bal.tab.time.list()].
#'         }
#'         \item{`weights`}{A vector, list, or `data.frame` containing weights for each unit or a string containing the names of the weights variables in `data`. Essentially the same input to `weights` in [bal.tab.data.frame()] or [bal.tab.time.list()].
#'         }
#'         \item{`distance`}{
#'             A vector, formula, or data frame containing distance values (e.g., propensity scores) or a character vector containing their names. If a formula or variable names are specified, `bal.tab()` will look in the argument to `data`, if specified. Essentially the same input to `distance` in [bal.tab.data.frame()].
#'         }
#'         \item{`formula.list`}{A list of vectors or `data.frame`s containing, for each time point, distance values (e.g., propensity scores) for each unit or a string containing the name of the distance variable in `data`. Essentially the same input to `distance.list` in [bal.tab.time.list()].
#'         }
#'         \item{`subclass`}{A vector containing subclass membership for each unit or a string containing the name of the subclass variable in `data`. Essentially the same input to `subclass` in [bal.tab.data.frame()].
#'         }
#'         \item{`match.strata`}{A vector containing matching stratum membership for each unit or a string containing the name of the matching stratum variable in `data`. Essentially the same input to `match.strata` in [bal.tab.data.frame()].
#'         }
#'         \item{`estimand`}{A `character` vector; whether the desired estimand is the "ATT", "ATC", or "ATE" for each set of weights. Essentially the same input to `estimand` in [bal.tab.data.frame()].
#'         }
#'         \item{`s.weights`}{A vector containing sampling weights for each unit or a string containing the name of the sampling weight variable in `data`. Essentially the same input to `s.weights` in [bal.tab.data.frame()] or [bal.tab.time.list()].
#'         }
#'         \item{`focal`}{The name of the focal treatment when multi-category treatments are used. Essentially the same input to `focal` in [bal.tab.data.frame()].
#'         }
#'         \item{`call`}{A `call` object containing the function call, usually generated by using [match.call()] inside the function that created `x`.
#'         }
#'     }
#' Any of these items can also be supplied directly to `bal.tab.default`, e.g., `bal.tab.default(x, formula = treat ~ x1 + x2)`. If supplied, it will override the object with the same role in `x`. In addition, any arguments to [bal.tab.formula()], [bal.tab.data.frame()], and [bal.tab.time.list()] are allowed and perform the same function.
#' 
#' At least some inputs containing information to create the treatment and covariates are required (e.g., `formula` and `data` or `covs` and `treat`). All other arguments are optional and have the same defaults as those in [bal.tab.data.frame()] or [bal.tab.time.list()]. If `treat.list`, `covs.list`, or `formula.list` are supplied in `x` or as an argument to `bal.tab.default()`, the function will proceed considering a longitudinal treatment. Otherwise, it will proceed considering a point treatment.
#' 
#' `bal.tab.default()`, like other `bal.tab()` methods, is just a shortcut to supply arguments to `bal.tab.data.frame()` or `bal.tab.time.list()`. Therefore, any matters regarding argument priority or function are described in the documentation for these methods.
#' 
#' @seealso 
#' * [bal.tab.formula()] and [bal.tab.time.list()] for additional arguments to be supplied.
#' * [bal.tab()] for output and details of calculations.
#' * [`bal.tab.cluster()`][class-bal.tab.cluster] for more information on clustered data.
#' * [`bal.tab.imp()`][class-bal.tab.imp] for more information on multiply imputed data.
#' * [`bal.tab.multi()`][class-bal.tab.multi] for more information on multi-category treatments.
#' 
#' @examples
#' data("lalonde", package = "cobalt")
#' covs <- subset(lalonde,  select = -c(treat, re78))
#' 
#' ##Writing a function the produces output for direct
#' ##use in bal.tab.default
#' 
#' ate.weights <- function(treat, covs) {
#'     data <- data.frame(treat, covs)
#'     formula <- formula(data)
#'     ps <- glm(formula, data = data, 
#'               family = "binomial")$fitted.values
#'     weights <- treat/ps + (1-treat)/(1-ps)
#'     call <- match.call()
#'     out <- list(treat = treat,
#'                 covs = covs,
#'                 distance = ps,
#'                 weights = weights,
#'                 estimand = "ATE",
#'                 call = call)
#'     return(out)
#' }
#' 
#' out <- ate.weights(lalonde$treat, covs)
#' 
#' bal.tab(out, un = TRUE)

#' @exportS3Method bal.tab default
bal.tab.default <- function(x, stats, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, thresholds = NULL, weights = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, s.weights = NULL, abs = FALSE, subset = NULL, quick = TRUE,
                            ...) {
    tryCatch(args <- c(as.list(environment()), list(...))[-1], error = function(e) .err(conditionMessage(e)))
    
    #Adjustments to arguments
    args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
    args[vapply(args, is_null, logical(1L)) & names(args) %nin% names(match.call())[-1]] <- NULL
    
    #Initializing variables
    X <- do.call("x2base.default", c(list(obj = x), args),
                 quote = TRUE)
    
    args[names(args) %in% names(X)] <- NULL
    
    X <- assign.X.class(X)
    
    out <- do.call("base.bal.tab", c(list(X), args),
                   quote = TRUE)
    
    return(out)
}

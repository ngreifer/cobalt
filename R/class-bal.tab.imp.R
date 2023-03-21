#' Using `bal.tab()` with Multiply Imputed Data
#' @name class-bal.tab.imp
#' 
#' @description
#' When using [bal.tab()] with multiply imputed data, the output will be different from the case with a single data set. Multiply imputed data can be used with all `bal.tab()` methods, and the `mimids` and `wimids` methods for \pkg{MatchThem} objects automatically incorporate multiply imputed data. This page outlines the outputs and options available with multiply imputed data.
#'     
#' There are two main components of the output of `bal.tab()` with multiply imputed data: the within-imputation balance summaries and the across-imputation balance summary. The within-imputation balance summaries display balance for units within each imputed data set separately. In general, this will not be very useful because interest rarely lies in the qualities of any individual imputed data set.
#'     
#' The across-imputation balance summary pools information across the within-imputation balance summaries to simplify balance assessment. It provides the average, smallest, and largest balance statistic for each covariate across all imputations. This allows you to see how bad the worst imbalance is and what balance looks like on average across the imputations. The summary behaves differently depending on whether `abs` is specified as `TRUE` or `FALSE`. When `abs = TRUE`, the across-imputation balance summary will display the mean absolute balance statistics and the maximum absolute balance statistics. When `abs = FALSE`, the across-imputation balance summary will display the minimum, mean, and maximum of the balance statistic in its original form.
#' 
#' @section Allowable arguments:
#' 
#' There are four arguments for each `bal.tab()` method that can handle multiply imputed data: `imp`, `which.imp`, `imp.summary`, and `imp.fun`.
#' 
#' \describe{
#'     \item{`imp`}{A vector of imputation membership. This can be factor, character, or numeric vector. This argument is required to let `bal.tab()` know that the data is multiply imputed unless \pkg{MatchThem} objects are used. If a `data` argument is specified, this can also be the name of a variable in `data` that contains imputation membership. If the `data` argument is a `mids` object, the output of a call to `mice()`, `imp` does not need to be specified and will automatically be extracted from the `mids` object.}
#'     \item{`which.imp`}{This is a display option that does not affect computation. If `.all`, all imputations in `imp` will be displayed. If `.none` (the default), no imputations will be displayed. Otherwise, can be a vector of imputation indices for which to display balance.}
#'     \item{`imp.summary`}{This is a display option that does not affect computation. If `TRUE`, the balance summary across imputations will be displayed. The default is `TRUE`, and if `which.imp` is `.none`, it will automatically be set to `TRUE`.}
#'     \item{`imp.fun`}{This is a display option that does not affect computation. Can be "min", "mean", or "max" and corresponds to which function is used in the across-imputation summary to combine results across imputations. For example, if `imp.fun = "mean"` the mean balance statistic across imputations will be displayed. The default when `abs = FALSE` in the `bal.tab()` call is to display all three. The default when `abs = FALSE` in the `bal.tab()` call is to display just the mean and max balance statistic.
#'     }
#' }
#' 
#' @section Output:
#' 
#' The output is a `bal.tab.imp` object, which inherits from `bal.tab`. It has the following elements:
#'         
#' * `Imputation.Balance`: For each imputation, a regular `bal.tab` object containing a balance table, a sample size summary, and other balance assessment tools, depending on which options are specified.
#' * `Balance.Across.Imputations`: The balance summary across imputations. This will include the combination of each balance statistic for each covariate across all imputations according to the value of `imp.fun`.
#' * `Observations`: A table of sample sizes or effective sample sizes averaged across imputations before and after adjustment.
#'     
#' As with other methods, multiple weights can be specified, and values for all weights will appear in all tables.
#' 
#' 
#' @seealso
#' * [bal.tab()]
#' * [bal.tab.data.frame()]
#' * [print.bal.tab()]
#' * `vignette("segmented-data")` for examples
#' 
NULL

base.bal.tab.imp <- function(X,
                             which.imp = NA,
                             imp.summary = getOption("cobalt_imp.summary"),
                             imp.fun = getOption("cobalt_imp.fun", NULL),
                             ...) {
    A <- list(...)
    
    #Preparations
    
    if (is_null(A[["quick"]])) A[["quick"]] <- TRUE
    
    imp <- factor(X$imp)
    
    if (missing(which.imp)) {
        which.imp <- NA
    }
    
    if (is_null(imp.summary)) {
        imp.summary <- is_not_null(which.imp) && 
            (anyNA(which.imp) || !is.numeric(which.imp) || 
                 (is.numeric(which.imp) && !any(which.imp %in% seq_len(nlevels(imp)))))
    }
    
    all.agg.funs <- c("min", "mean", "max")
    agg.fun <- tolower(as.character(if_null_then(imp.fun, A[["agg.fun"]], all.agg.funs)))
    agg.fun <- match_arg(agg.fun, all.agg.funs, several.ok = TRUE)
    
    X$covs <- do.call("get.C2", c(X, A[names(A) %nin% names(X)]), quote = TRUE)
    
    #Setup output object
    out.names <- c("Imputation.Balance", 
                   "Balance.Across.Imputations", 
                   "Observations", 
                   "call")
    out <- make_list(out.names)
    
    #Get list of bal.tabs for each imputation
    
    out[["Imputation.Balance"]] <- lapply(levels(imp), function(i) {
        X_i <- assign.X.class(subset_X(X, imp == i)) 
        X_i$call <- NULL
        tryCatch({
            do.call("base.bal.tab", c(list(X_i), A[names(A) %nin% names(X_i)]), quote = TRUE)
        },
        error = function(e) {
            .err(sprintf("in imputation %s: %s", i, conditionMessage(e)))
        })
    })
    
    names(out[["Imputation.Balance"]]) <- levels(imp)
    
    #Create summary of lists
    
    if (imp.summary || !A$quick) {
        out[["Balance.Across.Imputations"]] <- balance.summary(out[["Imputation.Balance"]], 
                                                               agg.funs = agg.fun)
        
        if (length(agg.fun) == 1) {
            out <- c(out, threshold.summary(compute = attr(out[["Imputation.Balance"]][[1]][["Balance"]], "compute"),
                                            thresholds = attr(out[["Imputation.Balance"]][[1]][["Balance"]], "thresholds"),
                                            no.adj = !attr(out[["Imputation.Balance"]][[1]], "print.options")$disp.adj,
                                            balance.table = out[["Balance.Across.Imputations"]],
                                            weight.names = attr(out[["Imputation.Balance"]][[1]], "print.options")$weight.names,
                                            agg.fun = agg.fun))
        }
        
        observations <- lapply(out[["Imputation.Balance"]], function(x) x[["Observations"]])
        
        out[["Observations"]] <- samplesize.across.imps(observations)
    }
    
    out[["call"]] <- X$call
    attr(out, "print.options") <- c(attr(out[["Imputation.Balance"]][[1]], "print.options"),
                                    list(which.imp = which.imp,
                                         imp.summary = imp.summary,
                                         imp.fun = agg.fun))
    class(out) <- c("bal.tab.imp", "bal.tab")
    
    return(out)
}

#' Using `bal.tab()` with Longitudinal Treatments
#' @name class-bal.tab.msm
#' 
#' @description
#'     When using [bal.tab()] with longitudinal treatments, the output will be different from the case with point treatments, and there are some options that are common across all `bal.tab()` methods for dealing with longitudinal data. This page outlines the outputs and options in this case.
#'     
#'     There are two main components of the output of `bal.tab()` with longitudinal treatments: the time-point-specific balance summary and across-time-points balance summary. The time-point-specific balance summaries are standard point treatment balance summaries at each time point.
#'     
#'     The across-time-points balance summary is, for each variable, the greatest imbalance across all time-point-specific balance summaries. If the greatest observed imbalance is tolerable, then all other imbalances for that variable will be tolerable too, so focusing on reducing the greatest imbalance is sufficient for reducing imbalance overall. The balance summary will not be computed if multi-category treatments or multiply imputed data are used.
#' 
#' @section Allowable arguments:
#' 
#' There are two additional arguments for each `bal.tab()` method that can handle longitudinal treatments: `which.time` and `msm.summary`.
#' 
#' \describe{
#'     \item{`which.time`}{This is a display option that does not affect computation. If `.all` (the default), all time points will be displayed. If `.none`, no time points will be displayed. Otherwise, can be a vector of treatment names or indices for which to display balance.}
#'     \item{`msm.summary`}{This is a display option that does not affect computation. If `TRUE`, the balance summary across time points will be displayed. The default is `TRUE`, and if `which.time` is `.none`, it will automatically be set to `TRUE`.}
#' }
#' 
#' @section Output: 
#' The output is a `bal.tab.msm` object, which inherits from `bal.tab`. It has the following elements:
#'         
#' * `Time.Balance`: For each time point, a regular `bal.tab` object containing a balance table, a sample size summary, and other balance assessment tools, depending on which options are specified.
#' * `Balance.Across.Times`: The balance summary across time points. This will include the maximum balance statistic(s) for each covariate across all time points.
#' * `Observations`: A table of sample sizes or effective sample sizes for each time point before and after adjustment.
#'     
#' As with other methods, multiple weights can be specified, and values for all weights will appear in all tables.
#' 
#' @note The balance tables presented here are not the same as those recommended by Jackson (2016) and computed in his R package, \CRANpkg{confoundr}, as these do not take into account treatment history. The balance statistics presented here should be used with caution and may not reflect balance in an accurate way.
#' 
#' @references 
#' Jackson, J. W. (2016). Diagnostics for Confounding of Time-varying and Other Joint Exposures: *Epidemiology*, 27(6), 859–869. \doi{10.1097/EDE.0000000000000547}
#' 
#' @seealso
#' * [bal.tab()]
#' * [bal.tab.time.list()]
#' * [print.bal.tab()]
#' * `vignette("longitudinal-treat")` for examples
#' 
NULL

base.bal.tab.msm <- function(X,
                             which.time,
                             msm.summary = getOption("cobalt_msm.summary"),
                             ...) {
    A <- list(...)
    
    #Preparations
    
    if (is_null(A[["quick"]])) A[["quick"]] <- TRUE
    
    treat.types <- vapply(X$treat.list, get.treat.type, character(1L))
    
    if (missing(which.time)) {
        if (all_the_same(treat.types) && "multinomial" %nin% treat.types && is_null(X$imp)) which.time <- NA
        else which.time <- NULL
    }
    
    if (is_null(msm.summary)) {
        msm.summary <- is_not_null(which.time) && (anyNA(which.time) ||
                                                       !(is.character(which.time) || is.numeric(which.time)) ||
                                                       (is.numeric(which.time) && !any(which.time %in% seq_along(X$treat.list))) ||
                                                       (is.character(which.time) && !any(which.time %in% names(X$treat.list))))
    }
    
    #Setup output object
    out <- list()
    
    #Get list of bal.tabs for each time period
    out[["Time.Balance"]] <- lapply(seq_along(X$covs.list), function(ti) {
        X_ti <- X

        X_ti$covs <- X_ti$covs.list[[ti]]
        X_ti$treat <- X_ti$treat.list[[ti]]
        X_ti$addl <- X_ti$addl.list[[ti]]
        X_ti$distance <- X_ti$distance.list[[ti]]
        
        X_ti[c("covs.list", "treat.list", "addl.list", "distance.list", "call")] <- NULL

        X_ti <- .assign_X_class(X_ti)
        
        X_ti$s.d.denom <- {
            if (attr(X_ti, "X.class") == "cont") "all"
            else "pooled"
        }
        
        do.call("base.bal.tab", c(list(X_ti), A[names(A) %nin% names(X_ti)]), quote = TRUE)
    })
    
    names(out[["Time.Balance"]]) <- {
        if (length(names(X$treat.list)) == length(X$treat.list)) names(X$treat.list)
        else seq_along(X$treat.list)
    }
    
    if ((!A$quick || msm.summary) && is_null(X$imp) && all_the_same(treat.types) &&
        !any(treat.types == "multinomial")) {
        out[["Balance.Across.Times"]] <- balance.summary(out[["Time.Balance"]],
                                                         agg.funs = "max",
                                                         include.times = TRUE)
        
        out <- c(out,
                 threshold.summary(compute = attr(out[["Time.Balance"]][[1]][["Balance"]], "compute"),
                                   thresholds = attr(out[["Time.Balance"]][[1]][["Balance"]], "thresholds"),
                                   no.adj = !attr(out[["Time.Balance"]][[1]], "print.options")$disp.adj,
                                   balance.table = out[["Balance.Across.Times"]],
                                   weight.names = attr(out[["Time.Balance"]][[1]], "print.options")$weight.names,
                                   agg.fun = "max"))
        
        out[["Observations"]] <- grab(out[["Time.Balance"]], "Observations")
    }
    
    out[["call"]] <- X$call
    
    attr(out, "print.options") <- c(attr(out[["Time.Balance"]][[1]], "print.options"),
                                    list(which.time = which.time,
                                         msm.summary = msm.summary))
    
    class(out) <- c("bal.tab.msm", "bal.tab")
    
    out
}

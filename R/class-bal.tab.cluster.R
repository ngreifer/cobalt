#' Using `bal.tab()` with Clustered Data
#' @name class-bal.tab.cluster
#' 
#' @description
#' When using [bal.tab()] with clustered data, the output will be different from the case with single-level data, and there are some options that are common across all `bal.tab()` methods. This page outlines the outputs and options in this case.
#'     
#' There are two main components of the output of `bal.tab()` with clustered data: the within-cluster balance summaries and the across-cluster balance summary. The within-cluster balance summaries display balance for units within each cluster separately.
#'     
#' The across-cluster balance summary pools information across the within-cluster balance summaries to simplify balance assessment. It provides a combination (e.g., mean or maximum) of each balance statistic for each covariate across all clusters. This allows you to see how bad the worst imbalance is and what balance looks like on average. The balance summary will not be computed if longitudinal treatments, multi-category treatments, or multiply imputed data are used.
#' 
#' @section Allowable arguments:
#' 
#' There are four arguments for each `bal.tab()` method that can handle clustered data: `cluster`, `which.cluster`, `cluster.summary`, and `cluster.fun`.
#' 
#' \describe{
#'     \item{`cluster`}{A vector of cluster membership. This can be factor, character, or numeric vector. This argument is required to let `bal.tab()` know that the data is clustered. If a `data` argument is specified, this can also be the name of a variable in `data` that contains cluster membership.}
#'     \item{`which.cluster`}{This is a display option that does not affect computation. If `.all` (the default), all clusters in `cluster` will be displayed. If `.none`, no clusters will be displayed. Otherwise, can be a vector of cluster names or numerical indices for which to display balance. Indices correspond to the alphabetical order of cluster names (or the order of cluster levels if a factor).}
#'     \item{`cluster.summary`}{This is a display option that does not affect computation. If `TRUE`, the balance summary across clusters will be displayed. The default is `TRUE`, and if `which.cluster` is `.none`, it will automatically be set to `TRUE`.}
#'     \item{`cluster.fun`}{This is a display option that does not affect computation. Can be "min", "mean", or "max" and corresponds to which function is used in the across-cluster summary to combine results across clusters. For example, if `cluster.fun = "mean"` the mean balance statistic across clusters will be displayed. The default when `abs = FALSE` in the `bal.tab()` call is to display all three. The default when `abs = FALSE` in the `bal.tab()` call is to display just the mean and max balance statistic.
#'     }
#' }
#' 
#' @section Output:
#' 
#' The output is a `bal.tab.cluster` object, which inherits from `bal.tab`. It has the following elements:
#'         
#' * `Cluster.Balance`: For each cluster, a regular `bal.tab` object containing a balance table, a sample size summary, and other balance assessment tools, depending on which options are specified.
#' * `Cluster.Summary`: The balance summary across clusters. This will include the combination of each balance statistic for each covariate across all clusters according to the value of `cluster.fun`.
#' * `Observations`: A table of sample sizes or effective sample sizes for each cluster before and after adjustment.
#'     
#' As with other methods, multiple weights can be specified, and values for all weights will appear in all tables.
#' 
#' @seealso
#' * [bal.tab()]
#' * [bal.tab.data.frame()]
#' * [print.bal.tab()]
#' * `vignette("segmented-data")` for examples
#' 
NULL

base.bal.tab.cluster <- function(X,
                                 which.cluster,
                                 cluster.summary = getOption("cobalt_cluster.summary"),
                                 cluster.fun = getOption("cobalt_cluster.fun", NULL),
                                 ...) {
    A <- list(...)
    
    #Preparations
    
    if (is_null(A[["quick"]])) A[["quick"]] <- TRUE
    if (is_null(A[["abs"]])) A[["abs"]] <- FALSE
    
    X$cluster <- factor(X$cluster)
    
    .cluster_check(X$cluster, X$treat)
    
    #Process cluster.summary
    if (missing(which.cluster)) {
        which.cluster <- NULL
    }
    
    if (is_null(cluster.summary)) {
        cluster.summary <- is_not_null(which.cluster) && anyNA(which.cluster)
    }
    
    all.agg.funs <- c("min", "mean", "max")
    agg.fun <- tolower(as.character(if_null_then(cluster.fun, A[["agg.fun"]], all.agg.funs)))
    agg.fun <- match_arg(agg.fun, all.agg.funs, several.ok = TRUE)
    
    X$covs <- do.call(".get_C2", c(X, A[names(A) %nin% names(X)]), quote = TRUE)
    
    var_types <- attr(X$covs, "var_types")
    
    if (get.treat.type(X$treat) != "continuous") {
        if (is_null(A$continuous)) A$continuous <- getOption("cobalt_continuous", "std")
        if (is_null(A$binary)) A$binary <- getOption("cobalt_binary", "raw")
    }
    else {
        if (is_null(A$continuous)) A$continuous <- getOption("cobalt_continuous", "std")
        if (is_null(A$binary)) A$binary <- getOption("cobalt_binary", "std")
    }
    
    if (get.treat.type(X$treat) != "continuous" &&
        "mean.diffs" %in% X$stats &&
        ((A$binary == "std" && any(var_types == "Binary")) ||
         (A$continuous == "std" && any(var_types != "Binary")))) {
        X$s.d.denom <- .get_s.d.denom(X$s.d.denom,
                                      estimand = X$estimand,
                                      weights = X$weights, 
                                      subclass = X$subclass,
                                      treat = X$treat,
                                      focal = X$focal)
    }
    else if (get.treat.type(X$treat) == "continuous" &&
             any(c("correlations", "spearman.correlations") %in% X$stats) &&
             ((A$binary == "std" && any(var_types == "Binary")) ||
              (A$continuous == "std" && any(var_types != "Binary")))) {
        X$s.d.denom <- .get_s.d.denom.cont(X$s.d.denom,
                                           weights = X$weights,
                                           subclass = X$subclass)
    }
    
    #Setup output object
    out <- list()
    
    #Get list of bal.tabs for each imputation
    out[["Cluster.Balance"]] <- lapply(levels(X$cluster), function(cl) {
        X_cl <- .assign_X_class(subset_X(X, X$cluster == cl)) 
        X_cl$call <- NULL
        
        tryCatch({
            do.call("base.bal.tab", c(list(X_cl), A[setdiff(names(A), names(X_cl))]), quote = TRUE)
        },
        error = function(e) {
            .err(sprintf("in cluster %s: %s", add_quotes(cl), conditionMessage(e)))
        })
    })
    
    names(out[["Cluster.Balance"]]) <- levels(X$cluster)
    
    #Create summary of lists
    
    if ((cluster.summary || !A$quick) && is_null(X$covs.list) && get.treat.type(X$treat) != "multinomial" && is_null(X$imp)) {
        out[["Balance.Across.Clusters"]] <- balance.summary(out[["Cluster.Balance"]], 
                                                            agg.funs = if_null_then(agg.fun, c("min", "mean", "max")))
        
        if (length(agg.fun) == 1) {
            out <- c(out, threshold.summary(compute = attr(out[["Cluster.Balance"]][[1]][["Balance"]], "compute"),
                                            thresholds = attr(out[["Cluster.Balance"]][[1]][["Balance"]], "thresholds"),
                                            no.adj = !attr(out[["Cluster.Balance"]][[1]], "print.options")$disp.adj,
                                            balance.table = out[["Balance.Across.Clusters"]],
                                            weight.names = attr(out[["Cluster.Balance"]][[1]], "print.options")$weight.names,
                                            agg.fun = agg.fun))
        }
        
        observations <- grab(out[["Cluster.Balance"]], "Observations")
        
        out[["Observations"]] <- samplesize.across.clusters(observations)
    }
    
    
    out[["call"]] <- X$call
    
    attr(out, "print.options") <- c(attr(out[["Cluster.Balance"]][[1]], "print.options"),
                                    list(which.cluster = which.cluster,
                                         cluster.summary = cluster.summary,
                                         cluster.fun = agg.fun))
    
    class(out) <- c("bal.tab.cluster", "bal.tab")
    
    out
}

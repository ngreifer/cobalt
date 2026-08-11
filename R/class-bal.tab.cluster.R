#' Using `bal.tab()` with Clustered Data
#' @name class-bal.tab.cluster
#' 
#' @description
#' When using [bal.tab()] with clustered or subgrouped data, the output will be different from the case with single-level data, and there are some options that are common across all `bal.tab()` methods. This page outlines the outputs and options in this case.
#'     
#' There are two main components of the output of `bal.tab()` with clustered data: the within-cluster balance summaries and the across-cluster balance summary. The within-cluster balance summaries display balance for units within each cluster separately.
#'     
#' The across-cluster balance summary pools information across the within-cluster balance summaries to simplify balance assessment. It provides a combination (e.g., mean or maximum) of each balance statistic for each covariate across all clusters. This allows you to see how bad the worst imbalance is and what balance looks like on average. The balance summary will not be computed if longitudinal treatments, multi-category treatments, or multiply imputed data are used.
#' 
#' In order to use the `thresholds` argument with `bal.tab()` with clustered data and the balance summary across clustered displayed, `cluster.fun` must be supplied and set to a single string, which is not the default.
#' 
#' @section Allowable arguments:
#' 
#' There are four arguments for each `bal.tab()` method that can handle clustered data: `cluster`, `which.cluster`, `cluster.summary`, and `cluster.fun`.
#' 
#' \describe{
#'     \item{`cluster`}{A vector of cluster membership. This can be factor, character, or numeric vector. This argument is required to let `bal.tab()` know that the data is clustered. If a `data` argument is specified, this can also be the name of a variable in `data` that contains cluster membership.}
#'     \item{`which.cluster`}{This is a display option that does not affect computation. If `.all` (the default), all clusters in `cluster` will be displayed. If `.none`, no clusters will be displayed. Otherwise, can be a vector of cluster names or numerical indices for which to display balance. Indices correspond to the alphabetical order of cluster names (or the order of cluster levels if a factor).}
#'     \item{`cluster.summary`}{This is a display option that does not affect computation. If `TRUE`, the balance summary across clusters will be displayed. The default is `TRUE`, and if `which.cluster` is `.none`, it will automatically be set to `TRUE`.}
#'     \item{`cluster.fun`}{This is a display option that does not affect computation. Can be "min", "mean", or "max" and corresponds to which function is used in the across-cluster summary to combine results across clusters. For example, if `cluster.fun = "mean"` the mean balance statistic across clusters will be displayed. The default when `abs = FALSE` in the `bal.tab()` call is to display all three. The default when `abs = TRUE` in the `bal.tab()` call is to display just the mean and maximum absolute balance statistic.
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
  
  X[["cluster"]] <- factor(X[["cluster"]])
  
  #A longitudinal `X` carries `treat.list` rather than `treat`, and `.cluster_check()`
  #accepts either one treatment or a list of them. This used to read `X$treat`, which
  #reached the list only by partial matching.
  .cluster_check(X[["cluster"]], X[["treat"]] %or% X[["treat.list"]])
  
  #Process cluster.summary
  if (missing(which.cluster)) {
    which.cluster <- NULL
  }
  
  if (is_null(cluster.summary)) {
    cluster.summary <- is_not_null(which.cluster) && anyNA(which.cluster)
  }
  
  all.agg.funs <- c("min", "mean", "max")
  agg.fun <- as.character(cluster.fun %or% A[["agg.fun"]] %or% all.agg.funs)
  agg.fun <- arg::match_arg(agg.fun, all.agg.funs, several.ok = TRUE)
  
  #With longitudinal treatments, `X` has `covs.list`/`treat.list` rather than
  #`covs`/`treat`, and `base.bal.tab.msm()` derives the covariates, treatment, and
  #`s.d.denom` for each time point itself, so this preparation is skipped.
  if (is_null(X[["covs.list"]])) {
    prep <- .bal.tab_prepare(X, A)
    X <- prep[["X"]]
    A <- prep[["A"]]

    #A wrapper keeps the user's `s.d.denom` when nothing is standardized, so
    #that each per-stratum child does not re-resolve it independently.
    X[["s.d.denom"]] <- .resolve_s.d.denom(X, prep[["var_types"]], A[["continuous"]], A[["binary"]]) %or%
      X[["s.d.denom"]]
  }

  #Setup output object
  out <- list()

  #Get list of bal.tabs for each cluster
  out[["Cluster.Balance"]] <- lapply(levels(X[["cluster"]]), function(cl) {
    #Subsetting is inside `tryCatch()` so that errors it raises (e.g., a cluster
    #in which the treatment takes only one value) are labeled with the cluster.
    tryCatch({
      X_cl <- subset_X(X, X[["cluster"]] == cl) |>
        .assign_X_class()

      X_cl[["call"]] <- NULL

      do.call("base.bal.tab", c(list(X_cl), A[setdiff(names(A), names(X_cl))]), quote = TRUE)
    },
    error = function(e) {
      arg::err("in cluster {.str {cl}}: {conditionMessage(e)}")
    })
  })
  
  names(out[["Cluster.Balance"]]) <- levels(X[["cluster"]])
  
  #Create summary of lists
  
  #A subclassified child has no single `Balance` table to summarize, so it is
  #excluded here as multiply imputed data already is.
  if ((cluster.summary || !A[["quick"]]) && is_null(X[["covs.list"]]) &&
      get.treat.type(X[["treat"]]) != "multinomial" && is_null(X[["imp"]]) &&
      is_null(X[["subclass"]])) {
    summ <- .bal.tab_summarize(out[["Cluster.Balance"]], "Balance.Across.Clusters",
                               agg.funs = agg.fun,
                               obs.fun = function(cl) {
                                 grab(cl, "Observations") |>
                                   samplesize_across_clusters()
                               })

    out[names(summ)] <- summ
  }
  
  
  out[["call"]] <- X[["call"]]
  
  attr(out, "print.options") <- c(.attr(out[["Cluster.Balance"]][[1L]], "print.options"),
                                  list(which.cluster = which.cluster,
                                       cluster.summary = cluster.summary,
                                       cluster.fun = agg.fun))
  
  set_class(out, c("bal.tab.cluster", "bal.tab"))
}

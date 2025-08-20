#' Using `bal.tab()` with Subclassified Data
#' @name class-bal.tab.subclass
#' 
#' @description
#' When using [bal.tab()] with subclassified data, i.e., data split into subclasses where balance may hold, the output will be different from the standard, non-subclassified case, and there is an additional option for controlling display. This page outlines the outputs and options in this case.
#'     
#' There are two main components of the output of `bal.tab()` with subclassified data: the balance within subclasses and the balance summary across subclasses. The within-subclass balance displays essentially are standard balance displays for each subclass, except that only "adjusted" values are available, because the subclassification itself is the adjustment.
#'     
#' The balance summary is, for each variable, like a weighted average of the balance statistics across subclasses. This is computed internally by assigning each individual a weight based on their subclass and treatment group membership and then computing weighted balance statistics as usual with these weights. This summary is the same one would get if subclasses were supplied to the `match.strata` argument rather than to `subclass`. Because the means and mean differences are additive, their computed values will be weighted averages of the subclass-specific values, but for other statistics, the computed values will not be. 
#'     
#' @section Allowable arguments:
#' 
#' There are three arguments for `bal.tab()` that relate to subclasses: `subclass`, `which.subclass`, and `subclass.summary`.
#' 
#' \describe{
#'     \item{`subclass`}{For the `data.frame` and formula methods of `bal.tab()`, a vector of subclass membership or the name of the variable in `data` containing subclass membership. When using subclassification with a function compatible with \pkg{cobalt}, such as `matchit()` in \pkg{MatchIt}, this argument can be omitted because the subclasses are in the output object.}
#'     \item{`which.subclass`}{This is a display option that does not affect computation. If `.all`, all subclasses in `subclass` will be displayed. If `.none` (the default), no subclasses will be displayed. Otherwise, can be a vector of subclass indices for which to display balance.}
#'     \item{`subclass.summary`}{This is a display option that does not affect computation. If `TRUE`, the balance summary across subclasses will be displayed. The default is `TRUE`, and if `which.subclass` is `.none`, it will automatically be set to `TRUE`.}
#' }
#' 
#' @section Output:
#' The output is a `bal.tab.subclass` object, which inherits from `bal.tab`. It has the following elements:
#'         
#' * `Subclass.Balance`: A list of data frames containing balance information for each covariate in each subclass.
#' * `Balance.Across.Subclass`: A data frame containing balance statistics for each covariate aggregated across subclasses and for the original sample (i.e., unadjusted). See [bal.tab()] for details on what this includes.
#' * `Observations`: A table of sample sizes in each subclass and overall.
#' 
#' 
#' @seealso
#' * [bal.tab()]
#' * [bal.tab.data.frame()]
#' * [print.bal.tab()]
#' 
NULL

base.bal.tab.subclass <- function(X,
                                  type,
                                  int = FALSE,
                                  poly = 1,
                                  continuous,
                                  binary,
                                  imbalanced.only = getOption("cobalt_imbalanced.only", FALSE),
                                  un = getOption("cobalt_un", FALSE),
                                  disp = NULL,
                                  which.subclass = NA,
                                  subclass.summary = getOption("cobalt_subclass.summary"),
                                  disp.bal.tab = getOption("cobalt_disp.bal.tab", TRUE),
                                  disp.call = getOption("cobalt_disp.call", FALSE),
                                  abs = FALSE,
                                  quick = TRUE,
                                  ...) {
  #Preparations
  A <- clear_null(list(...))
  A$subset <- NULL
  
  if (type == "bin" && get.treat.type(X$treat) != "binary") {
    .err("the treatment must be a binary variable")
  }
  
  if (missing(continuous)) continuous <- getOption("cobalt_continuous", "std")
  if (missing(binary)) binary <- getOption("cobalt_binary", switch(type, bin = "raw", "std"))
  
  X$subclass <- factor(X$subclass)
  
  if (missing(which.subclass)) {
    which.subclass <- {
      if (isTRUE(A[["disp.subclass"]])) seq_len(nlevels(X$subclass))
      else NA_integer_
    }
  }
  
  if (is_null(A[["disp.subclass"]])) {
    A[["disp.subclass"]] <- !anyNA(which.subclass)
  }
  
  if (is_null(subclass.summary)) {
    subclass.summary <- is_not_null(which.subclass) && 
      (anyNA(which.subclass) || !is.numeric(which.subclass) || 
         (is.numeric(which.subclass) && !any(which.subclass %in% seq_len(nlevels(X$subclass)))))
  }
  
  no.adj <- FALSE
  
  if (is_null(X$s.weights)) {
    X$s.weights <- rep_with(1, X$treat)
  }
  
  disp <- process_disp(disp, ...)
  
  #Actions
  
  out <- list()
  
  C <- do.call(".get_C2", c(X, A[setdiff(names(A), names(X))],
                            list(int = int, poly = poly)),
               quote = TRUE)
  co.names <- attr(C, "co.names")
  
  var_types <- attr(C, "var_types")
  
  if (get.treat.type(X$treat) != "continuous" &&
      "mean.diffs" %in% X$stats &&
      ((binary == "std" && any(var_types == "Binary")) ||
       (continuous == "std" && !all(var_types == "Binary")))) {
    X$s.d.denom <- .get_s.d.denom(X$s.d.denom,
                                  estimand = X$estimand,
                                  weights = X$weights, 
                                  subclass = X$subclass,
                                  treat = X$treat,
                                  focal = X$focal)
  }
  else if (get.treat.type(X$treat) == "continuous" &&
           any(c("correlations", "spearman.correlations", "distance.correlations") %in% X$stats) &&
           ((binary == "std" && any(var_types == "Binary")) ||
            (continuous == "std" && !all(var_types == "Binary")))) {
    X$s.d.denom <- .get_s.d.denom.cont(X$s.d.denom,
                                       weights = X$weights,
                                       subclass = X$subclass)
  }
  
  out[["Subclass.Balance"]] <- do.call("balance_table_subclass", 
                                       c(list(C,
                                              type = type, 
                                              weights = NULL, 
                                              treat = X$treat, 
                                              subclass = X$subclass,
                                              continuous = continuous,
                                              binary = binary, 
                                              s.d.denom = X$s.d.denom[1L], 
                                              thresholds = X$thresholds,
                                              disp = disp,
                                              stats = X$stats, 
                                              s.weights = X$s.weights, 
                                              abs = abs, 
                                              quick = quick), A), quote = TRUE)
  
  if (subclass.summary || !quick) {
    out[["Balance.Across.Subclass"]] <- {
      if (type == "bin") {
        do.call("balance_table", 
                c(list(C, 
                       type = type, 
                       weights = data.frame(Adj = strata2weights(X$subclass, X$treat,
                                                                 X$estimand, X$focal)), 
                       treat = X$treat, 
                       s.d.denom = X$s.d.denom[1L], 
                       s.weights = X$s.weights, 
                       continuous = continuous, 
                       binary = binary, 
                       thresholds = X$thresholds,
                       un = un, 
                       disp = disp,
                       stats = X$stats, 
                       abs = abs, 
                       no.adj = FALSE, quick = quick, 
                       var_types = attr(C, "var_types"),
                       s.d.denom.list = X$s.d.denom.list), A), quote = TRUE)
      }
      else if (type == "cont") {
        do.call("balance_table_across_subclass_cont", 
                c(list(do.call("balance_table", c(list(C, 
                                                       type = type, 
                                                       weights = NULL,
                                                       treat = X$treat, 
                                                       s.d.denom = X$s.d.denom[1L], 
                                                       s.weights = X$s.weights, 
                                                       continuous = continuous, 
                                                       binary = binary, 
                                                       thresholds = X$thresholds,
                                                       un = un, 
                                                       disp = disp,
                                                       stats = X$stats, 
                                                       abs = abs, 
                                                       no.adj = TRUE, 
                                                       quick = quick, 
                                                       var_types = attr(C, "var_types"),
                                                       s.d.denom.list = X$s.d.denom.list), A), quote = TRUE), 
                       balance.table.subclass.list = out[["Subclass.Balance"]], 
                       subclass.obs = out[["Observations"]], 
                       r.threshold = X$thresholds[["correlations"]]), A), quote = TRUE)
      }
    }
  }
  
  #Reassign disp... and ...threshold based on balance table output
  compute <- attr(out[["Subclass.Balance"]], "compute")
  thresholds <- attr(out[["Subclass.Balance"]], "thresholds")
  disp <- attr(out[["Subclass.Balance"]], "disp")
  
  for (s in compute) {
    if (is_not_null(thresholds[[s]])) {
      out[[paste.("Balanced", s, "Subclass")]] <- setNames(do.call("data.frame", lapply(out[["Subclass.Balance"]], function(x) .baltal(x[[STATS[[s]]$Threshold]]))),
                                                           paste("Subclass", levels(X$subclass)))
      max.imbal.list <- lapply(out[["Subclass.Balance"]], function(x) {
        .max_imbal(x[x[["Type"]] != "Distance", , drop = FALSE], 
                   col.name = paste.(STATS[[s]]$bal.tab_column_prefix, "Adj"), 
                   thresh.col.name = STATS[[s]]$Threshold, 
                   abs_stat = STATS[[s]]$abs)
      })
      out[[paste.("Max.Imbalance", s, "Subclass")]] <- as.data.frame(do.call("rbind", max.imbal.list), 
                                                                     row.names = paste("Subclass", levels(X$subclass)))
    }
  }
  
  out[["Observations"]] <- samplesize(treat = X$treat, 
                                      type = type,
                                      weights = NULL, 
                                      subclass = X$subclass,
                                      s.weights = X$s.weights, 
                                      method = X$method, 
                                      discarded = X$discarded)
  
  out[["call"]] <- X$call
  attr(out, "print.options") <- list(thresholds = thresholds,
                                     imbalanced.only = imbalanced.only,
                                     un = un, 
                                     disp = disp,
                                     compute = compute, 
                                     disp.adj = !no.adj, 
                                     which.subclass = which.subclass,
                                     disp.subclass = A[["disp.subclass"]],
                                     subclass.summary = subclass.summary,
                                     disp.bal.tab = disp.bal.tab, 
                                     disp.call = disp.call,
                                     abs = abs,
                                     continuous = continuous,
                                     binary = binary,
                                     quick = quick,
                                     treat_names = treat_vals(X$treat),
                                     type = type,
                                     co.names = co.names)
  
  set_class(out, c("bal.tab.subclass", "bal.tab"))
}

base.bal.tab.subclass.binary <- function(X, ...) {
  base.bal.tab.subclass(X, type = "bin", ...)
}

base.bal.tab.subclass.cont <- function(X, ...) {
  .err("subclasses are not yet compatible with continuous treatments")
  # base.bal.tab.subclass(X, type = "cont", ...)
}

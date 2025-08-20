#' Using `bal.tab()` with Multi-Category Treatments
#' @name class-bal.tab.multi
#' 
#' @description
#' When using [bal.tab()] with multi-category treatments, the output will be different from the case with binary or continuous treatments, and there are some options that are common across all `bal.tab()` methods. This page outlines the outputs and options in this case.
#'     
#' There are two main components of the output of `bal.tab()` with multi-category treatments: the two-group treatment comparisons and the balance summary. The two-group treatment comparisons are standard binary treatment comparison either for pairs of groups (e.g., for treatments A, B, and C, "A vs. B", "A vs. C", and "B vs. C") or each group against all the groups (i.e., the entire sample).
#'     
#' The balance summary is, for each variable, the greatest imbalance across all two-group comparisons. So, for variable X1, if "A vs. B" had a standardized mean difference of 0.52, "A vs. C" had a standardized mean difference of .17,  and "B vs. C" had a standardized mean difference of .35, the balance summary would have 0.52 for the value of the standardized mean difference for X1. The same goes for other variables and other measures of balance. If the greatest observed imbalance is tolerable, then all other imbalances for that variable will be tolerable too, so focusing on reducing the greatest imbalance is sufficient for reducing imbalance overall. (Note that when `s.d.denom = "pooled"`, i.e., when the estimand is the ATE, the pooled standard deviation in the denominator will be the average of the standard deviations across all treatment groups, not just those used in the pairwise comparison.) The balance summary will not be computed if multiply imputed data are used.
#' 
#' @section Allowable arguments:
#' 
#' There are four arguments for each `bal.tab()` method that can handle multi-category treatments: `pairwise`, `focal`, `which.treat`, and `multi.summary`.
#' 
#' \describe{
#'     \item{`pairwise`}{
#'         Whether to compute the two-group comparisons pairwise or not. If `TRUE`, `bal.tab()` will compute comparisons for each pair of treatments. This can be valuable if treatments are to be compared with one another (which is often the case). If `FALSE`, `bal.tab()` will compute balance for each treatment group against the full unadjusted sample when `focal` is `NULL` and for each non-focal group against the focal group otherwise.
#'     }
#'     \item{`focal`}{
#'         When one group is to be compared to multiple control groups in an ATT analysis, the group considered "treated" is the focal group. By specifying the name or index of the treatment condition considered focal, `bal.tab()` will only compute and display pairwise balance for treatment comparisons that include the focal group when `pairwise = FALSE`.
#'     }
#'     \item{`which.treat`}{
#'         This is a display option that does not affect computation. When displaying the `bal.tab` output, which treatments should be displayed? If a vector of length 1 is entered, all comparisons involving that treatment group will be displayed. If a vector of length 2 or more is entered, all comparisons involving treatments that both appear in the input will be displayed. For example, inputting `"A"` will display "A vs. B" and "A vs. C", while entering `c("A", "B")` will only display "A vs. B". `.none` indicates no treatment comparisons will be displayed, and `.all` indicates all treatment comparisons will be displayed. `.none` is the default.
#'     }
#'     \item{`multi.summary`}{
#'         If `TRUE`, the balance summary across all comparisons will be computed and displayed. This includes one row for each covariate with maximum balance statistic across all pairwise comparisons. Note that, if variance ratios or KS statistics are requested in addition to mean differences, the displayed values may not come from the same pairwise comparisons; that is, the greatest standardized mean difference and the greatest variance ratio may not come from the same comparison. The default is `TRUE`, and if `which.treat` is `.none`, it will automatically be set to `TRUE`.
#'     }
#' }
#' 
#' @section Output:
#' The output is a `bal.tab.multi` object, which inherits from `bal.tab`. It has the following elements:
#'         
#' * `Pair.Balance`:For each pair of treatment groups, a regular `bal.tab` object containing a balance table, a sample size summary, and other balance assessment tools, depending on which options are specified. If `pairwise` is `FALSE`, the comparisons will be between each group and the groups combined (labeled "All") when `focal` is `NULL` and between each non-focal group and the focal group otherwise.
#' * `Balance.Across.Pairs`: The balance summary across two-group comparisons. This will include the greatest (i.e., maximum) absolute balance statistics(s) for each covariate across all comparisons computed. Thresholds can be requested for each balance measure as with binary treatments.
#' * `Observations`: A table of sample sizes or effective sample sizes for each treatment group before and after adjustment.
#'     
#' As with other methods, multiple weights can be specified, and values for all weights will appear in all tables.
#' 
#' @note
#' In versions 4.3.1 and earlier, setting `pairwise = FALSE` would compare each group to the full adjusted sample. Now, each group is compared to the full *un*adjusted sample (unadjusted except for `s.weights`, if supplied).
#'     
#' In versions 4.3.1 and earlier, `pairwise` was ignored with non-`NULL` `focal` and was automatically set to `FALSE`. `pairwise` can be specified and its default is now `TRUE`, so balance between all treatment groups will be computed by default rather than only between each non-group and the focal group. To recover previous functionality, set `pairwise = FALSE` with non-`NULL` `focal`.
#' 
#' @seealso
#' * [bal.tab()]
#' * [bal.tab.data.frame()]
#' * [print.bal.tab()]
#' * `vignette("segmented-data")` for examples
#'
NULL

base.bal.tab.multi <- function(X,
                               pairwise = TRUE,
                               which.treat,
                               multi.summary = getOption("cobalt_multi.summary"),
                               ...) {
  A <- list(...)
  
  #Preparations
  if (is_not_null(X$weights))  {
    check_if_zero_weights(X$weights, X$treat)
    if (ncol(X$weights) == 1L) names(X$weights) <- "Adj"
  }
  if (is_null(A[["quick"]])) A[["quick"]] <- TRUE
  
  if (missing(which.treat)) {
    if (is_null(X$imp)) which.treat <- NA
    else which.treat <- NULL
  }
  
  if (is_null(multi.summary)) {
    multi.summary <- is_not_null(which.treat) && anyNA(which.treat)
  }
  
  #Treat is a factor variable
  if (is_null(X$focal)) {
    treat.combinations <- {
      if (pairwise) utils::combn(treat_names(X$treat), 2L, simplify = FALSE)
      else lapply(treat_names(X$treat), c, "All")
    }
  }
  else {
    if (length(X$focal) > 1L) {
      .err("`focal` must be a vector of length 1 containing the name or index of the focal treatment group")
    }
    
    if (is.numeric(X$focal)) {
      X$focal <- levels(X$treat)[X$focal]
    }
    
    if (!is.character(X$focal)) {
      .err("`focal` must be the name or index of the focal treatment group")
    }
    
    treat.combinations <- {
      if (pairwise) utils::combn(treat_names(X$treat), 2L, simplify = FALSE)
      else lapply(setdiff(treat_names(X$treat), X$focal), c, X$focal)
    }
  }
  
  X$covs <- do.call(".get_C2", c(X, A[setdiff(names(A), names(X))]), quote = TRUE)
  
  var_types <- attr(X$covs, "var_types")
  
  if (is_null(A$continuous)) A$continuous <- getOption("cobalt_continuous", "std")
  if (is_null(A$binary)) A$binary <- getOption("cobalt_binary", "raw")
  
  
  if ("mean.diffs" %in% X$stats &&
      ((A$binary == "std" && any(var_types == "Binary")) ||
       (A$continuous == "std" && !all(var_types == "Binary")))) {
    X$s.d.denom <- .get_s.d.denom(X$s.d.denom,
                                  estimand = X$estimand,
                                  weights = X$weights, 
                                  subclass = X$subclass,
                                  treat = X$treat,
                                  focal = X$focal)
    
    bin.vars <- var_types == "Binary"
    
    if (is_null(X$weights)) {
      X$s.d.denom.list <- list(.compute_s.d.denom(X$covs, X$treat, s.d.denom = X$s.d.denom,
                                                  s.weights = X$s.weights, bin.vars = bin.vars))
    }
    else {
      X$s.d.denom.list <- setNames(lapply(seq_along(X$s.d.denom),
                                          function(i) .compute_s.d.denom(X$covs, X$treat,
                                                                         s.d.denom = X$s.d.denom[i], s.weights = X$s.weights, 
                                                                         bin.vars = bin.vars, weighted.weights = X$weights[[i]])),
                                   names(X$s.d.denom))
    }
  }
  
  #Setup output object
  out <- list()
  
  if (pairwise || is_not_null(X$focal)) {
    balance.tables <- lapply(treat.combinations, function(t) {
      X_t <- subset_X(X, X$treat %in% treat_vals(X$treat)[t])
      # X_t$treat <- process_treat(X_t$treat)
      X_t <- .assign_X_class(X_t)
      X_t$call <- NULL
      do.call("base.bal.tab", c(list(X_t), A[setdiff(names(A), names(X_t))]), quote = TRUE)
    })
  }
  else {
    if (any(treat_vals(X$treat) == "All")) {
      .err('"All" cannot be the name of a treatment level. Please rename your treatments')
    }
    balance.tables <- lapply(treat.combinations, function(t) {
      n <- length(X$treat)
      X_t <- X
      X_t$call <- NULL
      X_t <- subset_X(X_t, c(seq_len(n), which(X$treat == treat_vals(X$treat)[t[1L]])))
      X_t$treat <- factor(rep(0:1, times = c(n, sum(X$treat == treat_vals(X$treat)[t[1L]]))),
                          levels = c(0, 1), labels = c("All", t[1L])) |>
        process_treat()
      
      if (is_not_null(X_t$weights)) {
        X_t$weights[X_t$treat == "All", ] <- 1 #Uncomment to compare each group to unweighted dist.
      }
      
      X_t <- .assign_X_class(X_t)
      
      do.call("base.bal.tab", c(list(X_t), A[setdiff(names(A), names(X_t))]), quote = TRUE)
    })
  }
  
  for (i in seq_along(balance.tables)) {
    names(balance.tables)[i] <- paste(rev(treat.combinations[[i]]), collapse = " vs. ")
  }
  
  out[["Pair.Balance"]] <- balance.tables
  
  if ((multi.summary || !A$quick) && is_null(X$imp)) {
    out[["Balance.Across.Pairs"]] <- balance_summary(balance.tables, 
                                                     agg.funs = "max")
    
    out <- c(out,
             threshold_summary(compute = attr(out[["Pair.Balance"]][[1L]][["Balance"]], "compute"),
                               thresholds = attr(out[["Pair.Balance"]][[1L]][["Balance"]], "thresholds"),
                               no.adj = !attr(out[["Pair.Balance"]][[1L]], "print.options")$disp.adj,
                               balance.table = out[["Balance.Across.Pairs"]],
                               weight.names = attr(out[["Pair.Balance"]][[1L]], "print.options")$weight.names,
                               agg.fun = "max"))
    
    out[["Observations"]] <- samplesize_multi(balance.tables, treat_names(X$treat), X$focal)
  }
  
  out[["call"]] <- X$call
  
  attr(out, "print.options") <- c(attr(out[["Pair.Balance"]][[1L]], "print.options"),
                                  list(treat_vals_multi = treat_vals(X$treat),
                                       which.treat = which.treat,
                                       multi.summary = multi.summary,
                                       pairwise = pairwise))
  
  attr(out, "print.options")[["treat_names"]] <- NULL
  
  set_class(out, c("bal.tab.multi", "bal.tab"))
}

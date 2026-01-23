#' @title Print Results of a Call to `bal.tab()`
#' 
#' @description Prints `bal.tab()` output in a clean way. Provides options for printing.
#' 
#' @param x a `bal.tab` object; the output of a call to [bal.tab()].
#' @param imbalanced.only `logical`; whether to display only the covariates that failed to meet at least one of balance thresholds. Depends only on whether threshold were initial set in the call to `bal.tab()` and not on any arguments to `print()` (except `disp.bal.tab`).
#' @param un `logical`; whether to display balance values for the unadjusted sample. Ignored (and set to `TRUE`) if no conditioning was performed.
#' @param disp.bal.tab `logical`; whether to display the table of balance statistics. If `FALSE`, only other values (e.g., the call, sample sizes, balance tallies, and maximum imbalances) will be presented.
#' @param disp.call `logical`; whether to display the function call for the input object, if any.
#' @param stats `character`; which statistic(s) should be reported. For binary or multi-category treatments, the options are "mean.diffs" for mean differences (standardized or not according the selected `bal.tab()` options), "variance.ratios" for variance ratios, and "ks.statistics" for Kolmogorov-Smirnov statistics. "mean.diffs" is the default. For continuous treatments, the only option is "correlations" for treatment-covariate correlations. Multiple options are allowed. Abbreviations allowed. Statistics that weren't requested in the original call to `bal.tab()` cannot be requested with `print()` unless `quick = FALSE` in the original call.
#' @param disp.thresholds `logical`; whether to display thresholds for each statistic for which thresholds were originally requested in the call to `bal.tab()`. Should be a named logical vector with names corresponding to the thresholds. For example, if thresholds for mean differences were requested in `bal.tab()`, set `disp.thresholds = c(m = FALSE)` to prevent them from being printed. If a statistic was prevented from being displayed by another argument to `print()`, the thresholds will not be displayed.
#' @param disp `character`; which distribution summary statistics to display. Allowable options include "means" and "sds". Statistics that weren't requested in the original call to `bal.tab()` cannot be requested with `print()` unless `quick = FALSE` in the original call.
#' @param which.subclass when used with subclassification, which subclass(es) to display. If `NULL`, all subclasses will be displayed. If `NA`, no subclasses will be displayed. Otherwise, can be a vector of subclass indices for which to display balance. To display the subclasses requested in the original call to `bal.tab()`, omit this argument. See [`class-bal.tab.subclass`] for details.
#' @param subclass.summary `logical`; when used with subclassification, whether to display the subclass balance summary table. If `which.subclass` is `NA`, `subclass.summary` will be set to `TRUE`. See [`class-bal.tab.subclass`] for details.
#' @param which.imp when used with multiply imputed data, which imputation(s) to display. If `NULL`, all imputations will be displayed. If `NA`, no imputations will be displayed. Otherwise, can be a vector of imputations numbers for which to display balance. To display the imputations requested in the original call to `bal.tab()`, omit this argument. See [`class-bal.tab.imp`] for details.
#' @param imp.summary `logical`; when used with multiply imputed data, whether to display the imputation summary table. If `which.imp` is `NA`, `imp.summary` will be set to `TRUE`. See [`class-bal.tab.imp`] for details.
#' @param imp.fun `character`; when used with multiply imputed data, a character vector of functions of balance statistics to display when displaying balance across imputations. Can be "mean", "min", or "max". More than one are allowed. See [`class-bal.tab.imp`] for details.
#' @param which.treat when used with multi-category treatments, which treatments to display. See [`bal.tab.multi()`][class-bal.tab.multi] for details.
#' @param multi.summary `logical`; when used with multi-category treatments, whether to display the balance summary table across pairwise comparisons. See [`bal.tab.multi()`][class-bal.tab.multi] for details.
#' @param which.time when used with longitudinal treatments, which time periods to display if longitudinal treatments are used. See [`class-bal.tab.msm`] for details.
#' @param msm.summary `logical`; when used with longitudinal treatments, whether to display the balance summary table across time periods. See [`class-bal.tab.msm`] for details.
#' @param which.cluster when used with clustered data, which cluster(s) to display. If `NULL`, all clusters will be displayed. If `NA`, no clusters will be displayed. Otherwise, can be a vector of cluster names or numerical indices for which to display balance. Indices correspond to the alphabetical order of cluster names. To display the clusters requested in the original call to `bal.tab()`, omit this argument. See [`class-bal.tab.cluster`] for details.
#' @param cluster.summary `logical`; when used with clustered data, whether to display the cluster summary table. If `which.cluster` is `NA`, `cluster.summary` will be set to `TRUE`. See [`class-bal.tab.cluster`] for details.
#' @param cluster.fun `character`; when used with clustered data, a character vector of functions of balance statistics to display when displaying balance across clusters. Can be "mean", "min", or "max". More than one are allowed. See [`class-bal.tab.cluster`] for details.
#' @param digits the number of digits to display.
#' @param ... further arguments passed to or from other methods.
#' 
#' @details Simply calling `bal.tab()` will print its results, but it can be useful to store the results into an object and print them again later, possibly with different print options specified. The `print()` function automatically dispatches the correct method for the `bal.tab` object given.
#' 
#' Any parameter used in `bal.tab()` for calculations, such as `int`, `addl`, or `distance`, cannot be used with `print()`; only those parameters listed above, those that solely determine printing options, can be used. To change computation options, a new call to `bal.tab()` must be performed.
#' 
#' Prior versions of `print()` had separate methods for each `bal.tab` class. Now they are dispatched internally.
#' 
#' @note Unless `quick = FALSE` in the original call to `bal.tab()` (which is not the default), some values may not be calculated, in which case using `print()` will not display these values even when requested. For example, if `stats = "m"` and `quick = TRUE` in the original call to `bal.tab()` (the default for both), setting `stats = "ks"` in `print()` will not print the KS statistics because they were not calculated.
#' 
#' @seealso 
#' [print()], [bal.tab()]
#' 
#' [`display-options`] for further information on some of these options.
#' 
#' @examplesIf rlang::is_installed("WeightIt")
#' data("lalonde", package = "cobalt")
#' library(WeightIt)
#' 
#' w.out <- weightit(treat ~ age + educ + married +
#'                     race + re74 + re75, 
#'                   data = lalonde)
#' 
#' b <- bal.tab(w.out, stats = c("m", "v", "ks"), 
#'              un = TRUE, v.threshold = 2)
#' 
#' print(b, un = FALSE, stats = c("m", "v"),
#'       disp.thresholds = c(v = FALSE))

#' @rdname print.bal.tab
#' @exportS3Method print bal.tab
print.bal.tab <- function(x, imbalanced.only, un, disp.bal.tab, disp.call,
                          stats, disp.thresholds, disp, 
                          which.subclass, subclass.summary,
                          which.imp, imp.summary, imp.fun, 
                          which.treat, multi.summary, 
                          which.time, msm.summary,
                          which.cluster, cluster.summary, cluster.fun, 
                          digits = max(3L, getOption("digits") - 3), ...) {
  
  #Replace .all and .none with NULL and NA respectively
  .call <- match.call(expand.dots = TRUE)
  .alls <- vapply(seq_along(.call), function(z) identical(.call[[z]], quote(.all)), logical(1L))
  .nones <- vapply(seq_along(.call), function(z) identical(.call[[z]], quote(.none)), logical(1L))
  if (any(c(.alls, .nones))) {
    .call[.alls] <- expression(NULL)
    .call[.nones] <- expression(NA)
    return(eval.parent(.call))
  }
  
  A <- try_chk(c(as.list(environment()), list(...))[-1L])
  
  A[vapply(A, rlang::is_missing, logical(1L))] <- NULL
  
  unpack_p.ops <- function(b) {
    out <- do.call("print_process", c(list(b), A), quote = TRUE)
    
    if (inherits(b, c("bal.tab.bin", "bal.tab.cont"))) {
      return(out)
    }
    
    b_ <- b[[which(endsWith(names(b), ".Balance"))]][[1L]]
    if (!inherits(b_, "bal.tab")) {
      return(out)
    }
    
    c(out, unpack_p.ops(b_))
  }
  
  p.ops <- unpack_p.ops(x)
  
  #Prevent exponential notation printing
  rlang::with_options({
    bal.tab_print(x, p.ops)
  }, scipen = 999)
}

bal.tab_print <- function(x, p.ops) {
  UseMethod("bal.tab_print")
}
#' @exportS3Method NULL
bal.tab_print.bal.tab <- function(x, p.ops) {
  
  call <- if (p.ops$disp.call) x$call else NULL
  
  balance <- x$Balance
  
  thresholds <- setdiff(names(p.ops$thresholds), p.ops$drop.thresholds)
  baltal <- setNames(x[paste.("Balanced", thresholds)], thresholds)
  maximbal <- setNames(x[paste.("Max.Imbalance", thresholds)], thresholds)
  
  nn <- x$Observations
  
  if (is_not_null(call)) {
    cat(.ul("Call") %+% "\n " %+% paste(deparse(call), collapse = "\n") %+% "\n\n")
  }
  
  if (p.ops$disp.bal.tab) {
    keep.row <- {
      if (p.ops$imbalanced.only)
        rowSums(apply(balance[grepl(".Threshold", names(balance), fixed = TRUE)], 2L,
                      function(x) !is.na(x) & startsWith(x, "Not Balanced"))) > 0
      else
        rep.int(TRUE, nrow(balance))
    }
    
    keep.col <- setNames(as.logical(c(TRUE, 
                                      rep(unlist(lapply(p.ops$compute[p.ops$compute %nin% all_STATS()], function(s) {
                                        p.ops$un && s %in% p.ops$disp
                                      })), switch(p.ops$type, bin = 2L, cont = 1L)),
                                      unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS()[!get_from_STATS("adj_only")]], function(s) {
                                        c(p.ops$un && s %in% p.ops$disp,
                                          if (p.ops$un && !p.ops$disp.adj && is_not_null(p.ops$thresholds[[s]])) s %in% thresholds)
                                      })),
                                      rep(c(rep(unlist(lapply(p.ops$compute[p.ops$compute %nin% all_STATS()], function(s) {
                                        p.ops$disp.adj && s %in% p.ops$disp
                                      })), switch(p.ops$type, bin = 2L, cont = 1L)),
                                      unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS()], function(s) {
                                        c(p.ops$disp.adj && s %in% p.ops$disp,
                                          if (p.ops$disp.adj && is_not_null(p.ops$thresholds[[s]])) s %in% thresholds)
                                      }))
                                      ), 
                                      p.ops$nweights + !p.ops$disp.adj))),
                         names(balance))
    
    cat(.ul("Balance Measures") %+% "\n")
    
    if (is_null(keep.row)) cat(.it("No covariates to display.") %+% "\n")
    else if (!any(keep.row)) cat(.it("All covariates are balanced.") %+% "\n")
    else .print_data_frame(round_df_char(balance[keep.row, keep.col, drop = FALSE],
                                         p.ops$digits, na_vals = "."))
    cat("\n")
  }
  
  for (s in p.ops$compute) {
    if (is_not_null(baltal[[s]])) {
      cat(.ul(sprintf("Balance tally for %s", STATS[[s]]$balance_tally_for)) %+% "\n")
      .print_data_frame(baltal[[s]])
      cat("\n")
    }
    if (is_not_null(maximbal[[s]])) {
      cat(.ul(sprintf("Variable with the greatest %s", STATS[[s]]$variable_with_the_greatest)) %+% "\n")
      .print_data_frame(round_df_char(maximbal[[s]], p.ops$digits, na_vals = "."), row.names = FALSE)
      cat("\n")
    }
  }
  
  if (is_not_null(nn)) {
    
    drop.nn <- rowSums(nn) == 0
    ss.type <- .attr(nn, "ss.type")[!drop.nn]
    nn <- nn[!drop.nn, , drop = FALSE]
    
    if (all(c("All (ESS)", "All (Unweighted)") %in% rownames(nn)) && 
        all(check_if_zero(nn["All (ESS)", ] - nn["All (Unweighted)", ]))) {
      nn <- nn[rownames(nn) != "All (Unweighted)", , drop = FALSE]
      rownames(nn)[rownames(nn) == "All (ESS)"] <- "All"
    }
    
    if (all(c("Matched (ESS)", "Matched (Unweighted)") %in% rownames(nn)) && 
        all(check_if_zero(nn["Matched (ESS)", ] - nn["Matched (Unweighted)", ]))) {
      nn <- nn[rownames(nn) != "Matched (Unweighted)", , drop = FALSE]
      rownames(nn)[rownames(nn) == "Matched (ESS)"] <- "Matched"
    }
    
    cat(.ul(.attr(nn, "tag")) %+% "\n")
    
    print.warning <- FALSE
    
    if (length(ss.type) > 1L && nunique.gt(ss.type[-1L], 1L)) {
      ess <- ifelse(ss.type == "ess", "*", "")
      nn <- setNames(cbind(nn, ess), c(names(nn), ""))
      print.warning <- TRUE
    }
    
    .print_data_frame(round_df_char(nn, digits = min(2L, p.ops$digits), pad = " "))
    
    if (print.warning) {
      cat(.it("* indicates effective sample size"))
    }
  }
  
  invisible(x)
}
#' @exportS3Method NULL
bal.tab_print.bal.tab.cluster <- function(x, p.ops) {
  
  call <- if (p.ops$disp.call) x$call else NULL
  c.balance <- x$Cluster.Balance
  c.balance.summary <- x$Balance.Across.Clusters
  
  thresholds <- setdiff(names(p.ops$thresholds), p.ops$drop.thresholds)
  baltal <- setNames(x[paste.("Balanced", thresholds)], thresholds)
  maximbal <- setNames(x[paste.("Max.Imbalance", thresholds)], thresholds)
  
  nn <- x$Observations
  
  #Printing
  if (is_not_null(call)) {
    cat(.ul("Call") %+% "\n " %+% paste(deparse(call), collapse = "\n") %+% "\n\n")
  }
  
  if (is_not_null(p.ops$which.cluster)) {
    cat(.ul("Balance by cluster") %+% "\n")
    for (i in p.ops$which.cluster) {
      cat("\n - - - " %+% .it("Cluster: " %+% names(c.balance)[i]) %+% " - - - \n")
      print(c.balance[[i]])
      # bal.tab_print(c.balance[[i]], p.ops)
    }
    cat(strrep(" -", round(nchar(sprintf("\n - - - Cluster: %s - - - ", names(c.balance)[i])) / 2)), "\n\n")
  }
  
  if (isTRUE(as.logical(p.ops$cluster.summary)) && is_not_null(c.balance.summary)) {
    s.keep.col <- setNames(as.logical(c(TRUE, 
                                        unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS()[!get_from_STATS("adj_only")]], function(s) {
                                          c(unlist(lapply(p.ops$computed.cluster.funs, function(af) {
                                            p.ops$un && s %in% p.ops$disp && af %in% p.ops$cluster.fun
                                          })), 
                                          if (p.ops$un && !p.ops$disp.adj && length(p.ops$cluster.fun) == 1L &&
                                              is_not_null(p.ops$thresholds[[s]]))
                                            s %in% thresholds)
                                        })),
                                        rep(
                                          unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS()], function(s) {
                                            c(unlist(lapply(p.ops$computed.cluster.funs, function(af) {
                                              p.ops$disp.adj && s %in% p.ops$disp && af %in% p.ops$cluster.fun
                                            })), 
                                            if (p.ops$disp.adj && length(p.ops$cluster.fun) == 1L &&
                                                is_not_null(p.ops$thresholds[[s]]))
                                              s %in% thresholds)
                                          })),
                                          p.ops$nweights + !p.ops$disp.adj)
    )), names(c.balance.summary))
    
    if (p.ops$disp.bal.tab) {
      cat(.ul("Balance summary across all clusters") %+% "\n")
      c.balance.summary[, s.keep.col, drop = FALSE] |>
        round_df_char(p.ops$digits, na_vals = ".") |>
        .print_data_frame()
      cat("\n")
    }
    
    for (s in p.ops$compute) {
      if (is_not_null(baltal[[s]])) {
        cat(.ul(sprintf("Balance tally for %s", STATS[[s]]$balance_tally_for)) %+% "\n")
        .print_data_frame(baltal[[s]])
        cat("\n")
      }
      if (is_not_null(maximbal[[s]])) {
        cat(.ul(sprintf("Variable with the greatest %s", STATS[[s]]$variable_with_the_greatest)) %+% "\n")
        maximbal[[s]] |>
          round_df_char(p.ops$digits, na_vals = ".") |>
          .print_data_frame(row.names = FALSE)
        cat("\n")
      }
    }
    
    if (is_not_null(nn)) {
      drop.nn <- rowSums(nn) == 0
      ss.type <- .attr(nn, "ss.type")[!drop.nn]
      nn <- nn[!drop.nn, , drop = FALSE]
      if (all(c("All (ESS)", "All (Unweighted)") %in% rownames(nn)) && 
          all(check_if_zero(nn["All (ESS)", ] - nn["All (Unweighted)", ]))) {
        nn <- nn[rownames(nn) != "All (Unweighted)", , drop = FALSE]
        rownames(nn)[rownames(nn) == "All (ESS)"] <- "All"
      }
      
      if (all(c("Matched (ESS)", "Matched (Unweighted)") %in% rownames(nn)) && 
          all(check_if_zero(nn["Matched (ESS)", ] - nn["Matched (Unweighted)", ]))) {
        nn <- nn[rownames(nn) != "Matched (Unweighted)", , drop = FALSE]
        rownames(nn)[rownames(nn) == "Matched (ESS)"] <- "Matched"
      }
      
      cat(.ul(.attr(nn, "tag")) %+% "\n")
      
      print.warning <- FALSE
      
      if (length(ss.type) > 1L && nunique.gt(ss.type[-1L], 1L)) {
        ess <- ifelse(ss.type == "ess", "*", "")
        nn <- setNames(cbind(nn, ess), c(names(nn), ""))
        print.warning <- TRUE
      }
      
      .print_data_frame(round_df_char(nn, digits = min(2L, p.ops$digits), pad = " "))
      
      if (print.warning) {
        cat(.it("* indicates effective sample size"))
      }
    }
  }
  
  invisible(x)
}
#' @exportS3Method NULL
bal.tab_print.bal.tab.imp <- function(x, p.ops) {
  
  call <- if (p.ops$disp.call) x$call else NULL
  i.balance <- x[["Imputation.Balance"]]
  i.balance.summary <- x[["Balance.Across.Imputations"]]
  
  thresholds <- setdiff(names(p.ops$thresholds), p.ops$drop.thresholds)
  baltal <- setNames(x[paste.("Balanced", thresholds)], thresholds)
  maximbal <- setNames(x[paste.("Max.Imbalance", thresholds)], thresholds)
  
  nn <- x$Observations
  
  #Printing output
  if (is_not_null(call)) {
    cat(.ul("Call") %+% "\n " %+% paste(deparse(call), collapse = "\n") %+% "\n\n")
  }
  
  if (is_not_null(p.ops$which.imp)) {
    cat(.ul("Balance by imputation") %+% "\n")
    for (i in p.ops$which.imp) {
      cat("\n - - - " %+% .it("Imputation " %+% names(i.balance)[i]) %+% " - - - \n")
      print(i.balance[[i]])
      # bal.tab_print(i.balance[[i]], p.ops)
    }
    cat(strrep(" -", round(nchar(sprintf("\n - - - Imputation: %s - - - ", names(i.balance)[i])) / 2)), "\n\n")
  }
  
  if (isTRUE(as.logical(p.ops$imp.summary)) && is_not_null(i.balance.summary)) {
    s.keep.col <- as.logical(c(TRUE, 
                               unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS()[!get_from_STATS("adj_only")]], function(s) {
                                 c(unlist(lapply(p.ops$computed.imp.funs, function(af) {
                                   p.ops$un && s %in% p.ops$disp && af %in% p.ops$imp.fun
                                 })), 
                                 if (p.ops$un && !p.ops$disp.adj && length(p.ops$imp.fun) == 1L && is_not_null(p.ops$thresholds[[s]])) s %in% thresholds)
                               })),
                               rep(
                                 unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS()], function(s) {
                                   c(unlist(lapply(p.ops$computed.imp.funs, function(af) {
                                     p.ops$disp.adj && s %in% p.ops$disp && af %in% p.ops$imp.fun
                                   })), 
                                   if (p.ops$disp.adj && length(p.ops$imp.fun) == 1L && is_not_null(p.ops$thresholds[[s]])) s %in% thresholds)
                                 })),
                                 p.ops$nweights + !p.ops$disp.adj)
    ))
    
    if (p.ops$disp.bal.tab) {
      cat(.ul("Balance summary across all imputations") %+% "\n")
      i.balance.summary[, s.keep.col, drop = FALSE] |>
        round_df_char(p.ops$digits, na_vals = ".") |>
        .print_data_frame()
      cat("\n")
    }
    
    for (s in p.ops$compute) {
      if (is_not_null(baltal[[s]])) {
        cat(.ul(sprintf("Balance tally for %s", STATS[[s]]$balance_tally_for)) %+% "\n")
        .print_data_frame(baltal[[s]])
        cat("\n")
      }
      if (is_not_null(maximbal[[s]])) {
        cat(.ul(sprintf("Variable with the greatest %s", STATS[[s]]$variable_with_the_greatest)) %+% "\n")
        maximbal[[s]] |>
          round_df_char(p.ops$digits, na_vals = ".") |>
          .print_data_frame(row.names = FALSE)
        cat("\n")
      }
    }
    
    if (is_not_null(nn)) {
      drop.nn <- rowSums(nn) == 0
      ss.type <- .attr(nn, "ss.type")[!drop.nn]
      nn <- nn[!drop.nn, , drop = FALSE]
      if (all(c("All (ESS)", "All (Unweighted)") %in% rownames(nn)) && 
          all(check_if_zero(nn["All (ESS)", ] - nn["All (Unweighted)", ]))) {
        nn <- nn[rownames(nn) != "All (Unweighted)", , drop = FALSE]
        rownames(nn)[rownames(nn) == "All (ESS)"] <- "All"
      }
      
      if (all(c("Matched (ESS)", "Matched (Unweighted)") %in% rownames(nn)) && 
          all(check_if_zero(nn["Matched (ESS)", ] - nn["Matched (Unweighted)", ]))) {
        nn <- nn[rownames(nn) != "Matched (Unweighted)", , drop = FALSE]
        rownames(nn)[rownames(nn) == "Matched (ESS)"] <- "Matched"
      }
      
      cat(.ul(.attr(nn, "tag")) %+% "\n")
      
      print.warning <- FALSE
      
      if (length(ss.type) > 1L && nunique.gt(ss.type[-1L], 1L)) {
        ess <- ifelse(ss.type == "ess", "*", "")
        nn <- setNames(cbind(nn, ess), c(names(nn), ""))
        print.warning <- TRUE
      }
      
      .print_data_frame(round_df_char(nn, digits = min(2L, p.ops$digits), pad = " "))
      
      if (print.warning) {
        cat(.it("* indicates effective sample size"))
      }
    }
  }
  
  invisible(x)
}
#' @exportS3Method NULL
bal.tab_print.bal.tab.multi <- function(x, p.ops) {
  
  call <- if (p.ops$disp.call) x$call else NULL
  m.balance <- x[["Pair.Balance"]]
  m.balance.summary <- x[["Balance.Across.Pairs"]]
  
  thresholds <- setdiff(names(p.ops$thresholds), p.ops$drop.thresholds)
  baltal <- setNames(x[paste.("Balanced", thresholds)], thresholds)
  maximbal <- setNames(x[paste.("Max.Imbalance", thresholds)], thresholds)
  
  nn <- x$Observations
  
  #Printing output
  if (is_not_null(call)) {
    cat(.ul("Call") %+% "\n " %+% paste(deparse(call), collapse = "\n") %+% "\n\n")
  }
  
  if (is_not_null(p.ops$disp.treat.pairs)) {
    headings <- setNames(character(length(p.ops$disp.treat.pairs)), p.ops$disp.treat.pairs)
    if (p.ops$pairwise) cat(.ul("Balance by treatment pair") %+% "\n")
    else cat(.ul("Balance by treatment group") %+% "\n")
    
    for (i in p.ops$disp.treat.pairs) {
      headings[i] <- "\n - - - " %+% .it(.attr(m.balance[[i]], "print.options")$treat_names[1L] %+% " (0) vs. " %+%
                                           .attr(m.balance[[i]], "print.options")$treat_names[2L] %+% " (1)") %+% " - - - \n"
      cat(headings[i])
      print(m.balance[[i]])
      # bal.tab_print(m.balance[[i]], p.ops)
    }
    cat(strrep(" -", round(max(nchar(headings)) / 2)), " \n\n")
  }
  
  if (isTRUE(as.logical(p.ops$multi.summary)) && is_not_null(m.balance.summary)) {
    keep.row <- {
      if (p.ops$imbalanced.only)
        rowSums(apply(m.balance.summary[grepl(".Threshold", names(m.balance.summary), fixed = TRUE)], 2L,
                      function(x) !is.na(x) & startsWith(x, "Not Balanced"))) > 0
      else
        rep.int(TRUE, nrow(m.balance.summary))
    }
    
    computed.agg.funs <- "max"
    s.keep.col <- as.logical(c(TRUE, 
                               unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS("bin")[!get_from_STATS("adj_only")]], function(s) {
                                 c(unlist(lapply(computed.agg.funs, function(af) {
                                   p.ops$un && s %in% p.ops$disp && af %in% "max"
                                 })), 
                                 if (p.ops$un && !p.ops$disp.adj && is_not_null(p.ops$thresholds[[s]])) s %in% thresholds)
                               })),
                               rep(
                                 unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS("bin")], function(s) {
                                   c(unlist(lapply(computed.agg.funs, function(af) {
                                     p.ops$disp.adj && s %in% p.ops$disp && af %in% "max"
                                   })), 
                                   if (p.ops$disp.adj && is_not_null(p.ops$thresholds[[s]])) s %in% thresholds)
                                 })),
                                 p.ops$nweights + !p.ops$disp.adj)
    ))
    names(s.keep.col) <- names(m.balance.summary)
    
    if (p.ops$disp.bal.tab) {
      cat(.ul("Balance summary across all treatment pairs") %+% "\n")
      
      if (is_null(keep.row)) cat(.it("No covariates to display.") %+% "\n")
      else if (!any(keep.row)) cat(.it("All covariates are balanced.") %+% "\n")
      else .print_data_frame(round_df_char(m.balance.summary[keep.row, s.keep.col, drop = FALSE],
                                           p.ops$digits, na_vals = "."))
      cat("\n")
    }
    
    for (s in p.ops$compute) {
      if (is_not_null(baltal[[s]])) {
        cat(.ul(sprintf("Balance tally for %s", STATS[[s]]$balance_tally_for)) %+% "\n")
        .print_data_frame(baltal[[s]])
        cat("\n")
      }
      if (is_not_null(maximbal[[s]])) {
        cat(.ul(sprintf("Variable with the greatest %s", STATS[[s]]$variable_with_the_greatest)) %+% "\n")
        .print_data_frame(round_df_char(maximbal[[s]], p.ops$digits, na_vals = "."), row.names = FALSE)
        cat("\n")
      }
    }
    
    if (is_not_null(nn)) {
      tag <- .attr(nn, "tag")
      drop.nn <- rowSums(nn) == 0
      ss.type <- .attr(nn, "ss.type")[!drop.nn]
      nn <- nn[!drop.nn, , drop = FALSE]
      if (all(c("All (ESS)", "All (Unweighted)") %in% rownames(nn)) && 
          all(check_if_zero(nn["All (ESS)", ] - nn["All (Unweighted)", ]))) {
        nn <- nn[rownames(nn) != "All (Unweighted)", , drop = FALSE]
        rownames(nn)[rownames(nn) == "All (ESS)"] <- "All"
      }
      if (all(c("Matched (ESS)", "Matched (Unweighted)") %in% rownames(nn)) && 
          all(check_if_zero(nn["Matched (ESS)", ] - nn["Matched (Unweighted)", ]))) {
        nn <- nn[rownames(nn) != "Matched (Unweighted)", , drop = FALSE]
        rownames(nn)[rownames(nn) == "Matched (ESS)"] <- "Matched"
      }
      
      cat(.ul(tag) %+% "\n")
      
      print.warning <- FALSE
      
      if (length(ss.type) > 1L && nunique.gt(ss.type[-1L], 1L)) {
        ess <- ifelse(ss.type == "ess", "*", "")
        nn <- setNames(cbind(nn, ess), c(names(nn), ""))
        print.warning <- TRUE
      }
      
      .print_data_frame(round_df_char(nn, digits = min(2L, p.ops$digits), pad = " "))
      
      if (print.warning) {
        cat(.it("* indicates effective sample size"))
      }
    }
  }
  
  invisible(x)
  
}
#' @exportS3Method NULL
bal.tab_print.bal.tab.msm <- function(x, p.ops) {
  
  call <- if (p.ops$disp.call) x$call else NULL
  msm.balance <- x[["Time.Balance"]]
  msm.balance.summary <- x[["Balance.Across.Times"]]
  
  thresholds <- setdiff(names(p.ops$thresholds), p.ops$drop.thresholds)
  baltal <- setNames(x[paste.("Balanced", thresholds)], thresholds)
  maximbal <- setNames(x[paste.("Max.Imbalance", thresholds)], thresholds)
  
  nn <- x$Observations
  
  #Printing output
  if (is_not_null(call)) {
    cat(.ul("Call") %+% "\n " %+% paste(deparse(call), collapse = "\n") %+% "\n\n")
  }
  
  if (is_not_null(p.ops$which.time)) {
    cat(.ul("Balance by Time Point") %+% "\n")
    for (i in p.ops$which.time) {
      cat("\n - - - " %+% .it("Time: " %+% as.character(i)) %+% " - - - \n")
      print(msm.balance[[i]])
    }
    cat(strrep(" -", round(nchar(sprintf("\n - - - Time: %s - - - ", i)) / 2)), "\n\n")
  }
  
  if (isTRUE(as.logical(p.ops$msm.summary)) && is_not_null(msm.balance.summary)) {
    keep.row <- {
      if (p.ops$imbalanced.only)
        rowSums(apply(msm.balance.summary[grepl(".Threshold", names(msm.balance.summary), fixed = TRUE)], 2L,
                      function(x) !is.na(x) & startsWith(x, "Not Balanced"))) > 0
      else
        rep.int(TRUE, nrow(msm.balance.summary))
    }
    
    computed.agg.funs <- "max"
    s.keep.col <- as.logical(c(TRUE, 
                               TRUE,
                               unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS()[!get_from_STATS("adj_only")]], function(s) {
                                 c(unlist(lapply(computed.agg.funs, function(af) {
                                   p.ops$un && s %in% p.ops$disp && af %in% "max"
                                 })), 
                                 if (p.ops$un && !p.ops$disp.adj && is_not_null(p.ops$thresholds[[s]])) s %in% thresholds)
                               })),
                               rep(
                                 unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS()], function(s) {
                                   c(unlist(lapply(computed.agg.funs, function(af) {
                                     p.ops$disp.adj && s %in% p.ops$disp && af %in% "max"
                                   })), 
                                   if (p.ops$disp.adj && is_not_null(p.ops$thresholds[[s]])) s %in% thresholds)
                                 })),
                                 p.ops$nweights + !p.ops$disp.adj)
    ))
    
    
    if (p.ops$disp.bal.tab) {
      cat(.ul("Balance summary across all time points") %+% "\n")
      
      if (is_null(keep.row)) cat(.it("No covariates to display.") %+% "\n")
      else if (!any(keep.row)) cat(.it("All covariates are balanced.") %+% "\n")
      else .print_data_frame(round_df_char(msm.balance.summary[keep.row, s.keep.col, drop = FALSE],
                                           p.ops$digits, na_vals = "."))
      cat("\n")
    }
    
    for (s in p.ops$compute) {
      if (is_not_null(baltal[[s]])) {
        cat(.ul(sprintf("Balance tally for %s", STATS[[s]]$balance_tally_for)) %+% "\n")
        .print_data_frame(baltal[[s]])
        cat("\n")
      }
      if (is_not_null(maximbal[[s]])) {
        cat(.ul(sprintf("Variable with the greatest %s", STATS[[s]]$variable_with_the_greatest)) %+% "\n")
        .print_data_frame(round_df_char(maximbal[[s]], p.ops$digits, na_vals = "."), row.names = FALSE)
        cat("\n")
      }
    }
    
    if (is_not_null(nn)) {
      print.warning <- FALSE
      cat(.ul(.attr(nn[[1L]], "tag")) %+% "\n")
      
      for (ti in seq_along(nn)) {
        cat(" - " %+% .it("Time " %+% as.character(ti)) %+% "\n")
        drop.nn <- rowSums(nn[[ti]]) == 0
        ss.type <- .attr(nn[[ti]], "ss.type")[!drop.nn]
        nn[[ti]] <- nn[[ti]][!drop.nn, , drop = FALSE]
        if (all(c("All (ESS)", "All (Unweighted)") %in% rownames(nn)) && 
            all(check_if_zero(nn["All (ESS)", ] - nn["All (Unweighted)", ]))) {
          nn <- nn[rownames(nn) != "All (Unweighted)", , drop = FALSE]
          rownames(nn)[rownames(nn) == "All (ESS)"] <- "All"
        }
        
        if (all(c("Matched (ESS)", "Matched (Unweighted)") %in% rownames(nn[[ti]])) && 
            all(check_if_zero(nn[[ti]]["Matched (ESS)", ] - nn[[ti]]["Matched (Unweighted)", ]))) {
          nn[[ti]] <- nn[[ti]][rownames(nn[[ti]]) != "Matched (Unweighted)", , drop = FALSE]
          rownames(nn[[ti]])[rownames(nn[[ti]]) == "Matched (ESS)"] <- "Matched"
        }
        
        if (length(ss.type) > 1L && nunique.gt(ss.type[-1L], 1L)) {
          ess <- ifelse(ss.type == "ess", "*", "")
          nn[[ti]] <- setNames(cbind(nn[[ti]], ess), c(names(nn[[ti]]), ""))
          print.warning <- TRUE
        }
        
        .print_data_frame(round_df_char(nn[[ti]], digits = min(2L, p.ops$digits), pad = " "))
      }
      
      if (print.warning) {
        cat(.it("* indicates effective sample size"))
      }
    }
  }
  
  invisible(x)
}

#' @exportS3Method NULL
bal.tab_print.bal.tab.subclass <- function(x, p.ops) {
  
  call <- if (p.ops$disp.call) x$call else NULL
  s.balance <- x$Subclass.Balance
  b.a.subclass <- x$Balance.Across.Subclass
  s.nn <- x$Observations
  
  thresholds <- setdiff(names(p.ops$thresholds), p.ops$drop.thresholds)
  baltal <- setNames(x[paste.("Balanced", thresholds, "Subclass")], thresholds)
  maximbal <- setNames(x[paste.("Max.Imbalance", thresholds, "Subclass")], thresholds)
  
  if (is_not_null(call)) {
    cat(.ul("Call") %+% "\n " %+% paste(deparse(call), collapse = "\n") %+% "\n\n")
  }
  
  #Print subclass balance
  if (p.ops$disp.bal.tab && is_not_null(p.ops$which.subclass)) {
    s.keep.col <- setNames(c(TRUE,
                             rep(unlist(lapply(p.ops$compute[p.ops$compute %nin% all_STATS()], function(s) {
                               s %in% p.ops$disp
                             })), switch(p.ops$type, bin = 2L, cont = 1L)),
                             unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS()], function(s) {
                               c(s %in% p.ops$disp,
                                 if (is_not_null(p.ops$thresholds[[s]])) s %in% thresholds)
                             }))),
                           names(s.balance[[1L]]))
    
    cat(.ul("Balance by subclass"))
    for (i in p.ops$which.subclass) {
      s.keep.row <- {
        if (p.ops$imbalanced.only)
          rowSums(apply(s.balance[[i]][grepl(".Threshold", names(s.balance), fixed = TRUE)], 2L,
                        function(x) !is.na(x) & startsWith(x, "Not Balanced"))) > 0
        else
          rep.int(TRUE, nrow(s.balance[[i]]))
      }
      
      cat("\n - - - " %+% .it("Subclass " %+% as.character(i)) %+% " - - - \n")
      
      if (is_null(s.keep.row)) cat(.it("No covariates to display.") %+% "\n")
      else if (!any(s.keep.row)) cat(.it("All covariates are balanced.") %+% "\n")
      else .print_data_frame(round_df_char(s.balance[[i]][s.keep.row, s.keep.col, drop = FALSE], p.ops$digits, na_vals = "."))
    }
    cat("\n")
  }
  
  #Print balance across subclasses
  if (p.ops$subclass.summary && is_not_null(b.a.subclass)) {
    
    if (p.ops$disp.bal.tab) {
      a.s.keep.row <- {
        if (p.ops$imbalanced.only)
          rowSums(apply(b.a.subclass[grepl(".Threshold", names(b.a.subclass), fixed = TRUE)], 2L,
                        function(x) !is.na(x) & startsWith(x, "Not Balanced"))) > 0
        else
          rep.int(TRUE, nrow(b.a.subclass))
      }
      
      a.s.keep.col <- setNames(as.logical(c(TRUE, 
                                            rep(unlist(lapply(p.ops$compute[p.ops$compute %nin% all_STATS()], function(s) {
                                              p.ops$un && s %in% p.ops$disp
                                            })), switch(p.ops$type, bin = 2L, cont = 1L)),
                                            unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS()[!get_from_STATS("adj_only")]], function(s) {
                                              c(p.ops$un && s %in% p.ops$disp,
                                                if (p.ops$un && !p.ops$disp.adj && is_not_null(p.ops$thresholds[[s]])) s %in% thresholds)
                                            })),
                                            rep(c(rep(unlist(lapply(p.ops$compute[p.ops$compute %nin% all_STATS()], function(s) {
                                              p.ops$disp.adj && s %in% p.ops$disp
                                            })), 2L),
                                            unlist(lapply(p.ops$compute[p.ops$compute %in% all_STATS()], function(s) {
                                              c(p.ops$disp.adj && s %in% p.ops$disp,
                                                if (p.ops$disp.adj && is_not_null(p.ops$thresholds[[s]])) s %in% thresholds)
                                            }))
                                            ), p.ops$disp.adj)
      )),
      names(b.a.subclass))
      
      cat(.ul("Balance measures across subclasses") %+% "\n")
      
      if (is_null(a.s.keep.row)) cat(.it("No covariates to display.") %+% "\n")
      else if (!any(a.s.keep.row)) cat(.it("All covariates are balanced.") %+% "\n")
      else .print_data_frame(round_df_char(b.a.subclass[a.s.keep.row, a.s.keep.col, drop = FALSE],
                                           p.ops$digits, na_vals = "."))
      
      cat("\n")
    }
    
    for (s in p.ops$compute) {
      if (is_not_null(baltal[[s]])) {
        cat(.ul(sprintf("Balance tally for %s across subclasses", STATS[[s]]$balance_tally_for)) %+% "\n")
        .print_data_frame(baltal[[s]])
        cat("\n")
      }
      if (is_not_null(maximbal[[s]])) {
        cat(.ul(sprintf("Variable with the greatest %s across subclasses", STATS[[s]]$variable_with_the_greatest)) %+% "\n")
        .print_data_frame(round_df_char(maximbal[[s]], p.ops$digits, na_vals = "."), row.names = FALSE)
        cat("\n")
      }
    }
    
    if (is_not_null(s.nn)) {
      cat(.ul(.attr(s.nn, "tag")) %+% "\n")
      .print_data_frame(round_df_char(s.nn, digits = min(2L, p.ops$digits), pad = " "))
    }
  }
  
  invisible(x)
}

#Process arguments
print_process <- function(x, ...) {
  UseMethod("print_process")
}
#' @exportS3Method NULL
print_process.bal.tab.cluster <- function(x, which.cluster, cluster.summary, cluster.fun, ...) {
  c.balance <- x$Cluster.Balance
  p.ops <- .attr(x, "print.options")
  
  if (!missing(cluster.summary)) {
    .chk_flag(cluster.summary)
    if (p.ops$quick && !p.ops$cluster.summary && cluster.summary) {
      .wrn("{.arg cluster.summary} cannot be set to {.val {TRUE}} if {.code quick = TRUE} in the original call to {.fun bal.tab}")
    }
    else {
      p.ops$cluster.summary <- cluster.summary
    }
  }
  
  if (!missing(which.cluster)) {
    which.cluster_deparse <- deparse1(substitute(which.cluster))
    
    if (which.cluster_deparse == ".none") which.cluster <- NA
    else if (which.cluster_deparse == ".all") which.cluster <- NULL
    
    p.ops$which.cluster <- which.cluster
  }
  
  computed.cluster.funs <- {
    if (!p.ops$quick || is_null(p.ops$cluster.fun)) c("min", "mean", "max")
    else p.ops$cluster.fun
  }
  
  if (missing(cluster.fun) || is_null(cluster.fun)) {
    cluster.fun <- {
      if (p.ops$abs) c("mean", "max")
      else c("min", "mean", "max")
    }
  }
  else if (!is.character(cluster.fun) || !all(cluster.fun %pin% computed.cluster.funs)) {
    .err("{.arg cluster.fun} must be {.or {.val {computed.cluster.funs}}}")
  }
  
  cluster.fun <- match_arg(tolower(cluster.fun), computed.cluster.funs,
                           several.ok = TRUE)
  
  #Checks and Adjustments
  if (is_null(p.ops$which.cluster)) {
    which.cluster <- seq_along(c.balance)
  }
  else if (anyNA(p.ops$which.cluster)) {
    which.cluster <- integer()
  }
  else if (is.numeric(p.ops$which.cluster)) {
    which.cluster <- intersect(seq_along(c.balance), p.ops$which.cluster)
    if (is_null(which.cluster)) {
      .wrn("no indices in {.arg which.cluster} are cluster indices. Displaying all clusters instead")
      which.cluster <- seq_along(c.balance)
    }
  }
  else if (is.character(p.ops$which.cluster)) {
    which.cluster <- seq_along(c.balance)[names(c.balance) %in% p.ops$which.cluster]
    if (is_null(which.cluster)) {
      .wrn("no names in {.arg which.cluster} are cluster names. Displaying all clusters instead")
      which.cluster <- seq_along(c.balance)
    }
  }
  else {
    .wrn("the argument to {.arg which.cluster} must be {.val {quote(.all)}}, {.val {quote(.none)}}, or a vector of cluster indices or cluster names. Displaying all clusters instead")
    which.cluster <- seq_along(c.balance)
  }
  
  list(cluster.summary = p.ops$cluster.summary,
       cluster.fun = cluster.fun,
       which.cluster = which.cluster,
       computed.cluster.funs = computed.cluster.funs)
}
#' @exportS3Method NULL
print_process.bal.tab.imp <- function(x, which.imp, imp.summary, imp.fun, ...) {
  i.balance <- x[["Imputation.Balance"]]
  p.ops <- .attr(x, "print.options")
  
  if (!missing(imp.summary)) {
    .chk_flag(imp.summary)
    if (p.ops$quick && !p.ops$imp.summary && imp.summary) {
      .wrn("{.arg imp.summary} cannot be set to {.val {TRUE}} if {.code quick = TRUE} in the original call to {.fun bal.tab}")
    }
    else {
      p.ops$imp.summary <- imp.summary
    }
  }
  
  if (!missing(which.imp)) {
    which.imp_deparse <- deparse1(substitute(which.imp))
    
    if (which.imp_deparse == ".none") which.imp <- NA
    else if (which.imp_deparse == ".all") which.imp <- NULL
    
    p.ops$which.imp <- which.imp
  }
  
  computed.imp.funs <- {
    if (!p.ops$quick || is_null(p.ops$imp.fun)) c("min", "mean", "max")
    else p.ops$imp.fun
  }
  
  if (missing(imp.fun) || is_null(imp.fun)) {
    imp.fun <- {
      if (p.ops$abs) c("mean", "max")
      else c("min", "mean", "max")
    }
  }
  else if (!is.character(imp.fun) || !all(imp.fun %pin% computed.imp.funs)) {
    .err("{.arg imp.fun} must be {.or {.val {computed.imp.funs}}}")
  }
  
  imp.fun <- match_arg(tolower(imp.fun), computed.imp.funs, several.ok = TRUE)
  
  #Checks and Adjustments
  if (is_null(p.ops$which.imp)) {
    which.imp <- seq_along(i.balance)
  }
  else if (anyNA(p.ops$which.imp)) {
    which.imp <- integer()
  }
  else if (is.numeric(p.ops$which.imp)) {
    which.imp <- intersect(seq_along(i.balance), p.ops$which.imp)
    if (is_null(which.imp)) {
      .wrn("no numbers in {.arg which.imp} are imputation numbers. No imputations will be displayed")
      which.imp <- integer()
    }
  }
  else {
    .wrn("the argument to {.arg which.imp} must be {.val {quote(.all)}}, {.val {quote(.none)}}, or a vector of imputation numbers")
    which.imp <- integer()
  }
  
  list(imp.summary = p.ops$imp.summary,
       imp.fun = imp.fun,
       which.imp = which.imp,
       computed.imp.funs = computed.imp.funs)
}
#' @exportS3Method NULL
print_process.bal.tab.multi <- function(x, which.treat, multi.summary, ...) {
  
  m.balance <- x[["Pair.Balance"]]
  
  p.ops <- .attr(x, "print.options")
  
  if (!missing(multi.summary)) {
    .chk_flag(multi.summary)
    if (p.ops$quick && !p.ops$multi.summary && multi.summary) {
      .wrn("{.arg multi.summary} cannot be set to {.val {TRUE}} if {.code quick = TRUE} in the original call to {.fun bal.tab}")
    }
    else {
      p.ops$multi.summary <- multi.summary
    }
  }
  
  if (!missing(which.treat)) {
    which.treat_deparse <- deparse1(substitute(which.treat))
    
    if (which.treat_deparse == ".none") which.treat <- NA
    else if (which.treat_deparse == ".all") which.treat <- NULL
    
    p.ops$which.treat <- which.treat
  }
  
  #Checks and Adjustments
  if (is_null(p.ops$which.treat)) {
    which.treat <- p.ops$treat_vals_multi
  }
  else if (anyNA(p.ops$which.treat)) {
    which.treat <- character()
  }
  else if (!is.character(p.ops$which.treat) && !is.numeric(p.ops$which.treat)) {
    .wrn("the argument to {.arg which.treat} must be {.val {quote(.all)}}, {.val {quote(.none)}}, or a vector of treatment names or indices. No treatment pairs will be displayed")
    which.treat <- character()
  }
  else {
    if (length(p.ops$treat_vals_multi) == 2L) {
      p.ops$which.treat <- as.character(p.ops$which.treat)
    }
    
    if (is.numeric(p.ops$which.treat)) {
      which.treat <- p.ops$treat_vals_multi[seq_along(p.ops$treat_vals_multi) %in% p.ops$which.treat]
      if (is_null(which.treat)) {
        .wrn("no numbers in {.arg which.treat} correspond to treatment values. No treatment pairs will be displayed")
        which.treat <- character()
      }
    }
    else if (is.character(p.ops$which.treat)) {
      which.treat <- p.ops$treat_vals_multi[p.ops$treat_vals_multi %in% p.ops$which.treat]
      if (is_null(which.treat)) {
        .wrn("no names in {.arg which.treat} correspond to treatment values. No treatment pairs will be displayed")
        which.treat <- character()
      }
    }
  }
  
  if (is_null(which.treat)) {
    disp.treat.pairs <- character()
  }
  else if (p.ops$pairwise) {
    if (length(which.treat) == 1L) {
      disp.treat.pairs <- names(m.balance)[vapply(m.balance, function(z) {
        treat_names <- .attr(z, "print.options")$treat_names
        any(p.ops$treat_vals_multi[treat_names] == which.treat)
      }, logical(1L))]
    }
    else {
      disp.treat.pairs <- names(m.balance)[vapply(m.balance, function(z) {
        treat_names <- .attr(z, "print.options")$treat_names
        all(p.ops$treat_vals_multi[treat_names] %in% which.treat)
      }, logical(1L))]
    }
  }
  else {
    if (length(which.treat) == 1L) {
      disp.treat.pairs <- names(m.balance)[vapply(m.balance, function(z) {
        treat_names <- .attr(z, "print.options")$treat_names
        any(p.ops$treat_vals_multi[setdiff(treat_names, "All")] == which.treat)
      }, logical(1L))]
    }
    else {
      disp.treat.pairs <- names(m.balance)[vapply(m.balance, function(z) {
        treat_names <- .attr(z, "print.options")$treat_names
        all(p.ops$treat_vals_multi[setdiff(treat_names, "All")] %in% which.treat)
      }, logical(1L))]
    }
  }
  
  list(disp.treat.pairs = disp.treat.pairs,
       multi.summary = p.ops$multi.summary,
       pairwise = p.ops$pairwise)
}
#' @exportS3Method NULL
print_process.bal.tab.msm <- function(x, which.time, msm.summary, ...) {
  
  msm.balance <- x[["Time.Balance"]]
  
  p.ops <- .attr(x, "print.options")
  
  if (!missing(msm.summary)) {
    .chk_flag(msm.summary)
    if (p.ops$quick && !p.ops$msm.summary && msm.summary) {
      .wrn("{.arg msm.summary} cannot be set to {.val {TRUE}} if {.code quick = TRUE} in the original call to {.fun bal.tab}")
    }
    else {
      p.ops$msm.summary <- msm.summary
    }
  }
  
  if (!missing(which.time)) {
    which.time_deparse <- deparse1(substitute(which.time))
    
    if (which.time_deparse == ".none") which.time <- NA
    else if (which.time_deparse == ".all") which.time <- NULL
    
    p.ops$which.time <- which.time
  }
  
  #Checks and Adjustments
  if (is_null(p.ops$which.time)) {
    which.time <- seq_along(msm.balance)
  }
  else if (anyNA(p.ops$which.time)) {
    which.time <- integer()
  }
  else if (is.numeric(p.ops$which.time)) {
    which.time <- seq_along(msm.balance)[seq_along(msm.balance) %in% p.ops$which.time]
    if (is_null(which.time)) {
      .wrn("no numbers in {.arg which.time} are treatment time points. No time points will be displayed")
      which.time <- integer()
    }
  }
  else if (is.character(p.ops$which.time)) {
    which.time <- seq_along(msm.balance)[names(msm.balance) %in% p.ops$which.time]
    if (is_null(which.time)) {
      .wrn("no names in {.arg which.time} are treatment names. No time points will be displayed")
      which.time <- integer()
    }
  }
  else {
    .wrn("the argument to {.arg which.time} must be {.val {quote(.all)}}, {.val {quote(.none)}}, or a vector of time point numbers. No time points will be displayed")
    which.time <- integer()
  }
  
  list(msm.summary = p.ops$msm.summary,
       which.time = which.time)
}
#' @exportS3Method NULL
print_process.bal.tab <- function(x, imbalanced.only, un, disp.bal.tab, disp.call, stats,
                                  disp.thresholds, disp, digits = max(3, getOption("digits") - 3),
                                  ...) {
  p.ops <- .attr(x, "print.options")
  
  drop.thresholds <- c()
  
  #Adjustments to print options
  if (!missing(un) && p.ops$disp.adj) {
    .chk_flag(un)
    if (p.ops$quick && !p.ops$un && un) {
      .wrn("{.arg un} cannot be set to {.val {TRUE}} if {.code quick = TRUE} in the original call to {.fun bal.tab}")
    }
    else {
      p.ops$un <- un
    }
  }
  
  if (!missing(disp)) {
    .chk_character(disp)
    allowable.disp <- c("means", "sds", all_STATS(p.ops$type))
    
    if (!all(disp %in% allowable.disp)) {
      .err("{.val {setdiff(disp, allowable.disp)}} {?is/are} not allowed in {.arg disp}")
    }
    
    if (all(disp %in% p.ops$compute)) {
      .wrn("{.arg disp} cannot include {.or {.val {setdiff(disp, p.ops$compute)}}} if {.code quick = TRUE} in the original call to {.fun bal.tab}")
    }
    else{
      p.ops$disp <- disp
    }
  }
  
  if (is_not_null(...get("disp.means"))) {
    .chk_flag(...get("disp.means"), "disp.means")
    
    if ("means" %in% p.ops$compute || !...get("disp.means")) {
      p.ops$disp <- unique(c(p.ops$disp, "means"[...get("disp.means")]))
    }
    else {
      .wrn("{.arg disp.means} cannot be set to {.val {TRUE}} if {.code quick = TRUE} in the original call to {.fun bal.tab}")
    }
  }
  
  if (is_not_null(...get("disp.sds"))) {
    .chk_flag(...get("disp.sds"), "disp.sds")
    
    if ("sds" %nin% p.ops$compute && ...get("disp.sds")) {
      .wrn("{.arg disp.sds} cannot be set to {.val {TRUE}} if {.code quick = TRUE} in the original call to {.fun bal.tab}")
    }
    else {
      p.ops$disp <- unique(c(p.ops$disp, "sds"[...get("disp.sds")]))
    }
  }
  
  if (!missing(stats)) {
    .chk_character(stats)
    stats <- match_arg(stats, all_STATS(p.ops$type), several.ok = TRUE)
    stats_in_p.ops <- stats %in% p.ops$compute
    
    if (!all(stats_in_p.ops)) {
      .err("{.arg stats} cannot contain {.or {.val {stats[!stats_in_p.ops]}}} when {?it/they} {?was/were} not requested in the original call to {.fun bal.tab}")
    }
    
    p.ops$disp <- unique(c(p.ops$disp[p.ops$disp %nin% all_STATS()], stats))
  }
  
  for (s in all_STATS(p.ops$type)) {
    if (is_not_null(...get(STATS[[s]]$disp_stat))) {
      .chk_flag(...get(STATS[[s]]$disp_stat), STATS[[s]]$disp_stat)
      if (s %nin% p.ops$compute && isTRUE(...get(STATS[[s]]$disp_stat))) {
        .wrn("{.arg {STATS[[s]]$disp_stat}} cannot be set to {.val {TRUE}} if {.code quick = TRUE} in the original call to {.fun bal.tab}")
      }
      else {
        p.ops$disp <- unique(c(p.ops$disp, s))
      }
    }
  }
  
  for (s in p.ops$compute[p.ops$compute %in% all_STATS(p.ops$type)]) {
    if (STATS[[s]]$threshold %in% ...names()) {
      temp.thresh <- ...get(STATS[[s]]$threshold)
      if (is_not_null(temp.thresh) &&
          (!is.numeric(temp.thresh) || length(temp.thresh) != 1L ||
           is_null(p.ops[["thresholds"]][[s]]) ||
           p.ops[["thresholds"]][[s]] != temp.thresh)) {
        .err("{.arg {STATS[[s]]$threshold}} must be {.val {list(NULL)}} or left unspecified")
      }
      
      if (is_null(temp.thresh)) {
        drop.thresholds <- c(drop.thresholds, s)
      }
    }
    
    if (s %nin% p.ops$disp) {
      drop.thresholds <- c(drop.thresholds, s)
    }
  }
  
  if (!missing(disp.thresholds)) {
    .chk_logical(disp.thresholds)
    .chk_not_any_na(disp.thresholds)
    
    if (is_null(names(disp.thresholds))) {
      if (length(disp.thresholds) > length(p.ops[["thresholds"]])) {
        .err("more entries were given to {.arg disp.thresholds} than there are thresholds in the {.cls bal.tab} object")
      }
      
      if (length(disp.thresholds) == 1L) {
        disp.thresholds <- rep_with(disp.thresholds, p.ops[["thresholds"]])
      }
      
      names(disp.thresholds) <- names(p.ops[["thresholds"]])[seq_along(disp.thresholds)]
    }
    
    if (!all(names(disp.thresholds) %pin% names(p.ops[["thresholds"]]))) {
      .wrn('{.val {names(disp.thresholds)[!names(disp.thresholds) %pin% names(p.ops[["thresholds"]])]}} {?is/are} not available in thresholds and will be ignored')
      disp.thresholds <- disp.thresholds[names(disp.thresholds) %pin% names(p.ops[["thresholds"]])]
    }
    
    names(disp.thresholds) <- match_arg(names(disp.thresholds), names(p.ops[["thresholds"]]), several.ok = TRUE)
    
    for (i in names(disp.thresholds)) {
      if (!disp.thresholds[i]) {
        drop.thresholds <- c(drop.thresholds, i)
      }
    }
  }
  
  if (!missing(disp.bal.tab)) {
    .chk_flag(disp.bal.tab)
    p.ops$disp.bal.tab <- disp.bal.tab
  }
  
  if (p.ops$disp.bal.tab) {
    if (!missing(imbalanced.only)) {
      .chk_flag(imbalanced.only)
      p.ops$imbalanced.only <- imbalanced.only
    }
    
    if (p.ops$imbalanced.only && is_null(p.ops$thresholds)) {
      .wrn("a threshold must be specified if {.code imbalanced.only = TRUE}. Displaying all covariates")
      p.ops$imbalanced.only <- FALSE
    }
  }
  else {
    p.ops$imbalanced.only <- FALSE
  }
  
  if (!missing(disp.call)) {
    .chk_flag(disp.call)
    if (disp.call && is_null(x$call)) {
      .wrn("{.arg disp.call} cannot be set to {.val {TRUE}} if the input object does not have a {.field call} component")
    }
    else {
      p.ops$disp.call <- disp.call
    }
  }
  
  list(un = p.ops$un,
       disp = p.ops$disp,
       compute = p.ops$compute,
       drop.thresholds = drop.thresholds,
       disp.bal.tab = p.ops$disp.bal.tab,
       imbalanced.only = p.ops$imbalanced.only,
       digits = digits,
       disp.adj = p.ops$disp.adj,
       thresholds = p.ops$thresholds,
       type = p.ops$type,
       nweights = p.ops$nweights,
       disp.call = p.ops$disp.call)
}
#' @exportS3Method NULL
print_process.bal.tab.subclass <- function(x, imbalanced.only, un, disp.bal.tab, disp.call, stats,
                                           disp.thresholds, disp, digits = max(3, getOption("digits") - 3),
                                           which.subclass, subclass.summary, ...) {
  s.balance <- x$Subclass.Balance
  p.ops <- .attr(x, "print.options")
  
  drop.thresholds <- c()
  
  #Adjustments to print options
  if (!missing(un) && p.ops$disp.adj) {
    .chk_flag(un)
    if (p.ops$quick && !p.ops$un && un) {
      .wrn("{.arg un} cannot be set to {.val {TRUE}} if {.code quick = TRUE} in the original call to {.fun bal.tab}")
    }
    else {
      p.ops$un <- un
    }
  }
  
  if (!missing(disp)) {
    .chk_character(disp)
    allowable.disp <- c("means", "sds", all_STATS(p.ops$type))
    
    if (!all(disp %in% allowable.disp)) {
      .err("{.val {setdiff(disp, allowable.disp)}} {?is/are} not allowed in {.arg disp}")
    }
    
    if (all(disp %in% p.ops$compute)) {
      p.ops$disp <- disp
    }
    else {
      .wrn("{.arg disp} cannot include {.or {.val {setdiff(disp, p.ops$compute)}}} if {.code quick = TRUE} in the original call to {.fun bal.tab}")
    }
  }
  
  if (is_not_null(...get("disp.means"))) {
    .chk_flag(...get("disp.means"), "disp.means")
    
    if ("means" %in% p.ops$compute || !...get("disp.means")) {
      p.ops$disp <- unique(c(p.ops$disp, "means"[...get("disp.means")]))
    }
    else {
      .wrn("{.arg disp.means} cannot be set to {.val {TRUE}} if {.code quick = TRUE} in the original call to {.fun bal.tab}")
    }
  }
  
  if (is_not_null(...get("disp.sds"))) {
    .chk_flag(...get("disp.sds"), "disp.sds")
    
    if ("sds" %nin% p.ops$compute && ...get("disp.sds")) {
      .wrn("{.arg disp.sds} cannot be set to {.val {TRUE}} if {.code quick = TRUE} in the original call to {.fun bal.tab}")
    }
    else {
      p.ops$disp <- unique(c(p.ops$disp, "sds"[...get("disp.sds")]))
    }
  }
  
  if (!missing(stats)) {
    .chk_character(stats)
    stats <- match_arg(stats, all_STATS(p.ops$type), several.ok = TRUE)
    stats_in_p.ops <- stats %in% p.ops$compute
    
    if (!all(stats_in_p.ops)) {
      .err("{.arg stats} cannot contain {.or {.val {stats[!stats_in_p.ops]}}} when {?it/they} {?was/were} not requested in the original call to {.fun bal.tab}")
    }
    
    p.ops$disp <- unique(c(p.ops$disp[p.ops$disp %nin% all_STATS()], stats))
  }
  
  for (s in all_STATS(p.ops$type)) {
    if (is_not_null(...get(STATS[[s]]$disp_stat))) {
      .chk_flag(...get(STATS[[s]]$disp_stat), STATS[[s]]$disp_stat)
      
      if (s %nin% p.ops$compute && isTRUE(...get(STATS[[s]]$disp_stat))) {
        .wrn("{.arg {STATS[[s]]$disp_stat}} cannot be set to {.val {TRUE}} if {.code quick = TRUE} in the original call to {.fun bal.tab}")
      }
      else {
        p.ops$disp <- unique(c(p.ops$disp, s))
      }
    }
  }
  
  for (s in p.ops$compute[p.ops$compute %in% all_STATS(p.ops$type)]) {
    if (STATS[[s]]$threshold %in% ...names()) {
      temp.thresh <- ...get(STATS[[s]]$threshold)
      
      if (is_not_null(temp.thresh) &&
          (!is.numeric(temp.thresh) || length(temp.thresh) != 1L ||
           is_null(p.ops[["thresholds"]][[s]]) ||
           p.ops[["thresholds"]][[s]] != temp.thresh)) {
        .err("{.arg {STATS[[s]]$threshold}} must be {.val {list(NULL)}} or left unspecified")
      }
      
      if (is_null(temp.thresh)) {
        drop.thresholds <- c(drop.thresholds, s)
      }
    }
    
    if (s %nin% p.ops$disp) {
      drop.thresholds <- c(drop.thresholds, s)
    }
  }
  
  if (!missing(disp.thresholds)) {
    .chk_logical(disp.thresholds)
    .chk_not_any_na(disp.thresholds)
    
    if (is_null(names(disp.thresholds))) {
      if (length(disp.thresholds) > length(p.ops[["thresholds"]])) {
        .err("more entries were given to {.arg disp.thresholds} than there are thresholds in the {.cls bal.tab} object")
      }
      
      if (length(disp.thresholds) == 1L) {
        disp.thresholds <- rep_with(disp.thresholds, p.ops[["thresholds"]])
      }
      
      names(disp.thresholds) <- names(p.ops[["thresholds"]])[seq_along(disp.thresholds)]
    }
    
    if (!all(names(disp.thresholds) %pin% names(p.ops[["thresholds"]]))) {
      .wrn('{names(disp.thresholds)[!names(disp.thresholds) %pin% names(p.ops[["thresholds"]])]} {?is/are} not available in thresholds and will be ignored')
      disp.thresholds <- disp.thresholds[names(disp.thresholds) %pin% names(p.ops[["thresholds"]])]
    }
    
    names(disp.thresholds) <- match_arg(names(disp.thresholds), names(p.ops[["thresholds"]]), several.ok = TRUE)
    for (x in names(disp.thresholds)) {
      if (!disp.thresholds[x]) {
        drop.thresholds <- c(drop.thresholds, x)
      }
    }
  }
  
  if (!missing(disp.bal.tab)) {
    .chk_flag(disp.bal.tab)
    p.ops$disp.bal.tab <- disp.bal.tab
  }
  
  if (p.ops$disp.bal.tab) {
    if (!missing(imbalanced.only)) {
      .chk_flag(imbalanced.only)
      p.ops$imbalanced.only <- imbalanced.only
    }
    
    if (p.ops$imbalanced.only && is_null(p.ops$thresholds)) {
      .wrn("a threshold must be specified if {.code imbalanced.only = TRUE}. Displaying all covariates")
      p.ops$imbalanced.only <- FALSE
    }
  }
  else {
    p.ops$imbalanced.only <- FALSE
  }
  
  if (!missing(disp.call)) {
    .chk_flag(disp.call)
    
    if (disp.call && is_null(x$call)) {
      .wrn("{.arg disp.call} cannot be set to {.val {TRUE}} if the input object does not have a {.field call} component")
    }
    else {
      p.ops$disp.call <- disp.call
    }
  }
  
  if (!missing(subclass.summary)) {
    .chk_flag(subclass.summary)
    if (p.ops$quick && !p.ops$subclass.summary && subclass.summary) {
      .wrn("`subclass.summary` cannot be set to `TRUE` if `quick = TRUE` in the original call to `bal.tab()`")
    }
    else {
      p.ops$subclass.summary <- subclass.summary
    }
  }
  
  if (!missing(which.subclass)) {
    p.ops$which.subclass <- which.subclass
  }
  else if (is_not_null(...get("disp.subclass"))) {
    p.ops$which.subclass <- if (isTRUE(...get("disp.subclass"))) NULL else NA
  }
  
  #Checks and Adjustments
  if (is_null(p.ops$which.subclass)) {
    which.subclass <- seq_along(s.balance)
  }
  else if (anyNA(p.ops$which.subclass)) {
    which.subclass <- integer()
  }
  else if (is.numeric(p.ops$which.subclass)) {
    which.subclass <- intersect(seq_along(s.balance), p.ops$which.subclass)
    if (is_null(which.subclass)) {
      .wrn("no values supplied {.arg which.subclass} are subclass indices. No subclasses will be displayed")
      which.subclass <- NA
    }
  }
  else {
    .wrn("the argument to {.arg which.subclass} must be {.val {quote(.all)}}, {.val {quote(.none)}}, or a vector of subclass indices. No subclasses will be displayed")
    which.subclass <- NA
  }
  
  list(un = p.ops$un,
       disp = p.ops$disp,
       compute = p.ops$compute,
       drop.thresholds = drop.thresholds,
       disp.bal.tab = p.ops$disp.bal.tab,
       imbalanced.only = p.ops$imbalanced.only,
       digits = digits,
       disp.adj = p.ops$disp.adj,
       thresholds = p.ops$thresholds,
       type = p.ops$type,
       subclass.summary = p.ops$subclass.summary,
       which.subclass = which.subclass,
       disp.call = p.ops$disp.call)
}

#Alternative to print.data.frame() that only prints non-length 0 data.frames
.print_data_frame <- function(x, ...) {
  if (is_not_null(x) && NROW(x) > 0L && NCOL(x) > 0L) {
    print.data.frame(x, ...)
  }
  else {
    invisible(x)
  }
}

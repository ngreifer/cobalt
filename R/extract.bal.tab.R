#' @title Extract a Balance Table for Further Use
#' @name extract.bal.tab
#'
#' @description
#' `as.data.frame()` extracts the balance statistics computed by [bal.tab()] as a tidy data frame, one row per covariate, sample, and statistic. `format()` returns the balance table exactly as [print.bal.tab()] displays it, as a data frame of formatted strings, ready to be passed to a table-rendering function.
#'
#' Together they are meant to remove the need to pick apart a `bal.tab` object by hand when reporting balance in a document.
#'
#' @param x a `bal.tab` object; the output of a call to [bal.tab()].
#' @param row.names,optional ignored; present for consistency with the [as.data.frame()] generic.
#' @param ... arguments passed to [print.bal.tab()] to control which statistics, samples, and covariates are included, e.g. `stats`, `disp`, `un`, `imbalanced.only`, or `disp.thresholds`. The `.all` and `.none` shorthands are accepted as they are by `print()`. Arguments that were not computed in the original call to `bal.tab()` cannot be requested here, for the same reason they cannot be requested in `print()`.
#' @param wide `logical`; for `as.data.frame()`, whether to return the table in the layout `print()` uses, with one column per sample and statistic, rather than the default tidy layout. Default is `FALSE`.
#' @param digits for `format()`, the number of significant digits to display. Default is the same as for `print()`.
#' @param component for `format()`, which table to return: `"balance"` (the default) for the balance table, or `"observations"` for the sample size table.
#'
#' @returns
#' `as.data.frame()` returns a data frame with one row per covariate, sample, and statistic, and these columns:
#' \describe{
#'   \item{`variable`}{the covariate name, as it appears in the row names of the balance table.}
#'   \item{`type`}{the covariate type: `"Binary"`, `"Contin."`, or `"Distance"`.}
#'   \item{`sample`}{`"Unadjusted"`, or the name of the set of weights.}
#'   \item{`stat`}{the name of the statistic, using the same names as `bal.tab()`'s `stats` argument (e.g., `"mean.diffs"`), with `"mean"` and `"sd"` for the distribution summary statistics.}
#'   \item{`group`}{for `"mean"` and `"sd"`, the name of the treatment group the value describes, as [bal.tab()] names it: the two group names for a binary treatment, the level for a multi-category one, `"Uncensored"` or `"Full"` for a censoring indicator, and `"All"` for a continuous treatment, which has no groups, or for the full sample when `pairwise = FALSE`. `NA` for a statistic that contrasts two groups, which belongs to neither of them.}
#'   \item{`estimate`}{the value of the statistic.}
#'   \item{`threshold`}{the balance verdict, e.g. `"Balanced, <0.1"`, when a threshold was requested for that statistic; `NA` otherwise.}
#'   \item{`threshold.value`}{the numeric threshold; `NA` when none was requested.}
#' }
#' When the data are segmented -- by cluster, imputation, treatment pair, time point, or subclass -- one further column per level of segmentation identifies it. Segmentation is always a column, never a nested list, so the result is a single rectangle whatever the shape of the input.
#'
#' A multi-category treatment is reported one pair of groups at a time, but a mean or a standard deviation belongs to a group rather than to a comparison, and is the same in every pair that group appears in. Such a row therefore appears once, with `pair` set to `NA`; only the statistics that contrast two groups carry a `pair`. The same applies to the full sample's own means when `pairwise = FALSE`, which would otherwise be repeated against every group.
#'
#' With `wide = TRUE`, the columns are those `print()` displays, with the covariate names moved from the row names into a `variable` column and any segmentation columns retained.
#'
#' `format()` returns a data frame of character vectors with the covariate names as row names, formatted exactly as `print()` displays them: rounded to `digits`, padded to a common number of decimal places, with `NA` shown as `"."`.
#'
#' @details
#' Both functions accept every argument `print()` accepts, and resolve them the same way, so the same warnings are raised when a requested value was not computed because `quick = TRUE` in the original call to `bal.tab()`.
#'
#' `as.data.frame()` returns the balance statistics themselves, from each innermost balance table. It does not return the summaries across clusters, imputations, treatment pairs, or time points, which are aggregates of those statistics; `format()` returns the summary when that is what `print()` displays.
#'
#' @seealso
#' * [bal.tab()]
#' * [print.bal.tab()]
#' * [love.plot()] for a graphical alternative
#'
#' @examples
#' data("lalonde", package = "cobalt")
#'
#' b <- bal.tab(treat ~ age + educ + race + re74, data = lalonde,
#'              s.d.denom = "pooled", stats = c("m", "ks"),
#'              thresholds = c(m = .1), un = TRUE)
#'
#' #Tidy: one row per covariate, sample, and statistic
#' head(as.data.frame(b))
#'
#' #The layout print() shows
#' as.data.frame(b, wide = TRUE)
#'
#' #Ready for knitr::kable() or any other table renderer
#' format(b)
#'
#' format(b, component = "observations")
#'
#' #print()'s arguments work here too
#' as.data.frame(b, stats = "ks", un = FALSE)
#'
NULL

#What a leaf's groups are called, indexed by the label its balance table's columns carry.
#`print()` labels a binary treatment's columns positionally -- `M.0`, `M.1` -- which says
#nothing about which group is which, and a multi-category object reuses those same two
#positions for a different pair of groups in every table. Returns NULL for a treatment
#with no groups at all.
.leaf_groups <- function(x) {
  p.ops <- .attr(x, "print.options")

  treat_names <- p.ops[["treat_names"]]

  if (is_null(treat_names)) {
    return(NULL)
  }

  setNames(unname(treat_names), p.ops[["group.labels"]] %or% c("0", "1"))
}

#What a per-group estimate's group is called. A continuous treatment has no groups, so its
#means and standard deviations describe the sample as a whole -- the same thing
#`pairwise = FALSE` calls `All`.
.ALL_GROUP <- "All"

.group_name <- function(label, groups) {
  if (is.na(label) || is_null(groups)) {
    return(.ALL_GROUP)
  }

  out <- unname(groups[label])

  if (is.na(out)) label else out
}

#The innermost balance tables, each with the segmentation that identifies it and the names
#of the groups its columns are about. A subclassified object is the one shape whose leaves
#are bare data frames rather than `bal.tab` objects, so it is reached explicitly.
.bal.tab_leaves <- function(x, keys = list(), groups = NULL) {
  #Each level of a segmented object may name its groups differently -- every pair of a
  #multi-category treatment does -- so the innermost naming wins.
  groups <- .leaf_groups(x) %or% groups

  if (inherits(x, c("bal.tab.bin", "bal.tab.cont"))) {
    return(list(list(keys = keys,
                     table = x[["Balance"]],
                     groups = groups)))
  }

  if (inherits(x, "bal.tab.subclass")) {
    return(lapply(names(x[["Subclass.Balance"]]), function(i) {
      list(keys = c(keys, list(subclass = i)),
           table = x[["Subclass.Balance"]][[i]],
           groups = groups)
    }))
  }

  level <- switch(class(x)[1L],
                  "bal.tab.cluster" = "cluster",
                  "bal.tab.imp" = "imp",
                  "bal.tab.multi" = "pair",
                  "bal.tab.msm" = "time",
                  "nesting")

  children <- x[[which(endsWith(names(x), ".Balance"))]]

  #Time points are named after the treatment at each one, which need not be distinct,
  #so fall back to the index -- which is how `print()` labels them.
  nms <- {
    if (anyDuplicated(names(children)) > 0L) NULL
    else names(children)
  } %or% as.character(seq_along(children))

  lapply(seq_along(children), function(i) {
    .bal.tab_leaves(children[[i]], c(keys, setNames(list(nms[i]), level)), groups)
  }) |>
    unlist(recursive = FALSE)
}

#' @rdname extract.bal.tab
#' @exportS3Method as.data.frame bal.tab
as.data.frame.bal.tab <- function(x, row.names = NULL, optional = FALSE, ...,
                                  wide = FALSE) {
  .call <- .rewrite_all_none(match.call(expand.dots = TRUE))

  if (is_not_null(.call)) {
    return(eval.parent(.call))
  }

  arg::arg_flag(wide)

  A <- .display_args(environment(), list(...))
  A[c("row.names", "optional", "wide")] <- NULL

  p.ops <- .resolve_p.ops(x, A)
  thresholds <- setdiff(names(p.ops[["thresholds"]]), p.ops[["drop.thresholds"]])

  leaves <- .bal.tab_leaves(x)

  out <- lapply(leaves, function(leaf) {
    tab <- leaf[["table"]]

    #The spec describes every column the object's options could produce; a given
    #table holds a subset of them (a subclass block has no unadjusted columns).
    spec <- .p.ops_col_spec(p.ops)
    spec <- spec[spec[["name"]] %in% names(tab), , drop = FALSE]
    spec <- spec[.keep_bal_cols(spec, p.ops, thresholds), , drop = FALSE]

    keep.row <- .keep_bal_rows(tab, p.ops)

    d <- {
      if (wide) .bal.tab_wide(tab, spec, keep.row)
      else .bal.tab_long(tab, spec, keep.row, p.ops, leaf[["groups"]])
    }

    if (is_null(leaf[["keys"]]) || NROW(d) == 0L) {
      return(d)
    }

    cbind(as.data.frame(leaf[["keys"]], stringsAsFactors = FALSE), d,
          row.names = NULL)
  })

  out <- do.call("rbind", c(out, list(make.row.names = FALSE, stringsAsFactors = FALSE)))

  .drop_repeated_groups(out)
}

#A multi-category treatment is reported one pair at a time, but a mean or a standard
#deviation belongs to a group rather than to a comparison, and is the same in every pair
#that group appears in -- as is the full sample's when `pairwise = FALSE`. Reported once
#per pair it would say the same thing several times and imply it depended on the
#comparison, so it is reported once, with no pair attached.
.drop_repeated_groups <- function(d) {
  if (is_null(d) || !all(c("pair", "group", "stat") %in% names(d))) {
    return(d)
  }

  per.group <- !is.na(d[["group"]])

  if (!any(per.group)) {
    return(d)
  }

  is.na(d[["pair"]][per.group]) <- TRUE

  d <- unique(d)

  rownames(d) <- NULL

  d
}

#The rows `print()` would show, given `imbalanced.only`.
.keep_bal_rows <- function(tab, p.ops) {
  if (!isTRUE(p.ops[["imbalanced.only"]])) {
    return(rep.int(TRUE, NROW(tab)))
  }

  thr <- tab[grepl(".Threshold", names(tab), fixed = TRUE)]

  if (is_null(thr) || NCOL(thr) == 0L) {
    return(rep.int(TRUE, NROW(tab)))
  }

  rowSums(apply(thr, 2L, function(z) !is.na(z) & startsWith(z, "Not Balanced"))) > 0
}

#One row per covariate, sample, and statistic.
.bal.tab_long <- function(tab, spec, keep.row, p.ops, groups = NULL) {
  vals <- spec[spec[["quantity"]] %in% c("means", "sds", "stat"), , drop = FALSE]
  thr <- spec[spec[["quantity"]] == "threshold", , drop = FALSE]

  if (nrow(vals) == 0L || !any(keep.row)) {
    return(.empty_bal.tab_long())
  }

  d <- lapply(seq_len(nrow(vals)), function(k) {
    s <- vals[["stat"]][k]
    samp <- vals[["sample"]][k]

    #Only statistics carry thresholds; `s` is NA for a mean or an SD, and matching on
    #NA would select an NA column name rather than none.
    thr.col <- {
      if (is.na(s)) character(0L)
      else thr[["name"]][!is.na(thr[["stat"]]) & thr[["stat"]] == s &
                           thr[["sample"]] == samp]
    }

    data.frame(variable = rownames(tab)[keep.row],
               type = tab[["Type"]][keep.row],
               sample = if (samp == "Un") "Unadjusted" else samp,
               stat = switch(vals[["quantity"]][k],
                             "means" = "mean",
                             "sds" = "sd",
                             s),
               #A statistic that contrasts two groups belongs to neither of them.
               group = {
                 if (vals[["quantity"]][k] == "stat") NA_character_
                 else .group_name(vals[["group"]][k], groups)
               },
               estimate = tab[[vals[["name"]][k]]][keep.row],
               threshold = {
                 if (is_null(thr.col)) NA_character_
                 else tab[[thr.col[1L]]][keep.row]
               },
               threshold.value = {
                 if (is.na(s)) NA_real_
                 else as.numeric(p.ops[["thresholds"]][[s]] %or% NA_real_)
               },
               stringsAsFactors = FALSE)
  })

  d <- do.call("rbind", d)

  d <- d[order(match(d[["variable"]], rownames(tab))), , drop = FALSE]
  
  rownames(d) <- NULL
  
  d
}

.empty_bal.tab_long <- function() {
  data.frame(variable = character(0L),
             type = character(0L),
             sample = character(0L),
             stat = character(0L),
             group = character(0L),
             estimate = numeric(0L),
             threshold = character(0L),
             threshold.value = numeric(0L),
             stringsAsFactors = FALSE)
}

#The layout `print()` uses, with the covariate names promoted to a column.
.bal.tab_wide <- function(tab, spec, keep.row) {
  cols <- c("Type", spec[["name"]])

  data.frame(variable = rownames(tab)[keep.row],
             tab[keep.row, intersect(cols, names(tab)), drop = FALSE],
             row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
}

#' @rdname extract.bal.tab
#' @exportS3Method format bal.tab
format.bal.tab <- function(x, ..., digits = max(3L, getOption("digits") - 3L),
                           component = "balance") {
  .call <- .rewrite_all_none(match.call(expand.dots = TRUE))

  if (is_not_null(.call)) {
    return(eval.parent(.call))
  }

  component <- arg::match_arg(component, c("balance", "observations"))

  A <- .display_args(environment(), list(...))
  A[["component"]] <- NULL

  p.ops <- .resolve_p.ops(x, A)

  if (component == "observations") {
    nn <- x[["Observations"]]

    if (is_null(nn)) {
      arg::err("the {.cls bal.tab} object has no sample size table")
    }

    #A longitudinal object holds one table per time point.
    if (!is.data.frame(nn)) {
      return(lapply(nn, .format_observations, digits = p.ops[["digits"]] %or% digits))
    }

    return(.format_observations(nn, p.ops[["digits"]] %or% digits))
  }

  #The balance table `print()` shows at the top level: the summary across segments
  #when there is one, otherwise the object's own table.
  across <- names(x)[startsWith(names(x), "Balance.Across.")]
  
  tab <- {
    if (is_not_null(across) && is.data.frame(x[[across[1L]]])) x[[across[1L]]]
    else x[["Balance"]]
  }

  #Some shapes have no single top-level table -- a clustered, multiply imputed object
  #has one table per cluster and imputation and no summary across them. Give the
  #segmented tables in print()'s layout rather than nothing.
  if (is_null(tab)) {
    d <- as.data.frame(x, ..., wide = TRUE)

    num <- vapply(d, is.numeric, logical(1L))

    return(cbind(d[!num], round_df_char(d[num], p.ops[["digits"]] %or% digits,
                                        na_vals = ".")))
  }

  keep.col <- .format_keep_col(x, tab, p.ops)
  keep.row <- .keep_bal_rows(tab, p.ops)

  round_df_char(tab[keep.row, keep.col, drop = FALSE], p.ops[["digits"]] %or% digits,
                na_vals = ".")
}

#`print()` chooses the summary table's columns using the aggregating functions the
#object carries, which differ by shape.
.format_keep_col <- function(x, tab, p.ops) {
  thresholds <- setdiff(names(p.ops[["thresholds"]]), p.ops[["drop.thresholds"]])

  spec <- {
    if (inherits(x, "bal.tab.cluster"))
      .p.ops_col_spec(p.ops, p.ops[["computed.cluster.funs"]], p.ops[["requested.cluster.funs"]])
    else if (inherits(x, "bal.tab.imp"))
      .p.ops_col_spec(p.ops, p.ops[["computed.imp.funs"]], p.ops[["requested.imp.funs"]])
    else if (inherits(x, "bal.tab.multi"))
      .p.ops_col_spec(p.ops, "max")
    else if (inherits(x, "bal.tab.msm"))
      .p.ops_col_spec(p.ops, "max", include.times = TRUE)
    else if (inherits(x, "bal.tab.subclass"))
      .p.ops_col_spec(p.ops, samples = c("Un", "Adj"),
                      compute = if (isTRUE(p.ops[["quick"]])) NULL
                                else c("means", "sds",
                                       intersect(all_STATS(p.ops[["type"]]), p.ops[["stats"]])))
    else .p.ops_col_spec(p.ops)
  }

  agg.funs <- {
    if (inherits(x, "bal.tab.cluster")) p.ops[["cluster.fun"]]
    else if (inherits(x, "bal.tab.imp")) p.ops[["imp.fun"]]
    else if (inherits(x, c("bal.tab.multi", "bal.tab.msm"))) "max"
    else NULL
  }

  spec <- spec[spec[["name"]] %in% names(tab), , drop = FALSE]

  spec[["name"]][.keep_bal_cols(spec, p.ops, thresholds, agg.funs)]
}

#Matches what `.print_observations()` displays, without printing it.
.format_observations <- function(nn, digits) {
  drop.nn <- rowSums(nn) == 0
  nn <- nn[!drop.nn, , drop = FALSE]

  for (r in c("All", "Matched")) {
    rows <- paste0(r, c(" (ESS)", " (Unweighted)"))

    if (!all(rows %in% rownames(nn)) ||
        !all(check_if_zero(nn[rows[1L], ] - nn[rows[2L], ]))) {
      next
    }

    nn <- nn[rownames(nn) != rows[2L], , drop = FALSE]
    rownames(nn)[rownames(nn) == rows[1L]] <- r
  }

  round_df_char(nn, digits = min(2L, digits), pad = " ")
}

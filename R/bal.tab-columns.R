#The single source for the columns of a balance table.
#
#`balance_table()`, `balance_summary()`, and `balance_table_subclass()` lay their
#columns out the same way -- one block per sample, each block holding the requested
#moments and then each statistic followed by its threshold -- and each used to spell
#that layout out for itself. `print()` and `love.plot()` then re-derived the names
#again by position.
#
#Returns one row per column, in order, with enough metadata to select or relabel a
#column without parsing its name.
#`paste()` treats a zero-length argument as `""` rather than dropping it, so the
#absent parts of a column name have to be removed before pasting.
.paste_col <- function(...) {
  paste(unlist(list(...)), collapse = ".")
}

.bal_tab_col_spec <- function(type, compute, thresholds = list(), samples,
                              quantities = c("means", "sds"), agg.funs = NULL,
                              threshold.samples = samples,
                              threshold.agg.fun = NULL, include.times = FALSE) {
  moments <- c(M = "means", SD = "sds")
  moments <- moments[moments %in% intersect(quantities, compute)]

  stats <- intersect(compute, all_STATS(type))
  adj_only <- get_from_STATS("adj_only")

  #A continuous treatment has no groups, so its moment columns carry no group part.
  groups <- switch(type, "bin" = as.list(c("0", "1")), list(NULL))
  aggs <- {
    if (is_null(agg.funs)) list(NULL)
    else as.list(firstup(agg.funs))
  }

  #With at most one weight set the threshold column carries no sample suffix.
  bare.threshold <- length(setdiff(samples, "Un")) <= 1L

  rows <- c(list(c("Times", "times", NA, NA, NA, NA))[include.times],
            list(c("Type", "type", NA, NA, NA, NA)))

  for (sample in samples) {
    for (g in groups) {
      for (m in names(moments)) {
        rows <- c(rows, list(c(.paste_col(m, g, sample), moments[[m]], NA, sample,
                               NA, g %or% NA_character_)))
      }
    }

    for (s in stats) {
      #The unadjusted sample has no value for a statistic defined only after adjustment.
      if (sample == "Un" && adj_only[s]) next

      for (a in aggs) {
        rows <- c(rows, list(c(.paste_col(a, STATS[[s]]$bal.tab_column_prefix, sample),
                               "stat", s, sample, tolower(a) %or% NA_character_, NA)))
      }

      if (is_null(thresholds[[s]]) || sample %nin% threshold.samples) next

      rows <- c(rows, list(c(if (bare.threshold && sample != "Un") STATS[[s]]$Threshold
                             else paste.(STATS[[s]]$Threshold, sample),
                             "threshold", s, sample,
                             threshold.agg.fun %or% NA_character_, NA)))
    }
  }

  spec <- as.data.frame(do.call("rbind", rows), stringsAsFactors = FALSE)
  names(spec) <- c("name", "quantity", "stat", "sample", "agg.fun", "group")

  spec
}

#The threshold columns for one statistic, as rows of a spec.
.threshold_cols <- function(spec, s) {
  spec[spec[["quantity"]] == "threshold" & spec[["stat"]] == s, , drop = FALSE]
}

#The columns of a table belonging to a `bal.tab` being printed, derived from its
#print options. `agg.funs` names the aggregating functions a summary table carries;
#`requested.agg.funs` is what the original `bal.tab()` call asked for, which is what
#decided whether that table has threshold columns at all.
.p.ops_col_spec <- function(p.ops, agg.funs = NULL, requested.agg.funs = agg.funs,
                            samples = NULL, include.times = FALSE, compute = NULL) {
  no.adj <- !isTRUE(p.ops$disp.adj)

  #Subclassification is the adjustment, so a subclassified object records no weight
  #names but does have an adjusted sample.
  wn <- if (no.adj) "Adj" else p.ops$weight.names %or% "Adj"

  .bal_tab_col_spec(p.ops$type, compute %or% p.ops$compute, p.ops$thresholds,
                    samples = samples %or% c("Un", wn),
                    quantities = if (is_null(agg.funs)) c("means", "sds") else NULL,
                    agg.funs = agg.funs,
                    threshold.samples = {
                      if (length(requested.agg.funs) > 1L) character(0L)
                      else if (no.adj) "Un"
                      else wn
                    },
                    threshold.agg.fun = requested.agg.funs,
                    include.times = include.times)
}

#Which of a table's columns to display. Selecting by name rather than by position is
#what keeps this honest: the display rules and the column layout are decided in
#different places, and a positional vector silently misaligns when they disagree.
.keep_bal_cols <- function(spec, p.ops, thresholds, agg.funs = NULL) {
  vapply(seq_len(nrow(spec)), function(i) {
    quantity <- spec[["quantity"]][i]

    if (quantity %in% c("times", "type")) {
      return(TRUE)
    }

    #Each column belongs to the unadjusted or the adjusted sample, and is displayed
    #only if that sample is.
    if (!isTRUE(if (spec[["sample"]][i] == "Un") p.ops$un else p.ops$disp.adj)) {
      return(FALSE)
    }

    if (quantity == "threshold") {
      return(spec[["stat"]][i] %in% thresholds &&
               (is_null(agg.funs) || length(agg.funs) == 1L))
    }

    if (quantity == "stat") {
      return(spec[["stat"]][i] %in% p.ops$disp &&
               (is_null(agg.funs) || spec[["agg.fun"]][i] %in% agg.funs))
    }

    quantity %in% p.ops$disp
  }, logical(1L), USE.NAMES = FALSE) |>
    setNames(spec[["name"]])
}

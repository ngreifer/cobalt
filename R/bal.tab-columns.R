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
